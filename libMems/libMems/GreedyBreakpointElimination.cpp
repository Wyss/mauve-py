#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include "libMems/GreedyBreakpointElimination.h"
#include "libMems/ProgressiveAligner.h"
#include "libMems/Aligner.h"
#include "libMems/Islands.h"
#include "libMems/DNAFileSML.h"
#include "libMems/MuscleInterface.h"	// it's the default gapped aligner
#include "libGenome/gnRAWSource.h"
#include "libMems/gnAlignedSequences.h"
#include "libMems/CompactGappedAlignment.h"
#include "libMems/MatchProjectionAdapter.h"
#include "libMems/PairwiseMatchFinder.h"
#include "libMems/TreeUtilities.h"
#include "libMems/PairwiseMatchAdapter.h"

#include <boost/dynamic_bitset.hpp>
#include <boost/tuple/tuple.hpp>

#include <map>
#include <fstream>	// for debugging
#include <sstream>
#include <stack>
#include <algorithm>
#include <limits>
#include <iomanip>

using namespace std;
using namespace genome;

namespace mems {
// working in mems

bool penalize_repeats = false;

void printProgress( uint prev_prog, uint cur_prog, ostream& os )
{
	if( prev_prog != cur_prog )
	{
		if( cur_prog / 10 != prev_prog / 10 )
			os << endl;
		os << cur_prog << "%..";
		os.flush();
	}
}




void getPairwiseLCBs( 
	uint nI, 
	uint nJ, 
	uint dI, 
	uint dJ, 
	vector< TrackingMatch* >& tracking_matches, 
	vector< TrackingLCB<TrackingMatch*> >& t_lcbs,
	boost::multi_array< double, 3 >& tm_score_array,
	boost::multi_array< size_t, 3 >& tm_lcb_id_array )
{
	// make a set of projection matches
	vector< AbstractMatch* > pair_matches;
	for( size_t mI = 0; mI < tracking_matches.size(); ++mI )
	{
		if( tracking_matches[mI]->node_match->LeftEnd(nI) == NO_MATCH ||
			tracking_matches[mI]->node_match->LeftEnd(nJ) == NO_MATCH )
			continue;
		PairwiseMatchAdapter pma(tracking_matches[mI]->node_match, nI, nJ );
		pma.tm = tracking_matches[mI];
		if( pma.Orientation(0) == AbstractMatch::reverse )
			pma.Invert();
		pair_matches.push_back(pma.Copy());
	}
	// find LCBs...
	vector< gnSeqI > breakpoints;
	IdentifyBreakpoints( pair_matches, breakpoints );

	vector< vector< AbstractMatch* > > LCB_list;
	ComputeLCBs_v2( pair_matches, breakpoints, LCB_list );

	//
	// now compute scores on them
	//
	vector< double > lcb_scores(LCB_list.size());
	for( size_t lcbI = 0; lcbI < LCB_list.size(); ++lcbI )
	{
		double lcb_score = 0;
		for( size_t mI = 0; mI < LCB_list[lcbI].size(); ++mI )
		{
			PairwiseMatchAdapter* pma = (PairwiseMatchAdapter*)LCB_list[lcbI][mI];
			lcb_score += tm_score_array[pma->tm->match_id][dI][dJ];
		}
		lcb_scores[lcbI] = lcb_score;
	}

	// and build the pairwise adjacency list
	vector< LCB > adjacencies;
	computeLCBAdjacencies_v3( LCB_list, lcb_scores, adjacencies );

	t_lcbs.resize(adjacencies.size());
	for( size_t lcbI = 0; lcbI < adjacencies.size(); ++lcbI )
	{
		t_lcbs[lcbI] = adjacencies[lcbI];
		t_lcbs[lcbI].matches.resize(LCB_list[lcbI].size());
		for( size_t mI = 0; mI < LCB_list[lcbI].size(); ++mI )
			t_lcbs[lcbI].matches[mI] = ((PairwiseMatchAdapter*)LCB_list[lcbI][mI])->tm;
		// sort them by ptr
		sort( t_lcbs[lcbI].matches.begin(), t_lcbs[lcbI].matches.end() );

		// set the match LCB ids appropriately
		for( size_t mI = 0; mI < t_lcbs[lcbI].matches.size(); ++mI )
			tm_lcb_id_array[t_lcbs[lcbI].matches[mI]->match_id][dI][dJ] = lcbI;
	}

	// free the memory used by pairwise matches
	for( size_t mI = 0; mI < pair_matches.size(); ++mI )
		pair_matches[mI]->Free();
}

/** creates an appropriately sized matrix for mapping individual TrackingMatches to their containing LCBs */
void initTrackingMatchLCBTracking( 
	const std::vector< TrackingMatch >& tracking_matches, 
	size_t n1_count, 
	size_t n2_count, 
	boost::multi_array< size_t, 3 >& tm_lcb_id_array )
{
	tm_lcb_id_array.resize( boost::extents[tracking_matches.size()][n1_count][n2_count] );
	for( size_t mI = 0; mI < tracking_matches.size(); ++mI )
	{
		for( size_t nI = 0; nI < n1_count; ++nI )
			for( size_t nJ = 0; nJ < n2_count; ++nJ )
				tm_lcb_id_array[mI][nI][nJ] = LCB_UNASSIGNED;
	}
}


/** removes an LCB from an LCB list and coalesces surrounding LCBs.  Returns the number of LCBs removed 
 *  After LCBs are removed, the adjacency list should be processed with filterLCBs()
 *  @param	id_remaps	This is populated with a list of LCB ids that were deleted or coalesced and now have a new LCB id
 *                      for each coalesced LCB, an entry of the form <old id, new id> is added, deleted LCBs have
 *						entries of the form <deleted, -1>.  Entries appear in the order operations were performed
 *						and the function undoLcbRemoval() can undo these operations in reverse order
 */
template< class LcbVector >
uint RemoveLCBandCoalesce( size_t lcbI, uint seq_count, LcbVector& adjacencies, std::vector< double >& scores, std::vector< std::pair< uint, uint > >& id_remaps, std::vector< uint >& impact_list )
{
	uint removed_count = 0;
	vector< uint > imp_tmp(seq_count * (2 + seq_count * 4), LCB_UNASSIGNED);
	swap(impact_list, imp_tmp);
	size_t impactI = 0;
	id_remaps.clear();

	adjacencies[ lcbI ].lcb_id = -2;
	
	// update adjacencies
	uint seqI;
	uint left_adj;
	uint right_adj;
	for( seqI = 0; seqI < seq_count; seqI++ ){
		left_adj = adjacencies[ lcbI ].left_adjacency[ seqI ];
		right_adj = adjacencies[ lcbI ].right_adjacency[ seqI ];
		if( left_adj != -1 )
			adjacencies[ left_adj ].right_adjacency[ seqI ] = right_adj;
		if( right_adj != -1 && right_adj != adjacencies.size() )
			adjacencies[ right_adj ].left_adjacency[ seqI ] = left_adj;
	}

	// populate the impact list -- LCBs whose removal scores may change due to this one's removal
	for( seqI = 0; seqI < seq_count; seqI++ ){
		left_adj = adjacencies[ lcbI ].left_adjacency[ seqI ];
		right_adj = adjacencies[ lcbI ].right_adjacency[ seqI ];
		impact_list[impactI++] = left_adj;
		impact_list[impactI++] = right_adj;
		for( uint seqJ = 0; seqJ < seq_count; seqJ++ ){
			if( left_adj != -1 )
			{
				impact_list[impactI++] = adjacencies[ left_adj ].left_adjacency[ seqJ ];
				impact_list[impactI++] = adjacencies[ left_adj ].right_adjacency[ seqJ ];
			}
			if( right_adj != -1 )
			{
				impact_list[impactI++] = adjacencies[ right_adj ].left_adjacency[ seqJ ];
				impact_list[impactI++] = adjacencies[ right_adj ].right_adjacency[ seqJ ];
			}
		}
	}

	// just deleted an lcb...
	id_remaps.push_back( make_pair( lcbI, -1 ) );
	removed_count++;

	// check for collapse
	for( seqI = 0; seqI < seq_count; seqI++ ){
		left_adj = adjacencies[ lcbI ].left_adjacency[ seqI ];
		right_adj = adjacencies[ lcbI ].right_adjacency[ seqI ];
		// find the real slim shady
		while( left_adj != -1 && adjacencies[ left_adj ].lcb_id != left_adj )
			left_adj = adjacencies[ left_adj ].left_adjacency[ seqI ];
		while( right_adj != -1 && adjacencies[ right_adj ].lcb_id != right_adj )
			right_adj = adjacencies[ right_adj ].right_adjacency[ seqI ];
		if( left_adj == -1 || right_adj == -1 )
			continue;	// can't collapse with a non-existant LCB!
		if( adjacencies[ left_adj ].lcb_id != left_adj ||
			adjacencies[ right_adj ].lcb_id != right_adj )
			if( seqI > 0 )
				continue;	// already coalesced
			else
				cerr << "trouble on down street\n";

		// check whether the two LCBs are adjacent in each sequence
		boolean orientation = adjacencies[ left_adj ].left_end[ seqI ] > 0 ? true : false;
		uint seqJ;
		for( seqJ = 0; seqJ < seq_count; seqJ++ ){
			boolean j_orientation = adjacencies[ left_adj ].left_end[ seqJ ] > 0;
			if( j_orientation == orientation &&
				adjacencies[ left_adj ].right_adjacency[ seqJ ] != right_adj )
				break;
			if( j_orientation != orientation &&
				adjacencies[ left_adj ].left_adjacency[ seqJ ] != right_adj )
				break;
			// check that they are both in the same orientation
			if( adjacencies[ right_adj ].left_end[ seqJ ] > 0 != j_orientation )
				break;
		}

		if( seqJ != seq_count ||
			adjacencies[ left_adj ].to_be_deleted ||
			adjacencies[ right_adj ].to_be_deleted )
			continue;	// if these two aren't collinear, or one or both will get deleted, then don't coalesce
		

		// these two can be coalesced
		// do it.  do it now.
		id_remaps.push_back( make_pair( adjacencies[ right_adj ].lcb_id, left_adj ) );
		adjacencies[ right_adj ].lcb_id = left_adj;
		scores[ left_adj ] += scores[ right_adj ];
		adjacencies[ left_adj ].weight += adjacencies[ right_adj ].weight;

		// unlink right_adj from the adjacency list and
		// update left and right ends of left_adj
		for( seqJ = 0; seqJ < seq_count; seqJ++ ){
			boolean j_orientation = adjacencies[ left_adj ].left_end[ seqJ ] > 0;
			uint rr_adj = adjacencies[ right_adj ].right_adjacency[ seqJ ];
			uint rl_adj = adjacencies[ right_adj ].left_adjacency[ seqJ ];
			if( j_orientation == orientation ){
				adjacencies[ left_adj ].right_end[ seqJ ] = adjacencies[ right_adj ].right_end[ seqJ ];
				adjacencies[ left_adj ].right_adjacency[ seqJ ] = rr_adj;
				if( rr_adj != -1 )
					adjacencies[ rr_adj ].left_adjacency[ seqJ ] = left_adj;
			}else{
				adjacencies[ left_adj ].left_end[ seqJ ] = adjacencies[ right_adj ].left_end[ seqJ ];
				adjacencies[ left_adj ].left_adjacency[ seqJ ] = rl_adj;
				if( rl_adj != -1 )
					adjacencies[ rl_adj ].right_adjacency[ seqJ ] = left_adj;
			}
		}
		// just coalesced two LCBs...
		removed_count++;
	}
	// uniquify the impact list and get rid of empty entries
	std::sort( impact_list.begin(), impact_list.end() );
	vector< uint >::iterator imp_end = std::unique( impact_list.begin(), impact_list.end() );
	vector< uint >::iterator imp_preend = std::lower_bound( impact_list.begin(), imp_end, LCB_UNASSIGNED );
	impact_list.erase( imp_preend, impact_list.end() );

	return removed_count;
}


template< class LcbVector >
void undoLcbRemoval( uint seq_count, LcbVector& adjs, std::vector< std::pair< uint, uint > >& id_remaps )
{
	for( size_t rI = id_remaps.size(); rI > 0; --rI )
	{
		if( id_remaps[rI-1].second == -1 )
		{
			// this one was deleted
			// revert adjacencies
			uint lcbI = id_remaps[rI-1].first;
			for( uint seqI = 0; seqI < seq_count; seqI++ )
			{
				uint left_adj = adjs[ lcbI ].left_adjacency[ seqI ];
				uint right_adj = adjs[ lcbI ].right_adjacency[ seqI ];
				if( left_adj != -1 )
					adjs[ left_adj ].right_adjacency[ seqI ] = lcbI;
				if( right_adj != -1 && right_adj != adjs.size() )
					adjs[ right_adj ].left_adjacency[ seqI ] = lcbI;
			}
			adjs[lcbI].lcb_id = lcbI;	// reset the lcb id
			adjs[lcbI].to_be_deleted = false;	// no longer TBD
		}else{
			// this one was coalesced
			// uncoalesce it
			uint lcbI = id_remaps[rI-1].first;
			uint lcbJ = id_remaps[rI-1].second;
			adjs[lcbI].lcb_id = lcbI;
			adjs[lcbJ].weight -= adjs[lcbI].weight;
			// link lcbI back in
			// TODO: fix right end and left end coordinates
			for( uint seqI = 0; seqI < seq_count; ++seqI )
			{
				uint ladj = adjs[lcbI].left_adjacency[seqI];
				uint radj = adjs[lcbI].right_adjacency[seqI];
				if(  ladj == lcbJ )
				{
					adjs[lcbJ].right_adjacency[seqI] = lcbI;
					if( radj != -1 && radj != adjs.size())
						adjs[radj].left_adjacency[seqI] = lcbI;
				}else
				if(  radj == lcbJ )
				{
					adjs[lcbJ].left_adjacency[seqI] = lcbI;
					if( ladj != -1 && ladj != adjs.size())
						adjs[ladj].right_adjacency[seqI] = lcbI;
				}
			}
		}
	}
}

EvenFasterSumOfPairsBreakpointScorer::EvenFasterSumOfPairsBreakpointScorer( 
	double breakpoint_penalty,
	double minimum_breakpoint_penalty,
	boost::multi_array<double,2> bp_weight_matrix, 
	boost::multi_array<double,2> conservation_weight_matrix,
	vector< TrackingMatch* > tracking_match,
	PairwiseLCBMatrix& pairwise_adjacency_matrix,
	vector<node_id_t>& n1_descendants,
	vector<node_id_t>& n2_descendants,
	boost::multi_array< double, 3 >& tm_score_array,
	boost::multi_array< size_t, 3 >& tm_lcb_id_array,
	size_t seqI_begin,
	size_t seqI_end,
	size_t seqJ_begin,
	size_t seqJ_end
	) :
  bp_penalty( breakpoint_penalty ),
  min_breakpoint_penalty( minimum_breakpoint_penalty ),
  bp_weights( bp_weight_matrix ), 
  conservation_weights( conservation_weight_matrix ),
  tracking_matches( tracking_match ),
  pairwise_adjacencies( pairwise_adjacency_matrix ),
  n1_des(n1_descendants),
  n2_des(n2_descendants),
  tm_score_array(tm_score_array),
  tm_lcb_id_array(tm_lcb_id_array),
  seqI_count(pairwise_adjacencies.shape()[0]),
  seqJ_count(pairwise_adjacencies.shape()[1]),
  seqI_first(seqI_begin),
  seqI_last(seqI_end),
  seqJ_first(seqJ_begin),
  seqJ_last(seqJ_end),
  first_time(true)
{
	std::sort(tracking_matches.begin(), tracking_matches.end());
	pairwise_lcb_count.resize( boost::extents[pairwise_adjacencies.shape()[0]][pairwise_adjacencies.shape()[1]] );
	pairwise_lcb_score.resize( boost::extents[pairwise_adjacencies.shape()[0]][pairwise_adjacencies.shape()[1]] );;
	all_id_remaps.resize( boost::extents[pairwise_lcb_count.shape()[0]][pairwise_lcb_count.shape()[1]] );
	full_impact_list.resize( boost::extents[pairwise_lcb_count.shape()[0]][pairwise_lcb_count.shape()[1]] );
	my_del_lcbs.resize(100);	// buffer for use during lcb removal score computation
	for( size_t i = 0; i < 3; ++i )
	{
		internal_lcb_score_diff[i].resize( boost::extents[pairwise_adjacencies.shape()[0]][pairwise_adjacencies.shape()[1]] );
		internal_lcb_removed_count[i].resize( boost::extents[pairwise_adjacencies.shape()[0]][pairwise_adjacencies.shape()[1]] );
	}
	lsd_zeros.resize( internal_lcb_score_diff[0].num_elements(), 0 );
	lrc_zeros.resize( internal_lcb_removed_count[0].num_elements(), 0 );
	using_lsd = -1;
	size_t max_pair_adj_size = 0;
	for( size_t i = 0; i < seqI_count; ++i )
	{
		for( size_t j = 0; j < seqJ_count; ++j )
		{
			pairwise_lcb_count[i][j] = pairwise_adjacencies[i][j].size();
			pairwise_lcb_score[i][j] = 0;
			max_pair_adj_size = (std::max)(max_pair_adj_size, pairwise_adjacencies[i][j].size());
			for( size_t lcbI = 0; lcbI < pairwise_adjacencies[i][j].size(); ++lcbI )
				pairwise_lcb_score[i][j] += pairwise_adjacencies[i][j][lcbI].weight;
		}
	}
	bogus_scores.resize(max_pair_adj_size+10);
};


/**
 * Returns the number of possible moves a search algorithm may make from the current 
 * location in LCB search space.  In this case it's simply the total number of pairwise LCBs
 */
size_t EvenFasterSumOfPairsBreakpointScorer::getMoveCount()
{
	size_t move_count = 0;
	for( size_t i = seqI_first; i < seqI_last; ++i )
		for( size_t j = seqJ_first; j < seqJ_last; ++j )
			move_count += pairwise_adjacencies[i][j].size();
	return move_count;
}

/** returns the score of the current state */
double EvenFasterSumOfPairsBreakpointScorer::score()
{
	// score is the sum of all pairwise LCB scores,
	// minus the sum of all pairwise breakpoint penalties
	double score = 0;
	for( size_t seqI = seqI_first; seqI < seqI_last; ++seqI )
	{
		for( size_t seqJ = seqJ_first; seqJ < seqJ_last; ++seqJ )
		{
			const double pw_lcb_score = pairwise_lcb_score[seqI][seqJ];
			// add LCB scores
			score += pairwise_lcb_score[seqI][seqJ];
			// subtract breakpoint penalty
			// subtract 1 from number of LCBs so that a single circular LCB doesn't get penalized
			double cweights = 1 - conservation_weights[seqI][seqJ];
			double bweights = 1 - bp_weights[seqI][seqJ];
			double penalty = max( bp_penalty * cweights * cweights * cweights * cweights * bweights * bweights, min_breakpoint_penalty );
			if(first_time)
				cout << "Scoring with scaled breakpoint penalty: " << penalty << endl;
			first_time = false;
			score -= ( penalty * (pairwise_lcb_count[seqI][seqJ]-1));
			if( !(score > -1e200 && score < 1e200) )
			{
				genome::breakHere();
				cerr << "bp_weights[seqI][seqJ] " << bp_weights[seqI][seqJ] << endl;
				cerr << "conservation_weights[seqI][seqJ] " << conservation_weights[seqI][seqJ] << endl;
				cerr << "pairwise_lcb_count[seqI][seqJ] " << pairwise_lcb_count[seqI][seqJ] << endl;
				cerr << "pairwise_lcb_score[seqI][seqJ] " << pw_lcb_score << endl;
				cerr << "Invalid score!!\n";
			}
		}
	}
	return score;
}

/** scores a move */
double EvenFasterSumOfPairsBreakpointScorer::operator()( pair< double, size_t >& the_move  )
{
	size_t new_move_count;
	vector< pair< double, size_t > > new_move_list;
	using_lsd++;
	std::copy(lsd_zeros.begin(),lsd_zeros.end(),internal_lcb_score_diff[using_lsd].data());
	std::copy(lrc_zeros.begin(),lrc_zeros.end(),internal_lcb_removed_count[using_lsd].data());
	remove( the_move, false, internal_lcb_score_diff[using_lsd], internal_lcb_removed_count[using_lsd], false, new_move_list, new_move_count );
	applyScoreDifference( internal_lcb_score_diff[using_lsd], internal_lcb_removed_count[using_lsd] );
	double m_score = score();
	undoScoreDifference( internal_lcb_score_diff[using_lsd], internal_lcb_removed_count[using_lsd] );
	using_lsd--;
	return m_score;
}

bool EvenFasterSumOfPairsBreakpointScorer::isValid( pair< double, size_t >& the_move )
{
	using_lsd++;
	std::copy(lsd_zeros.begin(),lsd_zeros.end(),internal_lcb_score_diff[using_lsd].data());
	std::copy(lrc_zeros.begin(),lrc_zeros.end(),internal_lcb_removed_count[using_lsd].data());
	vector< pair< double, size_t > > new_move_list;
	size_t new_move_count;
	bool success = remove( the_move, false, internal_lcb_score_diff[using_lsd], internal_lcb_removed_count[using_lsd], false, new_move_list, new_move_count );
	using_lsd--;
	return success;
}

bool EvenFasterSumOfPairsBreakpointScorer::remove( pair< double, size_t >& the_move, vector< pair< double, size_t > >& new_move_list, size_t& new_move_count )
{
	using_lsd++;
	std::copy(lsd_zeros.begin(),lsd_zeros.end(),internal_lcb_score_diff[using_lsd].data());
	std::copy(lrc_zeros.begin(),lrc_zeros.end(),internal_lcb_removed_count[using_lsd].data());
	bool success = remove( the_move, true, internal_lcb_score_diff[using_lsd], internal_lcb_removed_count[using_lsd], true, new_move_list, new_move_count );
	if( success )
		applyScoreDifference( internal_lcb_score_diff[using_lsd], internal_lcb_removed_count[using_lsd] );
	using_lsd--;
	return success;
}

void EvenFasterSumOfPairsBreakpointScorer::applyScoreDifference( boost::multi_array< double, 2 >& lcb_score_diff, boost::multi_array< size_t, 2 >& lcb_removed_count )
{
	size_t nelems = pairwise_lcb_count.num_elements();
	for( size_t elemI = 0; elemI < nelems; elemI++ )
	{
		if( !(lcb_score_diff.data()[elemI] > -1e200 && lcb_score_diff.data()[elemI] < 1e200) )
		{
			genome::breakHere();
			cerr << "Invalid score!!\n";
		}
		pairwise_lcb_count.data()[elemI] -= lcb_removed_count.data()[elemI];
		pairwise_lcb_score.data()[elemI] -= lcb_score_diff.data()[elemI];
		if( !(pairwise_lcb_score.data()[elemI] > -1e200 && pairwise_lcb_score.data()[elemI] < 1e200) )
		{
			genome::breakHere();
			cerr << "Invalid score!!\n";
		}
	}
}

void EvenFasterSumOfPairsBreakpointScorer::undoScoreDifference( boost::multi_array< double, 2 >& lcb_score_diff, boost::multi_array< size_t, 2 >& lcb_removed_count )
{
	size_t nelems = pairwise_lcb_count.num_elements();
	for( size_t elemI = 0; elemI < nelems; elemI++ )
	{
		if( !(lcb_score_diff.data()[elemI] > -1e200 && lcb_score_diff.data()[elemI] < 1e200) )
		{
			genome::breakHere();
			cerr << "Invalid score!!\n";
		}
		pairwise_lcb_count.data()[elemI] += lcb_removed_count.data()[elemI];
		pairwise_lcb_score.data()[elemI] += lcb_score_diff.data()[elemI];
		if( !(pairwise_lcb_score.data()[elemI] > -1e200 && pairwise_lcb_score.data()[elemI] < 1e200) )
		{
			genome::breakHere();
			cerr << "Invalid score!!\n";
		}
	}
}

size_t EvenFasterSumOfPairsBreakpointScorer::getMaxNewMoveCount()
{
	return 20 * seqI_count * seqJ_count;
}

/** call to indicate that the given LCB has been removed 
  * returns false if the move was invalid
  */
bool EvenFasterSumOfPairsBreakpointScorer::remove( pair< double, size_t >& the_move, bool really_remove, boost::multi_array< double, 2 >& lcb_score_diff, boost::multi_array< size_t, 2 >& lcb_removed_count, bool score_new_moves, vector< pair< double, size_t > >& new_move_list, size_t& new_move_count )
{
	if( score_new_moves && !really_remove )
	{
		cerr << "Error: Incompatible options in the breakpoint scorer!!!\n";
		throw "oh shit!";
	}
	new_move_count = 0;
	// figure out which lcb we're being asked to delete
	size_t moveI = the_move.second;
	size_t move_count = 0;
	size_t move_base = 0;
	size_t seqI = 0;
	size_t seqJ = 0;
	for( seqI = seqI_first; seqI < seqI_last; ++seqI )
	{
		for( seqJ = seqJ_first; seqJ < seqJ_last; ++seqJ )
		{
			all_id_remaps[seqI][seqJ].clear();
			full_impact_list[seqI][seqJ].clear();
		}
	}

	for( seqI = seqI_first; seqI < seqI_last; ++seqI )
	{
		for( seqJ = seqJ_first; seqJ < seqJ_last; ++seqJ )
		{
			move_count += pairwise_adjacencies[seqI][seqJ].size();
			if( move_count > moveI )
				break;
			move_base = move_count;
		}
		if( move_count > moveI )
			break;
	}
	// score deletion of the LCB at (moveI - move_base) from the pairwise alignment of seqI and seqJ
	size_t del_lcb = moveI - move_base;
	if( pairwise_adjacencies[seqI][seqJ][del_lcb].lcb_id != del_lcb && really_remove )
	{
		if( pairwise_adjacencies[seqI][seqJ][del_lcb].lcb_id == LCB_UNASSIGNED )
			cerr << "bad movement, dirty dancing\n";
		return false;	// this is an invalid move -- already deleted or coalesced with another
	}
	if( pairwise_adjacencies[seqI][seqJ][del_lcb].lcb_id != del_lcb )
	{
		return false;	// this is an invalid move -- already deleted
	}
	
	vector< TrackingMatch* > matches(pairwise_adjacencies[seqI][seqJ][del_lcb].matches);
	double cur_score = score();

	if( really_remove )
	{
		deleted_tracking_matches.insert( deleted_tracking_matches.end(), matches.begin(), matches.end() );
	}

	for( size_t i = seqI_first; i < seqI_last; ++i )
	{
		for( size_t j = seqJ_first; j < seqJ_last; ++j )
		{
			lcb_score_diff[i][j] = 0;
			vector< TrackingLCB< TrackingMatch* > >& adjs = pairwise_adjacencies[i][j];
			// create a list of LCBs affected by deletion of this match
			// check whether any of them will have all of their matches removed
			if( lcb_ids.size() < matches.size() )
				lcb_ids.resize( matches.size() + 100 );
			for( size_t mI = 0; mI < matches.size(); ++mI )
				lcb_ids[mI] = tm_lcb_id_array[matches[mI]->match_id][i][j];
			size_t lcb_id_count = matches.size();
			std::sort(lcb_ids.begin(), lcb_ids.begin()+lcb_id_count);
			vector< size_t >::iterator last = std::unique(lcb_ids.begin(), lcb_ids.begin()+lcb_id_count);
			lcb_id_count = last - lcb_ids.begin();
			// delete the last one if its unassigned
			if( lcb_ids[lcb_id_count-1] == LCB_UNASSIGNED )
				lcb_id_count--;

			vector< pair< size_t, vector< TrackingMatch* > > > aff_lcbs(lcb_id_count);
			for( size_t lI = 0; lI < lcb_id_count; ++lI )
				aff_lcbs[lI].first = lcb_ids[lI];

			// organize the deleted matches
			for( size_t mI = 0; mI < matches.size(); ++mI )
			{
				size_t id = tm_lcb_id_array[matches[mI]->match_id][i][j];
				if( id == LCB_UNASSIGNED )
					continue;
				vector< pair< size_t, vector< TrackingMatch* > > >::iterator iter = std::lower_bound( aff_lcbs.begin(), aff_lcbs.end(), make_pair(id,vector< TrackingMatch* >() ) );
				iter->second.push_back( matches[mI] );
			}

			// actually delete the matches and keep a list of LCBs that get completely deleted
			size_t my_del_count = 0;
			for( size_t lI = 0; lI < aff_lcbs.size(); ++lI )
			{
				vector< TrackingMatch* >& cur_matches = adjs[lcb_ids[lI]].matches;
				size_t diff = cur_matches.size() - aff_lcbs[lI].second.size();
				if( diff == 0 )
				{
					if( my_del_count + 1 >= my_del_lcbs.size() )
						my_del_lcbs.resize(2*my_del_lcbs.size());
					my_del_lcbs[my_del_count++] = lcb_ids[lI];
					adjs[lcb_ids[lI]].to_be_deleted = true;
					lcb_score_diff[i][j] += adjs[lcb_ids[lI]].weight;
					if( really_remove )
					{
						adjs[lcb_ids[lI]].weight = 0;
						cur_matches.clear();
					}
					continue;
				}

				// update the LCB score
				double del_score_sum = 0;
				for( size_t mI = 0; mI < aff_lcbs[lI].second.size(); ++mI )
					del_score_sum += tm_score_array[aff_lcbs[lI].second[mI]->match_id][i][j];
				lcb_score_diff[i][j] += del_score_sum;
				full_impact_list[i][j].push_back( aff_lcbs[lI].first );

				if( really_remove )
				{
					adjs[lcb_ids[lI]].weight -= del_score_sum;
				
					// remove the deleted matches
					vector< TrackingMatch* > dest( diff );
					std::set_difference( cur_matches.begin(), cur_matches.end(), 
						aff_lcbs[lI].second.begin(), aff_lcbs[lI].second.end(), dest.begin() );
					swap( dest, cur_matches );
				}
			}

			lcb_removed_count[i][j] = 0;

			// now remove each LCB that needs to be deleted
			std::vector< std::pair< uint, uint > >& fid_remaps = all_id_remaps[i][j];
			std::vector< uint >& fimp_list = full_impact_list[i][j];
			for( size_t delI = 0; delI < my_del_count; ++delI )
			{
				if( adjs[my_del_lcbs[delI]].lcb_id != my_del_lcbs[delI] )
					continue;	// skip this one if it's already been deleted

				std::vector< std::pair< uint, uint > > id_remaps;
				std::vector< uint > impact_list;
				uint removed_count = RemoveLCBandCoalesce( my_del_lcbs[delI], 2, adjs, bogus_scores, id_remaps, impact_list );
				fid_remaps.insert( fid_remaps.end(), id_remaps.begin(), id_remaps.end() );
				fimp_list.insert( fimp_list.end(), impact_list.begin(), impact_list.end() );

				lcb_removed_count[i][j] += removed_count;
				// only do this part if we're really deleting
				if( really_remove )
				{
					// move all matches to the new LCB
					for( size_t rI = 0; rI < id_remaps.size(); ++rI )
					{
						if( id_remaps[rI].second == -1 )
							continue;	// deletion
						vector< TrackingMatch* >& src_matches = adjs[id_remaps[rI].first].matches;
						vector< TrackingMatch* >& dest_matches = adjs[id_remaps[rI].second].matches;
						for( size_t mI = 0; mI < src_matches.size(); ++mI )
							tm_lcb_id_array[src_matches[mI]->match_id][i][j] = id_remaps[rI].second;
						dest_matches.insert( dest_matches.end(), src_matches.begin(), src_matches.end() );
						std::sort( dest_matches.begin(), dest_matches.end() );
						src_matches.clear();
					}
				}
			}
		}
	}

	// will be undone later
	applyScoreDifference( lcb_score_diff, lcb_removed_count );
	double new_score = score();

	if( score_new_moves )
	{
		size_t mbase = 0;
		for( size_t i = seqI_first; i < seqI_last; ++i )
		{
			for( size_t j = seqJ_first; j < seqJ_last; ++j )
			{
				vector< TrackingLCB< TrackingMatch* > >& adjs = pairwise_adjacencies[i][j];
				std::vector< uint >& fimp_list = full_impact_list[i][j];
				sort( fimp_list.begin(), fimp_list.end() );
				vector< uint >::iterator iter = std::unique( fimp_list.begin(), fimp_list.end() );
				fimp_list.erase( iter, fimp_list.end() );
				for( size_t fI = 0; fI < fimp_list.size(); fI++ )
				{
					if( adjs[fimp_list[fI]].lcb_id != fimp_list[fI] )
					{
						new_move_list[new_move_count++] = make_pair( -(std::numeric_limits<double>::max)(), mbase + fimp_list[fI] );
						continue;	// this one got trashed
					}
					// score removal of this block
					pair< double, size_t > p( 0, mbase + fimp_list[fI] );
					double scorediff = (*this)(p) - new_score;
					p.first = scorediff;
					new_move_list[new_move_count++] = p;
				}
				mbase += adjs.size();
			}
		}
	}


	// if we're not really removing, undo all the removals
	if( !really_remove )
		for( size_t i = seqI_first; i < seqI_last; ++i )
			for( size_t j = seqJ_first; j < seqJ_last; ++j )
				undoLcbRemoval( 2, pairwise_adjacencies[i][j], all_id_remaps[i][j] );

	undoScoreDifference( lcb_score_diff, lcb_removed_count );

	// if the change in score doesn't match then this is an invalid move!!
	// allow for some numerical instability
	bool valid = true;
	if( new_score - cur_score < the_move.first - 0.00001 ||
		new_score - cur_score  > the_move.first + 0.00001 )
		valid = false;

	return valid;
}

vector< TrackingMatch* > EvenFasterSumOfPairsBreakpointScorer::getResults() 
{
	std::sort(deleted_tracking_matches.begin(), deleted_tracking_matches.end());
	vector< TrackingMatch* > result_matches(tracking_matches.size()-deleted_tracking_matches.size());
	std::set_difference( tracking_matches.begin(), tracking_matches.end(), deleted_tracking_matches.begin(), deleted_tracking_matches.end(), result_matches.begin() );
	return result_matches;
}

	bool EvenFasterSumOfPairsBreakpointScorer::validate()
{
	vector< TrackingMatch* > trams = getResults();	// need to apply any deletions...
	bool success = true;	// be optimistic!
	// make sure all the tracking matches point to the right LCBs
	for( size_t tmI = 0; tmI < trams.size(); tmI++ )
	{
		TrackingMatch* tm = trams[tmI];
		for( size_t i = 0; i < tm_lcb_id_array.shape()[1]; ++i )
			for( size_t j = 0; j < tm_lcb_id_array.shape()[2]; ++j )
			{
				// skip this match if it's not defined
				if( tm->node_match->LeftEnd(n1_des[i]) == NO_MATCH ||
					tm->node_match->LeftEnd(n2_des[j]) == NO_MATCH ||
					tm_lcb_id_array[tm->match_id][i][j] == LCB_UNASSIGNED)
					continue;
				// find the tracking match in this LCB
				size_t id = tm_lcb_id_array[tm->match_id][i][j];
				vector< TrackingMatch* >& matches = pairwise_adjacencies[i][j][id].matches;
				vector< TrackingMatch* >::iterator iter = std::lower_bound( matches.begin(), matches.end(), tm );
				if( iter == matches.end() || *iter != tm )
				{
					cerr << "Missing match!!\n";
					cerr << "lcb_id: " << id << endl;
					cerr << "match: " << tm << endl;
					genome::breakHere();
					success = false;
				}
			}
	}
	// make sure all the LCBs point to valid tracking matches
	for( size_t i = 0; i < pairwise_adjacencies.shape()[0]; ++i )
		for( size_t j = 0; j < pairwise_adjacencies.shape()[1]; ++j )
		{
			vector< TrackingLCB< TrackingMatch* > >& adjs = pairwise_adjacencies[i][j];
			for( size_t lcbI = 0; lcbI < adjs.size(); lcbI++ )
			{
				for( size_t mI = 0; mI < adjs[lcbI].matches.size(); ++mI )
				{
					vector< TrackingMatch* >::iterator iter = std::lower_bound( trams.begin(), trams.end(), adjs[lcbI].matches[mI] );
					if( *iter != adjs[lcbI].matches[mI] )
					{
						cerr << "Missing match:  in adjacencies but not tracking_matches!!\n";
						cerr << "lcb_id: " << tm_lcb_id_array[adjs[lcbI].matches[mI]->match_id][i][j] << endl;
						genome::breakHere();
						success = false;
					}
				}
			}
		}

	// make sure that the number of breakpoints matches up with what tracking_matches suggests
	vector< TrackingMatch* > final = trams;
	// convert back to an LCB list
	vector< AbstractMatch* > new_matches(final.size());
	for( size_t mI = 0; mI < final.size(); ++mI )
		new_matches[mI] = final[mI]->original_match;

	vector< gnSeqI > breakpoints;
	IdentifyBreakpoints( new_matches, breakpoints );
	vector< vector< AbstractMatch* > > LCB_list;
	IdentifyBreakpoints( new_matches, breakpoints );
	ComputeLCBs_v2( new_matches, breakpoints, LCB_list );
	cout << "breakpoints.size(): " << breakpoints.size() << "\tpairwise_lcb_count[0][0]: " << pairwise_lcb_count[0][0] << endl;
	if( breakpoints.size() != pairwise_lcb_count[0][0] )
		success = false;
	size_t adjI = 0;
	vector< TrackingLCB< TrackingMatch* > >& adjs = pairwise_adjacencies[0][0];
	for( size_t lcbI = 0; lcbI < LCB_list.size(); ++lcbI )
	{
		// make sure each LCB exists...
		while( adjI != -1 && adjI != adjs[adjI].lcb_id )
			adjI++;

		// compare matches...
		vector< AbstractMatch* > ms(adjs[adjI].matches.size()+LCB_list[lcbI].size(), (AbstractMatch*)NULL);
		std::sort( LCB_list[lcbI].begin(), LCB_list[lcbI].end() );
		vector< AbstractMatch* > asdf(adjs[adjI].matches.size());
		for( size_t mI = 0; mI < adjs[adjI].matches.size(); ++mI )
			asdf[mI] = adjs[adjI].matches[mI]->original_match;
		std::sort( asdf.begin(), asdf.end() );
		std::set_symmetric_difference( LCB_list[lcbI].begin(), LCB_list[lcbI].end(), asdf.begin(), asdf.end(), ms.begin() );
		// this should throw a fit if the sets aren't equal.
		if( ms[0] != NULL )
		{
			cerr << "In adjacencies:\n";
			for( size_t asdfI = 0; asdfI < asdf.size(); asdfI++ )
			{
				printMatch(asdf[asdfI], cerr);
				cerr << endl;
			}
			cerr << "\nIn LCB_list:\n";
			for( size_t mI = 0; mI < LCB_list[lcbI].size(); mI++ )
			{
				printMatch(LCB_list[lcbI][mI], cerr);
				cerr << endl;
			}
			cerr << "\nAll matches ssc1\n";
			SingleStartComparator<AbstractMatch> ssc1(1);
			std::sort(new_matches.begin(), new_matches.end(), ssc1);
			for( size_t mI = 0; mI < new_matches.size(); mI++ )
			{
				printMatch(new_matches[mI], cerr);
				cerr << endl;
			}

			cerr << "\nAll matches ssc0\n";
			SingleStartComparator<AbstractMatch> ssc0(0);
			std::sort(new_matches.begin(), new_matches.end(), ssc0);
			for( size_t mI = 0; mI < new_matches.size(); mI++ )
			{
				printMatch(new_matches[mI], cerr);
				cerr << endl;
			}
			genome::breakHere();
		}
		adjI++;
	}

	return success;
}



SimpleBreakpointScorer::SimpleBreakpointScorer( std::vector< LCB >& adjacencies, double breakpoint_penalty, bool collinear ) : 
  adjs( adjacencies ),
  bp_penalty( breakpoint_penalty ),
  collinear( collinear )
{
	scores = std::vector< double >(adjs.size(), 0);
	total_weight = 0;
	bp_count = adjs.size();
	for( size_t lcbI = 0; lcbI < adjs.size(); lcbI++ )
		total_weight += adjs[lcbI].weight;
}

size_t SimpleBreakpointScorer::getMoveCount() 
{
	return adjs.size();
}

double SimpleBreakpointScorer::score()
{
	double bp_score = (double)bp_count * bp_penalty;
	return total_weight - bp_score;
}

bool SimpleBreakpointScorer::isValid( size_t lcbI, double move_score )
{
	if( adjs[lcbI].lcb_id != lcbI )
		return false;
	return (*this)(lcbI) == move_score;
}

/** return the relative change in score if lcbI were to be removed */
double SimpleBreakpointScorer::operator()( size_t lcbI )
{
	double cur_score = score();
	std::vector< std::pair< uint, uint > > id_remaps;
	std::vector< uint > impact_list;
	uint bp_removed = RemoveLCBandCoalesce( lcbI, adjs[0].left_adjacency.size(), adjs, scores, id_remaps, impact_list );
	undoLcbRemoval( adjs[0].left_adjacency.size(), adjs, id_remaps );
	double bp_score = (double)(bp_count - bp_removed) * bp_penalty;
	double move_score = total_weight - adjs[lcbI].weight - bp_score;
	double score_diff = move_score - cur_score;
	if( collinear && bp_count - bp_removed > 0 && score_diff < 0 )
		return 1/(-score_diff);	// ensure that we continue removing blocks until only one is left
	return move_score - cur_score;
}

/** call to indicate that the given LCB has been removed */
void SimpleBreakpointScorer::remove( uint lcbI, vector< pair< double, size_t > >& new_moves )
{
	std::vector< std::pair< uint, uint > > id_remaps;
	std::vector< uint > impact_list;
	uint bp_removed = RemoveLCBandCoalesce( lcbI, adjs[0].left_adjacency.size(), adjs, scores, id_remaps, impact_list );
	total_weight -= adjs[lcbI].weight;
	bp_count -= bp_removed;
	for( size_t impI = 0; impI < impact_list.size(); impI++ )
	{
		if( adjs[impact_list[impI]].lcb_id != impact_list[impI] )
			continue;
		double scorediff = (*this)(impact_list[impI]);
		new_moves.push_back(make_pair(scorediff, impact_list[impI]));
	}
}


GreedyRemovalScorer::GreedyRemovalScorer( std::vector< LCB >& adjacencies, double minimum_weight ) : 
adjs( adjacencies ),
min_weight( minimum_weight )
{
	scores = std::vector< double >(adjs.size(), 0);
	total_weight = 0;
	for( size_t lcbI = 0; lcbI < adjs.size(); lcbI++ )
		total_weight += adjs[lcbI].weight - min_weight;
}

size_t GreedyRemovalScorer::getMoveCount() 
{
	return adjs.size();
}

double GreedyRemovalScorer::score()
{
	return total_weight;
}

bool GreedyRemovalScorer::isValid( size_t lcbI, double move_score )
{
	if( adjs[lcbI].lcb_id != lcbI )
		return false;
	return (*this)(lcbI) == move_score;
}

/** return the relative change in score if lcbI were to be removed */
double GreedyRemovalScorer::operator()( size_t lcbI )
{
	return -(adjs[lcbI].weight-min_weight);
}

/** call to indicate that the given LCB has been removed */
void GreedyRemovalScorer::remove( uint lcbI, vector< pair< double, size_t > >& new_moves )
{
	std::vector< std::pair< uint, uint > > id_remaps;
	std::vector< uint > impact_list;
	uint bp_removed = RemoveLCBandCoalesce( lcbI, adjs[0].left_adjacency.size(), adjs, scores, id_remaps, impact_list );
	total_weight -= (adjs[lcbI].weight-min_weight);
	for( size_t impI = 0; impI < impact_list.size(); impI++ )
	{
		if( adjs[impact_list[impI]].lcb_id != impact_list[impI] )
			continue;
		double scorediff = (*this)(impact_list[impI]);
		new_moves.push_back(make_pair(scorediff, impact_list[impI]));
	}
}




}	// namespace mems

