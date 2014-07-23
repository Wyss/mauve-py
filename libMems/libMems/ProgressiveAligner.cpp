/*******************************************************************************
 * $Id: progressiveAligner.cpp,v 1.47 2004/04/19 23:10:30 darling Exp $
 * BEWARE!!
 * This code was created in the likeness of the flying spaghetti monster
 *
 * dedicated to Loren...
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include "libMems/ProgressiveAligner.h"
#include "libMems/GreedyBreakpointElimination.h"
#include "libMems/Aligner.h"
#include "libMems/Islands.h"
#include "libMems/DNAFileSML.h"
#include "libMems/MuscleInterface.h"	// it's the default gapped aligner
#include "libMems/gnAlignedSequences.h"
#include "libMems/CompactGappedAlignment.h"
#include "libMems/MatchProjectionAdapter.h"
#include "libMems/PairwiseMatchFinder.h"
#include "libMems/TreeUtilities.h"
#include "libMems/PairwiseMatchAdapter.h"
#include "libMems/DistanceMatrix.h"

#include <boost/dynamic_bitset.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/graph/undirected_dfs.hpp>

#include <map>
#include <fstream>	// for debugging
#include <sstream>
#include <stack>
#include <algorithm>
#include <limits>
#include <iomanip>

#include "stdlib.h"

using namespace std;
using namespace genome;

namespace mems {


bool progress_msgs = false;

bool debug_me = false;
static int dbg_count = 0; 	 


double min_window_size = 200;
double max_window_size = 20000;  // don't feed MUSCLE anything bigger than this
double min_density = .5;
double max_density = .9;
size_t max_gap_length = 5000;
size_t lcb_hangover = 300;


void mergeUnalignedIntervals( uint seqI, vector< Interval* >& iv_list, vector< Interval* >& new_list );

/**
 * Test code to ensure that an individual LCB is truly collinear
 * @return	true if the LCB is good
 */
boolean my_validateLCB( MatchList& lcb ){
	vector< Match* >::iterator lcb_iter = lcb.begin();
	if( lcb.size() == 0 )
		return true;
	uint seq_count = (*lcb_iter)->SeqCount();
	uint seqI = 0;
	boolean complain = false;
	for(; seqI < seq_count; seqI++ ){
		lcb_iter = lcb.begin();
		int64 prev_coord = 0;
		for(; lcb_iter != lcb.end(); ++lcb_iter ){
			if( (*lcb_iter)->Start( seqI ) == NO_MATCH )
				continue;
			else if( prev_coord != 0 && (*lcb_iter)->Start( seqI ) < prev_coord ){
				complain = true;
			}
			prev_coord = (*lcb_iter)->Start( seqI );
		}
	}
	return !complain;
}

template< class BoostMatType >
void print2d_matrix( BoostMatType& mat, std::ostream& os )
{
	for( size_t i = 0; i < mat.shape()[0]; ++i )
	{
		for( size_t j = 0; j < mat.shape()[1]; ++j )
		{
			if( j > 0 )
				os << "\t";
			os << mat[i][j];
		}
		os << endl;
	}
}

double getDefaultBreakpointPenalty( std::vector< gnSequence* >& sequences )
{
	uint default_mer_size = MatchList::GetDefaultMerSize( sequences );
	double avg_seq_len = 0;
	for( size_t seqI = 0; seqI < sequences.size(); ++seqI )
		avg_seq_len += (double)sequences[seqI]->length();
	avg_seq_len /= (double)sequences.size();
	avg_seq_len = log( avg_seq_len ) / log( 2.0 );
	return avg_seq_len * 7000;	  // seems to work reasonably well?
}


double getDefaultBpDistEstimateMinScore( std::vector< gnSequence* >& sequences )
{
	// this value was empirically derived by a process that involved burning incense
	// and uttering arcane words
	return 3.0 * getDefaultBreakpointPenalty(sequences);
}



/*
 * A progressive alignment algorithm for genomes with rearrangements.
 * Start simple, add complexity later.
 * TODO: rewrite the algorithm outline
 */

ProgressiveAligner::ProgressiveAligner( uint seq_count ) :
Aligner( seq_count ),
breakpoint_penalty( -1 ),
min_breakpoint_penalty( 4000 ),
debug(false),
refine(true),
scoring_scheme(ExtantSumOfPairsScoring),
use_weight_scaling(true),
conservation_dist_scale(1),
bp_dist_scale(.9),
max_gapped_alignment_length(20000),
bp_dist_estimate_score(-1),
use_seed_families(false),
using_cache_db(true)
{
	gapped_alignment = true;
	max_window_size = max_gapped_alignment_length;
}

void ProgressiveAligner::SetMaxGappedAlignmentLength( size_t len )
{ 
	max_gapped_alignment_length = len; 
	max_window_size = max_gapped_alignment_length;
}

/** determine which extant sequences have been aligned at a given node */
void ProgressiveAligner::getAlignedChildren( node_id_t node, vector< node_id_t >& descendants )
{
	// do a depth first search along edges that have been aligned
	stack< node_id_t > node_stack;
	node_stack.push( node );
	vector< bool > visited( alignment_tree.size(), false );
	descendants.clear();
	while( node_stack.size() > 0 )
	{
		node_id_t cur_node = node_stack.top();
		if(progress_msgs) cout << "Evaluating aligned nodes linked to node " << cur_node << endl;
		node_stack.pop();
		visited[cur_node] = true;
		for( uint childI = 0; childI < alignment_tree[cur_node].children.size(); childI++ )
		{
			node_id_t child_id = alignment_tree[cur_node].children[childI];
			if( alignment_tree[cur_node].children_aligned[childI] && !visited[child_id])
				node_stack.push( child_id );
		}
		if( alignment_tree[ cur_node ].sequence != NULL )
			descendants.push_back( cur_node );
	}
}


/** determine which extant sequences have been aligned at a given node */
void ProgressiveAligner::getPath( node_id_t first_n, node_id_t last_n, vector< node_id_t >& path )
{
	// do a depth first search along edges that have been aligned
	stack< node_id_t > node_stack;
	node_stack.push( last_n );
	vector< bool > visited( alignment_tree.size(), false );
	while( node_stack.top() != first_n )
	{
		node_id_t cur_node = node_stack.top();
		size_t pre_size = node_stack.size();
		visited[cur_node] = true;
		for( uint childI = 0; childI < alignment_tree[cur_node].children.size(); childI++ )
		{
			node_id_t child_id = alignment_tree[cur_node].children[childI];
			if(!visited[child_id])
			{
				node_stack.push( child_id );
				break;
			}
		}
		if( pre_size != node_stack.size() )
			continue;
		for( uint parentI = 0; parentI < alignment_tree[cur_node].parents.size(); parentI++ )
		{
			node_id_t parent_id = alignment_tree[cur_node].parents[parentI];
			if(!visited[parent_id])
			{
				node_stack.push( parent_id );
				break;
			}
		}
		if( pre_size != node_stack.size() )
			continue;
		node_stack.pop();	// didn't make any progress
	}
	path = vector< node_id_t >( node_stack.size() );
	for( size_t pI = 0; pI < path.size(); pI++ )
	{
		path[pI] = node_stack.top();
		node_stack.pop();
	}
}






template<class MatchType>
void ProgressiveAligner::propagateDescendantBreakpoints( node_id_t node1, uint seqI, std::vector<MatchType*>& iv_list )
{
	SSC<MatchType> ilc(seqI);
	sort( iv_list.begin(), iv_list.end(), ilc );
	vector< SuperInterval >& ord = alignment_tree[ node1 ].ordering;
	vector<gnSeqI> bp_list;
	for( size_t sI = 0; sI < ord.size(); sI++ )
		bp_list.push_back( ord[sI].LeftEnd() );

	GenericMatchSeqManipulator<MatchType> ism( seqI );
	applyBreakpoints( bp_list, iv_list, ism );
}

// T should be a pointer type
template<class T, class Manipulator>
void applyAncestralBreakpoints( const vector< SuperInterval >& siv_list, vector<T>& ord, uint seqI, Manipulator& m )
{
	// make bp list
	vector<gnSeqI> bp_list(siv_list.size()*2, 0);
	size_t cur = 0;
	for( size_t i = 0; i < siv_list.size(); i++ )
	{
		if( siv_list[i].reference_iv.Start(seqI) == NO_MATCH )
			continue;
		bp_list[cur++] = siv_list[i].reference_iv.LeftEnd(seqI);
		bp_list[cur++] = siv_list[i].reference_iv.LeftEnd(seqI) + siv_list[i].reference_iv.Length(seqI);
	}
	bp_list.resize(cur);
	// sort the breakpoints and apply...
	sort( bp_list.begin(), bp_list.end() );
	applyBreakpoints( bp_list, ord, m );
}


// assuming breakpoints have been propagated in both directions
// there should now be a 1-to-1 correspondence between superintervals
// in the ancestor and descendants.
void ProgressiveAligner::linkSuperIntervals( node_id_t node1, uint seqI, node_id_t ancestor )
{
	// TODO: speed this up by implementing O(N) instead of O(N^2)
	vector<SuperInterval>& a_ord = alignment_tree[ancestor].ordering;
	vector<SuperInterval>& c_ord = alignment_tree[node1].ordering;
	// initialize all linkages to nothing
	for( size_t aI = 0; aI < a_ord.size(); aI++ )
		if( seqI == 0 )
			a_ord[aI].c1_siv = (std::numeric_limits<size_t>::max)();
		else
			a_ord[aI].c2_siv = (std::numeric_limits<size_t>::max)();
	for( size_t cI = 0; cI < c_ord.size(); cI++ )
		c_ord[cI].parent_siv = (std::numeric_limits<size_t>::max)();

	for( size_t aI = 0; aI < a_ord.size(); aI++ )
	{
		if( a_ord[aI].reference_iv.LeftEnd(seqI) == NO_MATCH )
			continue;
		size_t cI = 0;
		for( ; cI < c_ord.size(); cI++ )
		{
			if( absolut(a_ord[aI].reference_iv.Start(seqI)) != c_ord[cI].LeftEnd() )
				continue;
			if( a_ord[aI].reference_iv.Length(seqI) != c_ord[cI].Length() )
			{
				breakHere();
				cerr << "mapping length mismatch\n";
				cerr << "ancestor: " << ancestor << "\t node1: " << node1 << endl;
				cerr << "a_ord[" << aI << "].reference_iv.Length(" << seqI << "): " << a_ord[aI].reference_iv.Length(seqI) << endl;
				cerr << "a_ord[" << aI << "].reference_iv.LeftEnd(" << seqI << "): " << a_ord[aI].reference_iv.LeftEnd(seqI) << endl;
				cerr << "c_ord[" << cI << "].Length(): " << c_ord[cI].Length() << endl;
				cerr << "c_ord[" << cI << "].LeftEnd(): " << c_ord[cI].LeftEnd() << endl;
				cerr << "";
				cerr << "";
			}
			// link these
			if( seqI == 0 )
				a_ord[aI].c1_siv = cI;
			else
				a_ord[aI].c2_siv = cI;
			c_ord[cI].parent_siv = aI;
			break;
		}
		if( cI == c_ord.size() )
		{
			breakHere();
			cerr << "error no mapping\n";
		}
	}
}


void ProgressiveAligner::translateGappedCoordinates( vector<AbstractMatch*>& ml, uint seqI, node_id_t extant, node_id_t ancestor )
{
	// determine the path that must be traversed
	vector< node_id_t > trans_path;
	getPath( extant, ancestor, trans_path );

	// set seqI to forward orientation 
	for( size_t mI = 0; mI < ml.size(); mI++ )
		if( ml[mI]->Orientation(seqI) == AbstractMatch::reverse )
			ml[mI]->Invert();

	// for each node on the path, construct a complete coordinate translation
	for( size_t nI = 1; nI < trans_path.size(); nI++ )
	{
		// first sort matches on start pos and make them all forward oriented
		// then split them on superinterval boundaries and assign each to a superinterval
		// then convert each match's coordinates to be superinterval-local
		// then apply the coordinate translation with transposeCoordinates
		// then shift each match's coordinates to the global ancestral coordinate space
		SSC<AbstractMatch> ssc(seqI);
		sort(ml.begin(), ml.end(), ssc);

		// split on superinterval boundaries
		vector< SuperInterval >& siv_list = alignment_tree[trans_path[nI]].ordering;
		vector< vector< AbstractMatch* > > siv_matches = vector< vector< AbstractMatch* > >(siv_list.size());
		size_t cur_child = 0;
		if( alignment_tree[trans_path[nI]].children[0] == trans_path[nI-1] )
			cur_child = 0;
		else if( alignment_tree[trans_path[nI]].children[1] == trans_path[nI-1] )
			cur_child = 1;
		else 
		{
			breakHere();
			cerr << "forest fire\n";
		}

		AbstractMatchSeqManipulator amsm( seqI );
		applyAncestralBreakpoints(siv_list, ml, cur_child, amsm );
		
		// sort matches again because new ones were added at the end
		sort(ml.begin(), ml.end(), ssc);

		// assign each match to a siv, and convert coords to siv-local
		for( size_t mI = 0; mI < ml.size(); mI++ )
		{
			if( ml[mI]->LeftEnd(seqI) == 0 )
			{
				breakHere();
				cerr << "fefefe";
			}
			size_t sivI = 0;
			for( ; sivI < siv_list.size(); sivI++ )
			{
				if( siv_list[sivI].reference_iv.LeftEnd(cur_child) == NO_MATCH )
					continue;
				if( ml[mI]->LeftEnd(seqI) >= siv_list[sivI].reference_iv.LeftEnd(cur_child) &&
					ml[mI]->LeftEnd(seqI) < siv_list[sivI].reference_iv.LeftEnd(cur_child) + siv_list[sivI].reference_iv.Length(cur_child) )
					break;
			}
			if( sivI == siv_list.size() )
			{
				cerr << "nI is: "<< nI << endl;
				cerr << "trans_path: ";
				for( size_t ttI = 0; ttI < trans_path.size(); ttI++ )
					cerr << "  " << trans_path[ttI];
				cerr << endl;
				cerr << "problem seq: " << seqI << std::endl;
				cerr << "ml[" << mI << "]->Start(0) == " << ml[mI]->Start(0) << endl;
				cerr << "ml[" << mI << "]->Length(0) == " << ml[mI]->Length(1) << endl;
				cerr << "ml[" << mI << "]->Start(1) == " << ml[mI]->Start(0) << endl;
				cerr << "ml[" << mI << "]->Length(1) == " << ml[mI]->Length(1) << endl;
				cerr << "ml.size(): " << ml.size() << endl;
				for( sivI = 0; sivI < siv_list.size(); sivI++ )
				{
					cerr << "siv_list[" << sivI << "] left end 0: " << siv_list[sivI].reference_iv.LeftEnd(0)  << endl;
					if( siv_list[sivI].reference_iv.LeftEnd(0) != 0 )
						cerr << "siv_list[" << sivI << "] right end 0: " << siv_list[sivI].reference_iv.LeftEnd(0) + siv_list[sivI].reference_iv.Length(0) << endl;
					cerr << "siv_list[" << sivI << "] left end 1: " << siv_list[sivI].reference_iv.LeftEnd(1)  << endl;
					if( siv_list[sivI].reference_iv.LeftEnd(1) != 0 )
						cerr << "siv_list[" << sivI << "] right end 1: " << siv_list[sivI].reference_iv.LeftEnd(1) + siv_list[sivI].reference_iv.Length(1) << endl;
				}
				breakHere();
			}
			if( ml[mI]->LeftEnd(seqI) + ml[mI]->Length(seqI) > 
				siv_list[sivI].reference_iv.LeftEnd(cur_child) + siv_list[sivI].reference_iv.Length(cur_child) )
			{
				cerr << "doesn't fit\n";
				cerr << "ml[" << mI << "]->LeftEnd(" << seqI << "): " << ml[mI]->LeftEnd(seqI) << endl;
				cerr << "ml[" << mI << "]->RightEnd(" << seqI << "): " << ml[mI]->RightEnd(seqI) << endl;
				cerr << "siv_list[" << sivI << "] left end 0: " << siv_list[sivI].reference_iv.LeftEnd(0)  << endl;
				if( siv_list[sivI].reference_iv.LeftEnd(0) != 0 )
					cerr << "siv_list[" << sivI << "] right end 0: " << siv_list[sivI].reference_iv.LeftEnd(0) + siv_list[sivI].reference_iv.Length(0) << endl;
				cerr << "siv_list[" << sivI << "] left end 1: " << siv_list[sivI].reference_iv.LeftEnd(1)  << endl;
					if( siv_list[sivI].reference_iv.LeftEnd(1) != 0 )
						cerr << "siv_list[" << sivI << "] right end 1: " << siv_list[sivI].reference_iv.LeftEnd(1) + siv_list[sivI].reference_iv.Length(1) << endl;
				cerr << "ml.size(): " << ml.size() << endl;
				cerr << "siv_list.size(): " << siv_list.size() << endl;
				cerr << "trans_path:";
				for( size_t tI = 0; tI < trans_path.size(); tI++ )
					cerr << " " << trans_path[tI];
				cerr << endl;
				cerr << "trans_path[" << nI << "]: " << trans_path[nI] << endl;
				breakHere();
			}

			ml[mI]->SetLeftEnd( seqI, ml[mI]->LeftEnd(seqI) - siv_list[sivI].reference_iv.LeftEnd(cur_child) + 1 );
			// if this interval matches the reverse strand then we should effectively invert all matches
			if( siv_list[sivI].reference_iv.Start(cur_child) < 0 )
			{
				int64 new_lend = siv_list[sivI].reference_iv.Length(cur_child) - ml[mI]->LeftEnd(seqI);
				new_lend -= ml[mI]->Length( seqI ) - 2;
				new_lend *= ml[mI]->Orientation(seqI) == AbstractMatch::forward ? 1 : -1;
				ml[mI]->Invert();
				ml[mI]->SetStart( seqI, new_lend ); 
			}
			siv_matches[sivI].push_back( ml[mI] );
		}

		// apply the coordinate translation
		ml.clear();
		for( size_t sivI = 0; sivI < siv_matches.size(); sivI++ )
		{
			if( siv_matches[sivI].size() == 0 )
				continue;
			
			// get a CompactGappedAlignment<> for this interval
			CompactGappedAlignment<>* siv_cga = dynamic_cast<CompactGappedAlignment<>*>(siv_list[sivI].reference_iv.GetMatches()[0]);
			if( siv_list[sivI].reference_iv.GetMatches().size() > 1 )
				siv_cga = NULL;
			bool alloc_new_siv = false;
			CompactGappedAlignment<> tmp_cga;
			if( siv_cga == NULL )
			{
				alloc_new_siv = true;
				siv_cga = tmp_cga.Copy();
				CompactGappedAlignment<> dorkas(siv_list[sivI].reference_iv);
				*siv_cga = dorkas;
			}

			// now translate each match...
			for( size_t mI = 0; mI < siv_matches[sivI].size(); mI++ )
			{
				CompactGappedAlignment<>* match_cga = dynamic_cast<CompactGappedAlignment<>*>(siv_matches[sivI][mI]);
				bool alloc_new = false;
				if( match_cga == NULL )
				{
					match_cga = tmp_cga.Copy();
					*match_cga = CompactGappedAlignment<>(*(siv_matches[sivI][mI]));
					alloc_new = true;
				}
				siv_cga->translate( *match_cga, seqI, cur_child );

				if( alloc_new )
				{
					siv_matches[sivI][mI]->Free();
					siv_matches[sivI][mI] = match_cga;
				}
			}

			// shift coordinates back to global space
			for( size_t mI = 0; mI < siv_matches[sivI].size(); mI++ )
			{
				int64 cur_start = siv_matches[sivI][mI]->Start(seqI);
				if( cur_start > 0 )
					siv_matches[sivI][mI]->SetStart( seqI, cur_start + siv_list[sivI].LeftEnd() - 1 );
				else
					siv_matches[sivI][mI]->SetStart( seqI, cur_start - siv_list[sivI].LeftEnd() + 1);
				if( (siv_matches[sivI][mI]->LeftEnd(seqI) + siv_matches[sivI][mI]->Length(seqI) > siv_list.back().LeftEnd() + siv_list.back().Length() )
					 )
				{
					// is there something wrong with the translation table?
					cerr << "siv left is: " << siv_list[sivI].LeftEnd() << endl;
					cerr << "siv right is: " << siv_list[sivI].LeftEnd() + siv_list[sivI].Length() << endl;
					cerr << "match right is: " << siv_matches[sivI][mI]->LeftEnd(seqI) + siv_matches[sivI][mI]->Length(seqI) << endl;
					cerr << "superseq right is: " << siv_list.back().LeftEnd() + siv_list.back().Length() << endl;
					cerr << "";
					breakHere();
				}
				if( debug_aligner && siv_matches[sivI][mI]->Start(seqI) == 0 )
				{
					breakHere();
				}
			}
			if(alloc_new_siv)
				siv_cga->Free();
			ml.insert( ml.end(), siv_matches[sivI].begin(), siv_matches[sivI].end() );
		}
	}
	// restore forward orientation seqI
	for( size_t mI = 0; mI < ml.size(); mI++ )
		if( ml[mI]->Orientation(seqI) == AbstractMatch::reverse )
			ml[mI]->Invert();
}

class SuperIntervalPtrComp
{
public:
	bool operator()( const SuperInterval* a, const SuperInterval* b )
	{
		return (*a) < (*b);
	}
};

void ProgressiveAligner::recursiveApplyAncestralBreakpoints( node_id_t ancestor )
{
	stack<node_id_t> node_stack;
	node_stack.push(ancestor);
	while( node_stack.size() > 0 )
	{
		// pop the current node, apply ancestral breakpoints, recurse on children
		node_id_t cur_node = node_stack.top();
		node_stack.pop();
		SuperIntervalManipulator sim;
		if( progress_msgs ) cout << "cur node: " << cur_node << endl;
		for( size_t childI = 0; childI < alignment_tree[cur_node].children.size(); childI++ )
		{
			AlignmentTreeNode& atn = alignment_tree[alignment_tree[cur_node].children[childI]];
			if( progress_msgs ) cout << "childI " << childI << " aab\n";
			applyAncestralBreakpoints( alignment_tree[cur_node].ordering, atn.ordering, childI, sim );
			if( progress_msgs ) cout << "sort childI " << childI << "\n";
			vector<SuperInterval*> siv_ptr_list(atn.ordering.size());
			for( size_t sivI = 0; sivI < atn.ordering.size(); ++sivI )
				siv_ptr_list[sivI] = &(atn.ordering[sivI]);
			SuperIntervalPtrComp sipc;
			sort( siv_ptr_list.begin(), siv_ptr_list.end(), sipc );
			vector< SuperInterval > siv_list;
			for( size_t sivI = 0; sivI < siv_ptr_list.size(); ++sivI )
				siv_list.push_back(*siv_ptr_list[sivI]);
			swap(siv_list, atn.ordering);
			node_stack.push( alignment_tree[cur_node].children[childI] );
		}
		if( debug_aligner && alignment_tree[cur_node].children.size() > 0 )
			validateSuperIntervals(alignment_tree[cur_node].children[0], alignment_tree[cur_node].children[1], cur_node);
		if( progress_msgs ) cout << "linking node " << cur_node << "'s" << alignment_tree[cur_node].ordering.size() << " superintervals\n"; 
		for( size_t childI = 0; childI < alignment_tree[cur_node].children.size(); childI++ )
			linkSuperIntervals( alignment_tree[cur_node].children[childI], childI, cur_node );
	}
}


boolean getInterveningCoordinates( const AbstractMatch* iv, uint oseqI, Match* r_begin, Match* r_end, uint seqI, int64& gap_lend, int64& gap_rend ){
	// skip this sequence if it's undefined
	if( (r_end != NULL && r_end->Start( seqI ) == NO_MATCH) ||
		(r_begin != NULL && r_begin->Start( seqI ) == NO_MATCH) ){
		gap_lend = 0;
		gap_rend = 0;
		return true;
	}
			
	// determine the size of the gap
	gap_rend = r_end != NULL ? r_end->Start( seqI ) : iv->RightEnd( oseqI ) + 1;
	gap_lend = r_begin != NULL ? r_begin->End( seqI ) + 1 : iv->LeftEnd( oseqI );
	if( gap_rend < 0 || gap_lend < 0 ){
		gap_rend = r_begin != NULL ? -r_begin->Start( seqI ) : iv->RightEnd( oseqI ) + 1;
		gap_lend = r_end != NULL ? -r_end->Start( seqI ) + r_end->Length() : 1;
	}
	if( gap_rend <= 0 || gap_lend <= 0 ){
		// if either is still < 0 then there's a problem...
		genome::ErrorMsg( "Error constructing intervening coordinates" );
	}
	return true;
}


void ProgressiveAligner::pairwiseAnchorSearch( MatchList& r_list, Match* r_begin, Match* r_end, const AbstractMatch* iv, uint oseqI, uint oseqJ )
{
	uint seqI = 0;
	MatchList gap_list;
	vector< int64 > starts;
// 
//	Get the sequence in the intervening gaps between these two matches
//
	for( seqI = 0; seqI < 2; seqI++ )
	{
		int64 gap_end = 0;
		int64 gap_start = 0;
		getInterveningCoordinates( iv, (seqI == 0 ? oseqI : oseqJ), r_begin, r_end, seqI, gap_start, gap_end);
		int64 diff = gap_end - gap_start;
		diff = diff > 0 ? diff - 1 : 0;

		starts.push_back( gap_start );
		gnSequence* new_seq = NULL;
		if(diff > 0 && gap_start + diff - 1 <= r_list.seq_table[ seqI ]->length())
			new_seq = new gnSequence( r_list.seq_table[ seqI ]->ToString( diff, gap_start ) );
		else
			new_seq = new gnSequence();
		gap_list.seq_table.push_back( new_seq );
		gap_list.sml_table.push_back( new DNAMemorySML() );
	}

	gnSeqI avg_len = (gap_list.seq_table[0]->length() + gap_list.seq_table[1]->length())/2;
	uint search_seed_size = getDefaultSeedWeight( avg_len );
	gap_mh.get().Clear();
	
	uint seed_count = use_seed_families ? 3 : 1;
	for( size_t seedI = 0; seedI < seed_count; seedI++ )
	{
		//
		//	Create sorted mer lists for the intervening gap region
		//
		uint64 default_seed = getSeed( search_seed_size, seedI );
		if( search_seed_size < MIN_DNA_SEED_WEIGHT )
		{
			for( uint seqI = 0; seqI < gap_list.seq_table.size(); seqI++ )
				delete gap_list.seq_table[ seqI ];
			for( uint seqI = 0; seqI < gap_list.sml_table.size(); seqI++ )
				delete gap_list.sml_table[ seqI ];
			return;
		}
		for( uint seqI = 0; seqI < gap_list.seq_table.size(); seqI++ ){
			gap_list.sml_table[ seqI ]->Clear();
			gap_list.sml_table[ seqI ]->Create( *(gap_list.seq_table[ seqI ]), default_seed );
		}

		//
		//	Find all matches in the gap region
		//
		gap_mh.get().ClearSequences();
		if(seed_count>1)
		{
			MatchList cur_list = gap_list;
			gap_mh.get().FindMatches( cur_list );
			for( size_t mI = 0; mI < cur_list.size(); mI++ )
				cur_list[mI]->Free();
		}else
			gap_mh.get().FindMatches( gap_list );
	}
	if(seed_count>1)
		gap_mh.get().GetMatchList(gap_list);

	EliminateOverlaps_v2( gap_list );

	// for anchor accuracy, throw out any anchors that are shorter than the minimum
	// anchor length after EliminateOverlaps()
	gap_list.LengthFilter( MIN_ANCHOR_LENGTH + 3 );

	for( size_t gI = 0; gI < gap_list.size(); gI++ )
	{
		for( seqI = 0; seqI < 2; seqI++ )
		{
			int64 gap_rend = 0;
			int64 gap_lend = 0;
			getInterveningCoordinates( iv, (seqI == 0 ? oseqI : oseqJ), r_begin, r_end, seqI, gap_lend, gap_rend);
			gap_list[gI]->SetLeftEnd(seqI, gap_list[gI]->LeftEnd(seqI) + gap_lend - 1);
		}
	}
	r_list.insert(r_list.end(), gap_list.begin(), gap_list.end());

	// delete sequences and smls
	for( uint seqI = 0; seqI < gap_list.seq_table.size(); seqI++ )
		delete gap_list.seq_table[ seqI ];
	for( uint seqI = 0; seqI < gap_list.sml_table.size(); seqI++ )
		delete gap_list.sml_table[ seqI ];
}

template<class GappedAlignmentType>
void ProgressiveAligner::recurseOnPairs( const vector<node_id_t>& node1_seqs, const vector<node_id_t>& node2_seqs, const GappedAlignmentType& iv, Matrix<MatchList>& matches, Matrix< std::vector< search_cache_t > >& search_cache_db, Matrix< std::vector< search_cache_t > >& new_cache_db, boost::multi_array< vector< vector< int64 > >, 2 >& iv_regions )
{
	matches = Matrix<MatchList>(node1_seqs.size(),node2_seqs.size());

	std::vector< bitset_t > aln_matrix;
	iv.GetAlignment(aln_matrix);
	Match tmp(2);
	const size_t sizer = node1_seqs.size() * node2_seqs.size();
	std::vector< std::pair<size_t,size_t> > node_pairs(sizer);
	int nni = 0;
	for( size_t n1 = 0; n1 < node1_seqs.size(); n1++ )
		for( size_t n2 = 0; n2 < node2_seqs.size(); n2++ )
			node_pairs[nni++] = make_pair(n1,n2);

#pragma omp parallel for
	for(int ni = 0; ni < node_pairs.size(); ni++)
	{
		size_t n1 = node_pairs[ni].first;
		size_t n2 = node_pairs[ni].second;
		vector<node_id_t>::const_iterator n1_iter = node1_seqs.begin() + n1;
		vector<node_id_t>::const_iterator n2_iter = node2_seqs.begin() + n2;
		
		uint seqI = node_sequence_map[*n1_iter];
		uint seqJ = node_sequence_map[*n2_iter];
		MatchList& mlist = matches(n1, n2);
		std::vector< search_cache_t >& cache = search_cache_db(n1, n2);
		std::vector< search_cache_t >& new_cache = new_cache_db(n1, n2);
		mlist.seq_table.push_back( alignment_tree[*n1_iter].sequence );
		mlist.seq_table.push_back( alignment_tree[*n2_iter].sequence );

		if( iv.LeftEnd(seqI) == NO_MATCH )
		{
			if( iv.LeftEnd(seqJ) != NO_MATCH )
			{
				iv_regions[n1][n2][1].push_back(iv.LeftEnd(seqJ));
				iv_regions[n1][n2][1].push_back(iv.RightEnd(seqJ));
			}
			continue;	// no sense searching one isn't defined!
		}
		if(iv.LeftEnd(seqJ) == NO_MATCH )
		{
			if( iv.LeftEnd(seqI) != NO_MATCH )
			{
				iv_regions[n1][n2][0].push_back(iv.LeftEnd(seqI));
				iv_regions[n1][n2][0].push_back(iv.RightEnd(seqI));
			}
			continue;	// no sense searching one isn't defined!
		}

		gnSeqI charI = 0;
		gnSeqI charJ = 0;
		const size_t iv_aln_length = iv.AlignmentLength();

// first determine the outer aligned boundaries of the LCB and record them for
// later use
		pair< int64, int64 > pair_1l(0,0);
		pair< int64, int64 > pair_1r(0,0);
		pair< int64, int64 > pair_2l(0,0);
		pair< int64, int64 > pair_2r(0,0);
		for( uint colI = 0; colI <= iv_aln_length; colI++ )
		{
			if( colI == iv_aln_length || (aln_matrix[seqI].test(colI) && aln_matrix[seqJ].test(colI)) )
			{
				if( colI == 0 )
					break;	// nothing to see here, move along...
				if( iv.Orientation(seqI) == AbstractMatch::forward )
					pair_1l = make_pair( iv.LeftEnd(seqI), iv.LeftEnd(seqI)+charI );
				else
					pair_1r = make_pair( iv.RightEnd(seqI)-charI+1, iv.RightEnd(seqI)+1 );
				if( iv.Orientation(seqJ) == AbstractMatch::forward )
					pair_2l = make_pair( iv.LeftEnd(seqJ), iv.LeftEnd(seqJ)+charJ );
				else
					pair_2r = make_pair( iv.RightEnd(seqJ)-charJ+1, iv.RightEnd(seqJ)+1 );
				break;
			}
			if( colI < iv_aln_length && aln_matrix[seqI].test(colI) )
				++charI;
			if( colI < iv_aln_length && aln_matrix[seqJ].test(colI) )
				++charJ;
		}

		charI = 0;
		charJ = 0;
		for( uint colI = iv_aln_length; colI > 0 ; colI-- )
		{
			if( (aln_matrix[seqI].test(colI-1) && aln_matrix[seqJ].test(colI-1)) )
			{
				if( colI == iv_aln_length )
					break;	// nothing to see here, move along...
				if( iv.Orientation(seqI) == AbstractMatch::forward )
					pair_1r = make_pair( iv.RightEnd(seqI)-charI+1, iv.RightEnd(seqI)+1 );
				else
					pair_1l = make_pair( iv.LeftEnd(seqI), iv.LeftEnd(seqI)+charI );
				if( iv.Orientation(seqJ) == AbstractMatch::forward )
					pair_2r = make_pair( iv.RightEnd(seqJ)-charJ+1, iv.RightEnd(seqJ)+1 );
				else
					pair_2l = make_pair( iv.LeftEnd(seqJ), iv.LeftEnd(seqJ)+charJ );
				break;
			}
			if( aln_matrix[seqI].test(colI-1) )
				++charI;
			if( aln_matrix[seqJ].test(colI-1) )
				++charJ;
		}
		if( pair_1l.first < pair_1l.second )
		{
			iv_regions[n1][n2][0].push_back(pair_1l.first);
			iv_regions[n1][n2][0].push_back(pair_1l.second);
		}
		if( pair_1r.first < pair_1r.second )
		{
			if( pair_1l.first < pair_1l.second && pair_1r.first == pair_1l.second )
			{
				// just merge them into a single interval
				iv_regions[n1][n2][0].back() = pair_1r.second;
			}else{
				iv_regions[n1][n2][0].push_back(pair_1r.first);
				iv_regions[n1][n2][0].push_back(pair_1r.second);
				if( pair_1r.first <= pair_1l.second && pair_1r.second >= pair_1l.first )
				{
					cout << "Ohno.  Overlap in outside LCB search intervals\n";
					cout << "Left: " << pair_1l.first << '\t' << pair_1l.second << " right:  " << pair_1r.first << '\t' << pair_1r.second << endl;
					cout << "0 iv.Start(" << seqI << "): " << iv.Start(seqI) << '\t' << "iv.RightEnd(" << seqI << "): " << iv.RightEnd(seqI) << endl;
					if( pair_1l.first == 0 )
						genome::breakHere();
				}
			}
		}

		if( pair_2l.first < pair_2l.second )
		{
			iv_regions[n1][n2][1].push_back(pair_2l.first);
			iv_regions[n1][n2][1].push_back(pair_2l.second);
		}
		if( pair_2r.first < pair_2r.second )
		{
			if( pair_2l.first < pair_2l.second && pair_2r.first == pair_2l.second )
			{
				// just merge them into a single interval
				iv_regions[n1][n2][1].back() = pair_2r.second;
			}else{
				iv_regions[n1][n2][1].push_back(pair_2r.first);
				iv_regions[n1][n2][1].push_back(pair_2r.second);
				if( pair_2r.first <= pair_2l.second && pair_2r.second >= pair_2l.first )
				{
					cout << "Ohno.  Overlap in outside LCB search intervals\n";
					cout << "Left: " << pair_2l.first << '\t' << pair_2l.second << " right:  " << pair_2r.first << '\t' << pair_2r.second << endl;
					cout << "1 iv.Start(" << seqJ << "): " << iv.Start(seqJ) << '\t' << "iv.RightEnd(" << seqJ << "): " << iv.RightEnd(seqJ) << endl;
					cout << "charI " << charI << "\tcharJ" << charJ << endl;
					if( pair_2l.first == 0 )
						genome::breakHere();
				}
			}
		}

		charI = 0;
		charJ = 0;
		gnSeqI prev_charI = 0;
		gnSeqI prev_charJ = 0;
		bool in_gap = false;

		for( uint colI = 0; colI <= iv_aln_length; colI++ )
		{
			if( colI == iv_aln_length || 
				(aln_matrix[seqI].test(colI) && aln_matrix[seqJ].test(colI)) )
			{
				if( in_gap && 
					charI - prev_charI > min_recursive_gap_length &&
					charJ - prev_charJ > min_recursive_gap_length )
				{

					Match* l_match = NULL;
					l_match = tmp.Copy();
					if(iv.Orientation(seqI) == AbstractMatch::forward)
						l_match->SetLeftEnd(0, iv.LeftEnd(seqI)+prev_charI);
					else
					{
						l_match->SetLeftEnd(0, iv.RightEnd(seqI)-prev_charI);
						l_match->SetOrientation(0, AbstractMatch::reverse );
					}
					if(iv.Orientation(seqJ) == AbstractMatch::forward)
						l_match->SetLeftEnd(1, iv.LeftEnd(seqJ)+prev_charJ);
					else
					{
						l_match->SetLeftEnd(1, iv.RightEnd(seqJ)-prev_charJ);
						l_match->SetOrientation(1, AbstractMatch::reverse );
					}
					l_match->SetLength(0);
					Match* r_match = NULL;
					if( charJ != iv.RightEnd(seqJ) && charI != iv.RightEnd(seqI) )
					{
						r_match = tmp.Copy();
						if(iv.Orientation(seqI) == AbstractMatch::forward)
							r_match->SetLeftEnd(0, iv.LeftEnd(seqI)+charI);
						else
						{
							r_match->SetLeftEnd(0, iv.RightEnd(seqI)-charI);
							r_match->SetOrientation(0, AbstractMatch::reverse );
						}
						if(iv.Orientation(seqJ) == AbstractMatch::forward)
							r_match->SetLeftEnd(1, iv.LeftEnd(seqJ)+charJ);
						else
						{
							r_match->SetLeftEnd(1, iv.RightEnd(seqJ)-charJ);
							r_match->SetOrientation(1, AbstractMatch::reverse );
						}
						r_match->SetLength(0);
					}

					if( iv.Orientation(seqI) == AbstractMatch::reverse )
					{
						swap(l_match,r_match);
						if( l_match != NULL ) l_match->Invert();
						if( r_match != NULL ) r_match->Invert();
					}
					// check whether the current cache already has the searched region
					search_cache_t cacheval = make_pair( l_match, r_match );
					std::vector< search_cache_t >::iterator cache_entry = std::upper_bound( cache.begin(), cache.end(), cacheval, mems::cache_comparator );
					if( cache_entry == cache.end() || 
						(mems::cache_comparator( cacheval, *cache_entry ) || mems::cache_comparator( *cache_entry, cacheval )) )
					{
						// search this region
							pairwiseAnchorSearch(mlist, l_match, r_match, &iv, seqI, seqJ);
					}
					if(using_cache_db)
						new_cache.push_back( cacheval );
				}
				prev_charI = charI;
				prev_charJ = charJ;
				in_gap = false;
			}
			else
				in_gap = true;
			if( colI < iv.AlignmentLength() )
			{
				if( aln_matrix[seqI].test(colI) )
					++charI;
				if( aln_matrix[seqJ].test(colI) )
					++charJ;
			}
		}
	}
}

void ProgressiveAligner::getAncestralMatches( const vector< node_id_t > node1_seqs, const vector< node_id_t > node2_seqs, node_id_t node1, node_id_t node2, node_id_t ancestor, std::vector< AbstractMatch* >& ancestral_matches )
{
	// to save memory, always make node1_seqs the bigger vector
//	if( node1_seqs.size() < node2_seqs.size() )
//		swap( node1_seqs, node2_seqs );

	// for each pair of genomes, extract pairwise matches and translate up
	// eliminate overlaps
	for( uint seqI = 0; seqI < node1_seqs.size(); seqI++ )
	{
		uint ii = this->node_sequence_map[node1_seqs[seqI]];
		vector< AbstractMatch* > seqI_matches;

		for( uint seqJ = 0; seqJ < node2_seqs.size(); seqJ++ )
		{
			uint jj = this->node_sequence_map[node2_seqs[seqJ]];
			vector< AbstractMatch* > cur_matches;
			for( size_t mI = 0; mI < original_ml.size(); mI++ )
			{
				if( original_ml[mI]->LeftEnd(ii) == NO_MATCH )
					continue;
				if( original_ml[mI]->LeftEnd(jj) == NO_MATCH )
					continue;
				Match mm( 2 );
				Match* new_m = mm.Copy();
				new_m->SetStart( 0, original_ml[mI]->Start(ii));
				new_m->SetStart( 1, original_ml[mI]->Start(jj));
				new_m->SetLength(original_ml[mI]->Length());
				if( new_m->Start(0) < 0 )
					new_m->Invert();	// assign reference orientation to seq 0
				cur_matches.push_back( new_m );
			}
			// now translate cur_matches
			translateGappedCoordinates( cur_matches, 1, node2_seqs[seqJ], node2 );
			seqI_matches.insert( seqI_matches.end(), cur_matches.begin(), cur_matches.end() );
		}
		EliminateOverlaps_v2( seqI_matches );
		translateGappedCoordinates( seqI_matches, 0, node1_seqs[seqI], node1 );
		ancestral_matches.insert( ancestral_matches.end(), seqI_matches.begin(), seqI_matches.end() );
	}
	EliminateOverlaps_v2( ancestral_matches );
}


void ProgressiveAligner::getPairwiseMatches( const vector< node_id_t >& node1_seqs, const vector< node_id_t >& node2_seqs, Matrix<MatchList>& pairwise_matches )
{
	pairwise_matches = Matrix< MatchList >( node1_seqs.size(), node2_seqs.size() );

	// copy sequence tables
	for( uint seqI = 0; seqI < node1_seqs.size(); seqI++ )
	{
		for( uint seqJ = 0; seqJ < node2_seqs.size(); seqJ++ )
		{
			uint ii = this->node_sequence_map[node1_seqs[seqI]];
			uint jj = this->node_sequence_map[node2_seqs[seqJ]];
			pairwise_matches(seqI, seqJ).seq_table.push_back(original_ml.seq_table[ii]);
			pairwise_matches(seqI, seqJ).seq_table.push_back(original_ml.seq_table[jj]);
			pairwise_matches(seqI, seqJ).seq_filename.push_back(original_ml.seq_filename[ii]);
			pairwise_matches(seqI, seqJ).seq_filename.push_back(original_ml.seq_filename[jj]);
		}
	}

	// now copy pairwise matches
	for( size_t mI = 0; mI < original_ml.size(); mI++ )
	{
		for( uint seqI = 0; seqI < node1_seqs.size(); seqI++ )
		{
			uint ii = this->node_sequence_map[node1_seqs[seqI]];
			if( original_ml[mI]->LeftEnd(ii) == NO_MATCH )
				continue;
			for( uint seqJ = 0; seqJ < node2_seqs.size(); seqJ++ )
			{
				uint jj = this->node_sequence_map[node2_seqs[seqJ]];
				if( original_ml[mI]->LeftEnd(jj) == NO_MATCH )
					continue;
				Match mm( 2 );
				Match* new_m = mm.Copy();
				new_m->SetStart( 0, original_ml[mI]->Start(ii));
				new_m->SetStart( 1, original_ml[mI]->Start(jj));
				new_m->SetLength(original_ml[mI]->Length());
				if( new_m->Start(0) < 0 )
					new_m->Invert();	// assign reference orientation to seq 0
				pairwise_matches(seqI,seqJ).push_back( new_m );
			}
		}
	}
}


int IsDenseEnough( GappedAlignment* gal_iter )
{
	double total_len = 0;
	gnSeqI seqs = 0;
	for( uint seqI = 0; seqI < gal_iter->SeqCount(); seqI++ )
	{
		if( gal_iter->LeftEnd(seqI) == NO_MATCH )
			continue;
		total_len += gal_iter->Length(seqI);
	}
	double density = total_len / (gal_iter->AlignmentLength() * (double)gal_iter->Multiplicity());
	// density of 1 is ideal
	// the shorter the alignment, the closer we should be to 1 to allow splitting
	// use a linear threshold with (min_window_size,1) and (max_window_size,min_gappiness)
	// as endpoints of the threshold line
	
	// determine the density threshold for the given alignment length
	double threshold = ((max_density - min_density)/(min_window_size - max_window_size)) * ( (double)gal_iter->AlignmentLength() - max_window_size ) + min_density;
	if( density > max_density )	// don't bother aligning this, it's so dense we'll wait until iterative refinement.
		return 2;
	if( density > threshold )
		return 1;
	return 0;
}

void splitGappedAlignment( const GappedAlignment& ga, GappedAlignment& ga1, GappedAlignment& ga2, std::vector<size_t>& seqs1, std::vector<size_t>& seqs2 )
{
	const vector< string >& aln = GetAlignment( ga, std::vector<gnSequence*>(ga.SeqCount()) );
	ga1 = ga;
	ga2 = ga;
	for( size_t seqI = 0; seqI < seqs1.size(); seqI++ )
		ga2.SetLeftEnd(seqs1[seqI], NO_MATCH);
	for( size_t seqI = 0; seqI < seqs2.size(); seqI++ )
		ga1.SetLeftEnd(seqs2[seqI], NO_MATCH);
}

void removeLargeGapsPP( GappedAlignment& gal, list< GappedAlignment* >& gal_list, vector<bool>& gap_iv, const vector< size_t >& group1, const vector< size_t >& group2 )
{
	// scan through and remove any section where members of group1 aren't aligned to members of group2
	// for more than some number of nucleotides
	gap_iv.clear();
	gal_list.clear();
	const vector< string >& aln_matrix = GetAlignment(gal, vector<gnSequence*>(gal.SeqCount(),NULL));
	size_t gap_cols = 0;
	size_t last_aln_col = (std::numeric_limits<size_t>::max)();
	size_t col_base = 0;
	GappedAlignment* galp = gal.Copy();
	for( size_t colI = 0; colI < gal.AlignmentLength(); colI++ )
	{
		 size_t g1 = 0;
		 size_t g2 = 0;
		 for( ; g1 < group1.size(); ++g1 )
		 {
			 if( aln_matrix[group1[g1]][colI] != '-' )
				 break;
		 }
		 for( ; g2 < group2.size(); ++g2 )
		 {
			 if( aln_matrix[group2[g2]][colI] != '-' )
				 break;
		 }
		 if( g1 < group1.size() && g2 < group2.size() )
		 {
			 // it's an aligned col
			 if( gap_cols > max_gap_length )
			 {
				// crop out the middle gapped section
				gnSeqI split_point = 0;
				if( last_aln_col != (std::numeric_limits<size_t>::max)() )
				{
					split_point = last_aln_col + lcb_hangover - col_base;
					gal_list.push_back( galp );
					gap_iv.push_back(false);
					galp = (GappedAlignment*)galp->Split(split_point);	// set galp to the right side after splitting
					col_base += split_point;
				}
				split_point = colI - lcb_hangover - col_base;
				gal_list.push_back( galp );
				gap_iv.push_back(true);
				galp = (GappedAlignment*)galp->Split(split_point);	// set galp to the right side after splitting
				col_base += split_point;
			 }
			 last_aln_col = colI;
			 gap_cols = 0;
		 }else
			 ++gap_cols;
	}

	if( gap_cols > max_gap_length )
	{
		gnSeqI split_point = 0;
		if( last_aln_col != (std::numeric_limits<size_t>::max)() )
		{
			split_point = last_aln_col + lcb_hangover - col_base;
			gal_list.push_back( galp );
			gap_iv.push_back(false);
			galp = (GappedAlignment*)galp->Split(split_point);	// set galp to the right side after splitting
		}
		gap_iv.push_back(true);
	}else
		gap_iv.push_back(false);
	gal_list.push_back( galp );
}

void ProgressiveAligner::refineAlignment( GappedAlignment& gal, node_id_t ancestor, bool profile_aln, AlnProgressTracker& apt )
{
	// divide the gapped alignment up into windows of a given size and have
	// muscle refine the alignments
	// when anchors are dense use smaller windows to improve speed efficiency
	list< GappedAlignment* > gal_list;
	vector<bool> gap_iv;
	std::vector<node_id_t> nodes1;
	std::vector<node_id_t> nodes2;
	getAlignedChildren( alignment_tree[ancestor].children[0], nodes1 );
	getAlignedChildren( alignment_tree[ancestor].children[1], nodes2 );
	std::vector<size_t> seqs1( nodes1.size() );
	std::vector<size_t> seqs2( nodes2.size() );
	for( size_t nI = 0; nI < nodes1.size(); nI++ )
		seqs1[nI] = node_sequence_map[nodes1[nI]];
	for( size_t nI = 0; nI < nodes2.size(); nI++ )
		seqs2[nI] = node_sequence_map[nodes2[nI]];
//	if( profile_aln )
//	{
		removeLargeGapsPP( gal, gal_list, gap_iv, seqs1, seqs2 );
//	}else{
//		gal_list.push_back( gal.Copy() );
//		gap_iv.push_back(false);
//	}
	list< GappedAlignment* >::iterator gal_iter = gal_list.begin();
	vector<bool>::iterator gap_iter = gap_iv.begin();
	while(gal_iter != gal_list.end())
	{
		int density = IsDenseEnough( *gal_iter );
		if( (density == 0 && (*gal_iter)->AlignmentLength() > max_window_size / 3) ||
			(density == 1 && (*gal_iter)->AlignmentLength() > max_window_size ) ||
			(density == 2 && (*gal_iter)->AlignmentLength() > max_window_size * 3 )

//			  || ( (*gal_iter)->AlignmentLength() > min_window_size && density == 1 && profile_aln == true ) 
			  )
		{
			// split in half
			gnSeqI split_point = (*gal_iter)->AlignmentLength() / 2;
			list< GappedAlignment* >::iterator ins_iter = gal_iter;
			++ins_iter;
//			ins_iter = gal_list.insert(ins_iter, new GappedAlignment(**gal_iter) );
			ins_iter = gal_list.insert(ins_iter, (*gal_iter)->Copy());
			vector<bool>::iterator gap_ins_iter = gap_iter;
			size_t gap_off = gap_iter - gap_iv.begin();
			++gap_ins_iter;
			gap_iv.insert( gap_ins_iter, *gap_iter );
			gap_iter = gap_iv.begin() + gap_off;
			(*gal_iter)->CropEnd( split_point );
			(*ins_iter)->CropStart( (*ins_iter)->AlignmentLength() - split_point );
			continue;
		}

		++gal_iter;
		++gap_iter;
	}
	MuscleInterface& mi = MuscleInterface::getMuscleInterface();
	// now that the alignment is all split up use muscle to refine it
	gnSeqI new_len = 0;

	gap_iter = gap_iv.begin();

	const size_t gal_count = gal_list.size();
// this section can not be paralellized b/c it makes calls to muscle
#pragma omp critical
{
	for( int galI = 0; galI < gal_count; galI++ )
	{
		list<GappedAlignment*>::iterator my_g_iter = gal_list.begin();
		vector<bool>::iterator my_b_iter = gap_iv.begin();
		for(uint a = 0; a < galI; a++)
		{
			++my_g_iter;
			++my_b_iter;
		}
		apt.cur_leftend += (*my_g_iter)->AlignmentLength();
		if( profile_aln && !(*my_b_iter) )
		{
			GappedAlignment ga1;
			GappedAlignment ga2;
			splitGappedAlignment( **my_g_iter, ga1, ga2, seqs1, seqs2 );
			if( ga1.Multiplicity() > 0 && ga2.Multiplicity() > 0 )
			{
				mi.ProfileAlignFast( ga1, ga2, **my_g_iter, true );
			}
		}else if(!(*my_b_iter))
		{
			int density = IsDenseEnough( *my_g_iter );
			if( density == 0 )
				mi.RefineFast( **my_g_iter );
			else if( density == 1 )
				mi.RefineFast( **my_g_iter, 500 );
			else
				mi.RefineFast( **my_g_iter, 200 );
		}

		new_len += (*my_g_iter)->AlignmentLength();
		// print a progress message
		double cur_progress = ((double)apt.cur_leftend / (double)apt.total_len)*100.0;
		printProgress((uint)apt.prev_progress, (uint)cur_progress, cout);
		apt.prev_progress = cur_progress;
	}
	gal_iter = gal_list.end();
}

	// put humpty dumpty back together
	vector< string > aln_matrix( gal.SeqCount(), string( new_len, '-' ) );
	vector< string::size_type > pos( gal.SeqCount(), 0 );
	for( gal_iter = gal_list.begin(); gal_iter != gal_list.end(); ++gal_iter )
	{
		const vector< string >& tmp_mat = GetAlignment(**gal_iter, vector<gnSequence*>( gal.SeqCount() ) );
		for( uint seqI = 0; seqI < tmp_mat.size(); seqI++ )
		{
			if( gal.LeftEnd(seqI) == 0 )
				continue;
			aln_matrix[seqI].replace(pos[seqI], tmp_mat[seqI].size(), tmp_mat[seqI]);
			pos[seqI] += tmp_mat[seqI].size();
		}
		(*gal_iter)->Free();
//		delete (*gal_iter);
	}
	gal.SetAlignment(aln_matrix);
}

void ProgressiveAligner::doGappedAlignment( node_id_t ancestor, bool profile_aln )
{
	AlnProgressTracker apt;
	gnSeqI total_len = 0;
	for( size_t aI = 0; aI < alignment_tree[ancestor].ordering.size(); aI++ )
		total_len += alignment_tree[ancestor].ordering[aI].Length();
	apt.total_len = total_len;
	apt.prev_progress = 0;

	printProgress(-1, 0, cout);
	apt.cur_leftend = 1;

	for( size_t aI = 0; aI < alignment_tree[ancestor].ordering.size(); aI++ )
	{
		if( alignment_tree[ancestor].ordering[aI].reference_iv.Multiplicity() == 1 )
		{
			apt.cur_leftend += alignment_tree[ancestor].ordering[aI].reference_iv.AlignmentLength();
			continue;	// don't bother re-refining intervals that didn't get aligned here
		}

//		printMemUsage();
//		cout << "extract aln\n";
		GappedAlignment gal;
		extractAlignment(ancestor, aI, gal);
//		printMemUsage();
//		cout << "refine aln\n";
		if( gal.Multiplicity() > 1 )	// no point in refining intervals that are unaligned anyways
			refineAlignment( gal, ancestor, profile_aln, apt );
		else
			apt.cur_leftend += gal.AlignmentLength();
//		printMemUsage();
//		cout << "construct siv\n";
		ConstructSuperIntervalFromMSA(ancestor, aI, gal);
//		printMemUsage();

		// print a progress message
		double cur_progress = ((double)apt.cur_leftend / (double)apt.total_len)*100.0;
		printProgress((uint)apt.prev_progress, (uint)cur_progress, cout);
		apt.prev_progress = cur_progress;
	}
	printMemUsage();
	cout << "Fix left ends\n";
	FixLeftEnds(ancestor);
	printMemUsage();

	if( debug_aligner )
		validateSuperIntervals(alignment_tree[ancestor].children[0], alignment_tree[ancestor].children[1], ancestor);
	cout << "\ndone.\n";
}

void ProgressiveAligner::FixLeftEnds( node_id_t ancestor )
{
	// fixes all SuperInterval left-end coordinates for nodes below ancestor
	stack< node_id_t > node_stack;
	node_stack.push( ancestor );
	vector<bool> visited( alignment_tree.size(), false );
	while( node_stack.size() > 0 )
	{
		node_id_t cur_node = node_stack.top();
		// visit post-order
		if( !visited[cur_node] )
		{
			for( size_t childI = 0; childI < alignment_tree[cur_node].children.size(); childI++ )
				node_stack.push( alignment_tree[cur_node].children[childI] );
			visited[cur_node] = true;
			continue;
		}
		node_stack.pop();
		if( alignment_tree[cur_node].sequence != NULL )
			continue;	// don't do anything on leaf nodes

		vector< SuperInterval >& siv_list = alignment_tree[cur_node].ordering;
		gnSeqI left_end = 1;
		for( size_t sivI = 0; sivI < siv_list.size(); sivI++ )
		{
			siv_list[sivI].SetLeftEnd(left_end);
			siv_list[sivI].SetLength(siv_list[sivI].reference_iv.AlignmentLength());
			left_end += siv_list[sivI].reference_iv.AlignmentLength();
			CompactGappedAlignment<>* m_cga = dynamic_cast<CompactGappedAlignment<>*>(siv_list[sivI].reference_iv.GetMatches()[0]);
			
			// this one wasn't refined, just move it appropriately
			if( m_cga == NULL || siv_list[sivI].reference_iv.GetMatches().size() > 1 )
			{
				for( uint childI = 0; childI <= 1; childI++ )
				{
					size_t cur_siv = childI == 0 ? alignment_tree[cur_node].ordering[sivI].c1_siv : alignment_tree[cur_node].ordering[sivI].c2_siv;
					if( cur_siv == (std::numeric_limits<size_t>::max)() )
						continue;
					const SuperInterval& c_siv = alignment_tree[ alignment_tree[cur_node].children[childI] ].ordering[ cur_siv ];
					int64 diff = c_siv.LeftEnd() - siv_list[sivI].reference_iv.LeftEnd(childI);
					siv_list[sivI].reference_iv.SetLeftEnd(childI, c_siv.LeftEnd());
					const vector< AbstractMatch* >& matches = siv_list[sivI].reference_iv.GetMatches();
					for( size_t mI = 0; mI < matches.size(); mI++ )
					{
						if( matches[mI]->LeftEnd(childI) != NO_MATCH )
							matches[mI]->SetLeftEnd(childI, matches[mI]->LeftEnd(childI) + diff);
					}
				}

			}else{

				size_t c1_siv = alignment_tree[cur_node].ordering[sivI].c1_siv;
				if( c1_siv != (std::numeric_limits<size_t>::max)() )
				{
					const SuperInterval& c_siv = alignment_tree[ alignment_tree[cur_node].children[0] ].ordering[ c1_siv ];
					m_cga->SetLeftEnd(0, c_siv.LeftEnd());
					siv_list[sivI].reference_iv.SetLeftEnd(0, c_siv.LeftEnd());
					m_cga->SetLength(c_siv.Length(), 0);
					siv_list[sivI].reference_iv.SetLength(c_siv.Length(), 0);
					siv_list[sivI].reference_iv.SetOrientation(0, m_cga->Orientation(0));
				}
				size_t c2_siv = alignment_tree[cur_node].ordering[sivI].c2_siv;
				if( c2_siv != (std::numeric_limits<size_t>::max)() )
				{
					const SuperInterval& c_siv = alignment_tree[ alignment_tree[cur_node].children[1] ].ordering[ c2_siv ];
					m_cga->SetLeftEnd(1, c_siv.LeftEnd());
					siv_list[sivI].reference_iv.SetLeftEnd(1, c_siv.LeftEnd());
					m_cga->SetLength(c_siv.Length(), 1);
					siv_list[sivI].reference_iv.SetLength(c_siv.Length(), 1);
					siv_list[sivI].reference_iv.SetOrientation(1, m_cga->Orientation(1));
				}
			}
			if( debug_cga && m_cga && !m_cga->validate() )
//			if( m_cga && !m_cga->validate() )
				cerr << "oh junkedy\n";

		}
	}
}

/**
 * propagates an inversion of an ancestral SuperInterval to SuperIntervals in descendant nodes
 */
void propagateInvert( PhyloTree< AlignmentTreeNode >& alignment_tree, node_id_t ancestor, size_t ans_siv )
{
	stack< pair< node_id_t, size_t > > node_siv_stack;
	node_siv_stack.push( make_pair(ancestor, ans_siv) );
	while( node_siv_stack.size() > 0 )
	{
		pair< node_id_t, size_t > cur = node_siv_stack.top();
		node_siv_stack.pop();
		node_id_t cur_node = cur.first;
		if( alignment_tree[cur_node].ordering[cur.second].c1_siv != (std::numeric_limits<size_t>::max)() )
			node_siv_stack.push( make_pair( alignment_tree[cur_node].children[0], alignment_tree[cur_node].ordering[cur.second].c1_siv ) );
		if( alignment_tree[cur_node].ordering[cur.second].c2_siv != (std::numeric_limits<size_t>::max)() )
			node_siv_stack.push( make_pair( alignment_tree[cur_node].children[1], alignment_tree[cur_node].ordering[cur.second].c2_siv ) );
		if( cur_node == ancestor )
			continue;	// don't do anything at the ancestor
		if( alignment_tree[cur_node].sequence != NULL )
			continue;	// don't do anything on leaf nodes

		// reverse the homology structure at this node
		Interval& ref_iv = alignment_tree[cur_node].ordering[cur.second].reference_iv;
		vector< AbstractMatch* > matches;
		ref_iv.StealMatches( matches );
		AbstractMatch::orientation o0 = matches[0]->Orientation(0);
		AbstractMatch::orientation o1 = matches[0]->Orientation(1);
		matches[0]->Invert();
		if( o0 != AbstractMatch::undefined )
			matches[0]->SetOrientation(0,o0);
		if( o1 != AbstractMatch::undefined )
			matches[0]->SetOrientation(1,o1);
		ref_iv.SetMatches( matches );
		if( o0 != AbstractMatch::undefined )
		{
			ref_iv.SetOrientation(0,o0);
			ref_iv.SetLeftEnd(0,0);
		}
		if( o1 != AbstractMatch::undefined )
		{
			ref_iv.SetOrientation(1,o1);
			ref_iv.SetLeftEnd(1,0);
		}
	}
}


void ProgressiveAligner::ConstructSuperIntervalFromMSA( node_id_t ancestor, size_t ans_siv, GappedAlignment& gal )
{
	const vector< string >& aln_matrix = GetAlignment( gal, vector< gnSequence* >() );
	stack< pair< node_id_t, size_t > > node_siv_stack;
	node_siv_stack.push( make_pair(ancestor, ans_siv) );
	vector<bool> visited( alignment_tree.size(), false );
	while( node_siv_stack.size() > 0 )
	{
		pair< node_id_t, size_t > cur = node_siv_stack.top();
		node_id_t cur_node = cur.first;
		// visit post-order
		if( !visited[cur_node] )
		{
			if( alignment_tree[cur_node].ordering[cur.second].c1_siv != (std::numeric_limits<size_t>::max)() )
				node_siv_stack.push( make_pair( alignment_tree[cur_node].children[0], alignment_tree[cur_node].ordering[cur.second].c1_siv ) );
			if( alignment_tree[cur_node].ordering[cur.second].c2_siv != (std::numeric_limits<size_t>::max)() )
				node_siv_stack.push( make_pair( alignment_tree[cur_node].children[1], alignment_tree[cur_node].ordering[cur.second].c2_siv ) );
			visited[cur_node] = true;
			continue;
		}
		node_siv_stack.pop();
		if( alignment_tree[cur_node].sequence != NULL )
			continue;	// don't do anything on leaf nodes

		// build a super-interval
		vector< node_id_t > node1_seqs;	/**< the node id's of extant sequences below node 1 */
		vector< node_id_t > node2_seqs;	/**< the node id's of extant sequences below node 2 */
		getAlignedChildren( alignment_tree[cur_node].children[0], node1_seqs );
		getAlignedChildren( alignment_tree[cur_node].children[1], node2_seqs );
		vector< bitset_t > m_aln(2, bitset_t( aln_matrix[0].size(), false ) );
		gnSeqI seqI_len = 0;
		gnSeqI seqJ_len = 0;
		gnSeqI cur_col = 0;
		for( size_t colI = 0; colI < aln_matrix[0].size(); colI++ )
		{
			uint seqI = 0;
			uint seqJ = 0;
			for( ; seqI < node1_seqs.size(); ++seqI )
				if( aln_matrix[node_sequence_map[node1_seqs[seqI]]][colI] != '-' )
					break;
			for( ; seqJ < node2_seqs.size(); ++seqJ )
				if( aln_matrix[node_sequence_map[node2_seqs[seqJ]]][colI] != '-' )
					break;

			if( seqI == node1_seqs.size() && seqJ == node2_seqs.size() )
				continue;	// nothing in this column
			if( seqI != node1_seqs.size() )
			{
				seqI_len++;
				m_aln[0].set(cur_col);
			}
			if( seqJ != node2_seqs.size() )
			{
				seqJ_len++;
				m_aln[1].set(cur_col);
			}
			cur_col++;
		}
		m_aln[0].resize(cur_col);
		m_aln[1].resize(cur_col);
		CompactGappedAlignment<> tmp_cga(m_aln.size(), cur_col);
		CompactGappedAlignment<>* cga = tmp_cga.Copy();
		cga->SetLeftEnd(0, seqI_len > 0 ? 1 : 0);	// at this point we have no idea where the left end should really be
		cga->SetLeftEnd(1, seqJ_len > 0 ? 1 : 0);
		if( cga->LeftEnd(0) != NO_MATCH )
			cga->SetOrientation(0, alignment_tree[cur_node].ordering[cur.second].reference_iv.Orientation(0));
		if( cga->LeftEnd(1) != NO_MATCH )
			cga->SetOrientation(1, alignment_tree[cur_node].ordering[cur.second].reference_iv.Orientation(1));
		cga->SetLength(seqI_len,0);
		cga->SetLength(seqJ_len,1);
		cga->SetAlignment(m_aln);	// do this afterwords so that it can create the bitcount

		// the alignment may need to be reversed if the aligned parent is reverse
		size_t p_siv = alignment_tree[cur_node].ordering[cur.second].parent_siv;
		bool reverse_me = false;
		if( p_siv != (std::numeric_limits<size_t>::max)() )
		{
			size_t p_node = alignment_tree[cur_node].parents[0];
			int p_child = alignment_tree[p_node].children[0] == cur_node ? 0 : 1;
			if( alignment_tree[p_node].ordering[p_siv].reference_iv.Orientation(p_child) == AbstractMatch::reverse )
				reverse_me = true;
		}
		if( reverse_me )
		{
			cga->Invert();
			if( cga->LeftEnd(0) != NO_MATCH )
				cga->SetOrientation(0, alignment_tree[cur_node].ordering[cur.second].reference_iv.Orientation(0));
			if( cga->LeftEnd(1) != NO_MATCH )
				cga->SetOrientation(1, alignment_tree[cur_node].ordering[cur.second].reference_iv.Orientation(1));
			propagateInvert( alignment_tree, cur_node, cur.second );
		}

		alignment_tree[cur_node].ordering[cur.second].reference_iv = Interval();
		vector< AbstractMatch* > am_list(1, cga);
		alignment_tree[cur_node].ordering[cur.second].reference_iv.SetMatches( am_list );
		// set these to zero so they don't interfere with coordinate translation
		alignment_tree[cur_node].ordering[cur.second].reference_iv.SetLeftEnd(0, 0);
		alignment_tree[cur_node].ordering[cur.second].reference_iv.SetLeftEnd(1, 0);
	}
}

typedef boost::tuple<CompactGappedAlignment<>*, vector< bitset_t >*, AbstractMatch* > _sort_tracker_type;

template< class CompType >
class CgaBsComp
{
public:
	CgaBsComp( CompType& c ) : comp(c) {};
	bool operator()( const _sort_tracker_type& a, const _sort_tracker_type& b )
	{
		return comp( a.get<0>(), b.get<0>() );
	}
protected:
	CompType& comp;
};

template< typename MatchVector >
void multFilter( MatchVector& matches, uint mult = 2 )
{
	// apply a multiplicity filter
	size_t cur = 0;
	for( size_t mI = 0; mI < matches.size(); ++mI )
	{
		if( matches[mI]->Multiplicity() == mult )
			matches[cur++] = matches[mI];
		else
			matches[mI]->Free();
	}
	matches.erase(matches.begin()+cur, matches.end());
}

template< typename MatchVector >
void alignedNtCountFilter( MatchVector& matches, uint length )
{
	// require at least some number of aligned pairs in the anchor
	size_t cur = 0;
	for( size_t mI = 0; mI < matches.size(); ++mI )
	{
		size_t len_sum = 0;
		for( size_t seqI = 0; seqI < matches[mI]->SeqCount(); seqI++ )
			if(matches[mI]->LeftEnd(seqI) != NO_MATCH)
				len_sum += matches[mI]->Length(seqI);

		if( len_sum - length > matches[mI]->AlignmentLength() )
			matches[cur++] = matches[mI];
		else
			matches[mI]->Free();
	}
	matches.erase(matches.begin()+cur, matches.end());
}


bool debugging_cltm = false;
void ProgressiveAligner::constructLcbTrackingMatches( 
	node_id_t ancestral_node, 
	vector< AbstractMatch* >& ancestral_matches, 
	vector< LcbTrackingMatch< AbstractMatch* > >& tracking_matches 
	)
{
	node_id_t child_0 = alignment_tree[ancestral_node].children[0];
	node_id_t child_1 = alignment_tree[ancestral_node].children[1];
	// split up matches at descendant's breakpoints
	propagateDescendantBreakpoints( child_0, 0, ancestral_matches );
	propagateDescendantBreakpoints( child_1, 1, ancestral_matches );

	// store alignment bitvectors for each match...
	vector< bitset_t > bs_tmp(alignment_tree.size());
	vector< vector< bitset_t > > bs(ancestral_matches.size(), bs_tmp);
	vector< _sort_tracker_type > cga_list;
	// initialize alignment bitvectors
	for( size_t mI = 0; mI < ancestral_matches.size(); mI++ )
	{
		vector< bitset_t > aln( alignment_tree.size(), bitset_t(ancestral_matches[mI]->AlignmentLength() ) );
		swap( bs[mI], aln );
		ancestral_matches[mI]->GetAlignment(aln);
		swap( bs[mI][child_0], aln[0] );
		swap( bs[mI][child_1], aln[1] );
		CompactGappedAlignment<> c(alignment_tree.size(),0);
		c.SetLeftEnd(child_0, ancestral_matches[mI]->LeftEnd(0));
		c.SetOrientation(child_0, ancestral_matches[mI]->Orientation(0));
		c.SetLength(ancestral_matches[mI]->Length(0), child_0);
		c.SetLeftEnd(child_1, ancestral_matches[mI]->LeftEnd(1));
		c.SetOrientation(child_1, ancestral_matches[mI]->Orientation(1));
		c.SetLength(ancestral_matches[mI]->Length(1), child_1);
		cga_list.push_back(make_tuple(c.Copy(), &bs[mI], ancestral_matches[mI]));
	}

	stack<node_id_t> node_stack;
	node_stack.push(child_0);
	node_stack.push(child_1);
	while(node_stack.size() > 0)
	{
		node_id_t cur_node = node_stack.top();
		node_stack.pop();
		if( alignment_tree[cur_node].children.size() == 0 )
			continue;
		node_stack.push(alignment_tree[cur_node].children[0]);
		node_stack.push(alignment_tree[cur_node].children[1]);

		// do processing for cur_node...
		// 1. determine which interval in the current node each match falls into
		// 2. determine the offset of this match in that interval
		// 3. translate with that interval

		vector< SuperInterval >& siv_list = alignment_tree[cur_node].ordering;
		SingleStartComparator< CompactGappedAlignment<> > ssc(cur_node);
		CgaBsComp< SingleStartComparator< CompactGappedAlignment<> > > comp( ssc );
		sort(cga_list.begin(), cga_list.end(), comp);
		size_t mI = 0;
		size_t sivI = 0;
		while( mI < cga_list.size() && sivI < siv_list.size() )
		{
			CompactGappedAlignment<>* cur_match = cga_list[mI].get<0>();
			if( cur_match->Start(cur_node) == 0 )
			{
				mI++;
				continue;	// this one doesn't match in this lineage!!
			}
			if( cur_match->LeftEnd(cur_node) >= siv_list[sivI].LeftEnd() + siv_list[sivI].Length() )
			{
				sivI++;
				continue;
			}

			if( cur_match->LeftEnd(cur_node) + cur_match->Length(cur_node) > 
				siv_list[sivI].LeftEnd() + siv_list[sivI].Length() )
			{
				cerr << "doesn't fit\n";
				cerr << "cga_list[" << mI << "]->LeftEnd(" << cur_node << "): " << cur_match->LeftEnd(cur_node) << endl;
				cerr << "cga_list[" << mI << "]->RightEnd(" << cur_node << "): " << cur_match->RightEnd(cur_node) << endl;
				breakHere();
			}

			// extract the region of the siv matched by the current match
			CompactGappedAlignment<>* siv_cga = dynamic_cast<CompactGappedAlignment<>*>(siv_list[sivI].reference_iv.GetMatches()[0]);
			if( siv_list[sivI].reference_iv.GetMatches().size() > 1 )
				siv_cga = NULL;
			if( siv_cga == NULL )
			{
				CompactGappedAlignment<> tmp_cga;
				siv_cga = tmp_cga.Copy();
				*siv_cga = CompactGappedAlignment<>(siv_list[sivI].reference_iv);
				vector<AbstractMatch*> tmp_matches(1,siv_cga);
				siv_list[sivI].reference_iv.SetMatches(tmp_matches);
			}
			CompactGappedAlignment<> new_cga;
			siv_cga->copyRange(new_cga, cur_match->LeftEnd(cur_node) - siv_list[sivI].LeftEnd(), cur_match->Length(cur_node));
			if( cur_match->Orientation(cur_node) == AbstractMatch::reverse )
				new_cga.Invert();
			if( new_cga.Multiplicity() == 0 )
			{
				cerr << "impossible!  there's no match!\n";
				genome::breakHere();
			}
			// set the leftend in cga_list
			for( uint cur_child = 0; cur_child < 2; cur_child++ )
			{
				node_id_t sweet_child = alignment_tree[cur_node].children[cur_child];
				cur_match->SetLeftEnd(sweet_child, new_cga.LeftEnd(cur_child));
				if( new_cga.LeftEnd(cur_child) != NO_MATCH )
				{
					cur_match->SetOrientation(sweet_child, new_cga.Orientation(cur_child));
					cur_match->SetLength(new_cga.Length(cur_child), sweet_child);
				}
			}

			// prepare a cga for translation
			CompactGappedAlignment<> c(1,(*cga_list[mI].get<1>())[cur_node].size());
			c.SetLeftEnd(0,1);
			c.SetLength((*cga_list[mI].get<1>())[cur_node].count(),0);
			vector<bitset_t> bivouac(1, (*cga_list[mI].get<1>())[cur_node]);
			c.SetAlignment(bivouac);

			// now translate each child
			for( uint cur_child = 0; cur_child < 2; cur_child++ )
			{
				if( new_cga.Orientation(cur_child) == AbstractMatch::undefined )
					continue;
				CompactGappedAlignment<> cga_tmp = new_cga;
				cga_tmp.SetStart(cur_child, 1);
				c.translate(cga_tmp, cur_child, 0, false);
				// adjust for end-gaps
				bitset_t bs = (cga_tmp.GetAlignment())[cur_child];
				bs.resize(c.GetAlignment()[0].size(), false);
				bs <<= c.GetAlignment()[0].find_first();
				node_id_t sweet_child = alignment_tree[cur_node].children[cur_child];
				swap( (*cga_list[mI].get<1>())[sweet_child], bs );
				for( size_t testI = 0; testI < cga_tmp.SeqCount(); ++testI )
				{
					if( ((*cga_list[mI].get<1>())[testI].size() != 0 && (*cga_list[mI].get<1>())[testI].size() != (*cga_list[mI].get<1>())[sweet_child].size() ) )
					{
						cerr << "bj0rk3l\n";
						genome::breakHere();
					}
				}
			}

			debugging_cltm = false;
			mI++;	// advance to the next match
		}
	}
	tracking_matches.resize( cga_list.size() );
	// finally, construct CompactGappedAlignments out of the bitsets
	for( size_t bsI = 0; bsI < cga_list.size(); ++bsI )
	{
		cga_list[bsI].get<0>()->SetAlignment(*cga_list[bsI].get<1>());
		cga_list[bsI].get<0>()->validate();
		TrackingMatch& ltm = tracking_matches[bsI];
		ltm.node_match = cga_list[bsI].get<0>();
		ltm.original_match = cga_list[bsI].get<2>();
		ltm.match_id = bsI;

		bool found_extant = false;
		for( size_t i = 0; i < alignment_tree.size()-1; ++i )
		{
			size_t im = node_sequence_map[i];
			if( im == (std::numeric_limits<size_t>::max)() )
				continue;
			if( ltm.node_match->LeftEnd(i) != NO_MATCH )
				found_extant = true;
		}
		if( !found_extant )
		{
			cout << "orig aln len: " << ltm.original_match->AlignmentLength() << endl;
			cout << "orig lend 0: " << ltm.original_match->Start(0) << endl;
			cout << "orig lend 1: " << ltm.original_match->Start(1) << endl;
			cout << "orig length 0: " << ltm.original_match->Length(0) << endl;
			cout << "orig length 1: " << ltm.original_match->Length(1) << endl;

			cerr << "this is an ungrounded match!!!\n";
			genome::breakHere();
		}
	}
}

size_t countUnrefined( PhyloTree< AlignmentTreeNode >& alignment_tree, node_id_t ancestor )
{
	stack< node_id_t > node_stack;
	node_stack.push(ancestor);
	size_t unrefined_count = 0;
	while( node_stack.size() > 0 )
	{
		node_id_t cur_node = node_stack.top();
		node_stack.pop();
		if( alignment_tree[cur_node].children.size() > 0 )
			for( size_t childI = 0; childI < alignment_tree[cur_node].children.size(); ++childI )
				node_stack.push( alignment_tree[cur_node].children[childI] );
		if( !alignment_tree[cur_node].refined )
			unrefined_count++;
	}
	return unrefined_count;
}

void markAsRefined( PhyloTree< AlignmentTreeNode >& alignment_tree, node_id_t ancestor )
{
	stack< node_id_t > node_stack;
	node_stack.push(ancestor);
	size_t refined_count = 0;
	while( node_stack.size() > 0 )
	{
		node_id_t cur_node = node_stack.top();
		node_stack.pop();
		if( alignment_tree[cur_node].children.size() > 0 )
			for( size_t childI = 0; childI < alignment_tree[cur_node].children.size(); ++childI )
				node_stack.push( alignment_tree[cur_node].children[childI] );
		alignment_tree[cur_node].refined = true;
	}
	alignment_tree[ancestor].refined = false;
}



void ProgressiveAligner::pairwiseScoreTrackingMatches( 
						std::vector< TrackingMatch >& tracking_matches, 
						std::vector<node_id_t>& node1_descendants,
						std::vector<node_id_t>& node2_descendants,
						boost::multi_array< double, 3 >& tm_score_array)
{
	tm_score_array.resize( boost::extents[tracking_matches.size()][node1_descendants.size()][node2_descendants.size()] );
	for( size_t mI = 0; mI < tracking_matches.size(); ++mI )
	{
		TrackingMatch* cur_match = &tracking_matches[mI];
		AbstractMatch* node_match = cur_match->node_match;
		for( size_t nI = 0; nI < node1_descendants.size(); ++nI )
		{
			if( node_sequence_map[node1_descendants[nI]] == (std::numeric_limits<uint>::max)()  ||
				node_match->LeftEnd(node1_descendants[nI]) == NO_MATCH )
				continue;
			for( size_t nJ = 0; nJ < node2_descendants.size(); ++nJ )
			{
				if( node_sequence_map[node2_descendants[nJ]] == (std::numeric_limits<uint>::max)() ||
					node_match->LeftEnd(node2_descendants[nJ]) == NO_MATCH )
					continue;	// not extant or no match between this pair

				node_id_t cur_n1 = node1_descendants[nI];
				node_id_t cur_n2 = node2_descendants[nJ];
				size_t nsmI = node_sequence_map[cur_n1];
				size_t nsmJ = node_sequence_map[cur_n2];
				PairwiseMatchAdapter pma( node_match, cur_n1, cur_n2 );
				vector< AbstractMatch* > lcb_vect( 1, &pma );
				vector< gnSequence* > ex_seqs(2);
				ex_seqs[0] = alignment_tree[ cur_n1 ].sequence;
				ex_seqs[1] = alignment_tree[ cur_n2 ].sequence;

				tm_score_array[mI][nI][nJ] = GetPairwiseAnchorScore(lcb_vect, ex_seqs, subst_scoring, sol_list[nsmI], sol_list[nsmJ]);
			}
		}
	}
	computeAvgAncestralMatchScores(tracking_matches, node1_descendants, node2_descendants, tm_score_array);
}

void ProgressiveAligner::computeAvgAncestralMatchScores( 
						std::vector< TrackingMatch >& tracking_matches, 
						std::vector<node_id_t>& node1_descendants,
						std::vector<node_id_t>& node2_descendants,
						boost::multi_array< double, 3 >& tm_score_array)
{
	// now build up the consensus (ancestral) match scores and bp distances
	for( uint nodeI = 0; nodeI < node1_descendants.size(); nodeI++ )
	{
		for( uint nodeJ = 0; nodeJ < node2_descendants.size(); nodeJ++ )
		{
			node_id_t n1 = node1_descendants[nodeI];
			node_id_t n2 = node2_descendants[nodeJ];

			vector<node_id_t> n1_ext;
			vector<node_id_t> n2_ext;
			getAlignedChildren(n1, n1_ext);
			getAlignedChildren(n2, n2_ext);
			if( n1_ext.size() == 1 && n2_ext.size() == 1 )
				continue;	// this node has two extant nodes below it and was already scored

			// map the nodes in n1_ext to their indices in n1_descendants
			vector< node_id_t > n1_ext_map(n1_ext.size());
			for( size_t i = 0; i < n1_ext.size(); ++i )
			{
				vector< node_id_t >::iterator iter = std::find( node1_descendants.begin(), node1_descendants.end(), n1_ext[i] );
				n1_ext_map[i] = iter - node1_descendants.begin();
			}
			vector< node_id_t > n2_ext_map(n2_ext.size());
			for( size_t i = 0; i < n2_ext.size(); ++i )
			{
				vector< node_id_t >::iterator iter = std::find( node2_descendants.begin(), node2_descendants.end(), n2_ext[i] );
				n2_ext_map[i] = iter - node2_descendants.begin();
			}

			// compute scores for all matches at this node
			for( size_t mI = 0; mI < tracking_matches.size(); ++mI )
			{
				uint tally = 0;
				double score_sum = 0;
				for( size_t i = 0; i < n1_ext.size(); ++i )
				{
					if( tracking_matches[mI].node_match->LeftEnd(n1_ext[i]) == NO_MATCH )
						continue;
					for( size_t j = 0; j < n2_ext.size(); ++j )
					{
						if( tracking_matches[mI].node_match->LeftEnd(n2_ext[j]) == NO_MATCH )
							continue;
						++tally;
						score_sum += tm_score_array[mI][n1_ext_map[i]][n2_ext_map[j]];
					}
				}
				if( tally > 0 )
					tm_score_array[mI][nodeI][nodeJ] = score_sum / (double)tally;
			}
		}
	}
}


void ProgressiveAligner::computeInternalNodeDistances( 
						boost::multi_array<double, 2>& bp_dist_mat, 
						boost::multi_array<double, 2>& cons_dist_mat, 
						std::vector<node_id_t>& node1_descendants,
						std::vector<node_id_t>& node2_descendants)
{
	// bp distances for the current node.
	bp_dist_mat.resize(boost::extents[node1_descendants.size()][node2_descendants.size()]);
	cons_dist_mat.resize(boost::extents[node1_descendants.size()][node2_descendants.size()]);
	for( size_t nI = 0; nI < node1_descendants.size(); ++nI )
	{
		if( node_sequence_map[node1_descendants[nI]] == (std::numeric_limits<uint>::max)() )
			continue;
		for( size_t nJ = 0; nJ < node2_descendants.size(); ++nJ )
		{
			if( node_sequence_map[node2_descendants[nJ]] == (std::numeric_limits<uint>::max)() )
				continue;
			size_t i = node_sequence_map[node1_descendants[nI]];
			size_t j = node_sequence_map[node2_descendants[nJ]];
			bp_dist_mat[nI][nJ] = this->bp_distance[i][j];
			cons_dist_mat[nI][nJ] = this->conservation_distance[i][j];
		}
	}

	// now build up the consensus (ancestral) bp distances
	for( uint nodeI = 0; nodeI < node1_descendants.size(); nodeI++ )
	{
		for( uint nodeJ = 0; nodeJ < node2_descendants.size(); nodeJ++ )
		{
			node_id_t n1 = node1_descendants[nodeI];
			node_id_t n2 = node2_descendants[nodeJ];

			vector<node_id_t> n1_ext;
			vector<node_id_t> n2_ext;
			getAlignedChildren(n1, n1_ext);
			getAlignedChildren(n2, n2_ext);
			if( n1_ext.size() == 1 && n2_ext.size() == 1 )
				continue;	// this node has two extant nodes below it, so already has a dist

			// map the nodes in n1_ext to their indices in n1_descendants
			vector< node_id_t > n1_ext_map(n1_ext.size());
			for( size_t i = 0; i < n1_ext.size(); ++i )
			{
				vector< node_id_t >::iterator iter = std::find( node1_descendants.begin(), node1_descendants.end(), n1_ext[i] );
				n1_ext_map[i] = iter - node1_descendants.begin();
			}
			vector< node_id_t > n2_ext_map(n2_ext.size());
			for( size_t i = 0; i < n2_ext.size(); ++i )
			{
				vector< node_id_t >::iterator iter = std::find( node2_descendants.begin(), node2_descendants.end(), n2_ext[i] );
				n2_ext_map[i] = iter - node2_descendants.begin();
			}

			// compute average bp distance
			for( size_t i = 0; i < n1_ext.size(); ++i )
			{
				for( size_t j = 0; j < n2_ext.size(); ++j )
				{
					bp_dist_mat[nodeI][nodeJ] += bp_dist_mat[n1_ext_map[i]][n2_ext_map[j]];
					cons_dist_mat[nodeI][nodeJ] += cons_dist_mat[n1_ext_map[i]][n2_ext_map[j]];
				}
			}
			bp_dist_mat[nodeI][nodeJ] /= (double)(n1_ext.size() * n2_ext.size());
			cons_dist_mat[nodeI][nodeJ] /= (double)(n1_ext.size() * n2_ext.size());
		}
	}

}

double computeID( GappedAlignment& gal, size_t seqI, size_t seqJ )
{
	const vector< string >& aln_mat = GetAlignment( gal, vector< gnSequence* >(gal.SeqCount(), NULL ));
	double id = 0;
	double possible = 0;
	for( size_t colI = 0; colI < gal.AlignmentLength(); colI++ )
	{
		if( aln_mat[seqI][colI] == '-' || aln_mat[seqJ][colI] == '-' )
			continue;
		if( toupper(aln_mat[seqI][colI]) == toupper(aln_mat[seqJ][colI]))
			id++;
		possible++;
	}
	return id / possible;
}


//
//
// different option -- just pick a representative from leaf(A) and leaf(B) to translate
void ProgressiveAligner::getRepresentativeAncestralMatches( const vector< node_id_t > node1_seqs, const vector< node_id_t > node2_seqs, node_id_t node1, node_id_t node2, node_id_t ancestor, std::vector< AbstractMatch* >& ancestral_matches )
{
	// for each match, extract a representative match from any pair of genomes in node1_seqs and node2_seqs
	// translate up the resulting set of matches and eliminate overlaps
	vector< AbstractMatch* > cur_matches;
	boost::multi_array< vector< AbstractMatch* >, 2 > seq_matches( boost::extents[node1_seqs.size()][node2_seqs.size()] );
	for( size_t mI = 0; mI < original_ml.size(); mI++ )
	{
		for( uint seqI = 0; seqI < node1_seqs.size(); seqI++ )
		{
			uint ii = this->node_sequence_map[node1_seqs[seqI]];
			if( original_ml[mI]->LeftEnd(ii) == NO_MATCH )
				continue;

			for( uint seqJ = 0; seqJ < node2_seqs.size(); seqJ++ )
			{
				uint jj = this->node_sequence_map[node2_seqs[seqJ]];
				if( original_ml[mI]->LeftEnd(jj) == NO_MATCH )
					continue;
				Match mm( 2 );
				Match* new_m = mm.Copy();
				new_m->SetStart( 0, original_ml[mI]->Start(ii));
				new_m->SetStart( 1, original_ml[mI]->Start(jj));
				new_m->SetLength(original_ml[mI]->Length());
				if( new_m->Start(0) < 0 )
					new_m->Invert();	// assign reference orientation to seq 0
				seq_matches[seqI][seqJ].push_back( new_m );
				break;
			}
			break;
		}
	}
	for( uint seqI = 0; seqI < node1_seqs.size(); seqI++ )
		for( uint seqJ = 0; seqJ < node2_seqs.size(); seqJ++ )
		{
			translateGappedCoordinates( seq_matches[seqI][seqJ], 0, node1_seqs[seqI], node1 );
			translateGappedCoordinates( seq_matches[seqI][seqJ], 1, node2_seqs[seqJ], node2 );
			ancestral_matches.insert( ancestral_matches.end(), seq_matches[seqI][seqJ].begin(), seq_matches[seqI][seqJ].end() );
		}

	EliminateOverlaps_v2( ancestral_matches, true );
}

int cachecomp( const void* e1, const void* e2 )
{
	bool a = mems::cache_comparator(*(search_cache_t*)e1, *(search_cache_t*)e2);
	bool b = mems::cache_comparator(*(search_cache_t*)e2, *(search_cache_t*)e1);
	if(!a && !b)
		return 0;
	return a ? -1 : 1;
}

void ProgressiveAligner::alignProfileToProfile( node_id_t node1, node_id_t node2, node_id_t ancestor )
{
	// 1) find all pairwise matches
	// 2) convert to pairwise matches among the ancestral sequences
	//    - delete inconsistently aligned regions?
	// 3) perform greedy b.p. elimination on the pairwise matches
	// 4) extend LCBs
	// 5)  if total alignment weight hasn't changed, go to (8)
	// 6) search for additional matches between each match among extant sequences
	// 7) go back to 2
	// 8) perform a MUSCLE/Clustal alignment of each intervening region

	vector< node_id_t > node1_seqs;	/**< the node id's of extant sequences below node 1 */
	vector< node_id_t > node2_seqs;	/**< the node id's of extant sequences below node 2 */
	getAlignedChildren( node1, node1_seqs );
	getAlignedChildren( node2, node2_seqs );

	uint seqI, seqJ;
	gnSeqI prev_ancestral_seq_len = (std::numeric_limits<gnSeqI>::max)();

	printMemUsage();
	cout << "get ancestral matches\n";

	Matrix<MatchList> pairwise_matches( node1_seqs.size(), node2_seqs.size() );
//	getPairwiseMatches( node1_seqs, node2_seqs, pairwise_matches );
	vector< AbstractMatch* > anc_pairwise_matches;
	getRepresentativeAncestralMatches( node1_seqs, node2_seqs, node1, node2, ancestor, anc_pairwise_matches );
	printMemUsage();
	
	PhyloTree< AlignmentTreeNode > aln_tree_backup;

	/** A cache of regions that were searched in the previous round of recursion */
	Matrix< std::vector< search_cache_t > > search_cache_db(node1_seqs.size(), node2_seqs.size());
	double prev_anchoring_score = -(std::numeric_limits<double>::max)();
	double cur_anchoring_score = -(std::numeric_limits<double>::max)();

	while(true)
	{
		vector<AbstractMatch*> ancestral_matches;
		if( anc_pairwise_matches.size() > 0 )
		{
			ancestral_matches.insert( ancestral_matches.begin(), anc_pairwise_matches.begin(), anc_pairwise_matches.end() );
			anc_pairwise_matches.clear();
		}
		else
		{
			// part 2, construct pairwise matches to the ancestral sequence
			// A)  for each pairwise match, translate its
			//     coordinates to the ancestral genome
			//	   -- try to use translateCoordinates
			//     -- build a translation table for translateCoordinates

			for( seqI = 0; seqI < node1_seqs.size(); seqI++ )
			{
				for( seqJ = 0; seqJ < node2_seqs.size(); seqJ++ )
				{
					cout << node_sequence_map[node1_seqs[seqI]] << "," << node_sequence_map[node2_seqs[seqJ]] << " has " << pairwise_matches(seqI,seqJ).size() << " pairwise matches\n";
					cout.flush();

					vector< AbstractMatch* > am_list( pairwise_matches(seqI, seqJ).begin(), pairwise_matches(seqI, seqJ).end() );
					pairwise_matches(seqI, seqJ).clear();
					translateGappedCoordinates( am_list, 1, node2_seqs[seqJ], node2 );
					translateGappedCoordinates( am_list, 0, node1_seqs[seqI], node1 );
					ancestral_matches.insert( ancestral_matches.end(), am_list.begin(), am_list.end() );
				}
			}
		}
		// include any matches from a previous iteration of this loop
		for( size_t aI = 0; aI < alignment_tree[ancestor].ordering.size(); aI++ )
		{
			Interval& ref_iv = alignment_tree[ancestor].ordering[aI].reference_iv;
			if( ref_iv.Multiplicity() == 2 )
				for( size_t mI = 0; mI < ref_iv.GetMatches().size(); mI++ )
					if( ref_iv.GetMatches()[mI]->Multiplicity() > 1 )
						ancestral_matches.push_back( ref_iv.GetMatches()[mI]->Copy() );
		}

		// set seq 0 to forward ref. orientation
		for( size_t mI = 0; mI < ancestral_matches.size(); ++mI )
			if( ancestral_matches[mI]->Start(0) < 0 )
				ancestral_matches[mI]->Invert();

		// eliminate overlaps as they correspond to inconsistently or
		// multiply aligned regions
		EliminateOverlaps_v2( ancestral_matches );
		
		multFilter( ancestral_matches );

		vector< vector< AbstractMatch* > > LCB_list;
		vector< LCB > adjacencies;
		vector< gnSeqI > breakpoints;

		if( !collinear_genomes )
		{
			cout << "Performing Sum-of-pairs Greedy Breakpoint Elimination\n";
			cout.flush();
			// project the pairwise matches at this node to all-possible pairs matches at descendant nodes
			// keep a mapping of ancestral to extant matches so that when an ancestral match gets removed
			// the match among extant nodes also gets removed
			// how should candidate matches to remove be generated?
			// one possibility is to remove entire ancestral LCBs...  this may be problematic since ancestral
			// LCBs don't correspond to the pairwise LCBs thus an ancestral LCB could be removed with no useful
			// change in alignment score
			//
			//
			// translate the matches into LcbTrackingMatches
			printMemUsage();
			cout << "construct LCB tracking matches\n";
			vector< TrackingMatch > tracking_matches;
			boost::multi_array< size_t, 3 > tm_lcb_id_array;
			boost::multi_array< double, 3 > tm_score_array;
			constructLcbTrackingMatches( ancestor, ancestral_matches, tracking_matches );

			cout << "There are " << tracking_matches.size() << " tracking matches\n";
			size_t used_components = 0;
			for( size_t tmI = 0; tmI < tracking_matches.size(); ++tmI )
			{
				for( uint ssI = 0; ssI < tracking_matches[tmI].node_match->SeqCount(); ++ssI )
					if( tracking_matches[tmI].node_match->LeftEnd(ssI) != NO_MATCH )
						used_components++;
			}
			size_t total_components = tracking_matches.size() == 0 ? 0 : tracking_matches.size() * tracking_matches[0].node_match->SeqCount();
			cout << "There are " << used_components << " / " << total_components << " components used\n";

			vector<node_id_t> node1_descendants;
			vector<node_id_t> node2_descendants;
			if( scoring_scheme == ExtantSumOfPairsScoring )
			{
				node1_descendants = node1_seqs;
				node2_descendants = node2_seqs;
			}else{
				getDescendants(alignment_tree, node1, node1_descendants);
				getDescendants(alignment_tree, node2, node2_descendants);
			}

			//
			// score the matches
			//
			printMemUsage();
			cout << "init tracking match LCB tracking\n";
			initTrackingMatchLCBTracking( tracking_matches, node1_descendants.size(), node2_descendants.size(), tm_lcb_id_array );
			printMemUsage();
			cout << "pairwise score tracking matches\n";
			pairwiseScoreTrackingMatches( tracking_matches, node1_descendants, node2_descendants, tm_score_array );
			printMemUsage();

			// compute bp distances for the current node.
			// ancestral nodes take the average distance of extant nodes
			boost::multi_array<double, 2> bp_dist_mat;
			boost::multi_array<double, 2> cons_dist_mat;
			computeInternalNodeDistances( bp_dist_mat, cons_dist_mat, node1_descendants, node2_descendants);

			vector< TrackingMatch* > t_matches(tracking_matches.size());
			for( size_t mI = 0; mI < tracking_matches.size(); ++mI )
				t_matches[mI] = &tracking_matches[mI];

			// now sort these out into pairwise LCBs
			cout << "get pairwise LCBs\n";
			size_t pair_lcb_count = 0;
			PairwiseLCBMatrix pairwise_adj_mat(boost::extents[node1_descendants.size()][node2_descendants.size()]);
			for( uint nodeI = 0; nodeI < node1_descendants.size(); nodeI++ )
				for( uint nodeJ = 0; nodeJ < node2_descendants.size(); nodeJ++ )
				{
					getPairwiseLCBs( node1_descendants[nodeI], node2_descendants[nodeJ], nodeI, nodeJ, t_matches, pairwise_adj_mat[nodeI][nodeJ], tm_score_array, tm_lcb_id_array );
					pair_lcb_count += pairwise_adj_mat[nodeI][nodeJ].size();
				}
			cout << "there are " << pair_lcb_count << " pairwise LCBs\n";
			printMemUsage();

			sort( t_matches.begin(), t_matches.end() );

			// other possibility, choose pairwise LCBs to remove.  a score improvement is always guaranteed
			// compute LCBs among descendant nodes
			// this is a good idea.  it factors out ancestral breakpoint decisions entirely
			// need a data structure to track all pairwise LCBs that contain a given match
			// template <class MatchType>
			// class LcbTrackingMatch <MatchType> 
			// { 
			// public:
			//  MatchType node_match;
			//	boost::multi_array< size_t, 2 > lcb_id;
			// }
			// all pairwise LCBs would be evaluated for removal and the one that provides the greatest
			// overall score improvement gets removed.
			// upon removal, matches associated with that LCB would get removed, and any LCBs in other 
			// genomes would get removed if they no longer had any matches
			// to pull this off, the LCB struct needs to store the set of matches directly
			// 
			// but what about small cycles that appear only in 3 or more-way comparisons?  are these
			// important?  umm, yeah, but only if you believe in evolution.
			// 
			// so here's the dilly-oh: score against the ancestral ordering(s) *and* all pairwise orderings
			// for an ancestor.  ancestor contributes the sum of all descendants to the score and breakpoints
			// are penalized as the sum of /participating/ descendants.  a descendant is participating
			// if it has some matching region defined within the LCB and if removal of that matching region
			// eliminates a breakpoint in the pairwise comparison
			cout << "scaling bp penalty by conservation weight:\n";
			print2d_matrix(cons_dist_mat, cout);
			cout << "\n\nscaling bp penalty by bp weight: \n";
			print2d_matrix(bp_dist_mat, cout);
			cout << "\nGreedy BPE\n";
			vector< TrackingMatch* > final;
			if(scoring_scheme == AncestralScoring)
			{
				vector<node_id_t>::iterator d1_iter = std::find( node1_descendants.begin(), node1_descendants.end(), node1 );
				vector<node_id_t>::iterator d2_iter = std::find( node2_descendants.begin(), node2_descendants.end(), node2 );
				size_t d1_index = d1_iter - node1_descendants.begin();
				size_t d2_index = d2_iter - node2_descendants.begin();
				EvenFasterSumOfPairsBreakpointScorer spbs( breakpoint_penalty, min_breakpoint_penalty, bp_dist_mat, cons_dist_mat, 
					t_matches, pairwise_adj_mat, node1_descendants, node2_descendants, 
					tm_score_array, tm_lcb_id_array, d1_index, d1_index+1, d2_index, d2_index+1 );
				cur_anchoring_score = greedySearch( spbs );
				final = spbs.getResults();
			}else{
				EvenFasterSumOfPairsBreakpointScorer spbs( breakpoint_penalty, min_breakpoint_penalty, bp_dist_mat, cons_dist_mat, 
					t_matches, pairwise_adj_mat, node1_descendants, node2_descendants, 
					tm_score_array, tm_lcb_id_array, 0, node1_descendants.size(), 0, node2_descendants.size() );
				cur_anchoring_score = greedySearch( spbs );
				final = spbs.getResults();
			}
			cout << "done\n";

			// free memory used by pairwise projections
			for( size_t mI = 0; mI < tracking_matches.size(); ++mI )
				tracking_matches[mI].node_match->Free();

			ancestral_matches.clear();

			// free memory from deleted matches here
			std::sort(final.begin(), final.end());
			vector< TrackingMatch* > deleted_t_matches( t_matches.size(), NULL );
			std::set_difference( t_matches.begin(), t_matches.end(), final.begin(), final.end(), deleted_t_matches.begin() );
			for( size_t delI = 0; delI < deleted_t_matches.size(); ++delI )
			{
				if( deleted_t_matches[delI] == NULL )
					break;
				deleted_t_matches[delI]->original_match->Free();
			}

			// convert back to an LCB list
			vector< AbstractMatch* > new_matches(final.size());
			for( size_t mI = 0; mI < final.size(); ++mI )
				new_matches[mI] = final[mI]->original_match;

			IdentifyBreakpoints( new_matches, breakpoints );
			ComputeLCBs_v2( new_matches, breakpoints, LCB_list );

		} // end if !collinear
		else
		{	// if we are assuming all genomes are collinear, then we don't need the 
			// sophisticated pairwise breakpoint scoring and can get by with simple breakpoint
			// penalties
			IdentifyBreakpoints( ancestral_matches, breakpoints );
			ComputeLCBs_v2( ancestral_matches, breakpoints, LCB_list );

			vector< double > lcb_scores( LCB_list.size() );
			double score_sum = 100;	// anything > 0 would work.  this will be the breakpoint penalty
			for( size_t lcbI = 0; lcbI < LCB_list.size(); ++lcbI )
			{
				lcb_scores[lcbI] = SimpleGetLCBCoverage( LCB_list[lcbI] );
				score_sum += lcb_scores[lcbI];
			}

			computeLCBAdjacencies_v3( LCB_list, lcb_scores, adjacencies );

			// want to eliminate all breakpoints
			SimpleBreakpointScorer wbs( adjacencies, score_sum, true );
			cur_min_coverage = greedyBreakpointElimination_v4( adjacencies, lcb_scores, wbs, NULL, false );
			vector<AbstractMatch*> deleted_matches;
			filterMatches_v2( adjacencies, LCB_list, lcb_scores, deleted_matches );
			for( size_t delI = 0; delI < deleted_matches.size(); ++delI )
				deleted_matches[delI]->Free();
		}
		printMemUsage();

		ancestral_matches.clear();

		cout << "Arrived at " << LCB_list.size() << " intervals\n";
		// create an ancestral ordering
		vector< Interval* > pairwise_intervals;
		Interval tmp_iv;
		for( size_t lcbI = 0; lcbI < LCB_list.size(); lcbI++ )
		{
			pairwise_intervals.push_back( tmp_iv.Copy() );
			pairwise_intervals.back()->SetMatches( LCB_list[lcbI] );
		}
		LCB_list.clear();

		vector<gnSeqI> seq_lengths = vector<gnSeqI>(2,0);
		for( size_t aI = 0; aI < alignment_tree[node1].ordering.size(); ++aI )
			seq_lengths[0] += alignment_tree[node1].ordering[aI].Length();
		for( size_t aI = 0; aI < alignment_tree[node2].ordering.size(); ++aI )
			seq_lengths[1] += alignment_tree[node2].ordering[aI].Length();

		cout << "Adding unaligned intervals\n";
		addUnalignedIntervals_v2(pairwise_intervals, set<uint>(), seq_lengths);

		cout << "addUnalignedIntervals yields " << pairwise_intervals.size() << " intervals\n";

		bool borked = false;
		if(debug_aligner)
			borked = validatePairwiseIntervals(node1, node2, pairwise_intervals);

		// merge unaligned intervals
		cout << "Merging unaligned intervals\n";
		cout.flush();
		vector<Interval*> new_list1;
		vector<Interval*> merged_intervals;
		mergeUnalignedIntervals( 1, pairwise_intervals, new_list1 );
		mergeUnalignedIntervals( 0, new_list1, merged_intervals );
		cout << "Marbling gaps\n";
		cout.flush();
		for( size_t ivI = 0; ivI < merged_intervals.size(); ivI++ )
			merged_intervals[ivI]->Marble(50);

		cout << "Propagating descendant breakpoints\n";

		// split up intervals at descendant's breakpoints
		propagateDescendantBreakpoints( node1, 0, merged_intervals );
		propagateDescendantBreakpoints( node2, 1, merged_intervals );

		cout << "descendant 0(" << node1 << ") has " << alignment_tree[node1].ordering.size() << " intervals\n";
		cout << "descendant 1(" << node2 << ") has " << alignment_tree[node2].ordering.size() << " intervals\n";
		cout << "propagateDescendantBreakpoints yields " << merged_intervals.size() << " intervals\n";

		if(debug_aligner)
			borked = validatePairwiseIntervals(node1, node2, merged_intervals);
		cout << "Creating ancestral ordering\n";
		alignment_tree[ancestor].ordering.clear();
		createAncestralOrdering( merged_intervals, alignment_tree[ancestor].ordering );
		for( size_t ivI = 0; ivI < merged_intervals.size(); ivI++ )
			merged_intervals[ivI]->Free();
		merged_intervals.clear();	// free up some memory

		if(debug_aligner)
			validateSuperIntervals( node1, node2, ancestor );

		// if we're not making any progress then bail out...
		gnSeqI cur_ancestral_seq_len = 0;
		for( size_t aI = 0; aI < alignment_tree[ancestor].ordering.size(); aI++ )
			cur_ancestral_seq_len += alignment_tree[ancestor].ordering[aI].Length();

		if( !collinear_genomes )
			cout << "Previous anchoring score: " << prev_anchoring_score << ", new anchor score: " << cur_anchoring_score << endl;
		else
			cout << "Prev alignment len: " << prev_ancestral_seq_len << ", new alignment length: " << cur_ancestral_seq_len << endl;
		// if cur_seq_len has decreased then we're improving
		// if not, then we're done finding matches
		if( collinear_genomes && cur_ancestral_seq_len >= prev_ancestral_seq_len )
			break;

		// stop unless we've increased the anchoring score by at least 0.5%
		// the 0.5% is important for large alignments where many slow iterations might otherwise occur
		// that only increase the anchoring score by a tiny amount
		if( !collinear_genomes && cur_anchoring_score <= prev_anchoring_score + (genome::absolut(prev_anchoring_score)/200.0) )
			break;
		prev_anchoring_score = cur_anchoring_score;
		prev_ancestral_seq_len = cur_ancestral_seq_len;

		// accept the new alignment tree...
		cout << "Backing up alignment tree...\n";
		cout.flush();
		aln_tree_backup = alignment_tree;

		cout << "propagating ancestral breakpoints\n";
		cout.flush();
		recursiveApplyAncestralBreakpoints(ancestor);


		if( debug_me ) 	 
		{
			for( size_t aI = 0; aI < alignment_tree[ancestor].ordering.size(); aI++ ) 	 
			{
				GappedAlignment gal; 	 
				extractAlignment(ancestor, aI, gal); 	 

				bool check = false;
				for( size_t ii = 0; ii < gal.SeqCount(); ++ii )
				{
					if( gal.LeftEnd(ii) == 0 )
						continue;
					for( size_t jj = 0; jj < gal.SeqCount(); ++jj )
					{
						if( gal.LeftEnd(jj) == 0 )
							continue;
						check = check || computeID( gal, ii, jj ) < .5;
					}
				}
				if( check )
					cerr << "check iv " << aI << " dbg_count " << dbg_count << endl;
				else
					continue;

				const vector< string >& aln_mat = GetAlignment(gal, this->original_ml.seq_table); 	 
				gnSequence seq; 	 
				for( size_t seqI = 0; seqI < gal.SeqCount(); ++seqI ) 	 
					if( gal.LeftEnd(seqI) != NO_MATCH ) 	 
						seq += aln_mat[seqI]; 	 

				stringstream dbg_fname; 	 
				dbg_fname << "prof_dbg_iv_" << aI << ".dbg." << dbg_count++ << ".fas"; 	 
				ofstream debug_file( dbg_fname.str().c_str() ); 	 
				gnFASSource::Write( seq, debug_file, false ); 	 
				debug_file.close(); 	 
			} 	 
		}

		if(debug_aligner)
			validateSuperIntervals( node1, node2, ancestor );

		if(recursive)
		{
		// search for additional alignment anchors
		cout << "recursive anchor search\n";
		cout.flush();
		Matrix<MatchList> matches;
		Matrix< std::vector< search_cache_t > > new_cache_db(node1_seqs.size(), node2_seqs.size());
		// initialize storage for intervening regions
		boost::multi_array< std::vector< std::vector< int64 > >, 2 > iv_regions( boost::extents[node1_seqs.size()][node2_seqs.size()] );
		for( seqI = 0; seqI < node1_seqs.size(); seqI++ )
			for( seqJ = 0; seqJ < node2_seqs.size(); seqJ++ )
				iv_regions[seqI][seqJ].resize(2);
		vector< gnSequence* > bseqs( node1_seqs.size() + node2_seqs.size() );
		for( size_t aI = 0; aI < alignment_tree[ancestor].ordering.size(); aI++ )
		{
			CompactGappedAlignment<> cga;
			extractAlignment(ancestor, aI, cga);
			recurseOnPairs(node1_seqs, node2_seqs, cga, matches, search_cache_db, new_cache_db, iv_regions);

			// add any new matches to the pairwise_matches matrix
			for( seqI = 0; seqI < node1_seqs.size(); seqI++ )
				for( seqJ = 0; seqJ < node2_seqs.size(); seqJ++ )
					pairwise_matches(seqI, seqJ).insert( pairwise_matches(seqI, seqJ).end(), matches(seqI, seqJ).begin(), matches(seqI, seqJ).end() );

		}

		// add seqs
		for( seqI = 0; seqI < node1_seqs.size(); seqI++ )
			bseqs[seqI] = alignment_tree[ node1_seqs[seqI] ].sequence;
		for( seqJ = 0; seqJ < node2_seqs.size(); seqJ++ )
			bseqs[seqI+seqJ] =  alignment_tree[ node2_seqs[seqJ] ].sequence;

		MaskedMemHash nway_mh;
		// now search intervening regions
		for( seqI = 0; seqI < node1_seqs.size(); seqI++ )
			for( seqJ = 0; seqJ < node2_seqs.size(); seqJ++ )
			{
				std::sort( iv_regions[seqI][seqJ][0].begin(), iv_regions[seqI][seqJ][0].end() );
				std::sort( iv_regions[seqI][seqJ][1].begin(), iv_regions[seqI][seqJ][1].end() );
				MatchList new_matches;
				new_matches.seq_table.resize(2);
				new_matches.seq_table[0] = bseqs[seqI];
				new_matches.seq_table[1] = bseqs[node1_seqs.size() + seqJ];
				SearchLCBGaps( new_matches, iv_regions[seqI][seqJ], nway_mh );
				cout << seqI << "," << seqJ << " have " << new_matches.size() << " new matches outside LCBs\n";
				pairwise_matches(seqI, seqJ).insert( pairwise_matches(seqI, seqJ).end(), new_matches.begin(), new_matches.end() );
			}

		if(using_cache_db)
		{

		for( seqI = 0; seqI < node1_seqs.size(); seqI++ )
		{
			for( seqJ = 0; seqJ < node2_seqs.size(); seqJ++ )
			{
				for( size_t mI = 0; mI < search_cache_db(seqI,seqJ).size(); mI++ )
				{
					if( search_cache_db(seqI,seqJ)[mI].first != NULL )
						search_cache_db(seqI,seqJ)[mI].first->Free();
					if( search_cache_db(seqI,seqJ)[mI].second != NULL )
						search_cache_db(seqI,seqJ)[mI].second->Free();
				}
				search_cache_db(seqI,seqJ).clear();
				if(new_cache_db(seqI, seqJ).size() > 0)
				{
					// try sorting using C's qsort -- maybe there's something wrong with std::sort?
					search_cache_t* sc_array = new search_cache_t[new_cache_db(seqI,seqJ).size()];
					for( size_t i = 0; i < new_cache_db(seqI,seqJ).size(); i++ )
						sc_array[i] = new_cache_db(seqI,seqJ)[i];
					qsort(sc_array, new_cache_db(seqI,seqJ).size(), sizeof(AbstractMatch*), cachecomp);

					search_cache_db(seqI, seqJ).resize(new_cache_db(seqI,seqJ).size());
					for( size_t i = 0; i < new_cache_db(seqI,seqJ).size(); i++ )
						search_cache_db(seqI, seqJ)[i] = sc_array[i];
					delete[] sc_array;

					new_cache_db(seqI, seqJ).clear();
				}
				if( pairwise_matches(seqI,seqJ).size() > 0 )
					cout << seqI << "," << seqJ << " has an additional " << pairwise_matches(seqI,seqJ).size() << " matches\n";
			}
		}
		
		}
		}	// if recursive

		// restore backed up tree since we only want the final set of ancestral
		// breakpoints applied to the descendants
		cout << "Restoring backed up alignment tree...\n";
		cout.flush();
		swap( alignment_tree, aln_tree_backup );

	}	// end while(true)

	if( using_cache_db )
	{
	// delete the search cache
	for( seqI = 0; seqI < node1_seqs.size(); seqI++ )
		for( seqJ = 0; seqJ < node2_seqs.size(); seqJ++ )
			for( size_t mI = 0; mI < search_cache_db(seqI,seqJ).size(); mI++ )
			{
				if( search_cache_db(seqI,seqJ)[mI].first != NULL )
					search_cache_db(seqI,seqJ)[mI].first->Free();
				if( search_cache_db(seqI,seqJ)[mI].second != NULL )
					search_cache_db(seqI,seqJ)[mI].second->Free();
			}
	}

	printMemUsage();

	// aln_tree_backup has the highest scoring alignment_tree
	swap( alignment_tree, aln_tree_backup );
	cout << "propagating ancestral breakpoints\n";
	recursiveApplyAncestralBreakpoints(ancestor);

	printMemUsage();

	// step 8) construct a muscle alignment in each intervening region
	if( gapped_alignment )
	{
		cout << "performing a gapped alignment\n";
		doGappedAlignment(ancestor, true);
	}else
		cout << "skipping gapped alignment\n";
	if( refine )
	{
		size_t unrefined = countUnrefined( alignment_tree, ancestor );
		if( unrefined > 5 && ancestor != alignment_tree.root )
		{
			cout << "performing iterative refinement\n";
			doGappedAlignment(ancestor, false);
			markAsRefined( alignment_tree, ancestor );
		}
	}
	printMemUsage();


	if( debug_me ) 	 
	{
		for( size_t aI = 0; aI < alignment_tree[ancestor].ordering.size(); aI++ ) 	 
		{ 	 

			static int dbg_count = 0; 	 
			GappedAlignment gal; 	 
			extractAlignment(ancestor, aI, gal); 	 

			bool check = false;
			for( size_t ii = 0; ii < gal.SeqCount(); ++ii )
			{
				if( gal.LeftEnd(ii) == 0 )
					continue;
				for( size_t jj = 0; jj < gal.SeqCount(); ++jj )
				{
					if( gal.LeftEnd(jj) == 0 )
						continue;
					check = check || computeID( gal, ii, jj ) < .5;
				}
			}
			if( check )
				cerr << "check iv " << aI << " dbg_count " << dbg_count << endl;
			else
				continue;

			const vector< string >& aln_mat = GetAlignment(gal, this->original_ml.seq_table); 	 
			gnSequence seq; 	 
			for( size_t seqI = 0; seqI < gal.SeqCount(); ++seqI ) 	 
				if( gal.LeftEnd(seqI) != NO_MATCH ) 	 
					seq += aln_mat[seqI]; 	 

			stringstream dbg_fname; 	 
			dbg_fname << "prof_dbg_iv_" << aI << ".dbg." << dbg_count++ << ".fas"; 	 
			ofstream debug_file( dbg_fname.str().c_str() ); 	 
			gnFASSource::Write( seq, debug_file, false ); 	 
			debug_file.close(); 	 
		} 	 
	}


}


void addGuy( uint seqI, AbstractMatch::orientation orient, 
			std::vector< AbstractMatch* >& new_ivs, 
			vector<Interval*>& new_list )
{
	Interval tmp_iv;
	// set the orientation for any unaligned intervals
	if( orient == AbstractMatch::reverse )
	{
		for( size_t nI = 0; nI < new_ivs.size(); nI++ )
			if( new_ivs[nI]->LeftEnd(seqI) != NO_MATCH && new_ivs[nI]->Orientation(seqI) != orient)
				new_ivs[nI]->Invert();
	}
	// add this guy
	Interval* added_iv = tmp_iv.Copy();
	added_iv->SetMatches( new_ivs );
	new_list.push_back(added_iv);
}

void mergeUnalignedIntervals( uint seqI, vector< Interval* >& iv_list, vector< Interval* >& new_list )
{
	SSC<Interval> ivlcJ(seqI);
	sort( iv_list.begin(), iv_list.end(), ivlcJ );

	Interval tmp_iv;
	AbstractMatch::orientation orient = AbstractMatch::undefined;
	vector< AbstractMatch* > new_ivs;
	vector< Interval* > to_delete;
	for( size_t ordI = 0; ordI < iv_list.size(); ordI++ )
	{
		if( iv_list[ordI]->LeftEnd(seqI) == NO_MATCH )
		{
			new_list.push_back(iv_list[ordI]);
			iv_list[ordI] = NULL;
			continue;
		}

		if( orient == AbstractMatch::undefined && iv_list[ordI]->Multiplicity() == 2 )
		{
			orient = iv_list[ordI]->Orientation(seqI);
			vector< AbstractMatch* > matches;
			iv_list[ordI]->StealMatches( matches );
			if( orient == AbstractMatch::forward )
				new_ivs.insert( new_ivs.end(), matches.begin(), matches.end() );
			else
				new_ivs.insert( new_ivs.begin(), matches.begin(), matches.end() );

			// if it's the last one then add
			if( ordI + 1 == iv_list.size() )
				addGuy( seqI, orient, new_ivs, new_list );
			continue;
		}
		if( orient != AbstractMatch::undefined && iv_list[ordI]->Multiplicity() == 2 )
		{
			// add this guy...
			// set the orientation for any unaligned intervals
			addGuy( seqI, orient, new_ivs, new_list );

			// prepare a new one
			vector< AbstractMatch* > matches;
			orient = iv_list[ordI]->Orientation(seqI);
			iv_list[ordI]->StealMatches( matches );
			new_ivs.insert( new_ivs.end(), matches.begin(), matches.end() );
			// if it's the last one then add
			if( ordI + 1 == iv_list.size() )
				addGuy( seqI, orient, new_ivs, new_list );
			continue;
		}
		if( new_ivs.size() == 0 )
		{
			vector< AbstractMatch* > matches;
			iv_list[ordI]->StealMatches( matches );
			new_ivs.insert( new_ivs.end(), matches.begin(), matches.end() );
			continue;
		}
		// split this one in half (if its not the last one and there's something to split)...
		Interval* left_iv = iv_list[ordI]->Copy();
		to_delete.push_back( left_iv );	// make sure this gets deleted later
		bool cropped = (ordI + 1 < iv_list.size() && iv_list[ordI]->Length(seqI) > 1);
		if( cropped )
		{
			gnSeqI lendo = left_iv->AlignmentLength() / 2;
			left_iv->CropEnd( left_iv->AlignmentLength() - lendo );
			iv_list[ordI]->CropStart( lendo );
		}
		vector< AbstractMatch* > matches;
		left_iv->StealMatches( matches );
		if( orient == AbstractMatch::forward )
			new_ivs.insert( new_ivs.end(), matches.begin(), matches.end() );
		else
			new_ivs.insert( new_ivs.begin(), matches.begin(), matches.end() );

		addGuy( seqI, orient, new_ivs, new_list );
		// prepare for the next
		orient = AbstractMatch::undefined;
		if(cropped)
			ordI--;	// if we split a match, make sure we get the rest of this match on the next run through the loop
	}

	if( new_ivs.size() > 0 )
	{
		// uh-oh. there must not have been anything aligned
		addGuy( seqI, AbstractMatch::forward, new_ivs, new_list );
	}

	// free up any left_ivs that were allocated
	for( size_t delI = 0; delI < to_delete.size(); delI++ )
		to_delete[delI]->Free();

	// free up ivs left in iv_list
	for( size_t ivI = 0; ivI < iv_list.size(); ivI++ )
		if( iv_list[ivI] != NULL )
			iv_list[ivI]->Free();
	iv_list.clear();
}


/**
 * 
 */
void ProgressiveAligner::createAncestralOrdering( vector<Interval*>& interval_list, vector< SuperInterval >& ancestral_sequence )
{
	// construct an ancestral SuperSequence
	int64 left_end = 1;
	ancestral_sequence.resize( interval_list.size() );
	for( uint ivI = 0; ivI < interval_list.size(); ++ivI ){
		if(debug_aligner)
			interval_list[ivI]->ValidateMatches();
		vector<AbstractMatch*> matches;
		interval_list[ivI]->StealMatches(matches);
		ancestral_sequence[ivI].reference_iv.SetMatches(matches);
		ancestral_sequence[ivI].SetLeftEnd(left_end);
		ancestral_sequence[ivI].SetLength(ancestral_sequence[ivI].reference_iv.AlignmentLength());
		if(debug_aligner)
			ancestral_sequence[ivI].ValidateSelf();
		left_end += ancestral_sequence[ivI].Length();
	}
}

void markAligned( PhyloTree< AlignmentTreeNode >& alignment_tree, node_id_t subject_node, node_id_t neighbor )
{
	for( uint parentI = 0; parentI < alignment_tree[subject_node].parents.size(); parentI++ )
		if( alignment_tree[subject_node].parents[parentI] == neighbor )
			alignment_tree[subject_node].parents_aligned[parentI] = true;
	for( uint childI = 0; childI < alignment_tree[subject_node].children.size(); childI++ )
		if( alignment_tree[subject_node].children[childI] == neighbor )
			alignment_tree[subject_node].children_aligned[childI] = true;
}


bool
ProgressiveAligner::validateSuperIntervals(node_id_t node1, node_id_t node2, node_id_t ancestor)
{
		// validate the ancestor
	bool borked = false;
	vector< SuperInterval >& siv_list = alignment_tree[ancestor].ordering;
	gnSeqI n1_len = 0;
	gnSeqI n2_len = 0;
	gnSeqI my_len = 0;
	gnSeqI my_iv_len = 0;
	for( size_t sivI = 0; sivI < siv_list.size(); sivI++ )
	{
		if( siv_list[sivI].reference_iv.Start(0) != 0 )
			n1_len += siv_list[sivI].reference_iv.Length(0);
		if( siv_list[sivI].reference_iv.Start(1) != 0 )
			n2_len += siv_list[sivI].reference_iv.Length(1);
		my_len += siv_list[sivI].Length();
		my_iv_len += siv_list[sivI].reference_iv.AlignmentLength();
		siv_list[sivI].ValidateSelf();
	}
	gnSeqI real_n1len = 0;
	gnSeqI real_n2len = 0;

	vector< SuperInterval >& siv1_list = alignment_tree[node1].ordering;
	for( size_t sivI = 0; sivI < siv1_list.size(); sivI++ )
	{
		if( siv1_list[sivI].Length() == 0 )
			borked = true;
		real_n1len += siv1_list[sivI].Length();
		siv1_list[sivI].ValidateSelf();
	}

	vector< SuperInterval >& siv2_list = alignment_tree[node2].ordering;
	for( size_t sivI = 0; sivI < siv2_list.size(); sivI++ )
	{
		if( siv2_list[sivI].Length() == 0 )
			borked = true;
		real_n2len += siv2_list[sivI].Length();
		siv2_list[sivI].ValidateSelf();
	}

	if( real_n1len != n1_len || real_n2len != n2_len )
			borked = true;

	// check that each picks up where the last left off
	for( size_t sivI = 1; sivI < siv1_list.size(); sivI++ )
		if( siv1_list[sivI].LeftEnd() != siv1_list[sivI-1].LeftEnd() + siv1_list[sivI-1].Length() )
		{
			borked = true;
		}
	for( size_t sivI = 1; sivI < siv2_list.size(); sivI++ )
		if( siv2_list[sivI].LeftEnd() != siv2_list[sivI-1].LeftEnd() + siv2_list[sivI-1].Length() )
		{
			borked = true;
		}

	if( my_len != my_iv_len )
		borked = true;

	if( my_len < real_n1len || my_len < real_n2len )
		borked = true;

	if( borked )
	{
		breakHere();
		cerr << "child1 has " << siv1_list.size() << " ivs totalling " << real_n1len << " nt\n";
		cerr << "child2 has " << siv2_list.size() << " ivs totalling " << real_n2len << " nt\n";
		cerr << "parent has " << siv_list.size() << " ivs, n1_len: " << n1_len << " n2_len: " << n2_len << endl;
	}
	return borked;

}

bool ProgressiveAligner::validatePairwiseIntervals(node_id_t node1, node_id_t node2, std::vector<Interval*>& pair_iv)
{
		// validate the ancestor
	bool borked = false;
	gnSeqI n1_len = 0;
	gnSeqI n2_len = 0;
	for( size_t sivI = 0; sivI < pair_iv.size(); sivI++ )
	{
		if( pair_iv[sivI]->Start(0) != 0 )
			n1_len += pair_iv[sivI]->Length(0);
		if( pair_iv[sivI]->Start(1) != 0 )
			n2_len += pair_iv[sivI]->Length(1);

		vector< bitset_t > aln_mat;
		pair_iv[sivI]->GetAlignment(aln_mat);
		if( aln_mat[0].size() != pair_iv[sivI]->AlignmentLength() )
		{
			cerr << "broked\n";
		}
		pair_iv[sivI]->ValidateMatches();
	}
	gnSeqI real_n1len = 0;
	gnSeqI real_n2len = 0;

	vector< SuperInterval >& siv1_list = alignment_tree[node1].ordering;
	for( size_t sivI = 0; sivI < siv1_list.size(); sivI++ )
	{
		if( siv1_list[sivI].Length() == 0 )
			borked = true;
		real_n1len += siv1_list[sivI].Length();
	}

	vector< SuperInterval >& siv2_list = alignment_tree[node2].ordering;
	for( size_t sivI = 0; sivI < siv2_list.size(); sivI++ )
	{
		if( siv2_list[sivI].Length() == 0 )
			borked = true;
		real_n2len += siv2_list[sivI].Length();
	}

	if( real_n1len != n1_len || real_n2len != n2_len )
			borked = true;

	// check for overlapping intervals
	vector< Interval* > tmp_iv_list = pair_iv;
	for( uint seqI = 0; seqI < 2; seqI++ )
	{
		SSC<Interval> ssc(seqI);
		sort( tmp_iv_list.begin(), tmp_iv_list.end(), ssc );
		for( size_t ivI = 1; ivI < tmp_iv_list.size(); ivI++ )
		{
			if( tmp_iv_list[ivI-1]->LeftEnd(seqI) == NO_MATCH || tmp_iv_list[ivI]->LeftEnd(seqI) == NO_MATCH )
				continue;
			if( tmp_iv_list[ivI-1]->RightEnd(seqI) >= tmp_iv_list[ivI]->LeftEnd(seqI) )
			{
				cerr << "overlap:\n";
				cerr << "tmp_iv_list[ivI-1].RightEnd(seqI): " << tmp_iv_list[ivI-1]->RightEnd(seqI) << endl;
				cerr << "tmp_iv_list[ivI].LeftEnd(seqI): " << tmp_iv_list[ivI]->LeftEnd(seqI) << endl;
				breakHere();
			}
		}
	}

	if( borked )
	{
		cerr << "child1 has " << siv1_list.size() << " ivs totalling " << real_n1len << " nt\n";
		cerr << "child2 has " << siv2_list.size() << " ivs totalling " << real_n2len << " nt\n";
		cerr << "parent has " << pair_iv.size() << " ivs, n1_len: " << n1_len << " n2_len: " << n2_len << endl;
		if( n2_len < real_n2len )
		{
			SSC<Interval> sortie(1);
			sort( pair_iv.begin(), pair_iv.end(), sortie );
			size_t prev_iv = 9999999;
			for( size_t ivI = 0; ivI < pair_iv.size(); ++ivI)
			{
				if( pair_iv[ivI]->LeftEnd(1) == NO_MATCH )
					continue;

				if( prev_iv != 9999999 )
					cerr << "diff: " << pair_iv[ivI]->LeftEnd(1) - pair_iv[prev_iv]->RightEnd(1) << endl;
				cerr << "Interval " << ivI << " LeftEnd(1): " << pair_iv[ivI]->LeftEnd(1) << " RightEnd(1): " << pair_iv[ivI]->RightEnd(1) << std::endl;
				prev_iv = ivI;
			}
		}else if( n2_len > real_n2len )
		{
			SSC<Interval> sortie(1);
			sort( pair_iv.begin(), pair_iv.end(), sortie );
			for( size_t ivI = 0; ivI < pair_iv.size(); ++ivI)
			{
				if( pair_iv[ivI]->LeftEnd(1) < real_n2len )
					continue;
				cerr << "Interval " << ivI << " LeftEnd(1): " << pair_iv[ivI]->LeftEnd(1) << " RightEnd(1): " << pair_iv[ivI]->RightEnd(1) << std::endl;
			}
		}
		breakHere();
	}
	return borked;
}

void ProgressiveAligner::alignNodes( node_id_t node1, node_id_t node2, node_id_t ancestor )
{
	cout << "Aligning node " << node1 << " to " << node2 << " via " << ancestor << "!\n";
	// if node1 and node2 are not already children of ancestor then make it so...
	if( alignment_tree[node1].parents[0] != ancestor || 
		alignment_tree[node2].parents[0] != ancestor )
	{
		breakHere();
		cerr << "rotten\n";
	}
	
	alignProfileToProfile(node1, node2, ancestor);

	// mark edges as aligned
	markAligned( alignment_tree, node1, node2 );
	markAligned( alignment_tree, node2, node1 );
	markAligned( alignment_tree, node1, ancestor );
	markAligned( alignment_tree, node2, ancestor );
	markAligned( alignment_tree, ancestor, node1 );
	markAligned( alignment_tree, ancestor, node2 );
}

/**
 * finds the midpoint of a phylogenetic tree, returns the ids of the surrounding nodes in n1 and n2
 */
void findMidpoint( PhyloTree< AlignmentTreeNode >& alignment_tree, node_id_t& n1, node_id_t& n2 )
{
	// use boost's all pairs shortest path to find the longest path on the tree 
	// Then actually traverse the path to determine which edge
	// is halfway.
	double scaling_factor = 100000;
	using namespace boost;
	typedef adjacency_list<vecS, vecS, undirectedS, no_property,
	property< edge_weight_t, int, property< edge_color_t, default_color_type > > > Graph;
	const int V = alignment_tree.size();
	const std::size_t E = alignment_tree.size()-1;
	typedef std::pair < int, int >Edge;
	Edge* edge_array = new Edge[ alignment_tree.size() - 1 ];
	int* weights = new int[ alignment_tree.size() - 1 ];
	bitset_t child_found( alignment_tree.size(), false );
	size_t eI = 0;
	for( size_t vI = 0; vI < V; ++vI )
	{
		if( alignment_tree[vI].parents.size() != 0 )
		{
			edge_array[eI] = Edge( vI, alignment_tree[vI].parents[0] );
			// for some reason boost insists on using an int for weights.  need to figure that out
			weights[eI] = (int)(scaling_factor * genome::absolut(alignment_tree[vI].distance)) + 1;
			eI++;
		}
	}

#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
	// VC++ can't handle the iterator constructor
	Graph g(V);
	for (std::size_t j = 0; j < E; ++j)
	add_edge(edge_array[j].first, edge_array[j].second, g);
#else
	Graph g(edge_array, edge_array + E, V);
#endif

	property_map < Graph, edge_weight_t >::type w = get(edge_weight, g);
	int *wp = weights;

	graph_traits < Graph >::edge_iterator e, e_end;
	for (boost::tie(e, e_end) = edges(g); e != e_end; ++e)
		w[*e] = *wp++;

	boost::multi_array<int,2> D( boost::extents[V][V] );
	bool success = johnson_all_pairs_shortest_paths(g, D);
	if( !success )
	{
		cerr << "failed, is this really a tree?\n";
		return;
	}

	// find the most distant pair of nodes
	int max_dist = (std::numeric_limits<int>::min)();
	for (int i = 0; i < V; ++i) {
		for (int j = 0; j < V; ++j) {
			if( D[i][j] > max_dist )
			{
				max_dist = D[i][j];
				n1 = i;
				n2 = j;
			}
		}
	}

	typedef graph_traits<Graph>::vertex_descriptor vertex_t;
	std::vector < vertex_t > pred(num_vertices(g));
	std::vector < int > dist(num_vertices(g));
	pred[n1] = n1;

	undirected_dfs(g,
		root_vertex( vertex( n1, g ) ).
		visitor( make_dfs_visitor( make_pair(
			record_predecessors(&pred[0], on_tree_edge()),
			record_distances(&dist[0], on_tree_edge())
		))).
		edge_color_map(get(edge_color, g))
		);

	int cur_node = n2;
	int prev_node = n2;
	max_dist /= 2;
	while( cur_node != n1 && max_dist > 0 )
	{
		if( alignment_tree[cur_node].parents.size() > 0 && 
			alignment_tree[cur_node].parents[0] == pred[cur_node] )
		{
			max_dist -= (int)(scaling_factor * alignment_tree[cur_node].distance) + 1;
			prev_node = cur_node;
			cur_node = pred[cur_node];
		}else
		{
			prev_node = cur_node;
			cur_node = pred[cur_node];
			max_dist -= (int)(scaling_factor * alignment_tree[cur_node].distance) + 1;
		}
	}
	n1 = cur_node;
	n2 = prev_node;

	delete[] edge_array;
	delete[] weights;
}

void extendRootBranches( PhyloTree< AlignmentTreeNode >& alignment_tree )
{
	// find the max branch length and set the root branch lengths to twice that
	// swap children while we're at it
	node_id_t ancestor = alignment_tree.root;
	double max_blen = -(std::numeric_limits<double>::max)();
	for( size_t nI = 0; nI < alignment_tree.size(); ++nI )
	{
		if( alignment_tree[nI].distance > max_blen )
			max_blen = alignment_tree[nI].distance;
		if( alignment_tree[nI].children.size() > 0 &&
			alignment_tree[nI].children[0] > alignment_tree[nI].children[1] )
		{
			std::swap( alignment_tree[nI].children[0], alignment_tree[nI].children[1] );
		}
	}
	for( size_t cI = 0; cI < alignment_tree[ancestor].children.size(); ++cI )
		alignment_tree[alignment_tree[ancestor].children[cI]].distance = 2.0 * max_blen;
}

void chooseNextAlignmentPair( PhyloTree< AlignmentTreeNode >& alignment_tree, node_id_t& node1, node_id_t& node2, node_id_t& ancestor )
{

	// find the nearest alignable neighbor
	node1 = 0;
	node2 = 0;
	ancestor = 0;
	double nearest_distance = (numeric_limits<double>::max)();
	for( node_id_t nodeI = 0; nodeI < alignment_tree.size(); nodeI++ )
	{
		AlignmentTreeNode& cur_node = alignment_tree[ nodeI ];

		// skip this node if it's already been completely aligned
		// or is an extant sequence
		boolean completely_aligned = true;
		for( uint alignedI = 0; alignedI < cur_node.children_aligned.size(); alignedI++ )
			completely_aligned = completely_aligned && cur_node.children_aligned[alignedI];
		for( uint alignedI = 0; alignedI < cur_node.parents_aligned.size(); alignedI++ )
			completely_aligned = completely_aligned && cur_node.parents_aligned[alignedI];
		if( cur_node.sequence != NULL || completely_aligned )
			continue;
		

		vector< node_id_t > neighbor_id;
		vector< boolean > alignable;
		vector< double > distance;
		
		for( uint parentI = 0; parentI < cur_node.parents.size(); parentI++ )
		{
			neighbor_id.push_back( cur_node.parents[parentI] );
			vector< node_id_t >::iterator cur_neighbor = neighbor_id.end() - 1;
			if( *cur_neighbor == alignment_tree.root )
			{
				// need special handling for the root since the alignment
				// tree is supposed to be unrooted
				// add all of root's children except this one
			}
			distance.push_back( cur_node.distance );
			alignable.push_back( !cur_node.parents_aligned[parentI] && (alignment_tree[*cur_neighbor].ordering.size() != 0 || alignment_tree[*cur_neighbor].sequence != NULL) );
		}

		for( uint childI = 0; childI < cur_node.children.size(); childI++ )
		{
			neighbor_id.push_back( cur_node.children[childI] );
			vector< node_id_t >::iterator cur_neighbor = neighbor_id.end() - 1;
			distance.push_back( alignment_tree[*cur_neighbor].distance );
			alignable.push_back( !cur_node.children_aligned[childI] && (alignment_tree[*cur_neighbor].ordering.size() != 0 || alignment_tree[*cur_neighbor].sequence != NULL) );
		}

		if( cur_node.ordering.size() != 0 )
		{
			// this one already has at least two sequences aligned, if another
			// is alignable then check its distance
			for( int i = 0; i < neighbor_id.size(); i++ ){
				if( !alignable[i] )
					continue;
				if( distance[i] < nearest_distance )
				{
					nearest_distance = distance[i];
					node1 = nodeI;
					node2 = neighbor_id[i];
					ancestor = nodeI;
				}
			}
		}else{
			// find the nearest alignable pair
			for( int i = 0; i < neighbor_id.size(); i++ )
			{
				if( !alignable[i] )
					continue;
				for( int j = i+1; j < neighbor_id.size(); j++ )
				{
					if( !alignable[j] )
						continue;
					if( distance[i] + distance[j] < nearest_distance )
					{
						nearest_distance = distance[i] + distance[j];
						node1 = neighbor_id[i];
						node2 = neighbor_id[j];
						ancestor = nodeI;
					}
				}
			}
		}
	}
}

/** use a list of precomputed matches instead of computing them */
void ProgressiveAligner::setPairwiseMatches( MatchList& pair_ml )
{
	original_ml = pair_ml;
	pair_ml.clear();	// ProgressiveAligner owns the matches now...
}


node_id_t createAlignmentTreeRoot( PhyloTree< AlignmentTreeNode >& alignment_tree, node_id_t node1, node_id_t node2 )
{
		// create a new node and link it inline between node1 and node2
		AlignmentTreeNode atn;
		alignment_tree.push_back( atn );
		AlignmentTreeNode& old_root = alignment_tree[alignment_tree.root];
		AlignmentTreeNode& new_root = alignment_tree.back();

		if( find( alignment_tree[node1].children.begin(), alignment_tree[node1].children.end(), node2 ) !=
			alignment_tree[node1].children.end() )
		{
			new_root.children.push_back(node2);
			new_root.parents.push_back(node1);
			alignment_tree[node2].parents.push_back(alignment_tree.size()-1);
			alignment_tree[node1].children.push_back(alignment_tree.size()-1);
		}else{
			new_root.parents.push_back(node2);
			new_root.children.push_back(node1);
			alignment_tree[node2].children.push_back(alignment_tree.size()-1);
			alignment_tree[node1].parents.push_back(alignment_tree.size()-1);
		}

		// completely unlink node1 and node2 from each other
		findAndErase( alignment_tree[node1].children, node2 );
		findAndErase( alignment_tree[node2].children, node1 );
		findAndErase( alignment_tree[node1].parents, node2 );
		findAndErase( alignment_tree[node2].parents, node1 );


		// re-root the tree on the new node
		rerootTree( alignment_tree, alignment_tree.size()-1 );

		new_root.children_aligned = vector< boolean >( new_root.children.size(), false );
		old_root.children_aligned = vector< boolean >( old_root.children.size(), false );
		old_root.parents_aligned = vector< boolean >( old_root.parents.size(), false );
		new_root.sequence = NULL;

	return alignment_tree.root;
}

void ProgressiveAligner::extractAlignment( node_id_t ancestor, size_t super_iv, GappedAlignment& gal )
{
	CompactGappedAlignment<> cga;
	extractAlignment( ancestor, super_iv, cga );
	vector< string > aln;
	GetAlignment( cga, this->original_ml.seq_table, aln );
	gal = GappedAlignment(cga.SeqCount(), 0);
	for( size_t seqI = 0; seqI < cga.SeqCount(); ++seqI )
	{
		gal.SetStart(seqI, cga.Start(seqI));
		if( cga.Orientation(seqI) != AbstractMatch::undefined )
			gal.SetLength(cga.Length(seqI), seqI);
	}
	gal.SetAlignment(aln);

}

void ProgressiveAligner::extractAlignment( node_id_t ancestor, size_t super_iv, CompactGappedAlignment<>& cga )
{
	// determine the leaf node intervals below this super_iv
	vector< pair< node_id_t, size_t > > node_siv_list;
	stack< pair<node_id_t,size_t> > node_stack;
	node_stack.push(make_pair(ancestor,super_iv));
	while( node_stack.size() > 0 )
	{
		pair<node_id_t,size_t> cur = node_stack.top();
		node_id_t cur_node = cur.first;
		node_stack.pop();
		if( alignment_tree[cur_node].children.size() == 0 )
			node_siv_list.push_back( cur );
		for( size_t childI = 0; childI < alignment_tree[cur_node].children.size(); childI++ )
		{
			if( alignment_tree[cur_node].ordering[cur.second].reference_iv.LeftEnd(childI) == NO_MATCH )
				continue;
			size_t child_siv = childI == 0 ? alignment_tree[cur_node].ordering[cur.second].c1_siv : 
				alignment_tree[cur_node].ordering[cur.second].c2_siv;
			node_stack.push(make_pair(alignment_tree[cur_node].children[childI], child_siv) );
			node_id_t n = alignment_tree[cur_node].children[childI];
			if( alignment_tree[cur_node].ordering[cur.second].reference_iv.Length(childI) != alignment_tree[n].ordering[child_siv].Length() )
			{
				breakHere();
				cerr << "alignment_tree[cur_node].ordering[cur.second].reference_iv.Length(childI): " << alignment_tree[cur_node].ordering[cur.second].reference_iv.Length(childI) << endl;
				cerr << "rotten in the state of denmark...\n";
			}
		}
	}

	// armed with the list of pairs, extract each one...

	// for each interval at the root write out the alignment
	SuperInterval& a_iv = alignment_tree[ancestor].ordering[super_iv];
	cga = CompactGappedAlignment<>(seq_count, a_iv.Length());
	vector< bitset_t > aln_mats( seq_count );

	// use translateCoordinates to map out each sequence's original coordinates
	// to the alignment coordinates
	for( size_t pairI = 0; pairI < node_siv_list.size(); pairI++ )
	{
		node_id_t nodeI = node_siv_list[pairI].first;
		size_t seq_siv = node_siv_list[pairI].second;
		
		// translate seq_siv into ancestor alignment coordinates?  
		// we can abuse translateCoordinates and the Match data structure :
		//   - add a single "match" covering the entire sequence
		//   - translate it up to alignment root coordinates
		uint seqI = node_sequence_map[nodeI];
		Match mm(2);
		mm.SetStart(0, alignment_tree[nodeI].ordering[seq_siv].LeftEnd());
		mm.SetStart(1, alignment_tree[nodeI].ordering[seq_siv].LeftEnd());
		mm.SetLength( alignment_tree[nodeI].ordering[seq_siv].Length() );
		
		vector< AbstractMatch* > aml( 1, mm.Copy() );
		translateGappedCoordinates( aml, 0, nodeI, ancestor );

		if( aml.size() > 1 )
		{
			cerr << "huh?";
			genome::breakHere();
			SingleStartComparator<AbstractMatch> ssc( 0 );
			sort( aml.begin(), aml.end(), ssc );	// huh?
		}
		CompactGappedAlignment<>* trans_cga = dynamic_cast<CompactGappedAlignment<>*>(aml[0]);
		if( trans_cga == NULL )
		{
			CompactGappedAlignment<> tmp_cga;
			trans_cga = tmp_cga.Copy();
			*trans_cga = CompactGappedAlignment<>(*aml[0]);
		}

		if( trans_cga->LeftEnd(0) + trans_cga->Length(0) > a_iv.LeftEnd() + a_iv.Length() )
		{
			cerr << "trans_cga->Start(0): " << trans_cga->Start(0) << " trans_cga->Length(0): " << trans_cga->Length(0) << endl;
			cerr << "a_iv.LeftEnd(): " << a_iv.LeftEnd() << " a_iv.Length(): " << a_iv.Length() << endl;
			breakHere();
		}
		bool parity = trans_cga->Orientation(0) == trans_cga->Orientation(1);
		cga.SetLeftEnd(seqI, trans_cga->LeftEnd(1));
		AbstractMatch::orientation o = parity ? AbstractMatch::forward : AbstractMatch::reverse;
		cga.SetOrientation(seqI, o);
		const vector< bitset_t >& tmp = trans_cga->GetAlignment();
		aln_mats[seqI] = tmp[1];

		size_t offset = trans_cga->LeftEnd(0) - a_iv.LeftEnd();
		if( aln_mats[seqI].size() < a_iv.Length() )
		{
			// need to resize and shift appropriately
			aln_mats[seqI].resize( a_iv.Length() );
			aln_mats[seqI] <<= offset;	// this is backwards in boost::dynamic_bitset for some reason...
		}
		if( trans_cga->LeftEnd(0) < a_iv.LeftEnd() )
		{
			cerr << "trans_cga->LeftEnd(0): " << trans_cga->LeftEnd(0) << endl;
			cerr << "a_iv.LeftEnd(): " << a_iv.LeftEnd() << endl;
			breakHere();
		}

		// validate match lengths
		if( trans_cga->Length(1) != alignment_tree[nodeI].ordering[seq_siv].Length() )
		{
			cerr << "b0rked\n";
			breakHere();
		}
		// set the length and alignment appropriately
		cga.SetLength(trans_cga->Length(1), seqI);

		// free storage used by trans_cga
		trans_cga->Free();
	}
	for( uint seqI = 0; seqI < aln_mats.size(); seqI++ )
		if( aln_mats[seqI].size() == 0 )
			aln_mats[seqI].resize( a_iv.Length() );
	cga.SetAlignment(aln_mats);
}

unsigned getDefaultBreakpointMax( const std::vector< genome::gnSequence* >& seq_table )
{
	double avg_len = 0;
	for( size_t seqI = 0; seqI < seq_table.size(); ++seqI )
		avg_len += seq_table[seqI]->length();
	avg_len /= (double)(seq_table.size());
	// heavily rearranged, recently diverged genomes like yersinia have up to 15 rearrangements per megabase of sequence
	avg_len /= 1000000.0;	// convert to number of megabases
	avg_len *= 15.0;	// "lots" of rearrangement
	return (unsigned)avg_len;
}

// get a pairwise bp distance
void ProgressiveAligner::CreatePairwiseBPDistance( boost::multi_array<double, 2>& bp_distmat )
{
	uint seq_count = original_ml.seq_table.size();
	bp_distmat.resize(boost::extents[seq_count][seq_count]);
	for( size_t i = 0; i < seq_count; ++i )
		for( size_t j = 0; j < seq_count; ++j )
			bp_distmat[i][j] = 1;

#ifdef LCB_WEIGHT_LOSS_PLOT
	stringstream pair_bp_ofname;
	pair_bp_ofname << "pair_bp_log.txt";
	ofstream pair_bp_out( pair_bp_ofname.str().c_str() );
#endif

	vector< pair<uint, uint> > seq_pairs( (seq_count * (seq_count-1))/2 );
	int ii = 0;
	for( uint seqI = 0; seqI < seq_count; seqI++ )
		for( uint seqJ = seqI + 1; seqJ < seq_count; seqJ++ )
			seq_pairs[ii++] = make_pair(seqI,seqJ);

#pragma omp parallel for
	for(int i = 0; i < seq_pairs.size(); i++)
	{
		uint seqI = seq_pairs[i].first;
		uint seqJ = seq_pairs[i].second;
		vector<uint>::iterator n1 = find( node_sequence_map.begin(), node_sequence_map.end(), seqI );
		vector<uint>::iterator n2 = find( node_sequence_map.begin(), node_sequence_map.end(), seqJ );
		vector<node_id_t> n1_seqs( 1, n1-node_sequence_map.begin() );
		vector<node_id_t> n2_seqs( 1, n2-node_sequence_map.begin() );
		Matrix<MatchList> mml;
		getPairwiseMatches(n1_seqs, n2_seqs, mml);
		MatchList& ml = mml(0,0);

		// eliminate overlaps as they correspond to inconsistently or
		// multiply aligned regions
		EliminateOverlaps_v2( ml, true );
		ml.MultiplicityFilter(2);

		// do greedy b.p. elimination on the matches
		vector< MatchList > LCB_list;
		vector< LCB > adjacencies;
		vector< gnSeqI > breakpoints;
		IdentifyBreakpoints( ml, breakpoints );
		ComputeLCBs_v2( ml, breakpoints, LCB_list );
		vector< double > lcb_scores( LCB_list.size() );
		cout << "Pair " << seq_pairs[i].first << ", " << seq_pairs[i].second << " has " << LCB_list.size() << " initial LCBs\n";
		for( size_t lcbI = 0; lcbI < LCB_list.size(); ++lcbI )
			lcb_scores[lcbI] = GetPairwiseAnchorScore( LCB_list[lcbI], ml.seq_table, this->subst_scoring, sol_list[seqI], sol_list[seqJ] );

		computeLCBAdjacencies_v3( LCB_list, lcb_scores, adjacencies );

		// want to discard all low-weight LCBs
		// to arrive at a set of reliable LCBs
		double cons_id = 1 - this->conservation_distance[seqI][seqJ];
		double scaled_score = max( bp_dist_estimate_score * cons_id * cons_id * cons_id * cons_id, min_breakpoint_penalty);
		cout << "Using scaled bp penalty: " << scaled_score << endl;
		GreedyRemovalScorer wbs( adjacencies, scaled_score );
#ifdef LCB_WEIGHT_LOSS_PLOT
		cur_min_coverage = greedyBreakpointElimination_v4( adjacencies, lcb_scores, wbs, &pair_bp_out, seqI, seqJ );
		pair_bp_out.flush();
#else
		cur_min_coverage = greedyBreakpointElimination_v4( adjacencies, lcb_scores, wbs, NULL );
#endif
		MatchList deleted_matches;
		filterMatches_v2( adjacencies, LCB_list, lcb_scores, deleted_matches );
		cout << "Pair (" << seqI << "," << seqJ << ") has " << LCB_list.size() << " well-supported breakpoints\n";
		
		// now set the distance entry
		bp_distmat[seqI][seqJ] = LCB_list.size();
		bp_distmat[seqJ][seqI] = LCB_list.size();

		// free the matches
		for( size_t dI = 0; dI < ml.size(); dI++ )
			ml[dI]->Free();
	}
	// normalize to [0,1]
	double bp_max = 0;
	for( uint i = 0; i < bp_distmat.shape()[0]; ++i )
		for( uint j = 0; j < bp_distmat.shape()[1]; ++j )
		{
			if( bp_distmat[i][j] > bp_max )
				bp_max = bp_distmat[i][j];
		}

	double default_max = getDefaultBreakpointMax(original_ml.seq_table);
	bp_max = bp_max > default_max ? bp_max : default_max;

	for( uint i = 0; i < bp_distmat.shape()[0]; ++i )
		for( uint j = 0; j < bp_distmat.shape()[1]; ++j )
		{
			if( i != j )
				bp_distmat[i][j] /= bp_max;
			bp_distmat[i][j] *= bp_dist_scale;
		}
}

template< typename MatchListType >
void makeAlignmentTree( PhyloTree< AlignmentTreeNode >& alignment_tree, MatchListType& mlist, vector< uint >& node_sequence_map )
{
	// initialize all nodes to unaligned
	for( node_id_t nodeI = 0; nodeI < alignment_tree.size(); nodeI++ )
	{
		alignment_tree[nodeI].children_aligned = vector< boolean >( alignment_tree[nodeI].children.size(), false );
		alignment_tree[nodeI].parents_aligned = vector< boolean >( alignment_tree[nodeI].parents.size(), false );
		alignment_tree[nodeI].sequence = NULL;
		alignment_tree[nodeI].refined = false;
	}

	// set the sequence appropriately for extant sequences
	node_sequence_map = vector< uint >( alignment_tree.size(), -1 );
	for( uint seqI = 0; seqI < mlist.seq_table.size(); seqI++ )
	{
		stringstream seq_name;
		seq_name << "seq" << seqI + 1;
		node_id_t nodeI = 0;
		for( ; nodeI < alignment_tree.size(); nodeI++ )
		{
			if( seq_name.str() == alignment_tree[nodeI].name )
			{
				alignment_tree[nodeI].sequence = mlist.seq_table[seqI];
				Match mm(1);
				Match* m = mm.Copy();
				m->SetStart(0,1);
				m->SetLength(alignment_tree[nodeI].sequence->length(),0);
				vector<AbstractMatch*> tmp(1,m);
				Interval iv( tmp.begin(), tmp.end() );
				m->Free();
				SuperInterval si( iv );
				si.SetLeftEnd(1);
				si.SetLength(alignment_tree[nodeI].sequence->length());
				alignment_tree[nodeI].ordering.push_back( si );
				node_sequence_map[nodeI] = seqI;
				break;
			}
		}
		if( nodeI == alignment_tree.size() )
			throw "Phylogenetic tree names unrecognized.  Should follow seqN naming format\n";
	}
}

void DistanceMatrix( IntervalList& iv_list, NumericMatrix<double>& distmat )
{
	IdentityMatrix( iv_list, distmat );
	TransformDistanceIdentity(distmat);
}

/*
void makeSuperIntervals( IntervalList& iv_list, PhyloTree< TreeNode >& alignment_tree, vector< uint >& node_sequence_map )
{
	std::stack< node_id_t > node_stack;
	node_stack.push( alignment_tree.root );
	bitset_t visited( alignment_tree.size(), false );
	while( node_stack.size() > 0 )
	{
		node_id_t cur_node = node_stack.top();
		// visit post-order
		for( size_t cI = 0; cI < alignment_tree[cur_node].children.size(); ++cI )
		{
			if( !visited[alignment_tree[cur_node].children[cI]] )
				node_stack.push(alignment_tree[cur_node].children[cI]);
		}
		if( node_stack.top() != cur_node )
			continue;
		node_stack.pop();
		if( alignment_tree[cur_node].children.size() == 0 )
			continue;	// only process internal nodes

		// process this node
		// construct pairwise LCBs

		uint seqI = node_sequence_map[alignment_tree[cur_node].children[0]];
		uint seqJ = node_sequence_map[alignment_tree[cur_node].children[0]];
		vector< uint > projection( 2 );
		projection[0] = seqI;
		projection[1] = seqJ;

		vector< vector< MatchProjectionAdapter* > > LCB_list;
		vector< LCB > projected_adjs;
		projectIntervalList( iv_list, projection, LCB_list, projected_adjs );

		// create a superinterval for each adj 
//		alignment_tree[cur_node].ordering.resize(adjs.size());
//		for( size_t adjI = 0; adjI < adjs.size(); ++adjI )
//		{
//			SuperInterval& siv = alignment_tree[cur_node].ordering[adjI];
//			Match mleft(2);
//			mleft.SetStart(0,adjI);
//			mleft.SetStart(1,adjI);
//			mleft.SetLength(1);
//			siv.SetLeftEnd( adjI );
//			siv.SetLength(1);
//		}

	}
}
*/

void ProgressiveAligner::alignPP(IntervalList& prof1, IntervalList& prof2, IntervalList& interval_list )
{
	if( debug_aligner )
	{
		debug_interval = true;
		debug_cga = true;
	}

	seq_count = prof1.seq_table.size() + prof2.seq_table.size();

	if( this->breakpoint_penalty == -1 )
		this->breakpoint_penalty = getDefaultBreakpointPenalty( original_ml.seq_table );

	if( this->bp_dist_estimate_score == -1 )
		this->bp_dist_estimate_score = getDefaultBpDistEstimateMinScore( original_ml.seq_table );
	cout << "using default bp penalty: " << breakpoint_penalty << endl;
	cout << "using default bp estimate min score: " << bp_dist_estimate_score << endl;

	if( this->collinear_genomes )
		this->breakpoint_penalty = -1;

	if( collinear_genomes )
		cout << "\nAssuming collinear genomes...\n";
		
	EliminateOverlaps_v2( original_ml );
	// use existing pairwise matches
	MatchList mlist;
	mlist.clear();
	mlist = original_ml;
	cout << "Starting with " << mlist.size() << " multi-matches\n";

//
// Step 1) Compute guide trees for each profile and join them
//
	NumericMatrix< double > distance1;
	DistanceMatrix( prof1, distance1 );
	NumericMatrix< double > distance2;
	DistanceMatrix( prof2, distance2 );

	// Make a phylogenetic tree
	// use the identity matrix method and convert to a distance matrix
	MuscleInterface& ci = MuscleInterface::getMuscleInterface();	
	string guide_tree_fname1 = CreateTempFileName("guide_tree");
	registerFileToDelete( guide_tree_fname1 );
	ci.CreateTree( distance1, guide_tree_fname1 );
	string guide_tree_fname2 = CreateTempFileName("guide_tree");
	registerFileToDelete( guide_tree_fname2 );
	ci.CreateTree( distance2, guide_tree_fname2 );

	// read the trees
	ifstream tree_file1( guide_tree_fname1.c_str() );
	if( !tree_file1.is_open() )
		throw "Error opening guide tree file";
	PhyloTree< AlignmentTreeNode > tree1;
	tree1.readTree( tree_file1 );
	tree_file1.close();
	ifstream tree_file2( guide_tree_fname2.c_str() );
	if( !tree_file2.is_open() )
		throw "Error opening guide tree file";
	PhyloTree< AlignmentTreeNode > tree2;
	tree2.readTree( tree_file2 );
	tree_file2.close();


	// compute pairwise distances among all nodes
	NumericMatrix< double > distance;
	DistanceMatrix( mlist, distance );
	conservation_distance.resize(boost::extents[seq_count][seq_count]);
	for( uint seqI = 0; seqI < seq_count; ++seqI )
		for( uint seqJ = 0; seqJ < seq_count; ++seqJ )
			if( seqJ > seqI )
				conservation_distance[seqI][seqJ] = distance(seqI,seqJ);
			else
				conservation_distance[seqI][seqJ] = distance(seqJ,seqI);


	if( !collinear_genomes )
	{
		cout << "Calculating pairwise breakpoint distances\n";
		CreatePairwiseBPDistance(bp_distance);
		cout << "bp distance matrix:\n";
		print2d_matrix(bp_distance, cout);
		cout << endl;
	}

	// rescale the conservation distance
	double conservation_range = 1;
	double bp_range = 1;
	for( uint seqI = 0; seqI < seq_count; ++seqI )
		for( uint seqJ = 0; seqJ < seq_count; ++seqJ )
			conservation_distance[seqI][seqJ] = distance(seqI,seqJ) / conservation_range;

	if( !(collinear_genomes && seq_count > 20 ) )
	{
		cout << "genome content distance matrix:\n";
		print2d_matrix(conservation_distance, cout);
		cout << endl;
	}

//
// construct the alignment tree by joining trees from each profile
//
	vector< uint > nsmap1;
	vector< uint > nsmap2;
	makeAlignmentTree( tree1, prof1, nsmap1 );
//	prepareAlignmentTree(tree1);
	makeAlignmentTree( tree2, prof2, nsmap2 );
//	prepareAlignmentTree(tree2);

	alignment_tree.resize( tree1.size() + tree2.size() + 1 );
	// set the sequence appropriately for extant sequences
	node_sequence_map = vector< uint >( alignment_tree.size(), -1 );

	// initialize all nodes to unaligned
	for( node_id_t nodeI = 0; nodeI < alignment_tree.size()-1; nodeI++ )
	{
		if( nodeI < tree1.size() )
		{
			alignment_tree[nodeI].sequence = tree1[nodeI].sequence;
			alignment_tree[nodeI].children = tree1[nodeI].children;
			alignment_tree[nodeI].parents = tree1[nodeI].parents;
			alignment_tree[nodeI].ordering = tree1[nodeI].ordering;
			alignment_tree[nodeI].distance = tree1[nodeI].distance;
			alignment_tree[nodeI].name = tree1[nodeI].name;
			node_sequence_map[nodeI] = nsmap1[nodeI];
		}else{
			alignment_tree[nodeI].sequence = tree2[nodeI-tree1.size()].sequence;
			alignment_tree[nodeI].children = tree2[nodeI-tree1.size()].children;
			alignment_tree[nodeI].parents = tree2[nodeI-tree1.size()].parents;
			alignment_tree[nodeI].ordering = tree2[nodeI-tree1.size()].ordering;
			alignment_tree[nodeI].distance = tree2[nodeI-tree1.size()].distance;
			alignment_tree[nodeI].name = tree2[nodeI-tree1.size()].name;
			for( size_t cI = 0; cI < alignment_tree[nodeI].children.size(); cI++ )
				alignment_tree[nodeI].children[cI] += tree1.size();
			for( size_t pI = 0; pI < alignment_tree[nodeI].parents.size(); pI++ )
				alignment_tree[nodeI].parents[pI] += tree1.size();
			node_sequence_map[nodeI] = nsmap2[nodeI-tree1.size()];
			if( node_sequence_map[nodeI] != (std::numeric_limits<uint>::max)() )
				node_sequence_map[nodeI] += prof1.seq_table.size();
		}

		alignment_tree[nodeI].children_aligned = vector< boolean >( alignment_tree[nodeI].children.size(), true );
		alignment_tree[nodeI].parents_aligned = vector< boolean >( alignment_tree[nodeI].parents.size(), true );
		alignment_tree[nodeI].refined = true;
	}

	alignment_tree.back().children.push_back( tree1.size()-1 );
	alignment_tree.back().children.push_back( alignment_tree.size()-2 );
	alignment_tree.back().distance = 100;
	alignment_tree.back().children_aligned = vector< boolean >( alignment_tree.back().children.size(), true );
	alignment_tree.back().parents_aligned = vector< boolean >( alignment_tree.back().parents.size(), true );
	alignment_tree.back().refined = false;


	getAlignment( interval_list );

}

void ProgressiveAligner::getAlignment( IntervalList& interval_list )
{
	cout << "Aligning...\n";
	// pick each pair of sequences and align until none are left
	while(true)
	{
		node_id_t node1;
		node_id_t node2;
		node_id_t ancestor;
		chooseNextAlignmentPair( alignment_tree, node1, node2, ancestor );
		if( node1 == node2 )
			break;	// all pairs have been aligned

		// this is the last alignable pair in the unrooted tree
		// create a root from which the complete alignment can be extracted
		alignNodes( node1, node2, ancestor );
		if( ancestor == alignment_tree.root )
			break;  // all done
	}

	if( refine )
	{
		// perform iterative refinement
		cout << "Performing final pass iterative refinement\n";
		doGappedAlignment(alignment_tree.root, false);
	}

	// peel off the alignment from the root node
	cout << "root alignment has " << alignment_tree[alignment_tree.root].ordering.size() << " superintervals\n";
	vector< SuperInterval >& a_ivs = alignment_tree[alignment_tree.root].ordering;
	gnSeqI len = 0;
	for( size_t ivI = 0; ivI < a_ivs.size(); ivI++ )
	{
		len += a_ivs[ivI].Length();
	}
	cout << "root alignment length: " << len << endl;


	// for each interval at the root write out the alignment
	for( size_t ivI = 0; ivI < a_ivs.size(); ivI++ )
	{
		GappedAlignment ga(seq_count, a_ivs[ivI].Length());
		extractAlignment(alignment_tree.root, ivI, ga);
		vector<AbstractMatch*> tmp(1, &ga);
		interval_list.push_back( Interval(tmp.begin(), tmp.end()) );
	}
}

/**
 * 
 */

void ProgressiveAligner::align( vector< gnSequence* >& seq_table, IntervalList& interval_list ){
	if( debug_aligner )
	{
		debug_interval = true;
		debug_cga = true;
	}

	seq_count = seq_table.size();
	this->currently_recursing = false;
	interval_list.seq_table = seq_table;

	// find pairwise matches
	MatchList mlist;
	mlist.seq_table = seq_table;

	if( this->breakpoint_penalty == -1 )
		this->breakpoint_penalty = getDefaultBreakpointPenalty( seq_table );
	if( this->bp_dist_estimate_score == -1 )
		this->bp_dist_estimate_score = getDefaultBpDistEstimateMinScore( original_ml.seq_table );
	cout << "using default bp penalty: " << breakpoint_penalty << endl;
	cout << "using default bp estimate min score: " << bp_dist_estimate_score << endl;

	if( this->collinear_genomes )
		this->breakpoint_penalty = -1;

	if( collinear_genomes )
		cout << "\nAssuming collinear genomes...\n";

	mlist.clear();
	mlist = original_ml;
	cout << "Starting with " << mlist.size() << " multi-matches\n";
	cout << "Computing genome content distance matrix...\n";

//
// Step 2) Compute a phylogenetic guide tree using the pairwise matches
//
	NumericMatrix< double > distance;
	SingleCopyDistanceMatrix( mlist, mlist.seq_table, distance );
	cout << "\n\nGenome conservation distance matrix: " << endl;
	distance.print(cout);
	cout << endl;

	bool input_tree_specified = input_guide_tree_fname != "";
	bool output_tree_specified = output_guide_tree_fname != "";
	if( !input_tree_specified )
	{
		// Make a phylogenetic guide tree
		if( !output_tree_specified )
			output_guide_tree_fname = CreateTempFileName("guide_tree");
		input_guide_tree_fname = output_guide_tree_fname;
		cout << "Writing guide tree to " << output_guide_tree_fname << endl;
		MuscleInterface& mi = MuscleInterface::getMuscleInterface();
		mi.CreateTree( distance, output_guide_tree_fname );

	//	ci.SetDistanceMatrix( distance, output_guide_tree_fname );
	}

	conservation_distance.resize(boost::extents[seq_count][seq_count]);
	for( uint seqI = 0; seqI < seq_count; ++seqI )
		for( uint seqJ = 0; seqJ < seq_count; ++seqJ )
			if( seqJ > seqI )
			{
				conservation_distance[seqI][seqJ] = distance(seqI,seqJ);
				conservation_distance[seqJ][seqI] = distance(seqI,seqJ);
			}
			else
			{
				conservation_distance[seqI][seqJ] = distance(seqJ,seqI);
				conservation_distance[seqJ][seqI] = distance(seqJ,seqI);
			}

	cout << "reading tree...\n";
	// load the guide tree
	ifstream tree_file( input_guide_tree_fname.c_str() );
	if( !tree_file.is_open() )
		throw "Error opening guide tree file";
	alignment_tree.readTree( tree_file );
	tree_file.close();

	cout << "initializing alignment tree...\n";
	node_id_t node1;
	node_id_t node2;
	findMidpoint( alignment_tree, node1, node2 );
	moveRootToBranch( alignment_tree, node1, node2 );

	makeAlignmentTree( alignment_tree, mlist, node_sequence_map );
	// midpoint root the tree
//	findMidpoint( alignment_tree, node1, node2 );
//	node_id_t ancestor = 0;
//	if( seq_count > 2 )	// if only two sequences then the tree already has a root
//		ancestor = createAlignmentTreeRoot( alignment_tree, node1, node2 );

	// write out the rooted guide tree, but don't clobber the user's input tree
	if( !input_tree_specified || output_tree_specified )
	{
		ofstream out_tree_file( output_guide_tree_fname.c_str() );
		if( !out_tree_file.is_open() )
			throw "Error opening guide tree file for write";
		alignment_tree.writeTree( out_tree_file );
		out_tree_file.close();
	}

	// ensure the root is the last to get aligned and swap children to canonical order
	extendRootBranches(alignment_tree);


	if( !collinear_genomes )
	{
		// need sol lists for scoring
		vector<SeedOccurrenceList> blah(seq_count);
		swap( blah, sol_list );
//		sol_list = ;
		// temporarily create a weight 11 SML
/*		MatchList w11_mlist;
		w11_mlist.seq_filename = original_ml.seq_filename;
		w11_mlist.seq_table = original_ml.seq_table;
		cout << "Creating weight 11 SMLs for repeat detection\n";
		w11_mlist.CreateMemorySMLs( 11, NULL );
*/
		cout << "Constructing seed occurrence lists for repeat detection\n";
#pragma omp parallel for
		for( int seqI = 0; seqI < seq_count; seqI++ )
		{
			sol_list[seqI].construct(*(mlist.sml_table[seqI]));
//			delete w11_mlist.sml_table[seqI];
		}
//		w11_mlist.sml_table.clear();
	}
	if( !collinear_genomes && use_weight_scaling )
	{
		cout << "Calculating pairwise breakpoint distances\n";
		CreatePairwiseBPDistance(bp_distance);
	}

	// rescale the conservation distance
	if( use_weight_scaling )
	{
		for( uint seqI = 0; seqI < seq_count; ++seqI )
			for( uint seqJ = 0; seqJ < seq_count; ++seqJ )
				conservation_distance[seqI][seqJ] = distance(seqI,seqJ) * conservation_dist_scale;
	}else{
		bp_distance.resize(boost::extents[seq_count][seq_count]);
		for( uint seqI = 0; seqI < seq_count; ++seqI )
			for( uint seqJ = 0; seqJ < seq_count; ++seqJ )
			{
				conservation_distance[seqI][seqJ] = 0;
				bp_distance[seqI][seqJ] = 0;
			}
	}

	if( !collinear_genomes )
	{
		cout << "genome content distance matrix:\n";
		print2d_matrix(conservation_distance, cout);
		cout << endl;
		cout << "bp distance matrix:\n";
		print2d_matrix(bp_distance, cout);
		cout << endl;
	}

	getAlignment( interval_list );
}


// broken and unused function graveyard

}
