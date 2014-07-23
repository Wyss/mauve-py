#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef __GreedyBreakpointElimination_h__
#define __GreedyBreakpointElimination_h__

#include <libMems/AbstractMatch.h>
#include <iostream>
#include <boost/multi_array.hpp>
#include <libMems/PhyloTree.h>
#include <libMems/SubstitutionMatrix.h>
#include <libMems/SeedOccurrenceList.h>
#include <libMems/IntervalList.h>
#include <libMems/LCB.h>
#include <stack>

namespace mems {

extern bool penalize_repeats;

/**
 * A wrapper that maps a match among extant sequences to a match among ancestral and extant seqs
 */
template <class MatchType>
class LcbTrackingMatch
{ 
public:
	MatchType original_match;
	MatchType node_match;
	size_t match_id;	// used to index into global arrays of lcb_id and score
};
typedef LcbTrackingMatch< mems::AbstractMatch* > TrackingMatch;

/** 
 * This class is used to track relationships between LCBs during the LCB determination process.
 */
template <class MatchType>
class TrackingLCB
{
public:
	TrackingLCB(){}
	TrackingLCB( const TrackingLCB& l ){ *this = l; }
	/** Constructs a TrackingLCB from a pairwise LCB */
	TrackingLCB( const mems::LCB& l ){ *this = l; }
	TrackingLCB& operator=( const mems::LCB& l )
	{
		left_end[0] = l.left_end[0];
		left_end[1] = l.left_end[1];
		right_end[0] = l.right_end[0];
		right_end[1] = l.right_end[1];
		left_adjacency[0] = l.left_adjacency[0];
		left_adjacency[1] = l.left_adjacency[1];
		right_adjacency[0] = l.right_adjacency[0];
		right_adjacency[1] = l.right_adjacency[1];
		lcb_id = l.lcb_id;
		weight = l.weight;
		to_be_deleted = false;
		return *this;
	}
	int64 left_end[2];	/**< The left end position of the LCB in each sequence */
	int64 right_end[2];  /**< The right end position of the LCB in each sequence */
	uint left_adjacency[2];	/**< 'Pointers' (actually IDs) to the LCBs on the left in each sequence */
	uint right_adjacency[2];	/**< 'Pointers' (actually IDs) to the LCBs on the right in each sequence */
	double weight;		/**< The weight (or coverage) of this LCB */
	std::vector< MatchType > matches;
	int lcb_id;			/**< A numerical ID that can be assigned to this LCB */
	bool to_be_deleted;
};

/** indicates an LCB identifier hasn't been assigned or is unknown */
const uint LCB_UNASSIGNED = (std::numeric_limits<uint>::max)();

typedef boost::multi_array< std::vector< TrackingLCB< TrackingMatch* > >, 2 > PairwiseLCBMatrix;


/**
 * computes an anchoring score for the matches contained inside an LCB
 */
template< class MatchVector >
double GetPairwiseAnchorScore( 
		MatchVector& lcb, std::vector< genome::gnSequence* >& seq_table, 
		const mems::PairwiseScoringScheme& subst_scoring, mems::SeedOccurrenceList& sol_1, 
		mems::SeedOccurrenceList& sol_2, bool penalize_gaps = false );

class MoveScoreHeapComparator
{
public:
	bool operator()( const std::pair< double, size_t >& a, const std::pair< double, size_t >& b ) const
	{
		return a.first < b.first;	// want to order by > instead of <
	}
};

/**
 * Computes all pairwise LCBs from a set of tracking matches
 */
void getPairwiseLCBs( 
	uint nI, 
	uint nJ, 
	uint dI, 
	uint dJ, 
	std::vector< TrackingMatch* >& tracking_matches, 
	std::vector< TrackingLCB<TrackingMatch*> >& t_lcbs,
	boost::multi_array< double, 3 >& tm_score_array,
	boost::multi_array< size_t, 3 >& tm_lcb_id_array );

/** creates an appropriately sized matrix for mapping individual TrackingMatches to their containing LCBs */
void initTrackingMatchLCBTracking( 
  const std::vector< mems::TrackingMatch >& tracking_matches, 
	size_t n1_count, 
	size_t n2_count, 
	boost::multi_array< size_t, 3 >& tm_lcb_id_array );


/** removes an LCB from an LCB list and coalesces surrounding LCBs.  Returns the number of LCBs removed 
 *  After LCBs are removed, the adjacency list should be processed with filterLCBs()
 *  @param	id_remaps	This is populated with a list of LCB ids that were deleted or coalesced and now have a new LCB id
 *                      for each coalesced LCB, an entry of the form <old id, new id> is added, deleted LCBs have
 *						entries of the form <deleted, -1>.  Entries appear in the order operations were performed
 *						and the function undoLcbRemoval() can undo these operations in reverse order
 */
template< class LcbVector >
uint RemoveLCBandCoalesce( size_t lcbI, uint seq_count, 
						  LcbVector& adjacencies, 
						  std::vector< double >& scores, 
						  std::vector< std::pair< uint, uint > >& id_remaps, 
						  std::vector< uint >& impact_list );


void printMatch( mems::AbstractMatch* m, std::ostream& os );

inline
void printMatch( mems::AbstractMatch* m, std::ostream& os )
{
	for( size_t ii = 0; ii < m->SeqCount(); ++ii )
	{
		if( ii > 0 )
			os << '\t';
		os << "(" << m->Start(ii) << "," << m->RightEnd(ii) << ")";
	}
}

void printProgress( uint prev_prog, uint cur_prog, std::ostream& os );


template< typename PairType >
class LabelSort 
{
public:
	LabelSort( uint seqI ) : ssc( seqI ) {};
	bool operator()( const PairType& pt1, const PairType& pt2 )
	{
		return ssc( pt1.first, pt2.first );
	}
private:
	LabelSort();
	mems::SSC<mems::AbstractMatch> ssc;
};

template<class MatchVector>
void IdentifyBreakpoints( MatchVector& mlist, std::vector<gnSeqI>& breakpoints )
{
	if( mlist.size() == 0 )
		return;
	breakpoints = std::vector<gnSeqI>(1, mlist.size()-1);

	mems::SSC<mems::AbstractMatch> ssc(0);
	std::sort( mlist.begin(), mlist.end(), ssc );
	typedef typename MatchVector::value_type value_type;
	typedef std::pair< value_type, size_t > LabelPairType;
	std::vector< LabelPairType > label_list;
	typename MatchVector::iterator cur = mlist.begin();
	typename MatchVector::iterator end = mlist.end();
	size_t i = 0;
	for( ;cur != end; ++cur )
	{
		label_list.push_back( std::make_pair( *cur, i ) );
		++i;
	}

	uint seq_count = mlist[0]->SeqCount();
	// check for breakpoints in each sequence
	for( uint seqI = 1; seqI < seq_count; seqI++ )
	{
		LabelSort< LabelPairType > ls(seqI); 
		std::sort( label_list.begin(), label_list.end(), ls );

		typename std::vector< LabelPairType >::const_iterator prev = label_list.begin();
		typename std::vector< std::pair< typename MatchVector::value_type, size_t > >::const_iterator iter = label_list.begin();
		typename std::vector< std::pair< typename MatchVector::value_type, size_t > >::const_iterator lab_end = label_list.end();

		bool prev_orient = (*prev).first->Orientation(seqI) == (*prev).first->Orientation(0);
		if( !prev_orient )	// if we start in a different orientation than the ref seq there's a bp here
			breakpoints.push_back(prev->second);

		for( ++iter; iter != lab_end; ++iter )
		{
			bool cur_orient = (*iter).first->Orientation(seqI) == (*iter).first->Orientation(0);
			if( prev_orient == cur_orient &&
				( ( prev_orient && (*prev).second + 1 == (*iter).second) ||
				  ( !prev_orient && (*prev).second - 1 == (*iter).second) 
				)
			  )
			{
				prev_orient = cur_orient;
				++prev;
				continue;	// no breakpoint here
			}

			// always add the last match in a new block (scanning from left to right in seq 0)
			if( prev_orient )
				breakpoints.push_back( prev->second );
			if( !cur_orient )
				breakpoints.push_back( iter->second );

			prev_orient = cur_orient;
			++prev;
		}
		if( prev_orient )
			breakpoints.push_back( prev->second );
	}
	std::sort( breakpoints.begin(), breakpoints.end() );
	std::vector<gnSeqI>::iterator uni = std::unique( breakpoints.begin(), breakpoints.end() );
	breakpoints.erase( uni, breakpoints.end() );
}


template< class MatchVector >
void ComputeLCBs_v2( const MatchVector& meml, const std::vector<gnSeqI>& breakpoints, std::vector< MatchVector >& lcb_list )
{
	// there must be at least one end of a block defined
	if( breakpoints.size() < 1 )
		return;
		
	lcb_list.clear();
	
	// organize the LCBs into different MatchVector instances
	std::vector<gnSeqI>::const_iterator break_iter = breakpoints.begin();
	uint prev_break = 0;	// prev_break is the first match in the current block
	MatchVector lcb;
	for( ; break_iter != breakpoints.end(); ++break_iter ){
		// add the new MatchList to the set if it made the cut
		lcb_list.push_back( lcb );
		lcb_list.back().insert( lcb_list.back().end(), meml.begin() + prev_break, meml.begin() + *break_iter + 1 );
		prev_break = *break_iter + 1;
	}
}


template <class MatchVector>
void computeLCBAdjacencies_v3( const std::vector< MatchVector >& lcb_list, std::vector< double >& weights, std::vector< mems::LCB >& adjacencies )
{
	adjacencies.clear(); // start with no LCB adjacencies
	if( lcb_list.size() == 0 )
		return;	// there aren't any LCBs so there aren't any adjacencies!

	uint seq_count = lcb_list.front().front()->SeqCount();
	uint seqI;
	uint lcbI;
	for( lcbI = 0; lcbI < lcb_list.size(); ++lcbI ){
		mems::LCB lcb;
		std::vector<gnSeqI> left_end;
		std::vector<gnSeqI> length;
		std::vector<bool> orientation;
		FindBoundaries( lcb_list[lcbI], left_end, length, orientation );

		lcb.left_adjacency = std::vector<uint>( left_end.size(), -1 );
		lcb.right_adjacency = std::vector<uint>( left_end.size(), -1 );
		lcb.left_end = std::vector<int64>( left_end.size(), 0 );
		lcb.right_end = std::vector<int64>( left_end.size(), 0 );

		for( seqI = 0; seqI < seq_count; seqI++ ){
			// support "ragged edges" on the ends of LCBs
			if( left_end[seqI] == mems::NO_MATCH )
				continue;
			lcb.left_end[seqI] = left_end[seqI];
			lcb.right_end[seqI] = left_end[seqI] + length[seqI];
			if( !orientation[seqI] )
			{
				lcb.left_end[seqI] = -lcb.left_end[seqI];
				lcb.right_end[seqI] = -lcb.right_end[seqI];
			}
		}
		lcb.lcb_id = adjacencies.size();
		lcb.weight = weights[ lcbI ];
		adjacencies.push_back( lcb );
	}

	for( seqI = 0; seqI < seq_count; seqI++ ){
		mems::LCBLeftComparator llc( seqI );
		std::sort( adjacencies.begin(), adjacencies.end(), llc );
		for( lcbI = 1; lcbI + 1 < lcb_list.size(); lcbI++ ){
			adjacencies[ lcbI ].left_adjacency[ seqI ] = adjacencies[ lcbI - 1 ].lcb_id;
			adjacencies[ lcbI ].right_adjacency[ seqI ] = adjacencies[ lcbI + 1 ].lcb_id;
		}
		if( lcbI == lcb_list.size() )
			lcbI--;	// need to decrement when there is only a single LCB

		// set first and last lcb adjacencies to -1
		adjacencies[ 0 ].left_adjacency[ seqI ] = (uint)-1;
		adjacencies[ lcbI ].right_adjacency[ seqI ] = (uint)-1;
		if( lcbI > 0 ){
			adjacencies[ 0 ].right_adjacency[ seqI ] = adjacencies[ 1 ].lcb_id;
			adjacencies[ lcbI ].left_adjacency[ seqI ] = adjacencies[ lcbI - 1 ].lcb_id;
		}
	}
	mems::LCBIDComparator lic;
	std::sort( adjacencies.begin(), adjacencies.end(), lic );

}

/**
 *  Redesign to be more intuitive.  left_adjacency is always left, regardless of LCB orientation
 */
inline
void computeLCBAdjacencies_v3( mems::IntervalList& iv_list, std::vector< double >& weights, std::vector< mems::LCB >& adjacencies ){
	std::vector< std::vector< mems::Interval* > > nivs;
	for( size_t ivI = 0; ivI < iv_list.size(); ivI++ )
		nivs.push_back( std::vector< mems::Interval* >( 1, &iv_list[ivI] ) );
	computeLCBAdjacencies_v3( nivs, weights, adjacencies );
}

/**
 * Takes a set of filtered LCB adjacencies and an unfiltered set of matches as input
 * returns a filtered set of matches that reflects the LCBs found
 */
template< class MatchVector >
void filterMatches_v2( std::vector< mems::LCB >& adjacencies, std::vector< MatchVector >& lcb_list, std::vector< double >& weights, MatchVector& deleted_matches ){
	if( lcb_list.size() < 1 )
		return;
	MatchVector lcb_tmp = lcb_list[ 0 ];
	lcb_tmp.clear();
	std::vector< MatchVector > filtered_lcbs( lcb_list.size(), lcb_tmp );
	uint lcbI;
	for( lcbI = 0; lcbI < adjacencies.size(); lcbI++ ){
		if( adjacencies[ lcbI ].lcb_id == lcbI ){
			filtered_lcbs[ lcbI ].insert( filtered_lcbs[ lcbI ].end(), lcb_list[ lcbI ].begin(), lcb_list[ lcbI ].end() );
			continue;
		}
		if( adjacencies[ lcbI ].lcb_id == -1 ){
			std::cerr << "weird";
			continue; 	// this one was removed
		}
		if( adjacencies[ lcbI ].lcb_id == -2 )
		{
			deleted_matches.insert( deleted_matches.end(), lcb_list[lcbI].begin(), lcb_list[lcbI].end() );
			continue; 	// this one was removed
		}

		// this one points elsewhere
		// search and update the union/find structure for the target
		std::stack< uint > visited_lcbs;
		visited_lcbs.push( lcbI );
		uint cur_lcb = adjacencies[ lcbI ].lcb_id;
		while( adjacencies[ cur_lcb ].lcb_id != cur_lcb ){
			visited_lcbs.push( cur_lcb );
			cur_lcb = adjacencies[ cur_lcb ].lcb_id;
			if( cur_lcb == -1 || cur_lcb == -2 ){
//				std::cerr << "improper hoodidge\n";
				break;	// this one points to an LCB that got deleted
			}
		}
		while( visited_lcbs.size() > 0 ){
			adjacencies[ visited_lcbs.top() ].lcb_id = cur_lcb;
			visited_lcbs.pop();
		}
		// add this LCB's matches to the target LCB.
		if( cur_lcb != -1 && cur_lcb != -2 )
			filtered_lcbs[ cur_lcb ].insert( filtered_lcbs[ cur_lcb ].end(), lcb_list[ lcbI ].begin(), lcb_list[ lcbI ].end() );
		else
			deleted_matches.insert( deleted_matches.end(), lcb_list[lcbI].begin(), lcb_list[lcbI].end() );
	}


	lcb_list.clear();
	std::vector< double > new_weights;
	for( lcbI = 0; lcbI < filtered_lcbs.size(); lcbI++ ){
		if( filtered_lcbs[ lcbI ].size() > 0 ){
			lcb_list.push_back( filtered_lcbs[ lcbI ] );
			new_weights.push_back( weights[lcbI] );
		}
	}

	// sort the matches inside consolidated LCBs
	mems::MatchStartComparator<mems::AbstractMatch> msc( 0 );
	for( lcbI = 0; lcbI < lcb_list.size(); lcbI++ ){
		std::sort( lcb_list[ lcbI ].begin(), lcb_list[ lcbI ].end(), msc );
	}

	// calculate the LCB adjacencies
	weights = new_weights;
	computeLCBAdjacencies_v3( lcb_list, weights, adjacencies );

}

// predeclared to avoid need to include Islands.h
const score_t INV_SCORE = (std::numeric_limits<score_t>::max)();
void computeMatchScores( const std::string& seq1, const std::string& seq2, const PairwiseScoringScheme& scoring, std::vector<score_t>& scores );
void computeGapScores( const std::string& seq1, const std::string& seq2, const PairwiseScoringScheme& scoring, std::vector<score_t>& scores );


template< class MatchVector >
double GetPairwiseAnchorScore( MatchVector& lcb, 
							  std::vector< genome::gnSequence* >& seq_table, 
							  const mems::PairwiseScoringScheme& subst_scoring, 
							  mems::SeedOccurrenceList& sol_1, 
							  mems::SeedOccurrenceList& sol_2, 
							  bool penalize_gaps )
{
	double lcb_score = 0;
	typename MatchVector::iterator match_iter = lcb.begin();
	for( ; match_iter != lcb.end(); ++match_iter )
	{
		typedef typename MatchVector::value_type MatchPtrType;
		MatchPtrType m = *match_iter;
		std::vector< score_t > scores(m->AlignmentLength(), 0);
		std::vector< std::string > et;
		mems::GetAlignment(*m, seq_table, et);

		// get substitution/gap score
		mems::computeMatchScores( et[0], et[1], subst_scoring, scores );
		if( penalize_gaps )
			mems::computeGapScores( et[0], et[1], subst_scoring, scores );

		// scale match scores by uniqueness
		size_t merI = 0;
		size_t merJ = 0;
		double uni_count = 0;
		double uni_score = 0;
		const size_t m_aln_length = m->AlignmentLength();
		const int64 m_leftend_0 = m->LeftEnd(0);
		const int64 m_leftend_1 = m->LeftEnd(1);
		for( size_t colI = 0; colI < m_aln_length; ++colI )
		{
			if(et[0][colI] != '-' && et[1][colI] != '-' )
			{
				mems::SeedOccurrenceList::frequency_type uni1 = sol_1.getFrequency(m_leftend_0 + merI - 1);
				mems::SeedOccurrenceList::frequency_type uni2 = sol_2.getFrequency(m_leftend_1 + merJ - 1);
				mems::SeedOccurrenceList::frequency_type uniprod = uni1*uni2;
				uniprod = uniprod == 0 ? 1 : uniprod;
				// scale by the uniqueness product, which approximates the number of ways to match up non-unique k-mers
				// in the worst case of a very repetitive match, the score becomes the negative of the match score
				if( scores[colI] > 0 )
				{
					if(penalize_repeats)
						scores[colI] = (score_t)((double)scores[colI] * (2.0 / uniprod)) - scores[colI];
					else
						scores[colI] = (score_t)((mems::SeedOccurrenceList::frequency_type)scores[colI] / uniprod);
				}
			}
			if(et[0][colI] != '-')
				merI++;
			if(et[1][colI] != '-')
				merJ++;
		}


		double m_score = 0;
		for( size_t i = 0; i < scores.size(); ++i )
			if( scores[i] != INV_SCORE )
				m_score += scores[i];

		if( !( m_score > -1000000000 && m_score < 1000000000 ) )
		{
			std::cerr << "scoring error\n";
			genome::breakHere();
		}
		lcb_score += m_score;
	}
	

	return lcb_score;
}



class EvenFasterSumOfPairsBreakpointScorer
{
public:
	EvenFasterSumOfPairsBreakpointScorer( 
		double breakpoint_penalty,
		double minimum_breakpoint_penalty,
		boost::multi_array<double,2> bp_weight_matrix, 
		boost::multi_array<double,2> conservation_weight_matrix,
		std::vector< TrackingMatch* > tracking_match,
		mems::PairwiseLCBMatrix& pairwise_adjacency_matrix,
		std::vector<node_id_t>& n1_descendants,
		std::vector<node_id_t>& n2_descendants,
		boost::multi_array< double, 3 >& tm_score_array,
		boost::multi_array< size_t, 3 >& tm_lcb_id_array,
		size_t seqI_begin,
		size_t seqI_end,
		size_t seqJ_begin,
		size_t seqJ_end
		);

	/**
	 * Returns the number of possible moves a search algorithm may make from the current 
	 * location in LCB search space.  In this case it's simply the total number of pairwise LCBs
	 */
	size_t getMoveCount();

	/** returns the score of the current state */
	double score();

	/** scores a move */
	double operator()( std::pair< double, size_t >& the_move  );

	/** checks whether a particular move is a valid move */
	bool isValid( std::pair< double, size_t >& the_move );

	bool remove( std::pair< double, size_t >& the_move, std::vector< std::pair< double, size_t > >& new_move_list, size_t& new_move_count );

	/** applies a score difference */
	void applyScoreDifference( boost::multi_array< double, 2 >& lcb_score_diff, boost::multi_array< size_t, 2 >& lcb_removed_count );

	/** undoes a score difference, if it wasn't accepted for example */
	void undoScoreDifference( boost::multi_array< double, 2 >& lcb_score_diff, boost::multi_array< size_t, 2 >& lcb_removed_count );

	/** returns the maximum number of new moves generated by any LCB removal */
	size_t getMaxNewMoveCount();

	/** call to indicate that the given LCB has been removed 
	  * @param really_remove	set to false if the move should merely be checked for validity
	  * returns false if the move was invalid
	  */
	bool remove( std::pair< double, size_t >& the_move, bool really_remove, 
		boost::multi_array< double, 2 >& lcb_score_diff, boost::multi_array< size_t, 2 >& lcb_removed_count, 
		bool score_new_moves, std::vector< std::pair< double, size_t > >& new_move_list, size_t& new_move_count );

	/** returns the final set of TrackingMatch values which remain after applying greedy breakpoint elimination */
	std::vector< mems::TrackingMatch* > getResults();

	/** sanity checks all internal data structures */
	bool validate();

protected:
	double bp_penalty;
	boost::multi_array<double,2> bp_weights;
	boost::multi_array<double,2> conservation_weights;
	std::vector< mems::TrackingMatch* > tracking_matches;
	mems::PairwiseLCBMatrix pairwise_adjacencies;
	std::vector<node_id_t> n1_des;
	std::vector<node_id_t> n2_des;

	boost::multi_array< size_t, 2 > pairwise_lcb_count;
	boost::multi_array< double, 2 > pairwise_lcb_score;

	std::vector< TrackingMatch* > deleted_tracking_matches;

	double min_breakpoint_penalty;

private:
	// avoid continuous size lookup
	const size_t seqI_count;
	const size_t seqJ_count;

	// variables used during score computation
	boost::multi_array< std::vector< std::pair< uint, uint > >, 2 > all_id_remaps;
	boost::multi_array< std::vector< uint >, 2 > full_impact_list;
	boost::multi_array< double, 2 > internal_lcb_score_diff[3];
	boost::multi_array< size_t, 2 > internal_lcb_removed_count[3];
	int using_lsd;
	std::vector< double > lsd_zeros;
	std::vector< size_t > lrc_zeros;
	std::vector< double > bogus_scores;
	std::vector< size_t > my_del_lcbs;
	std::vector< size_t > lcb_ids;

	boost::multi_array< double, 3 >& tm_score_array;
	boost::multi_array< size_t, 3 >& tm_lcb_id_array;

	// limit to a range of sequences
	const size_t seqI_first;
	const size_t seqJ_first;
	const size_t seqI_last;
	const size_t seqJ_last;

	// for debugging
	bool first_time;
};


template< class BreakpointScorerType >
int64 greedyBreakpointElimination_v4( std::vector< mems::LCB >& adjacencies, std::vector< double >& scores, BreakpointScorerType& bp_scorer, std::ostream* status_out, size_t g1_tag = 0, size_t g2_tag = 0 );

template< class SearchScorer >
double greedySearch( SearchScorer& spbs );


/**
 * A breakpoint scorer that applies a fixed penalty for each breakpoint that exists in a set of
 * two or more sequences 
 */
class SimpleBreakpointScorer
{
public:
	SimpleBreakpointScorer( std::vector< LCB >& adjacencies, double breakpoint_penalty, bool collinear );

	size_t getMoveCount();

	double score();

	bool isValid( size_t lcbI, double move_score );

	/** return the relative change in score if lcbI were to be removed */
	double operator()( size_t lcbI );

	/** call to indicate that the given LCB has been removed */
	void remove( uint lcbI, std::vector< std::pair< double, size_t > >& new_moves );

private:
	std::vector< mems::LCB > adjs;
	double bp_penalty;
	std::vector< double > scores;
	double total_weight;
	size_t bp_count;
	bool collinear;
};


class GreedyRemovalScorer
{
public:
	GreedyRemovalScorer( std::vector< LCB >& adjacencies, double minimum_weight );

	size_t getMoveCount();

	double score();

	bool isValid( size_t lcbI, double move_score );

	/** return the relative change in score if lcbI were to be removed */
	double operator()( size_t lcbI );

	/** call to indicate that the given LCB has been removed */
	void remove( uint lcbI, std::vector< std::pair< double, size_t > >& new_moves );

private:
	std::vector< mems::LCB > adjs;
	double min_weight;
	std::vector< double > scores;
	double total_weight;
};




template< class BreakpointScorerType >
int64 greedyBreakpointElimination_v4( std::vector< mems::LCB >& adjacencies, 
			std::vector< double >& scores, BreakpointScorerType& bp_scorer, std::ostream* status_out, 
			size_t g1_tag, size_t g2_tag )
{
	// repeatedly remove the low weight LCBs until the minimum weight criteria is satisfied
	uint lcb_count = adjacencies.size();
	double total_initial_lcb_weight = 0;
	for( size_t wI = 0; wI < scores.size(); wI++ )
		total_initial_lcb_weight += scores[wI];
	double total_current_lcb_weight = total_initial_lcb_weight;

	if( adjacencies.size() == 0 )
		return 0;	// nothing can be done
	uint seq_count = adjacencies[0].left_end.size();
	
	double prev_score = bp_scorer.score();
	uint report_frequency = 10;
	uint moves_made = 0;

	size_t move_count = bp_scorer.getMoveCount();
	std::vector< std::pair< double, size_t > > move_heap( move_count * 2 );
	size_t heap_end = move_count;
	for( size_t moveI = 0; moveI < move_count; ++moveI )
	{
		move_heap[moveI].first = bp_scorer(moveI);
		move_heap[moveI].second = moveI;
	}

#ifdef LCB_WEIGHT_LOSS_PLOT
	std::vector< double >::iterator min_iter = std::min_element(scores.begin(), scores.end());
	double mins = *min_iter;
	if( status_out != NULL )
	{
		(*status_out) << g1_tag << '\t' << g2_tag << '\t' << lcb_count << '\t' << 1 - (total_current_lcb_weight / total_initial_lcb_weight) << '\t' << mins << endl;
	}
#endif

	// make a heap of moves ordered by score
	// repeatedly:
	// 1) pop the highest scoring move off the heap
	// 2) attempt to apply the move
	// 3) add any new moves to the heap
	// 4) stop when the highest scoring move no longer increases the score
	MoveScoreHeapComparator mshc;
	std::make_heap( move_heap.begin(), move_heap.end(), mshc );
	while( heap_end > 0 )
	{
		std::pop_heap( move_heap.begin(), move_heap.begin()+heap_end, mshc );
		heap_end--;
		std::pair< double, size_t > best_move = move_heap[ heap_end ];
#ifdef LCB_WEIGHT_LOSS_PLOT
		if( total_current_lcb_weight == scores[best_move.second] )
			break;	// don't remove the last LCB
#else
		if( (best_move.first < 0 ) ||
			total_current_lcb_weight == scores[best_move.second] )
			break;	// can't improve score
#endif

		std::vector< std::pair< double, size_t > > new_moves;
		bool success = bp_scorer.isValid(best_move.second, best_move.first);
		if( !success )
			continue;
		bp_scorer.remove(best_move.second, new_moves);

		
		for( size_t newI = 0; newI < new_moves.size(); newI++ )
		{
			if( heap_end < move_heap.size() )
			{
				heap_end++;
				move_heap[heap_end-1] = new_moves[newI];
				std::push_heap( move_heap.begin(), move_heap.begin()+heap_end, mshc );
			}else{
				// just push the rest on all at once
				size_t prev_size = move_heap.size();
				move_heap.insert( move_heap.end(), new_moves.begin()+newI, new_moves.end() );
				for( size_t newdI = 0; newdI < new_moves.size()-newI; newdI++ )
					std::push_heap( move_heap.begin(), move_heap.begin()+prev_size+newdI+1, mshc );
				heap_end = move_heap.size();
				break;
			}
		}

		total_current_lcb_weight -= scores[best_move.second];
		std::vector< std::pair< uint, uint > > id_remaps;
		std::vector< uint > impact_list;
		lcb_count -= RemoveLCBandCoalesce( best_move.second, adjacencies[0].left_end.size(), adjacencies, scores, id_remaps, impact_list );
#ifdef LCB_WEIGHT_LOSS_PLOT
		mins = scores[best_move.second];
		if( status_out != NULL )
		{
			(*status_out) << g1_tag << '\t' << g2_tag << '\t' << lcb_count << '\t' << 1 - (total_current_lcb_weight / total_initial_lcb_weight) << '\t' << mins << endl;
		}
#endif
		double cur_score = bp_scorer.score();
		prev_score = cur_score;
		moves_made++;
#ifndef LCB_WEIGHT_LOSS_PLOT
		if( status_out != NULL && moves_made % report_frequency == 0 )
			(*status_out) << "move: " << moves_made << " alignment score " << cur_score << std::endl;
#endif
	}

	return 0;
}

extern bool debug_aligner;

/** finds the best anchoring, returns the anchoring score */
template< class SearchScorer >
double greedySearch( SearchScorer& spbs )
{
	double prev_score = spbs.score();
	uint report_frequency = 10;
	uint moves_made = 0;
	if( debug_aligner )
		spbs.validate();
	size_t move_count = spbs.getMoveCount();
	std::vector< double > current_moves( spbs.getMoveCount() );
	// use double the size for the move heap to avoid an almost instant reallocation
	// when a new move gets pushed onto the heap
	size_t heap_end = spbs.getMoveCount();
	std::vector< std::pair< double, size_t > > move_heap( spbs.getMoveCount() * 2 );
	std::vector< std::pair< double, size_t > > new_moves( spbs.getMaxNewMoveCount() + 10 );
	for( size_t moveI = 0; moveI < move_count; ++moveI )
	{
		std::pair< double, size_t > p( 0, moveI );
		double scorediff = spbs(p) - prev_score;
		p.first = scorediff;
		move_heap[moveI] = p;
		current_moves[moveI] = p.first;
	}

	if( debug_aligner )
		spbs.validate();
	// make a heap of moves ordered by score
	// repeatedly:
	// 1) pop the highest scoring move off the heap
	// 2) attempt to apply the move
	// 3) add any new moves to the heap
	// 4) stop when the highest scoring move no longer increases the score
	MoveScoreHeapComparator mshc;
	std::make_heap( move_heap.begin(), move_heap.begin() + heap_end, mshc );
	double successful = 0;
	double invalids = 0;
	int progress = 0;
	int prev_progress = -1;
	while( heap_end > 0 )
	{
		std::pop_heap( move_heap.begin(), move_heap.begin()+heap_end, mshc );
		std::pair< double, size_t > best_move = move_heap[--heap_end];
		if( best_move.first < 0 )
			break;	// can't improve score

		if( best_move.first != current_moves[best_move.second] )
			continue;

		if( !spbs.isValid(best_move) )
		{
			invalids++;
			continue;
		}

		size_t new_move_count = 0;
		bool success = spbs.remove(best_move, new_moves, new_move_count);
		if( !success )
		{
			std::cerr << "numerical instability?  need to investigate this...\n";
//			genome::breakHere();
			invalids++;
			continue;
		}

		successful++;
		if( debug_aligner )
			spbs.validate();

		current_moves[ best_move.second ] = -(std::numeric_limits<double>::max)();
		for( size_t newI = 0; newI < new_move_count; newI++ )
			current_moves[ new_moves[newI].second ] = new_moves[newI].first;

		for( size_t newI = 0; newI < new_move_count; newI++ )
		{
			if( heap_end < move_heap.size() )
			{
				heap_end++;
				move_heap[heap_end-1] = new_moves[newI];
				std::push_heap( move_heap.begin(), move_heap.begin()+heap_end, mshc );
			}else{
				// just push the rest on all at once
				move_heap.resize( (std::min)((size_t)(heap_end * 1.6), heap_end + new_move_count) );
				std::copy( new_moves.begin() + newI, new_moves.begin() + new_move_count, move_heap.begin()+heap_end );
				for( size_t newdI = 0; newdI < new_move_count-newI; newdI++ )
					std::push_heap( move_heap.begin(), move_heap.begin()+heap_end+newdI+1, mshc );
				heap_end = move_heap.size();
				break;
			}
		}

		moves_made++;
		prev_progress = progress;
		progress = (100 * moves_made) / move_count;
		printProgress( prev_progress, progress, std::cout );
//		if( moves_made % report_frequency == 0 )
//			cout << "move: " << moves_made << " alignment score " << cur_score << " success ratio " << successful / invalids << endl;
	}

	return spbs.score();
}

struct AlnProgressTracker
{
	gnSeqI total_len;
	gnSeqI cur_leftend;
	double prev_progress;
};


}	// namespace mems

#endif // __greedyBreakpointElimination_h__

