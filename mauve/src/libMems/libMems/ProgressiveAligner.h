/*******************************************************************************
 * $Id: ProgressiveAligner.h,v 1.23 2004/04/19 23:10:13 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _ProgressiveAligner_h_
#define _ProgressiveAligner_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/SuperInterval.h"
#include "libMems/Aligner.h"
#include "libMems/PhyloTree.h"
#include "libMems/GreedyBreakpointElimination.h"
#include "libMems/CompactGappedAlignment.h"
#include "libMems/Islands.h"
#include <boost/type_traits/remove_pointer.hpp>
#include <boost/multi_array.hpp>
#include "libMems/SeedOccurrenceList.h"
#include "libMems/SubstitutionMatrix.h"
#include "libMems/MatchProjectionAdapter.h"

namespace mems
{

/** controls whether copious debugging tests and output gets written to screen */
extern bool debug_aligner;


/** A class that stores alignment-related information as a node in a phylogenetic tree */
class AlignmentTreeNode : public TreeNode
{
public:
	AlignmentTreeNode() : TreeNode(), refined(false) {};
	std::vector< SuperInterval > ordering;	/**< A total ordering on alignments of sequence contained by leafs below this node */
	std::vector< boolean > parents_aligned;		/**< have parents been aligned? */
	std::vector< boolean > children_aligned;	/**< have children been aligned? */
	genome::gnSequence* sequence;	/**< The sequence associated with this node, NULL for ancestral nodes */
	bool refined;	/**< true if iterative refinement has been applied to the alignment at this node */
};


double getDefaultBreakpointPenalty( std::vector< genome::gnSequence* >& sequences );

/**
 * Computes multiple genome alignments using a progressive alignment algorithm
 */
class ProgressiveAligner : public mems::Aligner
{
public:
	/** 
	 * Constructs an aligner for the specified number of sequences.
	 * @param seq_count 	The number of sequences that will be aligned with this Aligner
	 */
	ProgressiveAligner( uint seq_count );
	ProgressiveAligner( const ProgressiveAligner& al );
	ProgressiveAligner& operator=( const ProgressiveAligner& al );
	~ProgressiveAligner();

	/** sets the breakpoint penalty */
	void setBreakpointPenalty( double bp_penalty ){ breakpoint_penalty = bp_penalty; }
	/** sets the the minimum breakpoint penalty after scaling */
	void setMinimumBreakpointPenalty( double min_bp_penalty ){ min_breakpoint_penalty = min_bp_penalty; }
	/** assume all genomes are collinear when set to true */
	void setCollinear( boolean collinear ){ this->collinear_genomes = collinear; }
	/** use a list of precomputed matches instead of computing them */
	void setPairwiseMatches( mems::MatchList& pair_ml );
	/** use a precomputed guide tree stored in the given file */
	void setInputGuideTreeFileName( std::string& fname ){ this->input_guide_tree_fname = fname; }
	/** write the guide tree stored to the given file */
	void setOutputGuideTreeFileName( std::string& fname ){ this->output_guide_tree_fname = fname; }
	/** set the max length (in columns) of alignments passed to MUSCLE */
	void SetMaxGappedAlignmentLength( size_t len );
	/** set whether a cache database should be used to speed up recursive anchor search */ 
	void SetUseCacheDb( bool cbd ){ this->using_cache_db = cbd; }

	/** Set whether iterative refinement using MUSCLE should be performed (true/false) */
	void setRefinement( bool refine ){ this->refine = refine; }
	/** Set whether iterative refinement using MUSCLE should be performed (true/false) */
	void setGappedAlignment( bool do_gapped_alignment ){ this->gapped_alignment = do_gapped_alignment; }

	void setPairwiseScoringScheme( const mems::PairwiseScoringScheme& pss ){ this->subst_scoring = pss; }

	enum LcbScoringScheme
	{
		AncestralScoring,
		AncestralSumOfPairsScoring,
		ExtantSumOfPairsScoring
	};

	/** set LCB the scoring scheme */
	void setLcbScoringScheme( LcbScoringScheme scheme ){ scoring_scheme = scheme; }
	LcbScoringScheme getLcbScoringScheme(void){ return scoring_scheme; }

	void setUseSeedFamilies( bool use_seed_families ){ this->use_seed_families = use_seed_families; }
	bool getUseSeedFamilies(void){ return this->use_seed_families; }

	void setUseLcbWeightScaling( bool use_weight_scaling ){ this->use_weight_scaling = use_weight_scaling; }
	bool getUseLcbWeightScaling(void){ return this->use_weight_scaling; }

	void setBreakpointDistanceScale( double bp_dist_scale ){ this->bp_dist_scale = bp_dist_scale; }
	double getBreakpointDistanceScale(void){ return this->bp_dist_scale; }

	void setConservationDistanceScale( double conservation_dist_scale ){ this->conservation_dist_scale = conservation_dist_scale; }
	double getConservationDistanceScale(void){ return this->conservation_dist_scale; }

	void setBpDistEstimateMinScore( double min_score ){ this->bp_dist_estimate_score = min_score; }
	double getBpDistEstimateMinScore(void){ return this->bp_dist_estimate_score; }

	/** determine which extant sequences have been aligned at a given node */
	void getAlignedChildren( node_id_t node, std::vector< node_id_t >& descendants );

	/** chooses an ordering for aligned intervals at an ancestor node */
	void createAncestralOrdering( std::vector< mems::Interval* >& interval_list, std::vector< SuperInterval >& ancestral_sequence );

	/** constructs an alignment of node1 and node2 at their ancestor */
	void alignProfileToProfile( node_id_t node1, node_id_t node2, node_id_t ancestor );

	/** align the sequences at the designated pair of alignment tree nodes */
	void alignNodes( node_id_t node1, node_id_t node2, node_id_t ancestor );


	/** Given a set of sequences, construct and output an alignment as an IntervalList */
	void align( std::vector< genome::gnSequence* >& seq_table, mems::IntervalList& interval_list );

	void getPath( node_id_t first_n, node_id_t last_n, std::vector< node_id_t >& path );
	template<class MatchType>
	void propagateDescendantBreakpoints( node_id_t node1, uint seqI, std::vector< MatchType* >& iv_list );
	void linkSuperIntervals( node_id_t node1, uint seqI, node_id_t ancestor );
	void recursiveApplyAncestralBreakpoints( node_id_t ancestor );
	void extractAlignment( node_id_t ancestor, size_t super_iv, mems::GappedAlignment& gal );
	void extractAlignment( node_id_t ancestor, size_t super_iv, mems::CompactGappedAlignment<>& cga );

	void getPairwiseMatches( const std::vector< node_id_t >& node1_seqs, const std::vector< node_id_t >& node2_seqs, Matrix<mems::MatchList>& pairwise_matches );
	void getAncestralMatches( const std::vector< node_id_t > node1_seqs, const std::vector< node_id_t > node2_seqs, node_id_t node1, node_id_t node2, node_id_t ancestor, std::vector< mems::AbstractMatch* >& ancestral_matches );
	void getRepresentativeAncestralMatches( const std::vector< node_id_t > node1_seqs, const std::vector< node_id_t > node2_seqs, node_id_t node1, node_id_t node2, node_id_t ancestor, std::vector< mems::AbstractMatch* >& ancestral_matches );
	
	// functions for recursive anchor search
	
	template<class GappedAlignmentType>
	void recurseOnPairs( const std::vector<node_id_t>& node1_seqs, 
		const std::vector<node_id_t>& node2_seqs, const GappedAlignmentType& iv, 
		Matrix<mems::MatchList>& matches, Matrix< std::vector< mems::search_cache_t > >& search_cache_db, 
		Matrix< std::vector< mems::search_cache_t > >& new_cache_db,
		boost::multi_array< std::vector< std::vector< int64 > >, 2 >& iv_regions);
	void pairwiseAnchorSearch( mems::MatchList& r_list, mems::Match* r_begin, mems::Match* r_end, const mems::AbstractMatch* iv, uint oseqI, uint oseqJ );

	void translateGappedCoordinates( std::vector<mems::AbstractMatch*>& ml, uint seqI, node_id_t extant, node_id_t ancestor );

	void doGappedAlignment( node_id_t ancestor, bool profile_aln );
	void refineAlignment( mems::GappedAlignment& gal, node_id_t ancestor, bool profile_aln, AlnProgressTracker& apt );
	void FixLeftEnds( node_id_t ancestor );
	void ConstructSuperIntervalFromMSA( node_id_t ancestor, size_t ans_siv, mems::GappedAlignment& gal );

	// determines LCBs among each pair of genomes using a somewhat stringent homology 
	// criteria.  fills the distance matrix with the number of breakpoints between each pair
	void CreatePairwiseBPDistance( boost::multi_array<double, 2>& bp_distmat );

	void constructLcbTrackingMatches( node_id_t ancestral_node, std::vector< mems::AbstractMatch* >& ancestral_matches, std::vector< mems::LcbTrackingMatch< mems::AbstractMatch* > >& tracking_matches );

	void pairwiseScoreTrackingMatches( 
						std::vector< mems::TrackingMatch >& tracking_matches, 
						std::vector<node_id_t>& node1_descendants, 
						std::vector<node_id_t>& node2_descendants,
						boost::multi_array< double, 3 >& tm_score_array
						);

	void computeAvgAncestralMatchScores( 
						std::vector< TrackingMatch >& tracking_matches, 
						std::vector<node_id_t>& node1_descendants,
						std::vector<node_id_t>& node2_descendants,
						boost::multi_array< double, 3 >& tm_score_array
						);

	void computeInternalNodeDistances( 
						boost::multi_array<double, 2>& bp_dist_mat, 
						boost::multi_array<double, 2>& cons_dist_mat, 
						std::vector<node_id_t>& node1_descendants,
						std::vector<node_id_t>& node2_descendants);

	bool validateSuperIntervals(node_id_t node1, node_id_t node2, node_id_t ancestor);
	bool validatePairwiseIntervals(node_id_t node1, node_id_t node2, std::vector<mems::Interval*>& pair_iv);


	void alignPP(mems::IntervalList& prof1, mems::IntervalList& prof2, mems::IntervalList& interval_list );

protected:
	void getAlignment( mems::IntervalList& interval_list );

	mems::MatchList original_ml;	/**< The list of matches calculated among all sequences.  Also contains the full sequences and sorted mer lists */
	PhyloTree< AlignmentTreeNode > alignment_tree;
	std::vector< uint > node_sequence_map;
	double breakpoint_penalty;
	double min_breakpoint_penalty;
	std::string input_guide_tree_fname;
	std::string output_guide_tree_fname;
	boolean debug;
	boolean refine;
	bool using_cache_db;

	std::vector< SeedOccurrenceList > sol_list;
	boost::multi_array<double, 2> bp_distance;	/**< pairwise breakpoint distances.  dims will be [seq_count][seq_count] */
	boost::multi_array<double, 2> conservation_distance;	/**< pairwise genome conservation distances.  dims will be [seq_count][seq_count] */

	LcbScoringScheme scoring_scheme;
	bool use_weight_scaling;
	bool use_seed_families;

	double bp_dist_scale;
	double conservation_dist_scale;

	double bp_dist_estimate_score;	/**< the minimum LCB score to use when estimating BP distance.  should be conservative (high) */

	size_t max_gapped_alignment_length;

	mems::PairwiseScoringScheme subst_scoring;
};

extern bool debug_aligner;

	/** Select the next pair of nodes to align
	 *  The chosen pair will either be unaligned extant sequences or unaligned
	 *  ancestral sequences whose descendants have all been aligned.  The chosen pair has
	 *  the shortest path on the tree
	 *  When no sequences remain to be aligned, returns node1 == node2
	 */
void chooseNextAlignmentPair( PhyloTree< AlignmentTreeNode >& alignment_tree, node_id_t& node1, node_id_t& node2, node_id_t& ancestor );

void markAligned( PhyloTree< AlignmentTreeNode >& alignment_tree, node_id_t subject_node, node_id_t neighbor );

node_id_t createAlignmentTreeRoot( PhyloTree< AlignmentTreeNode >& alignment_tree, node_id_t node1, node_id_t node2 );

// homogenizes an alignment tree and ordering to prepare for alignment
void prepareAlignmentTree( PhyloTree< AlignmentTreeNode >& alignment_tree );

inline
ProgressiveAligner::~ProgressiveAligner()
{
	for( size_t mI = 0; mI < original_ml.size(); mI++ )
		original_ml[mI]->Free();
}

template<class T>
class AbsolutComparator
{
public:
	boolean operator()(const T& a, const T& b) const
	{
		return (genome::absolut(a) < genome::absolut(b));
	}
};



template <class MatchVector>
void processNewMatch( uint seqI, MatchVector& new_matches, typename MatchVector::value_type& new_match )
{
	new_match->SetStart( seqI, 0 );
	if( new_match->Multiplicity() > 1 && new_match->Length(seqI) > 0 )
		new_matches.push_back( new_match );
	else
	{
		new_match->Free();
		new_match = NULL;
	}
}
inline
bool checkConsistent(const AbstractMatch* a, const AbstractMatch* b)
{
	bool consistent_overlap = true;
	int64 o = (std::numeric_limits<int64>::max)();
	int64 inter = 0;
	uint seq_count = a->SeqCount();
	for( size_t seqI = 0; seqI < seq_count; seqI++ )
	{
		if(b->LeftEnd(seqI) == 0 || a->LeftEnd(seqI) == 0)
			continue;
		inter++;
		if(o == (std::numeric_limits<int64>::max)())
			o = b->Start(seqI) - a->Start(seqI);
		if(o != b->Start(seqI) - a->Start(seqI))
			consistent_overlap = false;
	}
	consistent_overlap = consistent_overlap && inter > 1;
	return consistent_overlap;
}

/**
 * Delete overlapping regions in favor of the larger match.
 * This code isn't perfect, it can delete too many base pairs in some cases
 * @param	ml	The vector of matches
 * @param	seq_ids	The indexes of sequences in which overlaps should be eliminated
 * @param	eliminate_both	Delete both of the overlapping matches, instead of leaving one remaining
 */
template <class MatchVector>
void EliminateOverlaps_v2( MatchVector& ml, const std::vector< uint >& seq_ids, bool eliminate_both = false ){
	if( ml.size() < 2 )
		return;
	uint seq_count = ml[0]->SeqCount();
	for( uint sidI = 0; sidI < seq_ids.size(); sidI++ ){
		uint seqI = seq_ids[ sidI ];
		mems::SingleStartComparator<mems::AbstractMatch> msc( seqI );
		std::sort( ml.begin(), ml.end(), msc );
		int64 matchI = 0;
		int64 nextI = 0;
		int64 deleted_count = 0;
		MatchVector new_matches;

		// scan forward to first defined match
		for(; matchI != ml.size(); matchI++ )
			if( ml[ matchI ]->Start( seqI ) != mems::NO_MATCH )
				break;

		for(; matchI < ml.size(); matchI++ ){
			if( ml[ matchI ] == NULL )
				continue;
			
			for( nextI = matchI + 1; nextI < ml.size(); nextI++ ){
				if( ml[ nextI ] == NULL )
					continue;

				boolean deleted_matchI = false;
				// check for overlaps
				int64 startI = ml[ matchI ]->Start( seqI );
				int64 lenI = ml[ matchI ]->Length( seqI );
				int64 startJ = ml[ nextI ]->Start( seqI );
				int64 diff =  genome::absolut( startJ ) - genome::absolut( startI ) - lenI;

				if( diff >= 0 )
					break;	// there are no more overlaps

				diff = -diff;
				typename MatchVector::value_type new_match;
				bool mem_iter_smaller = ( ml[ nextI ]->Multiplicity() > ml[ matchI ]->Multiplicity() ) ||
					( ml[ nextI ]->Multiplicity() == ml[ matchI ]->Multiplicity() && ml[ nextI ]->Length(seqI) > ml[ matchI ]->Length(seqI) );

				bool consistent_overlap = checkConsistent( ml[ matchI ], ml[ nextI ] );

				// delete bases from the smaller match
				if( (!consistent_overlap && eliminate_both) || mem_iter_smaller )
				{
					// mem_iter is smaller
					new_match = ml[matchI]->Copy();
					// erase base pairs from new_match
					if( diff >= lenI ){
//							cerr << "Deleting " << **mem_iter << " at the hands of\n" << **next_iter << endl;
						ml[ matchI ]->Free();
						ml[ matchI ] = NULL;
						matchI--;
						deleted_matchI = true;
						deleted_count++;
					}else{
						ml[ matchI ]->CropRight( diff, seqI );
						new_match->CropLeft( new_match->Length(seqI) - diff, seqI );
					}
					processNewMatch( seqI, new_matches, new_match );
				}
				if( (!consistent_overlap && eliminate_both) || !mem_iter_smaller )
				{
					// match_iter is smaller
					new_match = ml[nextI]->Copy();
					// erase base pairs from new_match
					if( diff >= ml[ nextI ]->Length(seqI) ){
//							cerr << "Deleting " << **next_iter << " at the hands of\n" << **mem_iter << endl;
						ml[ nextI ]->Free();
						ml[ nextI ] = NULL;
						deleted_count++;
					}else{
						ml[ nextI ]->CropLeft( diff, seqI );
						new_match->CropRight( new_match->Length(seqI) - diff, seqI );
					}
					processNewMatch( seqI, new_matches, new_match );
				}
				if( deleted_matchI )
					break;
			}
		}

		if( deleted_count > 0 ){
			size_t cur = 0;
			for( size_t mI = 0; mI < ml.size(); ++mI )
				if( ml[mI] != NULL )
					ml[cur++] = ml[mI];
			ml.erase( ml.begin() + cur, ml.end() );
		}
		ml.insert( ml.end(), new_matches.begin(), new_matches.end() );
		new_matches.clear();
	}
}

template <class MatchVector>
void EliminateOverlaps_v2( MatchVector& ml, bool eliminate_both = false )
{
	if( ml.size() < 2 )
		return;	// can't eliminate overlaps between fewer than 2 matches
	uint seq_count = ml[0]->SeqCount();
	std::vector< uint > seq_ids( seq_count );
	for( uint i = 0; i < seq_count; ++i )
		seq_ids[i] = i;
	EliminateOverlaps_v2( ml, seq_ids, eliminate_both );
};

template< class MatchVector >
uint64 SimpleGetLCBCoverage( MatchVector& lcb ){
	typename MatchVector::iterator match_iter = lcb.begin();
	uint64 coverage = 0;
	bool debug = true;
	for( ; match_iter != lcb.end(); ++match_iter ){
		double maxlen = 0;
		double minlen = 0;
		for( uint seqI = 0; seqI < (*match_iter)->SeqCount(); seqI++ )
		{
			if( (*match_iter)->LeftEnd(seqI) != mems::NO_MATCH )
			{
				maxlen += (double)(*match_iter)->Length(seqI);
				if( (*match_iter)->Length(seqI) > minlen )
					minlen = (double)(*match_iter)->Length(seqI);
			}
		}
		double score = exp( ((*match_iter)->AlignmentLength() - minlen) / (maxlen - minlen) );
		score *= maxlen;
		coverage += (uint64)score;
	}
	return coverage;
}

template< class MatchVectorType >
void addUnalignedIntervals_v2( MatchVectorType& iv_list, std::set< uint > seq_set, std::vector<gnSeqI> seq_lengths )
{
	std::vector< mems::LCB > adjacencies;
	uint lcbI;
	uint seqI;
	uint seq_count = seq_lengths.size();


	if( seq_set.size() == 0 )
	{
		// if an empty seq set was passed then assume all seqs
		// should be processed
		for( seqI = 0; seqI < seq_count; seqI++ )
			seq_set.insert( seqI );
	}
	std::vector< std::vector< typename MatchVectorType::value_type > > ymmv;
	for( size_t ivI = 0; ivI < iv_list.size(); ++ivI )
		ymmv.push_back( std::vector< typename MatchVectorType::value_type >( 1, iv_list[ivI] ) );

	std::vector< double > scores( iv_list.size(), 0 );
	computeLCBAdjacencies_v3( ymmv, scores, adjacencies );

	std::vector< int > rightmost;
	for( seqI = 0; seqI < seq_count; seqI++ ){
		rightmost.push_back( -1 );
	}

	for( lcbI = 0; lcbI <= adjacencies.size(); lcbI++ ){
		std::set< uint >::iterator seq_set_iterator = seq_set.begin();
		for( ; seq_set_iterator != seq_set.end(); seq_set_iterator++ ){
			seqI = *seq_set_iterator;
			// scan left
			int leftI;
			if( lcbI < adjacencies.size() ){
// left is always to the left!!
				leftI = adjacencies[ lcbI ].left_adjacency[ seqI ];
			}else
				leftI = rightmost[ seqI ];

			int rightI = lcbI < adjacencies.size() ? lcbI : -1;
// right is always to the right!!
			if( lcbI < adjacencies.size() )
				if( adjacencies[ lcbI ].right_adjacency[ seqI ] == -1 )
					rightmost[ seqI ] = lcbI;
			
			int64 left_start, right_start;
			mems::getGapBounds( seq_lengths, adjacencies, seqI, leftI, rightI, left_start, right_start );
			int64 gap_len =  genome::absolut( right_start ) - genome::absolut( left_start );
			if( gap_len > 0 ){
				mems::Match mm( seq_count );
				mems::Match* m = mm.Copy();
				for( uint seqJ = 0; seqJ < seq_count; seqJ++ ){
					m->SetStart( seqJ, 0 );
				}
				m->SetStart( seqI, left_start );
				m->SetLength( gap_len );
				mems::Interval iv;
				std::vector< mems::AbstractMatch* > tmpvec(1, m);
				iv.SetMatches( tmpvec );
				iv_list.push_back( iv.Copy() );
			}
		}
	}
}

inline
void projectIntervalList( mems::IntervalList& iv_list, std::vector< uint >& projection, std::vector< std::vector< mems::MatchProjectionAdapter* > >& LCB_list, std::vector< mems::LCB >& projected_adjs )
{
	std::vector< size_t > proj(projection.size());
	for( size_t i = 0; i < projection.size(); ++i )
		proj[i] = projection[i];
	std::vector< mems::MatchProjectionAdapter* > mpa_list;
	// construct pairwise Interval projections
	for( size_t corI = 0; corI < iv_list.size(); corI++ )
	{
		size_t projI = 0;
		for( ; projI < projection.size(); ++projI )
			if( iv_list[corI].LeftEnd(projection[projI]) == mems::NO_MATCH )
				break;
		if( projI != projection.size() )
			continue;
		mems::MatchProjectionAdapter mpa_tmp( &iv_list[corI], proj );
		mpa_list.push_back( mpa_tmp.Copy() );
		if( mpa_list.back()->Orientation(0) == mems::AbstractMatch::reverse )
			mpa_list.back()->Invert();
	}
	std::vector< gnSeqI > breakpoints;
	IdentifyBreakpoints( mpa_list, breakpoints );
	ComputeLCBs_v2( mpa_list, breakpoints, LCB_list );
	std::vector< double > lcb_scores( LCB_list.size(), 0 );
	computeLCBAdjacencies_v3( LCB_list, lcb_scores, projected_adjs );
}


template< class MatchType = mems::AbstractMatch >
class GenericMatchSeqManipulator
{
public:
	GenericMatchSeqManipulator( uint seq ) : m_seq(seq) {}
	gnSeqI LeftEnd(MatchType*& m) const{ return m->LeftEnd(m_seq); }
	gnSeqI Length(MatchType*& m) const{ return m->Length(m_seq); }
	void CropLeft(MatchType*& m, gnSeqI amount ) const{ m->CropLeft(amount, m_seq); }
	void CropRight(MatchType*& m, gnSeqI amount ) const{ m->CropRight(amount, m_seq); }
	template< typename ContainerType >
	void AddCopy(ContainerType& c, MatchType*& m) const{ c.push_back( m->Copy() ); }
private:
	uint m_seq;
};

typedef GenericMatchSeqManipulator<> AbstractMatchSeqManipulator;

class SuperIntervalManipulator
{
public:
	gnSeqI LeftEnd(const SuperInterval& siv) const{ return siv.LeftEnd(); }
	gnSeqI Length(const SuperInterval& siv) const{ return siv.Length(); }
	void CropLeft( SuperInterval& siv, gnSeqI amount ) const{ siv.CropLeft( amount );}
	void CropRight( SuperInterval& siv, gnSeqI amount ) const{ siv.CropRight( amount );}
	template< typename ContainerType >
	void AddCopy(ContainerType& c, const SuperInterval& siv) const{ c.push_back( siv ); }
};


// iv_list is a container class that contains pointers to intervals or 
// matches of some sort
// precondition: both bp_list and intervals *must* be sorted
template< class T, class Maniplator >
void applyBreakpoints( std::vector< gnSeqI >& bp_list, std::vector<T>& iv_list, Maniplator& manip )
{

	size_t iv_count = iv_list.size();
	size_t bpI = 0;
	size_t ivI = 0;
	while( ivI < iv_count && bpI < bp_list.size() )
	{
		if( manip.LeftEnd(iv_list[ivI]) == NO_MATCH )
		{
			++ivI;
			continue;	// undefined in seqI, so no breakpoint here
		}
		//  -(ivI)----
		//  -------|--
		if( manip.LeftEnd(iv_list[ivI]) + manip.Length(iv_list[ivI]) <= bp_list[bpI] )
		{
			++ivI;
			continue;
		}
		//  -----(ivI)-
		//  --|--------
		if( bp_list[bpI] <= manip.LeftEnd(iv_list[ivI]) )
		{
			++bpI;
			continue;
		}

		// if split_at isn't 0 then we need to split cur_iv
		// put the left side in the new list and crop cur_iv
		gnSeqI crop_amt = bp_list[bpI] - manip.LeftEnd(iv_list[ivI]);
		manip.AddCopy( iv_list, iv_list[ivI] );
		T& left_iv = iv_list.back();

		manip.CropLeft( iv_list[ivI], crop_amt );
		manip.CropRight( left_iv, manip.Length(left_iv)-crop_amt );
		// restore ordering
		size_t nextI = ivI + 1;
		while( nextI < iv_count && manip.LeftEnd( iv_list[nextI-1] ) > manip.LeftEnd( iv_list[nextI] ) )
		{
			std::swap( iv_list[nextI-1], iv_list[nextI] );
			nextI++;
		}

// assume that crop works correctly and that it's okay to pass matches with NO_MATCH		
/**/
		if( manip.Length( iv_list[ivI] ) == 0 )
		{
			std::cerr << "Big fat generic zero 1\n";
			genome::breakHere();
		}
		if( manip.Length( left_iv ) == 0 )
		{
			std::cerr << "Big fat generic zero 2\n";
			genome::breakHere();
		}
		if( manip.LeftEnd( iv_list[ivI] ) == 0 )
		{
			std::cerr << "uh oh\n";
			genome::breakHere();
		}
		if( manip.LeftEnd( left_iv ) == 0 )
		{
			std::cerr << "uh oh 2\n";
			genome::breakHere();
		}
/**/
	}
}


}

//namespace std {
//	void swap( PhyloTree<mems::AlignmentTreeNode>& a, PhyloTree<mems::AlignmentTreeNode>& b);
//}

#endif // _ProgressiveAligner_h_
