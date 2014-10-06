/*******************************************************************************
 * $Id: Aligner.h,v 1.23 2004/04/19 23:10:13 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _Aligner_h_
#define _Aligner_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnSequence.h"
#include "libMems/DNAMemorySML.h"
#include "libMems/GappedAligner.h"
#include "libMems/MatchList.h"
#include "libMems/Interval.h"
#include "libMems/IntervalList.h"
#include "libMems/MemHash.h"
#include "libMems/MaskedMemHash.h"
#include <map>
#include "libMems/NumericMatrix.h"
#include "libMems/GreedyBreakpointElimination.h"
#include <list>
#include "libMems/LCB.h"
#include "libMUSCLE/threadstorage.h"

namespace mems {

/**
 * A mem labeled with a number.
 * Used by LCB construction algorithm 
 */
class LabeledMem{
public:
	Match* mem;
	uint label;
};

/**
 * Compares Matches labeled with a number.
 * Used by LCB construction algorithm 
 */
class LabeledMemComparator {
public:
	LabeledMemComparator( uint seq ){
		m_seq = seq;
	}
	LabeledMemComparator( LabeledMemComparator& lmc ){
		m_seq = lmc.m_seq;
	}
	boolean operator()(const LabeledMem& a, const LabeledMem& b) const{
		
		int64 a_start = a.mem->Start( m_seq ), b_start = b.mem->Start( m_seq );
		if( a_start == NO_MATCH || b_start == NO_MATCH ){
			if( b_start != NO_MATCH )
				return true;
			return false;
		}
		if(a_start < 0)
			a_start = -a_start;
//			a_start = -a_start + a.mem->Length();
		if(b_start < 0)
			b_start = -b_start;
//			b_start = -b_start + b.mem->Length();
		int64 diff = a_start - b_start;
		return diff < 0;
	}
protected:
	uint m_seq;
private:
	LabeledMemComparator();
};

/**
 * A match with an associated list iterator.
 * Used by LCB construction algorithm 
 */
class PlacementMatch{
public:
	Match* mem;
	std::list< LabeledMem >::iterator iter;
};

/**
 * Compares Matches.
 * Used by LCB construction algorithm 
 */
class PlacementMatchComparator {
public:
	PlacementMatchComparator( uint seq ){
		m_seq = seq;
	}
	PlacementMatchComparator( PlacementMatchComparator& lmc ){
		m_seq = lmc.m_seq;
	}
	boolean operator()(const PlacementMatch& a, const PlacementMatch& b) const{
		
		int64 a_start = a.mem->Start( m_seq ), b_start = b.mem->Start( m_seq );
		if( a_start == NO_MATCH || b_start == NO_MATCH ){
			if( b_start != NO_MATCH )
				return true;
			return false;
		}
		if(a_start < 0)
			a_start = -a_start;
//			a_start = -a_start + a.mem->Length();
		if(b_start < 0)
			b_start = -b_start;
//			b_start = -b_start + b.mem->Length();

		int64 diff = a_start - b_start;
		return diff < 0;
	}
protected:
	uint m_seq;
private:
	PlacementMatchComparator();
};


/** a cache type to remember which intervals have already been searched */
typedef std::pair< mems::Match*, mems::Match* > search_cache_t;


/**
 * Used to find locally colinear blocks (LCBs) and do recursive
 * alignments on the blocks
 * To create an alignment one need only use the align method.
 * LCB lists are typically stored using the IntervalList class.  They can be
 * read and written in interval format using that class.  For input and output
 * of gapped alignments in other formats, see the gnAlignedSequences class.
 * Other methods in this class are available for experimentation.
 */
class Aligner {
public:
	/** 
	 * Constructs an aligner for the specified number of sequences.
	 * @param seq_count 	The number of sequences that will be aligned with this Aligner
	 */
	Aligner( uint seq_count );
	Aligner( const Aligner& al );
	Aligner& operator=( const Aligner& al );

	/**
	 * Performs an alignment.  Takes a MatchList as input and outputs a list of LCBs as an IntervalList.
	 * Several of the options can be used to filter out unlikely LCBs.  If the recursive option is
	 * specified, the regions between matches in each LCB are searched for further homology and a full
	 * gapped alignment is produced.
	 * @param mlist					The MatchList to use as input for the alignment process
	 * @param interval_list			The IntervalList that is created by the alignment process
	 * @param LCB_minimum_density	The minimum density that an LCB may have to be considered a valid block
	 *   						 	This should be a number between 0 and 1.
	 * @param LCB_minimum_range		A misnomer: really it's the minimum number of matching base pairs an LCB 
	 *								must contain to be considered an LCB. Coverage is defined as 
	 *								(length of match) * (# of matching sequences)
	 * @param recursive 			Option for performing a recursive alignment.  If this is set to
	 * 								true, all regions which have gaps will be searched for exact matches.
	 * @param extend_lcbs		If true, attempt to extend the boundaries of LCBs by searching for 
	 *                          additional matches between LCBs
	 * @param tree_filename		The name of the output file to write the phylogenetic guide tree into.  If
	 *                          an empty string is specified then a temporary file is created.  
	 * @throws AlignerError 	may be thrown if an error occurs
	 */
	void align( MatchList& mlist, IntervalList& interval_list, double LCB_minimum_density, double LCB_minimum_range, boolean recursive, boolean extend_lcbs, boolean gapped_alignment, std::string tree_filename = "" );
	
	void Recursion( MatchList& r_list, Match* r_begin, Match* r_end, boolean nway_only = false );
	void GetBestLCB( MatchList& r_list, MatchList& best_lcb );
	void DoSomethingCool( MatchList& mlist, Interval& iv );
	
	/**
	 * Set the minimum size of intervening region between two anchor matches that will
	 * be considered for recursive anchor determination.  When the gaps between two anchors
	 * are less than this cutoff value the region is handed off to the dynamic programming aligner
	 * e.g. ClustalW 
	 */
	void SetMinRecursionGapLength( gnSeqI min_r_gap );

	void SetMaxExtensionIterations( uint ext_iters ){ this->max_extension_iters = ext_iters; }

	void SearchWithinLCB( MatchList& mlist, std::vector< search_cache_t >& new_cache, bool leftmost = false, bool rightmost = false );
	void RecursiveAnchorSearch( MatchList& mlist, gnSeqI minimum_weight, std::vector< MatchList >& LCB_list, boolean entire_genome, std::ostream* status_out = NULL );

	void AlignLCB( MatchList& mlist, Interval& iv );
	void SetGappedAligner( GappedAligner& gal );
	/** forwards the request to whatever gapped aligner is being used */
	void SetMaxGappedAlignmentLength( gnSeqI len );

	/** Set output parameters for permutation matrices */
	void SetPermutationOutput( std::string& permutation_filename, int64 permutation_weight );
	void WritePermutation( std::vector< LCB >& adjacencies, std::string out_filename );

	void SetRecursive( bool value ){ this->recursive = value; }
protected:
	TLS<MemHash> gap_mh;			/**< Used during recursive alignment */
	MaskedMemHash nway_mh;	/**< Used during recursive alignment to find nway matches only */
	uint32 seq_count;		/**< The number of sequences this aligner is working with */
	boolean debug;			/**< Flag for debugging output */
	
	double LCB_minimum_density;
	double LCB_minimum_range;

	uint max_extension_iters;	/**< maximum number of attempts at LCB extension */
	
	int64 cur_min_coverage;	/**< Tracks the minimum weight of the least weight LCB */
	
	gnSeqI min_recursive_gap_length;	/**< Minimum size of gap regions that will be recursed on */

	void consistencyCheck( uint lcb_count, std::vector< LCB >& adjacencies, std::vector< MatchList >& lcb_list, std::vector< int64 >& weights );
	
	boolean recursive;		/**< Set to true if a recursive anchor search/gapped alignment should be performed */
	boolean extend_lcbs;	/**< Set to true if LCB extension should be attempted */
	boolean gapped_alignment;	/**< Set to true to complete a gapped alignment */
	boolean currently_recursing;	/**< True when the recursive search has begun */
	boolean collinear_genomes;	/**< Set to true if all genomes are assumed to be collinear */
	
	GappedAligner* gal;

	std::string permutation_filename;
	int64 permutation_weight;

	std::vector< search_cache_t > search_cache; /**< a list of recursive searches that have already been done */
};

/**
 * Thrown if some error occurs during alignment
 */
CREATE_EXCEPTION( AlignerError );

void transposeMatches( MatchList& mlist, uint seqI, const std::vector< int64 >& seq_regions );

/**
 * Deletes overlapping regions in a set of matches.  Always removes matching base pairs from the
 * match covering fewer bases.  Coverage is defined as (length of match) * (# of matching sequences)
 */
void EliminateOverlaps( MatchList& ml );

/**
 * Function to determine the breakpoints in a set of matches. 
 * Sorts the matches in mlist and returns the indices of breakpoints.
 * This function attempts (sometimes unsuccessfully) to determine subset LCBs.  If a set of
 * matches containing subset LCBs has been passed to it, the resulting breakpoint set may
 * be incorrect.  You have been warned.
 * @param mlist 		A list of matches to search for LCBs.
 * @param breakpoints 	The indices of matches in the sorted match list that are at LCB boundaries
 */
void AaronsLCB( MatchList& mlist, std::set<uint>& breakpoints );


void ComputeLCBs( MatchList& meml, std::set<uint>& breakpoints, std::vector<MatchList>& lcb_list, std::vector<int64>& weights );
void computeLCBAdjacencies_v2( std::vector<MatchList>& lcb_list, std::vector< int64 >& weights, std::vector< LCB >& adjacencies );
void computeLCBAdjacencies_v2( IntervalList& iv_list, std::vector< int64 >& weights, std::vector< LCB >& adjacencies );
void scanLeft( int& left_recurseI, std::vector< LCB >& adjacencies, int min_weight, int seqI );
void scanRight( int& right_recurseI, std::vector< LCB >& adjacencies, int min_weight, int seqI );
void GetLCBCoverage( MatchList& lcb, uint64& coverage );

int64 greedyBreakpointElimination( gnSeqI minimum_weight, std::vector< LCB >& adjacencies, std::vector< int64 >& weights, std::ostream* status_out = NULL );
void filterMatches( std::vector< LCB >& adjacencies, std::vector< MatchList >& lcb_list, std::vector< int64 >& weights );

void CreateGapSearchList( std::vector< LCB >& adjacencies, const std::vector< genome::gnSequence* >& seq_table, std::vector< std::vector< int64 > >& iv_regions, boolean entire_genome );
void SearchLCBGaps( MatchList& new_matches, const std::vector< std::vector< int64 > >& iv_regions, MaskedMemHash& nway_mh );

static const uint MIN_ANCHOR_LENGTH = 9;


/** used for search cache lookups */
class SearchCacheComparator
{
public:
	SearchCacheComparator() : msc(0){};
	bool operator()( const search_cache_t& a, const search_cache_t& b ) const
	{
		bool lt = true;
		if( a.first == NULL )
		{
			if( b.first == NULL )
				lt = false;
		}else if( b.first == NULL )
		{
			lt = false;
		}else if( !msc( a.first, b.first ) )
			lt = false;
		else if( a.second == NULL )
		{
			if( b.second == NULL )
				lt = false;
		}else if( b.second == NULL )
		{
			lt = false;
		}else if( !msc( a.second, b.second ) )
			lt = false;

		return lt;		
	}
protected:
	mems::MatchStartComparator<mems::Match> msc;
};

static SearchCacheComparator cache_comparator;


}

#endif // _Aligner_h_
