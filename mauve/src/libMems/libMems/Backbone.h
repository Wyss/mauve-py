/*******************************************************************************
 * $Id: Backbone.h,v 1.7 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef __Backbone_h__
#define __Backbone_h__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnSequence.h"
#include "libMems/SubstitutionMatrix.h"
#include "libMems/IntervalList.h"
#include "libMems/NumericMatrix.h"
#include "libMems/MatchList.h"
#include "libMems/GappedAlignment.h"
#include "libMems/CompactGappedAlignment.h"
#include "libMems/Aligner.h"
#include "libMems/Islands.h"
#include <boost/multi_array.hpp>

#include <sstream>
#include <vector>

namespace mems {

typedef mems::UngappedLocalAlignment< mems::HybridAbstractMatch<> > ULA;
typedef std::vector< std::vector< ULA* > > backbone_list_t;
// indexed by seqI, seqJ, ivI, hssI (left col, right col)
typedef boost::multi_array< std::vector< std::pair< size_t, size_t > >, 3 > pairwise_genome_hss_t;

class HssDetector;

/** compute the GC content of a set of sequences */
double computeGC( std::vector< genome::gnSequence* >& seq_table );

/**
 * collapse Intervals that are trivially collinear with each other
 */
void collapseCollinear( IntervalList& iv_list );

/**
 * sanity checks for alignment columns that contain only gaps
 */
void checkForAllGapColumns( IntervalList& iv_list );

/**
 * Applies pairwise transitive homology statistics to detect backbone in a single collinear alignment
 * Unaligns any regions found to be non-homologous, returns coordinates of the homologous segments in bb_list
 * @param	m			The input match in which homology detection will be applied
 * @param	seq_table	A sequence table with one gnSequence pointer per match component
 * @param	result		(output) A newly allocated CompactGappedAlignment that contains the resulting alignment of 
 *						homologous sequence.  It is the caller's responsibility to free the memory using AbstractMatch::Free()
 * @param	bb_list		(output) A list of homologous segments among each component of the output match
 * @param	left_homologous	Set to true if the detection code should assume that sequence beyond the left-most alignment
 *							column is homologous sequence
 * @param	right_homologous	Set to true if the detection code should assume that sequence beyond the right-most alignment
 *							column is homologous sequence
 */
void detectAndApplyBackbone( AbstractMatch* m, std::vector< genome::gnSequence* >& seq_table, CompactGappedAlignment<>*& result, backbone_list_t& bb_list, const Params& hmm_params, boolean left_homologous = false, boolean right_homologous = false );

/**
 * Applies pairwise transitive homology statistics to detect backbone in a genome alignment
 * Unaligns any regions found to be non-homologous, returns coordinates of the homologous segments in bb_list
 */
void detectAndApplyBackbone( IntervalList& iv_list, backbone_list_t& bb_list, const Params& hmm_params );

/**
 * Simply detects backbone using the particular algorithm implemented by HssDetector
 */
void detectBackbone( IntervalList& iv_list, backbone_list_t& bb_list, const HssDetector* detector );

/**
 * Writes a backbone column file.  This file type gets used by the Mauve GUI.
 */
void writeBackboneColumns( std::ostream& bb_out, backbone_list_t& bb_list );

/**
 * Writes a backbone sequence coordinate file.  This file type is easier to analyze with statistical packages.
 */
void writeBackboneSeqCoordinates( backbone_list_t& bb_list, IntervalList& iv_list, std::ostream& bb_out );

class HssDetector
{
public:
	typedef std::vector< CompactGappedAlignment<>* > MatchListType;
	virtual void operator() ( 
		const MatchListType& iv_list, 
		std::vector< genome::gnSequence* >& seq_table,  
		hss_array_t& hss_array ) const = 0;
};

class HomologyHmmDetector : public HssDetector
{
public:
	HomologyHmmDetector( const Params& hmm_params, bool left_homologous, bool right_homologous ) :
		p(hmm_params), left(left_homologous), right(right_homologous) {}
	virtual void operator() ( const MatchListType& iv_list, std::vector< genome::gnSequence* >& seq_table, hss_array_t& hss_array ) const
	{
		findHssHomologyHMM( iv_list, seq_table, hss_array, p, left, right );
	}
private:
	const Params& p;
	bool left; 
	bool right;
};

class BigGapsDetector : public HssDetector
{
public:
	BigGapsDetector( size_t big_gap_size ) : big(big_gap_size) {}
	virtual void operator() ( const MatchListType& iv_list, std::vector< genome::gnSequence* >& seq_table, hss_array_t& hss_array ) const
	{
		hss_array_t gap_array;
		findBigGaps( iv_list, seq_table, gap_array, big );
		// we want the cols that represent regions without big gaps...
		HssColsToIslandCols( iv_list, seq_table, gap_array, hss_array );
	}
private:
	size_t big;
};



typedef std::vector< std::pair< int64, int64 > > bb_seqentry_t;
typedef struct bb_entry_s
{
	bb_seqentry_t bb_seq;
	ULA bb_cols;
	size_t iv;
} bb_entry_t;

void addUniqueSegments( std::vector< bb_seqentry_t >& bb_seq_list, size_t min_length = 20 );
void mergeAdjacentSegments( std::vector< bb_seqentry_t >& bb_seq_list );

class BbSeqEntrySorter
{
public:
	BbSeqEntrySorter( size_t seqI ){ m_seq = seqI; }
	bool operator()( const bb_seqentry_t& a, const bb_seqentry_t& b )
	{
		return genome::absolut(a[m_seq].first) < genome::absolut(b[m_seq].first);
	}
	size_t m_seq;
};

inline
void printBbSeq( std::ostream& os, const bb_seqentry_t& bbseq )
{
	for( size_t i = 0; i < bbseq.size(); ++i )
	{
		if( i > 0 )
			os << '\t';
		os << "(" << bbseq[i].first << ", " << bbseq[i].second << ")";
	}
}

inline
void readBackboneSeqFile( std::istream& bbseq_input, std::vector< bb_seqentry_t >& backbone )
{
	std::string cur_line;
	std::getline( bbseq_input, cur_line );	// read off the header line
	while( std::getline( bbseq_input, cur_line ) )
	{
		bb_seqentry_t bb;
		std::stringstream line_str( cur_line );
		int64 lpos = 0;
		while( line_str >> lpos )
		{
			int64 rpos = 0;
			line_str >> rpos;
			bb.push_back( std::make_pair( lpos, rpos ) );
		}
		backbone.push_back(bb);
	}
}

inline
void writeBackboneSeqFile( std::ostream& bbseq_out, std::vector< bb_seqentry_t >& backbone )
{
	if(backbone.size()==0)
		return;	// can't write if there's no backbone!
	for( size_t seqI = 0; seqI < backbone[0].size(); seqI++ )
	{
		if( seqI > 0 )
			bbseq_out << '\t';
		stringstream ss;
		ss << "seq" << seqI;
		bbseq_out << ss.str() << "_leftend\t" << ss.str() << "_rightend";
	}
	bbseq_out << std::endl;
	for( size_t bbI = 0; bbI < backbone.size(); bbI++ )
	{
		for( size_t seqI = 0; seqI < backbone[bbI].size(); seqI++ )
		{
			if( seqI > 0 )
				bbseq_out << '\t';
			bbseq_out << backbone[bbI][seqI].first << '\t' << backbone[bbI][seqI].second;
		}
		bbseq_out << std::endl;
	}
}

inline
void readBackboneColsFile( std::istream& bbcol_input, std::vector< std::pair< size_t, ULA > >& bb_list )
{
	std::string cur_line;
	while( std::getline( bbcol_input, cur_line ) )
	{
		ULA tmp_ula;
		size_t ivI;
		std::stringstream ss( cur_line );
		ss >> ivI;
		size_t left_col;
		size_t len;
		ss >> left_col;
		ss >> len;
		gnSeqI bbseq;
		while( ss >> bbseq )
		{
			tmp_ula.SetStart( bbseq, left_col );
		}
		tmp_ula.SetLength( len );
		bb_list.push_back( std::make_pair( ivI, tmp_ula ) );
	}
}

void makeAllPairwiseGenomeHSS( IntervalList& iv_list, std::vector< CompactGappedAlignment<>* >& iv_ptrs, std::vector< CompactGappedAlignment<>* >& iv_orig_ptrs, pairwise_genome_hss_t& hss_cols, const HssDetector* detector );
void mergePairwiseHomologyPredictions( 	std::vector< CompactGappedAlignment<>* >& iv_orig_ptrs, pairwise_genome_hss_t& hss_cols, std::vector< std::vector< ULA* > >& ula_list );


}

#endif	// __Backbone_h__

