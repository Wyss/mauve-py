/*******************************************************************************
 * $Id: Islands.h,v 1.7 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef __Islands_h__
#define __Islands_h__

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
#include <boost/multi_array.hpp>
#include "libMems/HomologyHMM/homology.h"
#include "libMems/Scoring.h"

namespace mems {

/**
 * A class to represent an island in an alignment.  Islands are generally
 * large insertions of a region of sequence relative to
 * another sequence.
 */
class Island{
public:
	uint seqI;
	uint seqJ;
	int64 leftI;
	int64 leftJ;
	int64 rightI;
	int64 rightJ;
};

/**
 * Identifies gaps in the alignment between pairs of sequences that are longer than
 * some number of base pairs in length.  Prints islands to an output stream
 */
void simpleFindIslands( IntervalList& iv_list, uint island_size, std::ostream& island_out );
void findIslandsBetweenLCBs( IntervalList& iv_list, uint island_size, std::ostream& island_out );
void simpleFindIslands( IntervalList& iv_list, uint island_size, std::vector< Island >& island_list );

class HssCols{
public:
	uint seqI;
	uint seqJ;
	size_t left_col;
	size_t right_col;
};

typedef std::vector< HssCols > hss_list_t;
typedef boost::multi_array< hss_list_t, 3 > hss_array_t;

typedef HssCols IslandCols;	// use the same structure for island segs

template<typename MatchVector>
void hssColsToIslandCols( const MatchVector& iv_list, std::vector< genome::gnSequence* >& seq_table, std::vector< HssCols >& hss_list, std::vector< IslandCols >& island_col_list );

/**
 *  Find regions in each sequence that do not belong to any LCB, add them to their own
 * Interval (LCB) in the IntervalList.
 */
void addUnalignedIntervals( IntervalList& iv_list, std::set< uint > seq_set = std::set< uint >(), std::vector<gnSeqI> seq_lengths = std::vector<gnSeqI>() );

/**
 * Identifies stretches of alignment existing in all sequences that doesn't
 * contain a gap larger than a particular size.  Such regions are considered
 * the backbone of the alignment.
 */
void simpleFindBackbone( IntervalList& iv_list, uint backbone_size, uint max_gap_size, std::vector< GappedAlignment >& backbone_regions );

/**
 * writes out a list of backbone regions
 */
void outputBackbone( const std::vector< GappedAlignment >& backbone_regions, std::ostream& backbone_out );

void getGapBounds( std::vector<gnSeqI>& seq_lengths, std::vector< LCB >& adjacencies, uint seqJ, int leftI, int rightI, int64& left_start, int64& right_start );


static char charmap[128];
inline
char* getCharmap()
{
	static bool initialized = false;
	if(initialized)
		return charmap;
	memset(charmap, 0, 128);
	charmap['a'] = 0;
	charmap['c'] = 1;
	charmap['g'] = 2;
	charmap['t'] = 3;
	charmap['-'] = 4;
	charmap['A'] = 0;
	charmap['C'] = 1;
	charmap['G'] = 2;
	charmap['T'] = 3;
	charmap['-'] = 4;
	initialized = true;
	return charmap;
}
// a mapping from pairwise alignment columns to HomologyHMM emission codes
// row/column indices are as given by the charmap above (ACGT- == 01234).
static char colmap[5][5] = {
//    A   C   G   T   -
	{'1','3','4','5','7'},	// A
	{'3','2','6','4','7'},  // C
	{'4','6','2','3','7'},  // G
	{'5','4','3','1','7'},  // T
	{'7','7','7','7','\0'},  // -
};


inline
void findHssHomologyHMM( std::vector< std::string >& aln_table, hss_list_t& hss_list, uint seqI, uint seqJ, const Params& hmm_params,
						boolean left_homologous, boolean right_homologous )
{
	static char* charmap = getCharmap();

	// encode the alignment as column states
	std::string column_states(aln_table[0].size(),'q');
	vector< size_t > col_reference(column_states.size(), (std::numeric_limits<size_t>::max)() );
	size_t refI = 0;
	for( size_t colI = 0; colI < column_states.size(); colI++ )
	{
		char a = charmap[aln_table[seqI][colI]];
		char b = charmap[aln_table[seqJ][colI]];
		column_states[colI] = colmap[a][b];
		if(column_states[colI] != 0 )
			col_reference[refI++] = colI;
	}
	// filter out the gap/gap cols
	std::string::iterator sitr = std::remove(column_states.begin(), column_states.end(), 0);
	column_states.resize(sitr - column_states.begin());

	for( size_t colI = 2; colI < column_states.size(); colI++ )
	{
		if( column_states[colI] == '7' &&
			column_states[colI-1] == '7' &&
			(column_states[colI-2] == '7' || column_states[colI-2] == '8') )
			column_states[colI-1] = '8';
	}
	if( column_states.size() > 1 && column_states[0] == '7' && (column_states[1] == '7' || column_states[1] == '8'))
		column_states[0] = '8';
	if( column_states.size() > 1 && column_states[column_states.size()-1] == '7' && (column_states[column_states.size()-2] == '7'|| column_states[column_states.size()-2] == '8') )
		column_states[column_states.size()-1] = '8';
	// now feed it to the Homology prediction HMM
	string prediction;
	if( right_homologous && !left_homologous )
		std::reverse(column_states.begin(), column_states.end());

	run(column_states, prediction, hmm_params);

	if( right_homologous && !left_homologous )
		std::reverse(prediction.begin(), prediction.end());
	size_t prev_h = 0;
	size_t i = 1;
	for( ; i < prediction.size(); i++ )
	{
		if( prediction[i] == 'H' && prediction[i-1] == 'N' )
		{
			prev_h = i;
		}
		if( prediction[i] == 'N' && prediction[i-1] == 'H' )
		{
			HssCols hc;
			hc.seqI = seqI;
			hc.seqJ = seqJ;
			hc.left_col = col_reference[prev_h];
			hc.right_col = col_reference[i-1];
			hss_list.push_back(hc);
			prev_h = i;
		}
	}
	// get the last one
	if( prediction[i-1] == 'H' )
	{
		HssCols hc;
		hc.seqI = seqI;
		hc.seqJ = seqJ;
		hc.left_col = col_reference[prev_h];
		hc.right_col = col_reference[i-1];
		hss_list.push_back(hc);
	}
}


template< typename MatchVector >
void findHssHomologyHMM( const MatchVector& iv_list, std::vector< genome::gnSequence* >& seq_table,  hss_array_t& hss_array, const Params& hmm_params, boolean left_homologous, boolean right_homologous )
{
	typedef typename MatchVector::value_type MatchType;
	if( iv_list.size() == 0 )
		return;
	uint seq_count = seq_table.size();
	hss_array.resize( boost::extents[seq_count][seq_count][iv_list.size()] );
	for( uint iv_listI = 0; iv_listI < iv_list.size(); iv_listI++ ){
		const MatchType& iv = iv_list[ iv_listI ];
		std::vector< std::string > aln_table;
		GetAlignment( *iv, seq_table, aln_table );
		
		for( uint seqI = 0; seqI < seq_count; seqI++ ){
			uint seqJ;
			for( seqJ = seqI + 1; seqJ < seq_count; seqJ++ ){

				hss_list_t& hss_list = hss_array[seqI][seqJ][iv_listI];
				hss_list.clear();
				findHssHomologyHMM( aln_table, hss_list, seqI, seqJ, hmm_params, left_homologous, right_homologous );
			}
		}
	}
}


template< typename MatchVector >
void HssColsToIslandCols( const MatchVector& iv_list, std::vector< genome::gnSequence* >& seq_table, hss_array_t& hss_array, hss_array_t& island_col_array )
{

	typedef typename MatchVector::value_type MatchType;
	uint seq_count = seq_table.size();
	island_col_array.resize( boost::extents[seq_count][seq_count][iv_list.size()] );
	for( uint iv_listI = 0; iv_listI < iv_list.size(); iv_listI++ ){
		const MatchType& iv = iv_list[ iv_listI ];
		for( uint seqI = 0; seqI < seq_count; seqI++ ){
			uint seqJ;
			for( seqJ = seqI + 1; seqJ < seq_count; seqJ++ ){
				hss_list_t& hss_list = hss_array[seqI][seqJ][iv_listI];
				hss_list_t& island_col_list = island_col_array[seqI][seqJ][iv_listI];
				ComplementHss(iv_list[iv_listI]->AlignmentLength(),hss_list,island_col_list,seqI,seqJ);
			}
		}
	}
}
inline
void ComplementHss( const size_t alignment_length, hss_list_t& hss_list, hss_list_t& island_col_list, uint seqI=0, uint seqJ=0 )
{


	size_t left_col = 0;
	for( size_t hssI = 0; hssI < hss_list.size(); ++hssI )
	{
		if( left_col >= hss_list[hssI].left_col ) 
		{
			left_col = hss_list[hssI].right_col + 1;
			continue;	// handle the case where the HSS starts at col 0
		}
		// ending an island
		IslandCols isle;
		isle.seqI = seqI;
		isle.seqJ = seqJ;
		isle.left_col = left_col;
		isle.right_col = hss_list[hssI].left_col;
		island_col_list.push_back(isle);
		left_col = hss_list[hssI].right_col + 1;
	}

	if( left_col < alignment_length )
	{
		// add the last island
		IslandCols isle;
		isle.seqI = seqI;
		isle.seqJ = seqJ;
		isle.left_col = left_col;
		isle.right_col = alignment_length-1;
		island_col_list.push_back(isle);
	}
}

template< typename MatchVector >
void HssArrayToCga( const MatchVector& iv_list, std::vector< genome::gnSequence* >& seq_table, hss_array_t& hss_array, std::vector< CompactGappedAlignment<>* >& cga_list )
{
	typedef typename MatchVector::value_type MatchType;
	uint seq_count = seq_table.size();
	for( uint iv_listI = 0; iv_listI < iv_list.size(); iv_listI++ ){
		const MatchType& iv = iv_list[ iv_listI ];
		
		CompactGappedAlignment<>* iv_cga = dynamic_cast< CompactGappedAlignment<>* >(iv);
		bool allocated = false;
		if( iv_cga == NULL )
		{
			CompactGappedAlignment<> tmp_cga;
			iv_cga = tmp_cga.Copy();
			new (iv_cga) CompactGappedAlignment<>(*iv);
			allocated = true;
		}
		for( uint seqI = 0; seqI < seq_count; seqI++ ){
			for( uint seqJ = seqI + 1; seqJ < seq_count; seqJ++ ){
				hss_list_t& isle_list = hss_array[seqI][seqJ][iv_listI];
				for( size_t curI = 0; curI < isle_list.size(); ++curI )
				{
					// extract a cga
					CompactGappedAlignment<> tmp_cga;
					cga_list.push_back( tmp_cga.Copy() );
					iv_cga->copyRange( *(cga_list.back()), isle_list[curI].left_col, isle_list[curI].right_col - isle_list[curI].left_col + 1 );
					if( cga_list.back()->LeftEnd(0) == NO_MATCH )
					{
						// this one must have been covering an invalid region (gaps aligned to gaps)
						cga_list.back()->Free();
						cga_list.erase( cga_list.end()-1 );
					}
				}
			}
		}
		if( allocated )
			iv_cga->Free();
	}
}


template< class IntervalListType >
void addUnalignedRegions( IntervalListType& iv_list)
{
	std::vector< AbstractMatch* > new_ivs;
	std::vector< AbstractMatch* > iv_ptrs(iv_list.size());
	for( size_t i = 0; i < iv_list.size(); ++i )
		iv_ptrs[i] = &iv_list[i];
	for( size_t seqI = 0; seqI < iv_list.seq_table.size(); ++seqI )
	{
		SingleStartComparator< AbstractMatch > ssc( seqI );
		std::sort( iv_ptrs.begin(), iv_ptrs.end(), ssc );
		size_t ivI = 0;
		for( ; ivI < iv_ptrs.size(); ++ivI )
			if( iv_ptrs[ivI]->LeftEnd(seqI) != NO_MATCH )
				break;
		std::list< AbstractMatch* > iv_ptr_list;
		iv_ptr_list.insert( iv_ptr_list.end(), iv_ptrs.begin()+ivI, iv_ptrs.end() );
		AddGapMatches( iv_ptr_list, iv_ptr_list.begin(), iv_ptr_list.end(), seqI, 1, iv_list.seq_table[seqI]->length()+1, AbstractMatch::forward, iv_list.seq_table.size() );
		std::list< AbstractMatch* >::iterator iter = iv_ptr_list.begin();
		while( ivI != iv_ptrs.size() && iter != iv_ptr_list.end() )
		{
			if( iv_ptrs[ivI] == *iter )
				ivI++;
			else
				new_ivs.push_back( *iter );
			++iter;
		}
		while( iter != iv_ptr_list.end() )
		{
			new_ivs.push_back( *iter );
			++iter;
		}
	}
	// now add all the new intervals to iv_list
	size_t prev_size = iv_list.size();
	iv_list.resize( iv_list.size() + new_ivs.size() );
	for( size_t newI = 0; newI < new_ivs.size(); ++newI )
	{
		Interval iv( new_ivs.begin() + newI, new_ivs.begin() + newI + 1 );
		iv_list[prev_size + newI] = iv;
		new_ivs[newI]->Free();
	}
}


template< typename MatchVector >
void findBigGaps( const MatchVector& iv_list, std::vector< genome::gnSequence* >& seq_table,  hss_array_t& hss_array, size_t big_gap_size  )
{
	typedef typename MatchVector::value_type MatchType;
	if( iv_list.size() == 0 )
		return;
	uint seq_count = seq_table.size();
	hss_array.resize( boost::extents[seq_count][seq_count][iv_list.size()] );
	for( uint iv_listI = 0; iv_listI < iv_list.size(); iv_listI++ ){
		const MatchType& iv = iv_list[ iv_listI ];
		std::vector< std::string > aln_table;
		GetAlignment( *iv, seq_table, aln_table );
		
		for( uint seqI = 0; seqI < seq_count; seqI++ ){
			uint seqJ;
			for( seqJ = seqI + 1; seqJ < seq_count; seqJ++ )
			{
				if( iv->LeftEnd(seqI) == NO_MATCH || iv->LeftEnd(seqJ) == NO_MATCH )
					continue;

				hss_list_t& hss_list = hss_array[seqI][seqJ][iv_listI];
				hss_list.clear();
				size_t gap_count = 0;
				size_t gap_lend = 0;
				for( size_t cI = 0; cI < aln_table[seqI].size(); cI++ )
				{
					if( aln_table[seqI][cI] == '-' || aln_table[seqJ][cI] == '-' )
					{
						if( aln_table[seqI][cI] == '-' ^ aln_table[seqJ][cI] == '-' )
						{
							if( gap_count == 0 )
								gap_lend = cI;
							gap_count++;
						}
					}else if( gap_count >= big_gap_size )
					{
						HssCols hc;
						hc.seqI = seqI;
						hc.seqJ = seqJ;
						hc.left_col = gap_lend;
						hc.right_col = cI-1;
						hss_list.push_back( hc );
						gap_count = 0;
					}else
						gap_count = 0;
				}
			}
		}
	}
}


}

#endif // __Islands_h__
