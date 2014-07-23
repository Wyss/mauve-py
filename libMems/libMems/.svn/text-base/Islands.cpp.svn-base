/*******************************************************************************
 * $Id: Islands.cpp,v 1.12 2004/04/19 23:11:19 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing, copying, and modification
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/Islands.h"
#include "libMems/Aligner.h"
#include "libMems/GappedAlignment.h"

using namespace std;
using namespace genome;
namespace mems {

/**
 * Identifies gaps in the alignment between pairs of sequences that are longer than
 * some number of base pairs in length.  Prints islands to an output stream
 */
void simpleFindIslands( IntervalList& iv_list, uint island_size, ostream& island_out ){
	vector< Island > island_list;
	simpleFindIslands( iv_list, island_size, island_list );
	for( size_t isleI = 0; isleI < island_list.size(); isleI++ )
	{
		Island& i = island_list[isleI];
		island_out << i.seqI << '\t' << i.leftI << '\t' << i.rightI << '\t' 
				<< i.seqJ << '\t' << i.leftJ << '\t' << i.rightJ << endl;
	}
}


void simpleFindIslands( IntervalList& iv_list, uint island_size, vector< Island >& island_list ){
	if( iv_list.size() == 0 )
		return;
	for( uint iv_listI = 0; iv_listI < iv_list.size(); iv_listI++ ){
		Interval& iv = iv_list[ iv_listI ];
		gnAlignedSequences gnas;
		iv.GetAlignedSequences( gnas, iv_list.seq_table );
		uint seq_count = iv_list.seq_table.size();
		
		for( uint seqI = 0; seqI < seq_count; seqI++ ){
			uint seqJ;
			for( seqJ = seqI + 1; seqJ < seq_count; seqJ++ ){
				uint columnI = 0;
				gnSeqI curI = 0;
				gnSeqI curJ = 0;
				gnSeqI lastI = 0;
				gnSeqI lastJ = 0;
				for( columnI = 0; columnI < gnas.alignedSeqsSize(); columnI++ ){
					if( gnas.sequences[ seqI ][ columnI ] != '-' )
						curI++;
					if( gnas.sequences[ seqJ ][ columnI ] != '-' )
						curJ++;
					if( toupper( gnas.sequences[ seqI ][ columnI ] ) == 
						toupper( gnas.sequences[ seqJ ][ columnI ] ) &&
						gnas.sequences[ seqJ ][ columnI ] != '-' ){
						// check for an island that was big enough
						if( curI - lastI > island_size ||
							curJ - lastJ > island_size ){
							int64 leftI = iv.Start( seqI );
							int64 rightI = leftI < 0 ? leftI - curI : leftI + curI;
							leftI = leftI < 0 ? leftI - lastI : leftI + lastI;
							int64 leftJ = iv.Start( seqJ );
							int64 rightJ = leftJ < 0 ? leftJ - curJ : leftJ + curJ;
							leftJ = leftJ < 0 ? leftJ - lastJ : leftJ + lastJ;
							Island isle;
							isle.seqI = seqI;
							isle.seqJ = seqJ;
							isle.leftI = leftI;
							isle.leftJ = leftJ;
							isle.rightI = rightI;
							isle.rightJ = rightJ;
							island_list.push_back(isle);
						}
						
						lastI = curI;
						lastJ = curJ;
					}
				}
			}
		}
	}
}


/**
 * Identifies stretches of alignment existing in all sequences that doesn't
 * contain a gap larger than a particular size.  Such regions are considered
 * the backbone of the alignment.
 */
void simpleFindBackbone( IntervalList& iv_list, uint backbone_size, uint max_gap_size, vector< GappedAlignment >& backbone_regions ){
	if( iv_list.size() == 0 )
		return;
	for( uint iv_listI = 0; iv_listI < iv_list.size(); iv_listI++ ){
		Interval& iv = iv_list[ iv_listI ];
		gnAlignedSequences gnas;
		uint seqI;
		uint seq_count = iv_list.seq_table.size();
		vector< int64 > positions( seq_count );
		vector< int64 > starts( seq_count );
		vector< int64 > ends( seq_count );
		vector< uint > gap_size( seq_count, 0 );
		uint seqJ;
		gnSeqI bb_start_col = 0;
		gnSeqI bb_end_col = 0;
		GappedAlignment cur_backbone( seq_count, 0 );
		
		// initialize positions and starts
		for( seqI = 0; seqI < seq_count; seqI++ ){
			positions[ seqI ] = iv_list[ iv_listI ].Start( seqI );
			if( positions[ seqI ] < 0 )
				positions[ seqI ] -= iv_list[ iv_listI ].Length( seqI ) + 1;
		}
		starts = positions;
		ends = positions;

		iv.GetAlignedSequences( gnas, iv_list.seq_table );
		bool backbone = true;	// assume we are starting out with a complete alignment column
		uint columnI = 0;
		vector< int64 > prev_positions;
		for( ; columnI < gnas.alignedSeqsSize(); columnI++ ){
			bool no_gaps = true;
			prev_positions = positions;
			for( seqI = 0; seqI < seq_count; seqI++ ){
				char cur_char = gnas.sequences[ seqI ][ columnI ];
				if( cur_char != '-' && toupper(cur_char) != 'N' ){
					if( gap_size[ seqI ] > max_gap_size && backbone ){
						// end a stretch of backbone here only
						// if the backbone meets size requirements in each
						// sequence.
						for( seqJ = 0; seqJ < seq_count; seqJ++ ){
							if( ends[ seqJ ] - starts[ seqJ ] < backbone_size ){
								break;
							}
						}
						if( seqJ == seq_count ) {
							// it's a legitimate stretch of backbone
							backbone_regions.push_back( cur_backbone );
							uint bbI = backbone_regions.size() - 1;
							vector< string > aln_mat( seq_count );
							for( seqJ = 0; seqJ < seq_count; seqJ++ ){
								if( starts[ seqJ ] < 0 )
									backbone_regions[ bbI ].SetStart( seqJ, ends[ seqJ ] + 1);
								else
									backbone_regions[ bbI ].SetStart( seqJ, starts[ seqJ ] );
								backbone_regions[ bbI ].SetLength( ends[ seqJ ] - starts[ seqJ ], seqJ );
								aln_mat[ seqJ ] = gnas.sequences[ seqJ ].substr( bb_start_col, bb_end_col - bb_start_col + 1);
							}
							backbone_regions[ bbI ].SetAlignment(aln_mat);
							
						}
						// we either just finished backbone or a short area that didn't
						// qualify as backbone
						// look for a new backbone region
						backbone = false;
					}
					positions[ seqI ]++;
					gap_size[ seqI ] = 0;
				}else{
					gap_size[ seqI ]++;
					no_gaps = false;
				}
			}
			if( no_gaps ){
				bb_end_col = columnI;
				ends = positions;
				if( !backbone ){
					starts = prev_positions;
					bb_start_col = columnI;
					backbone = true;
				}
			}
		}

		// check for backbone one last time
		for( seqJ = 0; seqJ < seq_count; seqJ++ ){
			if( ends[ seqJ ] - starts[ seqJ ] < backbone_size ){
				break;
			}
		}
		if( seqJ == seq_count ) {
			// it's a legitimate stretch of backbone
			backbone_regions.push_back( cur_backbone );
			uint bbI = backbone_regions.size() - 1;
			vector< string > aln_mat( seq_count );
			for( seqJ = 0; seqJ < seq_count; seqJ++ ){
				if( starts[ seqJ ] < 0 )
					backbone_regions[ bbI ].SetStart( seqJ, ends[ seqJ ] + 1);
				else
					backbone_regions[ bbI ].SetStart( seqJ, starts[ seqJ ] );
				backbone_regions[ bbI ].SetLength( ends[ seqJ ] - starts[ seqJ ], seqJ );
				aln_mat[ seqJ ] = gnas.sequences[ seqJ ].substr( bb_start_col, bb_end_col - bb_start_col + 1);
			}
			backbone_regions[ bbI ].SetAlignment( aln_mat );
		}
	}
}


void outputBackbone( const vector< GappedAlignment >& backbone_regions, ostream& backbone_out ){
	for( uint bbI = 0; bbI < backbone_regions.size(); bbI++ ){
		for( uint seqJ = 0; seqJ < backbone_regions[ bbI ].SeqCount(); seqJ++ ){
			if( seqJ > 0 )
				backbone_out << '\t';
			int64 bb_rend = backbone_regions[ bbI ].Start( seqJ );
			if( backbone_regions[ bbI ].Start( seqJ ) < 0 )
				bb_rend -= (int64)backbone_regions[ bbI ].Length( seqJ );
			else
				bb_rend += (int64)backbone_regions[ bbI ].Length( seqJ );
			backbone_out << backbone_regions[ bbI ].Start( seqJ ) << '\t' <<  bb_rend;
		}
		backbone_out << endl;
	}
}


// always return the left end of the one to the left and the right of the one to the right

void getGapBounds( vector<gnSeqI>& seq_lengths, vector< LCB >& adjacencies, uint seqJ, int leftI, int rightI, int64& left_start, int64& right_start ){
	if( rightI != -1 )
		right_start = absolut( adjacencies[ rightI ].left_end[ seqJ ] );
	else
		right_start = seq_lengths[seqJ] + 1;
	
	if( leftI != -1 )
		left_start = absolut( adjacencies[ leftI ].right_end[ seqJ ] );
	else
		left_start = 1;
}


void addUnalignedIntervals( IntervalList& iv_list, set< uint > seq_set, vector<gnSeqI> seq_lengths ){
	vector< LCB > adjacencies;
	vector< int64 > weights;
	uint lcbI;
	uint seqI;
	if( seq_lengths.size() == 0 )
		for( seqI = 0; seqI < iv_list.seq_table.size(); seqI++ )
			seq_lengths.push_back(iv_list.seq_table[seqI]->length());

	uint seq_count = seq_lengths.size();


	if( seq_set.size() == 0 )
	{
		// if an empty seq set was passed then assume all seqs
		// should be processed
		for( seqI = 0; seqI < seq_count; seqI++ )
			seq_set.insert( seqI );
	}
	
	weights = vector< int64 >( iv_list.size(), 0 );
	computeLCBAdjacencies_v2( iv_list, weights, adjacencies );

	vector< int > rightmost;
	for( seqI = 0; seqI < seq_count; seqI++ ){
		rightmost.push_back( -1 );
	}
	for( lcbI = 0; lcbI <= adjacencies.size(); lcbI++ ){
		set< uint >::iterator seq_set_iterator = seq_set.begin();
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
			getGapBounds( seq_lengths, adjacencies, seqI, leftI, rightI, left_start, right_start );
			int64 gap_len =  absolut( right_start ) - absolut( left_start );
			if( gap_len > 0 ){
				Match mm( seq_count );
				Match* m = mm.Copy();
				for( uint seqJ = 0; seqJ < seq_count; seqJ++ ){
					m->SetStart( seqJ, 0 );
				}
				m->SetStart( seqI, left_start );
				m->SetLength( gap_len );
				vector<AbstractMatch*> tmp(1, m);
				iv_list.push_back( Interval(tmp.begin(), tmp.end()) );
				m->Free();
			}
		}
	}
}


void findIslandsBetweenLCBs( IntervalList& iv_list, uint island_size, ostream& island_out ){
	IntervalList iv_list_tmp = iv_list;
	addUnalignedIntervals( iv_list_tmp );
	uint seq_count = iv_list.seq_table.size();
	
	for( int ivI = iv_list.size(); ivI < iv_list_tmp.size(); ivI++ ){
		for( uint seqI = 0; seqI < seq_count; seqI++ ){
			if( iv_list_tmp[ ivI ].Length( seqI ) < island_size )
				continue;

			// this is an island, write the LCB island out
			gnSeqI left_end = absolut( iv_list_tmp[ ivI ].Start( seqI ) );
			gnSeqI right_end = left_end + iv_list_tmp[ ivI ].Length( seqI ) - 1;
			island_out << "LCB island:\t" << seqI << '\t' << left_end << '\t' << right_end << endl;
		}
	}
}

}
