/*******************************************************************************
 * $Id: scoreAlignment.cpp,v 1.14 2004/02/28 00:01:31 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/MatchList.h"
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include "libMems/IntervalList.h"
#include "libGenome/gnFilter.h"

using namespace std;
using namespace genome;
using namespace mems;

class IntervalCompare {
public:
	boolean operator()(const pair< gnSeqI, gnSeqI >& a, const pair< gnSeqI, gnSeqI >& b) const{
		if( ( a.first <= b.first && b.second <= a.second ) ||
			( b.first <= a.first && a.second <= b.second ) )
			return false;	// one contains the other, this must be a query, they are equal.
		if( a.first == b.first )
			return a.second < b.second;
		return a.first < b.first;
	}
};

class IntervalMap {
public:
	virtual void add( gnSeqI left, gnSeqI right ) = 0;
	virtual void find( gnSeqI point, vector< uint >& intervals ) const = 0;
};

class TreeIntervalMap : public IntervalMap {
public:
	virtual void add( gnSeqI left, gnSeqI right );
	virtual void find( gnSeqI point, vector< uint >& intervals ) const;
protected:
	map< pair< gnSeqI, gnSeqI >, uint, IntervalCompare > iv_map;
	
};

void TreeIntervalMap::add( gnSeqI left, gnSeqI right ) {
	pair< gnSeqI, gnSeqI > cur_pos;
	cur_pos.first = left;
	cur_pos.second = right;
	iv_map.insert( map< pair< gnSeqI, gnSeqI >, uint, IntervalCompare >::value_type( cur_pos, iv_map.size() ) );
}

void TreeIntervalMap::find( gnSeqI point, vector< uint >& intervals ) const{
	pair< gnSeqI, gnSeqI > cur_loc = pair< gnSeqI, gnSeqI >( point, point );
	map< pair< gnSeqI, gnSeqI >, uint, IntervalCompare >::const_iterator ivmap_iter = iv_map.lower_bound( cur_loc );
	map< pair< gnSeqI, gnSeqI >, uint, IntervalCompare >::const_iterator upper_iter = iv_map.upper_bound( cur_loc );
	while( ivmap_iter != upper_iter ){
		if( !iv_map.key_comp()( cur_loc, ivmap_iter->first ) &&
			!iv_map.key_comp()( ivmap_iter->first, cur_loc ) )
			intervals.push_back( ivmap_iter->second );

		ivmap_iter++;
	}
}

class VectorIntervalMap : public IntervalMap {
public:
	virtual void add( gnSeqI left, gnSeqI right );
	virtual void find( gnSeqI point, vector< uint >& intervals ) const;
protected:
	vector< pair< gnSeqI, gnSeqI > > iv_map;
};

void VectorIntervalMap::add( gnSeqI left, gnSeqI right ) {
	pair< gnSeqI, gnSeqI > cur_pos;
	cur_pos.first = left;
	cur_pos.second = right;
	iv_map.push_back( cur_pos );
}

void VectorIntervalMap::find( gnSeqI point, vector< uint >& intervals ) const{
	for( uint ivI = 0; ivI < iv_map.size(); ivI++ ){
		if( iv_map[ ivI ].first <= point && point <= iv_map[ ivI ].second )
			intervals.push_back( ivI );
	}
}

/**
 * program to score alignments
 * reads in a "correct" alignment and a calculated alignment
 * scores the calculated alignment based on the correct one
 */
int main( int argc, char* argv[] ){
	
	if( argc < 3 ){
		cout << "scoreAlignment <correct alignment> <calculated alignment> [evolved sequence file] [slagan]\n";
		return -1;
	}
	
	boolean debug_mismatches = false;	/**< turns on code to debug mismatches in evolved and aligned base pairs */
	boolean slagan_mode = false;	/**< Set to true if scoring SLAGAN alignments */
	string correct_fname = argv[ 1 ];
	string calculated_fname = argv[ 2 ];
	string evolved_fname;
	if( argc > 3 ){
		debug_mismatches = true;
		evolved_fname = argv[ 3 ];
	}
	if( argc > 4 ){
		string slagan = "slagan";
		if( slagan == argv[ 4 ] )
			slagan_mode = true;
	}
	ifstream correct_in;
	correct_in.open( correct_fname.c_str() );
	if( !correct_in.is_open() ){
		cerr << "Error opening " << correct_fname << endl;
		return -1;
	}
	ifstream calculated_in;
	calculated_in.open( calculated_fname.c_str() );
	if( !calculated_in.is_open() ){
		cerr << "Error opening " << calculated_fname << endl;
		return -1;
	}
try{
	IntervalList correct_ivs;
	IntervalList calculated_ivs;
	correct_ivs.ReadStandardAlignment( correct_in );
	correct_in.close();
	calculated_ivs.ReadStandardAlignment( calculated_in );
	calculated_in.close();
	gnSequence empty_seq;
	vector< gnSequence* > seq_table( correct_ivs[0].SeqCount(), &empty_seq );
	uint seq_count = seq_table.size();
	const gnFilter* comp_filter = gnFilter::DNAComplementFilter();
	
	gnSequence evolved_gnseqs;
	vector< string > evolved_seqs( seq_count );
	if( debug_mismatches ){
		evolved_gnseqs.LoadSource( evolved_fname );
		for( uint i = 0; i < seq_count; i++ ){
			evolved_seqs[ i ] = evolved_gnseqs.contig( i ).ToString();
		}
	}
	
	/** A map of locations of each interval to the interval's array index */
	vector< IntervalMap* > iv_map;
	uint seqI = 0;
	for( ; seqI < seq_count; seqI++ ){
		if( seqI > 0 && slagan_mode ){
			iv_map.push_back( new VectorIntervalMap() );
		}else{
			iv_map.push_back( new TreeIntervalMap() );
		}

		for( uint map_ivI = 0; map_ivI < calculated_ivs.size(); map_ivI++ ){
			pair< gnSeqI, gnSeqI > cur_pos;
			cur_pos.first = absolut( calculated_ivs[ map_ivI ].Start( seqI ) );
			cur_pos.second = cur_pos.first + calculated_ivs[ map_ivI ].Length( seqI ) - 1;
			iv_map[ seqI ]->add( cur_pos.first, cur_pos.second );
		}
	}
	
	// now compare these alignments somehow (use the evil megaloop)
	gnSeqI true_pos = 0;	/**< when a base is correctly aligned to an orthologous base */
	gnSeqI true_neg = 0;	/**< when a base is correctly aligned to a gap */
	gnSeqI false_pos = 0;	/**< when a base is wrongly aligned to another base */
	gnSeqI false_neg = 0;	/**< when a base is wrongly aligned to a gap */
	gnSeqI total = 0;
	gnSeqI unaligned_fn = 0;	/**< tally for errors due to unaligned regions */
	gnSeqI unaligned_tn = 0;

	gnSeqI bad_context = 0;
	gnSeqI multiple_intersection = 0;
	gnSeqI no_j = 0;
	
	for( uint cor_ivI = 0; cor_ivI < correct_ivs.size(); cor_ivI++ ){
		uint calc_ivI = 0;
		int64 calc_iv_lend = 0;
		int64 calc_iv_lendJ = 0;
		boolean parity_match = true;
		gnAlignedSequences cor_gnas;
		gnAlignedSequences calc_gnas;
		correct_ivs[ cor_ivI ].GetAlignedSequences( cor_gnas, seq_table );
		
		for( seqI = 0; seqI < seq_count; seqI++ ){
			int64 cor_iv_lend = correct_ivs[ cor_ivI ].Start( seqI );
			if( cor_iv_lend == NO_MATCH )
				continue;	// not defined in seqI, skip it
				
			for( uint seqJ = 0; seqJ < seq_count; seqJ++ ){
				if( seqI == seqJ )
					continue;

				int64 cor_iv_lendJ = correct_ivs[ cor_ivI ].Start( seqJ );

				/** base index for seqI in correct alignment */
				int64 baseI = cor_iv_lend < 0 ? -correct_ivs[ cor_ivI ].Length( seqI ) + 1 : 0;
				/** base index for seqJ in correct alignment */
				int64 baseJ = cor_iv_lendJ < 0 ? -correct_ivs[ cor_ivI ].Length( seqJ ) + 1 : 0;
				int64 calc_baseI;	/**< The current base pair in sequence I of the calculated alignment */
				int64 calc_baseJ;	/**< The current base pair in sequence J of the calculated alignment */
				int64 calc_colI = 0;	/**< The current column of the calculated alignment */
				// update calc_* variables with the current seqI/seqJ pair
				if( calc_ivI < calculated_ivs.size() && calc_iv_lend != 0 ){
					calc_iv_lend = calculated_ivs[ calc_ivI ].Start( seqI );
					calc_iv_lendJ = calculated_ivs[ calc_ivI ].Start( seqJ );
					calc_baseI = calculated_ivs[ calc_ivI ].Start( seqI );
					calc_baseJ = calculated_ivs[ calc_ivI ].Start( seqJ );
					if( ( calc_iv_lend > 0 && cor_iv_lend > 0 ) || ( calc_iv_lend < 0 && cor_iv_lend < 0 ) ){
						parity_match = true;
						calc_baseI += calc_baseI < 0 ? -calculated_ivs[ calc_ivI ].Length( seqI ) + 1 : 0;
						calc_baseJ += calc_baseJ < 0 ? -calculated_ivs[ calc_ivI ].Length( seqJ ) + 1 : 0;
					}else{
						parity_match = false;
						calc_baseI += calc_baseI > 0 ? calculated_ivs[ calc_ivI ].Length( seqI ) - 1 : 0;
						calc_baseJ += calc_baseJ > 0 ? calculated_ivs[ calc_ivI ].Length( seqJ ) - 1 : 0;
					}
					calc_colI = parity_match ? 0 : calc_gnas.alignedSeqsSize() - 1;
					// scan calc_colI to the first actual residue
					boolean saw_baseJ = false;
					while( true ){
						if( calc_colI < 0 || calc_colI >= calc_gnas.alignedSeqsSize() ){
							cerr << "Error locating residue in alignment, calculated alignment is corrupt\n";
							break;
						}
						if( calc_gnas.sequences[ seqI ][ calc_colI ] == '-' ){
							if( calc_gnas.sequences[ seqJ ][ calc_colI ] != '-' ){
								calc_baseJ += parity_match ? 1 : -1;
								saw_baseJ = true;
							}
							calc_colI += parity_match ? 1 : -1;
								
						}else
							break;
					}
					// if seqJ still contains a gap in calc_baseJ we haven't actually seen calc_baseJ yet
					if( !saw_baseJ && calc_gnas.sequences[ seqJ ][ calc_colI ] == '-' ){
						calc_baseJ += parity_match ? -1 : 1;
					}
				}
				
				for( gnSeqI colI = 0; colI < cor_gnas.alignedSeqsSize(); colI++ ){
					if( cor_gnas.sequences[ seqI ][ colI ] == '-' ){
						if( cor_gnas.sequences[ seqJ ][ colI ] != '-' )
							baseJ++;
						continue;
					}else if( seqJ < seqI && ( cor_gnas.sequences[ seqJ ][ colI ] != '-' )){
						// this one was already scored when seqI had the current value of seqJ
						baseI++;
						baseJ++;
						continue;
					}

					total++;	/** this aligned pair counts towards the totals */
					
					// calculate the actual base index in seqJ for the correct alignment
					int64 cor_baseJ = cor_iv_lendJ + baseJ;

					// check if the current correct alignment entry for seqI is in
					// the current interval of the calculated alignment
					// if not, scan through the calculated intervals until we find the right one
					// also check wether cor_baseJ fits (for the benefit of shuffle-lagan)
					if( calc_iv_lend == 0 || !(absolut( calc_iv_lend ) <= absolut( cor_iv_lend + baseI ) &&
						  absolut( cor_iv_lend + baseI ) < absolut( calc_iv_lend ) + calculated_ivs[ calc_ivI ].Length( seqI ) &&
						  absolut( calc_iv_lendJ ) <= absolut( cor_baseJ ) &&
						  absolut( cor_baseJ ) < absolut( calc_iv_lendJ ) + calculated_ivs[ calc_ivI ].Length( seqJ ) - 1 ) ){

						boolean possibly_incorrect = false;
						vector< uint > possible_ivsI, possible_ivsJ;
						iv_map[ seqI ]->find( absolut( cor_iv_lend + baseI ), possible_ivsI );
						iv_map[ seqJ ]->find( absolut( cor_baseJ ), possible_ivsJ );
						calc_ivI = calculated_ivs.size();
						if( possible_ivsI.size() == 0 )
							no_j++;
						// determine the intersection of possible_ivI and possible_ivJ
						vector< uint > intersection;
						uint pivI = 0;
						for( ; pivI < possible_ivsI.size(); pivI++ ){
							possibly_incorrect = true;
							uint pivJ = 0;
							for( ; pivJ < possible_ivsJ.size(); pivJ++ ){
								int64 s = absolut( calculated_ivs[ possible_ivsJ[ pivJ ] ].Start( seqJ ) );
								if( !(s <= cor_baseJ <= s + calculated_ivs[ possible_ivsJ[ pivJ ] ].Length( seqJ ) - 1 ) )
									cerr << "cor_baseJ doesn't fit!\n";
								if( possible_ivsI[ pivI ] == possible_ivsJ[ pivJ ] )
									intersection.push_back( pivI );
							}
						}
						if( intersection.size() > 0 ){
							calc_ivI = possible_ivsI[ intersection[ 0 ] ];
							calc_iv_lend = calculated_ivs[ calc_ivI ].Start( seqI );
							calc_iv_lendJ = calculated_ivs[ calc_ivI ].Start( seqJ );
						}
						if( intersection.size() > 1 ){
							multiple_intersection++;
						}

						// if we couldn't find baseI anywhere in the calculated alignment then treat
						// it as aligned to a gap, otherwise
						// update the gnAlignedSequences object for the new interval
						if( calc_ivI < calculated_ivs.size() ){
							calculated_ivs[ calc_ivI ].GetAlignedSequences( calc_gnas, seq_table );
							calc_baseI = calc_iv_lend;
							calc_baseJ = calc_iv_lendJ;
							if( ( calc_iv_lend > 0 && cor_iv_lend > 0 ) || ( calc_iv_lend < 0 && cor_iv_lend < 0 ) ){
								parity_match = true;
								calc_baseI += calc_baseI < 0 ? -calculated_ivs[ calc_ivI ].Length( seqI ) + 1 : 0;
								calc_baseJ += calc_baseJ < 0 ? -calculated_ivs[ calc_ivI ].Length( seqJ ) + 1 : 0;
							}else{
								parity_match = false;
								calc_baseI += calc_baseI > 0 ? calculated_ivs[ calc_ivI ].Length( seqI ) - 1 : 0;
								calc_baseJ += calc_baseJ > 0 ? calculated_ivs[ calc_ivI ].Length( seqJ ) - 1 : 0;
							}
							calc_colI = parity_match ? 0 : calc_gnas.alignedSeqsSize() - 1;
							boolean saw_baseJ = false;
							while( true ){
								if( calc_colI < 0 || calc_colI >= calc_gnas.alignedSeqsSize() ){
									cerr << "Error locating residue in alignment, calculated alignment is corrupt\n";
									break;
								}
								if( calc_gnas.sequences[ seqI ][ calc_colI ] == '-' ){
									if( calc_gnas.sequences[ seqJ ][ calc_colI ] != '-' ){
										calc_baseJ += parity_match ? 1 : -1;
										saw_baseJ = true;
									}
									calc_colI += parity_match ? 1 : -1;
										
								}else
									break;
							}
							// if seqJ still contains a gap in calc_baseJ we haven't actually seen calc_baseJ yet
							if( !saw_baseJ && calc_gnas.sequences[ seqJ ][ calc_colI ] == '-' ){
								calc_baseJ += parity_match ? -1 : 1;
							}

						}else{
							if( possibly_incorrect ){
								// aligned to the wrong context
								bad_context++;
								false_pos++;
								if( cor_gnas.sequences[ seqJ ][ colI ]  != '-' )
									baseJ++;
							}else if( cor_gnas.sequences[ seqJ ][ colI ]  != '-' ){
								// wrongly aligned to a gap
								unaligned_fn++;
								false_neg++;
								baseJ++;
							}else{
								// correctly aligned to a gap
								unaligned_tn++;
								true_neg++;
							}
							baseI++;
							calc_iv_lend = 0;	// reset calc_iv_lend
							continue;
						}
					}

					int64 diffI;
					if( parity_match )
						diffI = baseI + cor_iv_lend - calc_baseI;
					else
						diffI = baseI + cor_iv_lend + calc_baseI;
					
					gnSeqI cbI = 0, cbJ = 0;
					while( cbI < diffI ){
						gnSeqI next_colI = parity_match ? calc_colI + 1 : calc_colI - 1;
						if ( next_colI > 100000000 )
							cerr << "bug?\n";
						if( calc_gnas.sequences[ seqI ][ next_colI ] != '-' )
							cbI++;
						if( calc_gnas.sequences[ seqJ ][ next_colI ] != '-' )
							cbJ++;
						calc_colI += parity_match ? 1 : -1;
					}

					calc_baseI += parity_match ? cbI : -cbI;
					calc_baseJ += parity_match ? cbJ : -cbJ;
					// if cor_baseJ == calc_baseJ then this pair of sequences were correctly aligned!
					// classify the correctness of the aligned pair
					char cor_chI = cor_gnas.sequences[ seqI ][ colI ];
					char cor_chJ = cor_gnas.sequences[ seqJ ][ colI ];
					char calc_chI = calc_gnas.sequences[ seqI ][ calc_colI ];
					char calc_chJ = calc_gnas.sequences[ seqJ ][ calc_colI ];
					if( cor_iv_lend < 0 ){
						cor_chI = comp_filter->Filter( cor_chI );
					}
					if( cor_iv_lendJ < 0 ){
						cor_chJ = comp_filter->Filter( cor_chJ );
					}
					if( calc_iv_lend < 0 ){
						calc_chI = comp_filter->Filter( calc_chI );
					}
					if( calc_iv_lendJ < 0 ){
						calc_chJ = comp_filter->Filter( calc_chJ );
					}
					if( cor_chI != calc_chI && debug_mismatches ){
						if( evolved_seqs[ seqI ][ absolut( calc_baseI ) - 1 ] == cor_chI ){
							cerr << "The calculated alignment has incorrect base: " << calc_chI;
							cerr << " instead of " << evolved_seqs[ seqI ][ absolut( calc_baseI ) - 1 ]; 
							cerr << " at " << absolut( calc_baseI ) << " in sequence " << seqI << endl;
						}else{
							cerr << "The \"correct\" alignment has incorrect base: " << cor_chI;
							cerr << " instead of " << evolved_seqs[ seqI ][ absolut( calc_baseI ) - 1 ]; 
							cerr << " at " << absolut( calc_baseI ) << " in sequence " << seqI << endl;
						}
					}

					if( calc_chJ != '-' ){
						// make sure the calculated base actually matches the original sequence
						if( debug_mismatches && calc_chJ != evolved_seqs[ seqJ ][ absolut( calc_baseJ ) - 1 ] ){
							cerr << "The calculated alignment has incorrect base: " << calc_chJ;
							cerr << " instead of " << evolved_seqs[ seqJ ][ absolut( calc_baseJ ) - 1 ]; 
							cerr << " at " << absolut( calc_baseJ ) << " in sequence " << seqJ << endl;
						}
						if( cor_chJ != '-' &&
							( ( parity_match && cor_baseJ == calc_baseJ ) ||
							( !parity_match && cor_baseJ == -calc_baseJ ) ) ){
							true_pos++;
							// sanity check that the bases are really identical:
							if( cor_chI != calc_chI || cor_chJ != calc_chJ )
								cerr << "Calculated alignment contains a different base than the correct!\n";
						}else if( cor_chJ == '-' )
							false_neg++;
						else
							false_pos++;
					}else{
						if( cor_chJ == '-' )
							true_neg++;
						else
							false_pos++;
					}

					if( cor_gnas.sequences[ seqJ ][ colI ] != '-' )
						baseJ++;
					baseI++;
				}
			}
		}
	}

	cout << "Sensitivity: TP / TP + FN = " << (double)(true_pos) / (double)(true_pos + false_neg) << endl;
	cout << "Specificity: TN / TN + FP = " << (double)(true_neg) / (double)(true_neg + false_pos) << endl;
	cout << "TP + TN / total = " << (double)(true_pos + true_neg) / (double)(total) << endl;
	cout << "FP + FN / total = " << (double)(false_pos + false_neg) / (double)(total) << endl;
	cout << "unaligned error = " << (double)unaligned_fn / (double)total << endl;
	cout << "bad_context = " << (double)bad_context / (double)total << endl;
	cout << "multiple_intersection = " << (double)multiple_intersection / (double)total << endl;
	cout << "no_j = " << (double)no_j / (double)total << endl;
	
}catch( gnException& gne ){
	cerr << gne << endl;
}catch( exception& e ){
	cerr << e.what() << endl;
}

}


