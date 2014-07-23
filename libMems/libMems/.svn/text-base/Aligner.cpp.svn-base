/*******************************************************************************
 * $Id: Aligner.cpp,v 1.47 2004/04/19 23:10:30 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing, copying, and modification
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#include "libMems/Aligner.h"
#include "libMems/Islands.h"
#include "libMems/DNAFileSML.h"
#include "libMems/MuscleInterface.h"	// it's the default gapped aligner
#include "libGenome/gnRAWSource.h"
#include "libMems/DistanceMatrix.h"
#include "libMems/Files.h"

#include <map>
#include <fstream>	// for debugging
#include <sstream>
#include <stack>
#include <algorithm>
#include <limits>

using namespace std;
using namespace genome;
namespace mems {


boolean validateLCB( MatchList& lcb );
void validateRangeIntersections( vector< MatchList >& lcb_list  );
bool debug_shite = false;

/**
 * Test code to ensure that an individual LCB is truly collinear
 */
boolean validateLCB( MatchList& lcb ){
	vector< Match* >::iterator lcb_iter = lcb.begin();
	if( lcb.size() == 0 )
		return true;
	uint seq_count = (*lcb_iter)->SeqCount();
	uint seqI = 0;
	boolean complain = false;
	for(; seqI < seq_count; seqI++ ){
		lcb_iter = lcb.begin();
		int64 prev_coord = 0;
		for(; lcb_iter != lcb.end(); lcb_iter++ ){
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

/**
 * Delete overlapping regions in favor of the larger match.
 * This code isn't perfect, it can delete too many base pairs in some cases
 */
void EliminateOverlaps( MatchList& ml ){
	if( ml.size() < 2 )
		return;
	vector< Match* > result_matches;
	uint seq_count = ml[0]->SeqCount();
	for( uint seqI = 0; seqI < seq_count; seqI++ ){
		SingleStartComparator<AbstractMatch> msc( seqI );
		sort( ml.begin(), ml.end(), msc );
		int64 matchI = 0;
		int64 nextI = 0;
		int64 deleted_count = 0;
		vector< Match* > new_matches;

		// scan forward to first defined match
		for(; matchI != ml.size(); matchI++ )
			if( ml[ matchI ]->Start( seqI ) != NO_MATCH )
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
				int64 lenI = ml[ matchI ]->Length();
				int64 startJ = ml[ nextI ]->Start( seqI );
//				int64 diff =  absolut( startJ ) - absolut( startI ) - lenI;
				int64 diff =  absolut( startJ ) - absolut( startI ) - lenI;

				if( diff < 0 ){
					diff = -diff;
					Match* new_match;
					// delete bases from the smaller match
//					if( ml[ nextI ]->Length() * ml[ nextI ]->Multiplicity() >= 
//						lenI * ml[ matchI ]->Multiplicity() ){
					if( ( ml[ nextI ]->Multiplicity() > ml[ matchI ]->Multiplicity() ) ||
 						( ml[ nextI ]->Multiplicity() == ml[ matchI ]->Multiplicity() && ml[ nextI ]->Length() > ml[ matchI ]->Length() ) ){
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
							if( startI > 0 ){
								ml[ matchI ]->CropEnd( diff );
								new_match->CropStart( new_match->Length() - diff );
							}else{
								ml[ matchI ]->CropStart( diff );
								new_match->CropEnd( new_match->Length() - diff );
							}
						}
					}else{
						// match_iter is smaller
						new_match = ml[nextI]->Copy();
						// erase base pairs from new_match
						if( diff >= ml[ nextI ]->Length() ){
//							cerr << "Deleting " << **next_iter << " at the hands of\n" << **mem_iter << endl;
							ml[ nextI ]->Free();
							ml[ nextI ] = NULL;
							deleted_count++;
						}else{
							if( startJ > 0 ){
								ml[ nextI ]->CropStart( diff );
								new_match->CropEnd( new_match->Length() - diff );
							}else{
								ml[ nextI ]->CropEnd( diff );
								new_match->CropStart( new_match->Length() - diff );
							}
						}

					}
					new_match->SetStart( seqI, 0 );
					if( new_match->Multiplicity() > 1 && new_match->Length() > 0 )
						new_matches.push_back( new_match );
					else
					{
						new_match->Free();
						new_match = NULL;
					}
					if( deleted_matchI )
						break;
				}else
					break;	// there are no more overlaps
			}
//			if( nextI > 1 )
//				cerr << "There were " << nextI << " overlaps\n";
//			if( nextI > config_value_2 )
//				__asm(nop);
		}

		if( deleted_count > 0 ){
			result_matches.reserve( ml.size() - deleted_count );
			for( int64 copyI = 0; copyI < ml.size(); copyI++ ){
				if( ml[ copyI ] != NULL )
					result_matches.push_back( ml[ copyI ] );
			}
			ml.clear();
			ml.insert( ml.end(), result_matches.begin(), result_matches.end() );
		}
		ml.insert( ml.end(), new_matches.begin(), new_matches.end() );
		result_matches.clear();
		new_matches.clear();
	}
		
}


const gnSeqI default_min_r_gap_size = 200;
Aligner::Aligner( uint seq_count ) :
debug(false),
seq_count(seq_count),
min_recursive_gap_length(default_min_r_gap_size),
collinear_genomes(false),
gal(&(MuscleInterface::getMuscleInterface())),
permutation_weight(-1),
cur_min_coverage(-1),
max_extension_iters(4)
{}

Aligner::Aligner( const Aligner& al ) :
//gap_mh( al.gap_mh ),
nway_mh( al.nway_mh ),
seq_count( al.seq_count ),
debug( al.debug),
LCB_minimum_density( al.LCB_minimum_density),
LCB_minimum_range( al.LCB_minimum_range ),
cur_min_coverage( al.cur_min_coverage),
min_recursive_gap_length( al.min_recursive_gap_length ),
collinear_genomes( al.collinear_genomes ),
gal( al.gal ),
permutation_weight( al.permutation_weight ),
permutation_filename( al.permutation_filename ),
max_extension_iters( al.max_extension_iters )
{}

Aligner& Aligner::operator=( const Aligner& al )
{
	gap_mh = al.gap_mh;
	nway_mh = al.nway_mh;
	seq_count = al.seq_count;
	debug = al.debug;
	
	LCB_minimum_density = al.LCB_minimum_density;
	LCB_minimum_range = al.LCB_minimum_range;
	
	cur_min_coverage = al.cur_min_coverage;
	min_recursive_gap_length = al.min_recursive_gap_length;
	collinear_genomes = al.collinear_genomes;

	gal = al.gal;

	permutation_weight = al.permutation_weight;
	permutation_filename = al.permutation_filename;

	max_extension_iters = al.max_extension_iters;

	return *this;
}

void Aligner::SetMinRecursionGapLength( gnSeqI min_r_gap ) {
	min_recursive_gap_length = min_r_gap;
}

void Aligner::SetGappedAligner( GappedAligner& gal ){
	this->gal = &(gal);
}

void Aligner::SetMaxGappedAlignmentLength( gnSeqI len ){
	gal->SetMaxAlignmentLength( len );
}


/* returns true if all labels between start_label and end_label are contained in the no_match_labels set */
void scanLabels( set< uint >& no_match_labels, uint& start_label, boolean forward ){
	uint labelI;
	// scan no_match_labels for consecutive labels starting at start_label until one is missing
	if( forward ){
		for( labelI = start_label + 1; ; labelI++){
			set< uint >::iterator  label_iter = no_match_labels.find( labelI );
			if( label_iter == no_match_labels.end() ){
				start_label = labelI - 1;
				break;
			}
		}
	}else{
		for( labelI = start_label; labelI > 0; labelI--){
			set< uint >::iterator  label_iter = no_match_labels.find( labelI - 1 );
			if( label_iter == no_match_labels.end() ){
				start_label = labelI;
				break;
			}
		}
	}
}

boolean checkCollinearity( Match* m1, Match* m2 ){
	for( uint seqI = 0; seqI < m1->SeqCount(); seqI++ ){
		if( m1->Start( seqI ) == NO_MATCH ||
			m2->Start( seqI ) == NO_MATCH )
			continue;
		if((( m1->Start( seqI ) > 0 &&
			m2->Start( seqI ) > 0 ) ||
			(m1->Start( seqI ) < 0 &&
			m2->Start( seqI ) < 0 )) &&
			m1->Start( seqI ) <= m2->Start( seqI ) )
			continue;
		return false;
	}
	return true;
}

void scanFit( list< LabeledMem >& pair_list, list< LabeledMem >::iterator& list_iter, Match* new_match, uint sort_seq ){

	list< LabeledMem >::iterator cur_iter = list_iter;
	list< LabeledMem >::iterator last_iter = list_iter;
//	int64 initial_start = absolut( list_iter->mem->Start( sort_seq ) );
	int64 initial_start = absolut( list_iter->mem->Start( sort_seq ) );

	uint match_count = 0;
	for(; last_iter != pair_list.end(); ++last_iter ){
		if( last_iter->mem->Start( sort_seq ) == NO_MATCH ){
			++match_count;
			continue;
		}
//		if( absolut( last_iter->mem->Start( sort_seq ) ) < initial_start ||
//			absolut( last_iter->mem->Start( sort_seq ) ) > new_match->Start( sort_seq ) )
		if( absolut( last_iter->mem->Start( sort_seq ) ) < initial_start ||
			absolut( last_iter->mem->Start( sort_seq ) ) > new_match->Start( sort_seq ) )
			break;
		++match_count;
	}
	vector< vector< int > > score_vector;
	score_vector.reserve( new_match->SeqCount() - sort_seq - 1 );
	for( uint seqI = sort_seq + 1; seqI < new_match->SeqCount(); ++seqI ){
		vector< int > sv;
		score_vector.push_back( sv );
		score_vector[ score_vector.size() - 1 ].reserve( match_count );
	}
	uint matchI = 0;
	for(; cur_iter != last_iter; ++cur_iter ){
		
		for( uint seqI = sort_seq + 1; seqI < new_match->SeqCount(); ++seqI ){
			int64 p_start = cur_iter->mem->Start( seqI );
			int64 m_start = new_match->Start( seqI );
			p_start = p_start < 0 ? -p_start : p_start;
			m_start = m_start < 0 ? -m_start : m_start;
			if( m_start == NO_MATCH ){
				score_vector[ seqI - sort_seq - 1 ].push_back( 0 );
			}else if( p_start == NO_MATCH ){
				score_vector[ seqI - sort_seq - 1 ].push_back( 0 );
			}else if( p_start < m_start ){
				score_vector[ seqI - sort_seq - 1 ].push_back( 1 );
			}else
				score_vector[ seqI - sort_seq - 1 ].push_back( -1 );
		}
	}
	vector< int > scores;
	scores.reserve( match_count );
	for( matchI = match_count; matchI > 0; matchI-- )
		scores.push_back( 0 );
	for( uint seqI = 0; seqI < new_match->SeqCount() - sort_seq - 1; ++seqI ){
		boolean redefined = false;
		for( matchI = match_count; matchI > 0; matchI-- ){
			if( !redefined ){
				if( score_vector[ seqI ][ matchI - 1 ] >= 0 ){
					if( score_vector[ seqI ][ matchI - 1 ] == 1 )
						redefined = true;
					++scores[ matchI - 1 ];
				}
			}else{
				if( score_vector[ seqI ][ matchI - 1 ] == -1 )
					redefined = false;
			}
		}
	}
	// find the first highest scoring match
	cur_iter = list_iter;
	int max_score = 0;
	for( matchI = 0; matchI < match_count; ++matchI ){
		if( scores[ matchI ] > max_score ){
			max_score = scores[ matchI ];
			list_iter = cur_iter;
		}
		++cur_iter;
	}
}

/**
 * Aaron's subset LCB algorithm.  
 */
void AaronsLCB( MatchList& mlist, set<uint>& breakpoints ){
	breakpoints.clear(); // make sure this is empty
	if( mlist.size() == 0 )
		return;
	// can only look for breakpoints if there is more than one match!!
	if( mlist.size() == 1 ){
		breakpoints.insert( 0 );
		return;
	}
	uint seq_count = mlist[0]->SeqCount();

	SingleStartComparator<AbstractMatch> msc( 0 );
	sort( mlist.begin(), mlist.end(), msc );
	vector<Match*>::iterator mem_iter = mlist.begin();
	list<LabeledMem> pair_list;
	
	map<uint, Match*> debug_label_map;
	boolean debugging = false;
	
	
	list< PlacementMatch > placement_list;
	
	for(; mem_iter != mlist.end(); ++mem_iter ){
		if( (*mem_iter)->Start( 0 ) != NO_MATCH ){		
			// add this one to the list.
			LabeledMem lm;
			lm.mem = *mem_iter;
			lm.label = 0;
			pair_list.push_back( lm );
		}else{
			PlacementMatch pm;
			pm.mem = *mem_iter;
			pm.iter = pair_list.end();
			placement_list.push_back( pm );
		}
	}
	LabeledMemComparator lmc( 0 );
	pair_list.sort( lmc );
	list< LabeledMem >::iterator pair_iter = pair_list.begin();
	for(; pair_iter != pair_list.end(); ++pair_iter ){
		PlacementMatch pm;
		pm.mem = pair_iter->mem;
		pm.iter = pair_iter;
		placement_list.push_back( pm );
	}
	
	// place all the subset matches from each sequence in the correct place in pair_list.
	for( uint seqI = 1; seqI < seq_count; ++seqI ){
		PlacementMatchComparator pmc( seqI );
		placement_list.sort( pmc );
		list< PlacementMatch >::iterator placement_prev;
		list< PlacementMatch >::iterator placement_iter = placement_list.begin();
		if( placement_iter->iter == pair_list.end() &&
			placement_iter->mem->Start( seqI ) != NO_MATCH ){
			LabeledMem lm;
			lm.mem = placement_iter->mem;
			lm.label = 0;
			pair_list.insert( pair_list.begin(), lm );
			placement_iter->iter = pair_list.begin();
		}

		for( ++placement_iter; placement_iter != placement_list.end(); ++placement_iter ){
			placement_prev = placement_iter;
			placement_prev--;
			
			if( placement_iter->iter != pair_list.end() )
				continue;
			
			if( placement_iter->mem->Start( seqI ) == NO_MATCH )
				continue;
			
			list< LabeledMem >::iterator insert_iter = placement_prev->iter;
			if( insert_iter == pair_list.end() || placement_prev->mem->Start( seqI ) == NO_MATCH )
				insert_iter = pair_list.begin();
			else{
				if( insert_iter->mem->Start( seqI ) < 0 ){
					// invert if necessary and insert before
					if( placement_iter->mem->Start( seqI ) > 0 )
						placement_iter->mem->Invert();
					if( !checkCollinearity( placement_iter->mem, insert_iter->mem ) ){
						placement_iter->mem->Invert();
						scanFit( pair_list, insert_iter, placement_iter->mem, seqI );
						++insert_iter;
					}
				}else{
					// insert in the earliest place this match fits with surrounding matches
					scanFit( pair_list, insert_iter, placement_iter->mem, seqI );
					++insert_iter;
				}
			}
			
			LabeledMem lm;
			lm.mem = placement_iter->mem;
			lm.label = 0;
			pair_list.insert( insert_iter, lm );
			placement_iter->iter = insert_iter;
			placement_iter->iter--;
		}
	}
	boolean debug_labels = false;
	ofstream debug_label_file;
	if( debug_labels )
		debug_label_file.open( "label_debug.txt" );
	// number the LabeledMems in the pair_list
	uint cur_label = 0;
	mlist.clear();
	vector< LabeledMem > pair_vec;
	pair_vec.reserve( pair_list.size() );
	mlist.reserve( pair_list.size() );
	for( pair_iter = pair_list.begin(); pair_iter != pair_list.end(); ++pair_iter ){
		pair_iter->label = cur_label++;
		mlist.push_back( pair_iter->mem );
		pair_vec.push_back( *pair_iter );
		if( debug_labels ){
			debug_label_map.insert( map<uint, Match*>::value_type( pair_iter->label, pair_iter->mem ) );
			debug_label_file << pair_iter->label << '\t' << (*pair_iter->mem) << endl;
		}
	}
	if( debug_labels )
		debug_label_file.close();
	
	breakpoints.clear();
	pair_list.clear();
	vector< LabeledMem >::iterator pair_vec_iter;
	for( uint seqI = 1; seqI < seq_count; seqI++ ){
		// sort the list on the current genome
		LabeledMemComparator lmc( seqI );
		sort( pair_vec.begin(), pair_vec.end(), lmc );
		set< uint > no_match_labels;

		// debugging code
/*		stringstream debug_filename;
		debug_filename << "label_sort_" << seqI << ".txt";
		ofstream debug_file( debug_filename.str().c_str() );
		for( uint pairI = 0; pairI < pair_vec.size(); pairI++ ){
			debug_file << pair_vec[ pairI ].label << *pair_vec[ pairI ].mem << endl;
		}
		debug_file.close();
*/		// end debugging code
		
		pair_vec_iter = pair_vec.begin();
		uint block_start = pair_vec_iter->label;
		uint break_label = 0;
		for( ++pair_vec_iter; pair_vec_iter != pair_vec.end(); ++pair_vec_iter ){
			vector<LabeledMem>::iterator pair_prev = pair_vec_iter;
			pair_prev--;
			break_label = 0;
			uint scan_label = 0;
			if( pair_prev->mem->Start( seqI ) == NO_MATCH ){
				no_match_labels.insert( set< uint >::value_type( pair_prev->label ) );
				// get the correct block start
				if( pair_vec_iter->mem->Start( seqI ) < 0 ){
					block_start = pair_vec_iter->label;
					scanLabels( no_match_labels, block_start, true );
				}else if( pair_vec_iter->mem->Start( seqI ) > 0 ){
					block_start = pair_vec_iter->label;
					scanLabels( no_match_labels, block_start, false );
				}
				
				continue;
			}

			if( pair_prev->mem->Start( seqI ) < 0 ){
				// this block would break at its start
				break_label = block_start;
			}else{
				// this block would break at its end
				break_label = pair_prev->label;
				scanLabels( no_match_labels, break_label, true );
			}
			if( pair_vec_iter->mem->Start( seqI ) < 0 ){
				// scan forward to the beginning of new block
				scan_label = pair_vec_iter->label;
				scanLabels( no_match_labels, scan_label, true );
			}else{
				// scan back to the beginning of new block
				scan_label = pair_vec_iter->label;
				scanLabels( no_match_labels, scan_label, false );
			}

			if( pair_vec_iter->mem->Start( seqI ) < 0 &&
				pair_prev->mem->Start( seqI ) < 0 ){
				if( scan_label + 1 == pair_prev->label )
					continue;
				if( debugging ){
					map< uint, Match* >::const_iterator debug_iter = debug_label_map.find( pair_vec_iter->label );
					while( debug_iter->first <= pair_prev->label ){
						cout << debug_iter->first << '\t' << *(debug_iter->second) << endl;
						++debug_iter;
					}
				}
			}else
			if( pair_vec_iter->mem->Start( seqI ) > 0 &&
				pair_prev->mem->Start( seqI ) > 0 ){
				
				if( scan_label - 1 == pair_prev->label )
					continue;
				if( debugging ){
					map< uint, Match* >::const_iterator debug_iter = debug_label_map.find( pair_prev->label );
					while( debug_iter->first <= pair_vec_iter->label ){
						cout << debug_iter->first << '\t' << *(debug_iter->second) << endl;
						++debug_iter;
					}
				}
			}
			// check if the missing matches are in the set of non-matches

			// since it didn't meet any of the above
			// criteria it's a breakpoint.  insert the label of the end of the current block
			// note that if it's a reverse complement block, the end label is really the start label
			breakpoints.insert( break_label );
			block_start = scan_label;
		}

		// insert the correct block ending
		if( pair_vec_iter != pair_vec.begin() ){
			pair_vec_iter--;
			
			if( pair_vec_iter->mem->Start( seqI ) < 0 ){
				break_label = block_start;
			}else{
				break_label = pair_vec_iter->label;
				scanLabels( no_match_labels, break_label, true );
			}
			breakpoints.insert( break_label );
		}
	}
}

/** Set output parameters for permutation matrices */
void Aligner::SetPermutationOutput( std::string& permutation_filename, int64 permutation_weight )
{
	this->permutation_filename = permutation_filename;
	this->permutation_weight = permutation_weight;
}


void GetLCBCoverage( MatchList& lcb, uint64& coverage ){
	vector< Match* >::iterator match_iter = lcb.begin();
	coverage = 0;
	bool debug = true;
	for( ; match_iter != lcb.end(); ++match_iter ){
		coverage += (*match_iter)->Length() * (*match_iter)->Multiplicity();

		// if we have sequence information then
		// subtract the coverage for any position that contains an N
		if( lcb.seq_table.size() > 0 )
		{
			for( uint seqI = 0; seqI < (*match_iter)->SeqCount(); ++seqI )
			{
				gnSeqI lend = absolut((*match_iter)->Start(seqI));
				gnSeqI length = (*match_iter)->Length();
				if( lend == 0 )
					continue;
				string match_seq = lcb.seq_table[seqI]->ToString(length, lend);
				for( size_t s = 0; s < match_seq.size(); ++s )
					if( match_seq[s] == 'n' || match_seq[s] == 'N' )
						if( (*match_iter)->Start(seqI) > 0 )
							coverage--;
			}
		}
	}
}


void computeLCBAdjacencies_v2( vector<MatchList>& lcb_list, vector< int64 >& weights, vector< LCB >& adjacencies ){
	IntervalList iv_list;
	for( uint lcbI = 0; lcbI < lcb_list.size(); ++lcbI ){
		vector<AbstractMatch*> asdf;
		asdf.push_back( lcb_list[ lcbI ].front() );
		if( lcb_list[lcbI].size() > 1 )
			asdf.push_back( lcb_list[ lcbI ].back() );
		Interval iv( asdf.begin(), asdf.end() );
		iv_list.push_back( iv );
	}
	computeLCBAdjacencies_v2( iv_list, weights, adjacencies );
}

const uint NO_ADJACENCY = (std::numeric_limits<uint>::max)();

/**
 *  Redesign to be more intuitive.  left_adjacency is always left, regardless of LCB orientation
 */
void computeLCBAdjacencies_v2( IntervalList& iv_list, vector< int64 >& weights, vector< LCB >& adjacencies ){
	adjacencies.clear(); // start with no LCB adjacencies
	if( iv_list.size() == 0 )
		return;	// there aren't any LCBs so there aren't any adjacencies!

	uint seq_count = iv_list[0].SeqCount();
	uint seqI;
	uint lcbI;
	adjacencies.resize(iv_list.size());
	for( lcbI = 0; lcbI < iv_list.size(); ++lcbI ){
		LCB& lcb = adjacencies[lcbI];
		lcb.left_end.resize(seq_count);
		lcb.right_end.resize(seq_count);
		lcb.left_adjacency.resize(seq_count);
		lcb.right_adjacency.resize(seq_count);
		for( seqI = 0; seqI < seq_count; seqI++ ){
			// support "ragged edges" on the ends of LCBs
			int64 leftI = iv_list[lcbI].LeftEnd(seqI);
			int64 rightI = NO_MATCH;
			if( leftI != NO_MATCH )
			{
				leftI = iv_list[lcbI].Orientation(seqI) == AbstractMatch::forward ? leftI : -leftI;
				rightI = iv_list[lcbI].RightEnd(seqI)+1;
				rightI = iv_list[lcbI].Orientation(seqI) == AbstractMatch::forward ? rightI : -rightI;
			}

			lcb.left_end[seqI] = leftI;
			lcb.right_end[seqI] = rightI;
			lcb.left_adjacency[seqI] = NO_ADJACENCY;
			lcb.right_adjacency[seqI] = NO_ADJACENCY;
		}
		lcb.lcb_id = lcbI;
		lcb.weight = weights[ lcbI ];
		lcb.to_be_deleted = false;
	}

	for( seqI = 0; seqI < seq_count; seqI++ ){
		LCBLeftComparator llc( seqI );
		sort( adjacencies.begin(), adjacencies.end(), llc );
		for( lcbI = 1; lcbI + 1 < iv_list.size(); lcbI++ ){
			adjacencies[ lcbI ].left_adjacency[ seqI ] = adjacencies[ lcbI - 1 ].lcb_id;
			adjacencies[ lcbI ].right_adjacency[ seqI ] = adjacencies[ lcbI + 1 ].lcb_id;
		}
		if( lcbI == iv_list.size() )
			lcbI--;	// need to decrement when there is only a single LCB

		// set first and last lcb adjacencies to -1
		adjacencies[ 0 ].left_adjacency[ seqI ] = NO_ADJACENCY;
		adjacencies[ lcbI ].right_adjacency[ seqI ] = NO_ADJACENCY;
		if( lcbI > 0 ){
			adjacencies[ 0 ].right_adjacency[ seqI ] = adjacencies[ 1 ].lcb_id;
			adjacencies[ lcbI ].left_adjacency[ seqI ] = adjacencies[ lcbI - 1 ].lcb_id;
		}
	}
	LCBIDComparator lic;
	sort( adjacencies.begin(), adjacencies.end(), lic );
	
}


void scanLeft( int& left_recurseI, vector< LCB >& adjacencies, int min_weight, int seqI ){
	while( left_recurseI != -1 && adjacencies[ left_recurseI ].weight < min_weight )
		left_recurseI = adjacencies[ left_recurseI ].left_adjacency[ seqI ];
}
void scanRight( int& right_recurseI, vector< LCB >& adjacencies, int min_weight, int seqI ){
	while( right_recurseI != -1 && adjacencies[ right_recurseI ].weight < min_weight )
		right_recurseI = adjacencies[ right_recurseI ].right_adjacency[ seqI ];
}



/** iv_regions -- lists of intervening regions between LCBs in each sequence
  * start positions organized as iv_regions[ seqI ][ lcbI * 2 ]
  * end positions organized as iv_regions[ seqI ][ lcbI * 2 + 1 ] 
 */
void CreateGapSearchList( vector< LCB >& adjacencies, const vector< gnSequence* >& seq_table, vector< vector< int64 > >& iv_regions, boolean entire_genome ) 
{
	iv_regions.clear();
	if( adjacencies.size() == 0 )
		return;		// there aren't any intervening LCB regions!
	if( adjacencies.size() == 1 && !entire_genome )
		return; 	// there aren't any interveniing LCB regions in the local area
	boolean debug_lcb_extension = false;	/**< enables debugging output */
	const uint seq_count = seq_table.size();

	uint seqI = 0;
	int lcbI = 0;
	iv_regions = vector< vector< int64 > >( seq_count );

	// extract a gnSequence containing only the intervening regions
	for( seqI = 0; seqI < seq_count; seqI++ ){

		// find the first LCB in this sequence
		for( lcbI = 0; lcbI < adjacencies.size(); lcbI++ ){
			if( adjacencies[ lcbI ].left_adjacency[ seqI ] == -1 )
				break;
		}
		// start concatenating the intervening regions
		// scan right
		int right_recurseI = lcbI;
		lcbI = -1;
		if( !entire_genome && right_recurseI != -1 ){
			lcbI = right_recurseI;
			right_recurseI = adjacencies[ lcbI ].right_adjacency[ seqI ];
		}
		gnSeqI seq_len = 0;
		while( (lcbI != -1 || right_recurseI != -1 ) && right_recurseI < (int)adjacencies.size() ){
			int64 l_end = lcbI == -1 ? 1 : adjacencies[ lcbI ].right_end[ seqI ];
			int64 r_end = right_recurseI == -1 ? seq_table[ seqI ]->length() : adjacencies[ right_recurseI ].left_end[ seqI ];

			// break out if outside the last LCB and not searching the entire genome
			if( !entire_genome && right_recurseI == -1 )
				break;

			l_end = absolut( l_end );
			r_end = absolut( r_end );
			
			if( l_end > r_end && !( r_end + 1 == l_end && right_recurseI == -1 ) ){
				std::cerr << "Overlapping LCBs.  lcbI " << lcbI << " right_recurseI " << right_recurseI << endl;
				std::cerr << "lend: " << l_end << " rend: " << r_end << endl;
				l_end = r_end;
				
			}
			
			lcbI = right_recurseI;
			if( right_recurseI != -1 )
				right_recurseI = adjacencies[ right_recurseI ].right_adjacency[ seqI ];
			if( r_end + 1 == l_end && right_recurseI == -1 )
				continue;	// we're at the right end and there's nothing to add
			seq_len += r_end - l_end;
			iv_regions[ seqI ].push_back( l_end );
			iv_regions[ seqI ].push_back( r_end );
		}
		if( debug_lcb_extension )
			std::cerr << "seqI " << seqI << " seq_len: " << seq_len << endl;
	}

}

void SearchLCBGaps( MatchList& new_matches, const std::vector< std::vector< int64 > >& iv_regions, MaskedMemHash& nway_mh ) {
	if( iv_regions.size() == 0 )
		return;		// there aren't any intervening LCB regions!
	size_t sI = 0;
	for( ; sI < iv_regions.size(); sI++ )
		if( iv_regions[sI].size() > 0 )
			break;
	if( sI == iv_regions.size() )
		return;		// there aren't any intervening LCB regions!

	boolean debug_lcb_extension = false;	/**< enables debugging output */

	const uint seq_count = new_matches.seq_table.size();
	uint seqI = 0;
	int lcbI = 0;
	MatchList gap_list;
	gap_list.seq_table = vector< gnSequence* >( seq_count );	/**< intervening regions of sequences */
	gap_list.sml_table = vector< SortedMerList* >( seq_count );

	// extract a gnSequence containing only the intervening regions
	for( seqI = 0; seqI < seq_count; seqI++ ){
		gap_list.seq_table[ seqI ] = new gnSequence();
		gap_list.sml_table[ seqI ] = new DNAMemorySML();
		gnSeqI seq_len = 0;
		for( size_t ivI = 0; ivI < iv_regions[seqI].size(); ivI += 2 )
		{
			int64 l_end = iv_regions[seqI][ivI];
			int64 r_end = iv_regions[seqI][ivI+1];
			try{
			if( debug_lcb_extension )
				cerr << "Adding " << seqI << "\t" << l_end << "\t" << r_end << "\t(" << r_end - l_end << " bp)" << endl;
			gap_list.seq_table[ seqI ]->append( new_matches.seq_table[ seqI ]->ToString(r_end - l_end, l_end ) );
//			gap_list.seq_table[ seqI ]->append( new_matches.seq_table[ seqI ]->subseq( l_end, r_end - l_end ) );
			}catch(...){
				cout << "";
			}
			seq_len += r_end - l_end;
		}
		if( debug_lcb_extension )
			cerr << "seqI " << seqI << " seq_len: " << seq_len << endl;
	}
	//
	// search for MUMs in the intervening sequence regions
	//

	// calculate potential mer sizes for searches
	gnSeqI total_iv_length = 0;
	for( seqI = 0; seqI < seq_count; seqI++ ){
		total_iv_length += gap_list.seq_table[ seqI ]->length();
/*		cerr << "seqI: " << seqI << " length: " << gap_list.seq_table[ seqI ]->length();
		cerr << "\n";
*/
	}
	total_iv_length /= seq_count;

	uint search_mer_size = getDefaultSeedWeight( total_iv_length );
	if( search_mer_size < MIN_DNA_SEED_WEIGHT )
		return;		// The seed size is too small to be significant
	uint64 default_seed = getSeed( search_mer_size );
	
	//	Create sorted mer lists for the intervening gap region
	vector< boost::filesystem::path > delete_files;
	boolean create_succeeded = true;
	for( seqI = 0; seqI < seq_count; seqI++ ){
		gap_list.sml_table[ seqI ]->Clear();
		try{
			if( debug_lcb_extension )
				cerr << "Creating memory SML for seqI " << seqI << endl;
			gap_list.sml_table[ seqI ]->Create( *(gap_list.seq_table[ seqI ]), default_seed );
		}catch(...){
			create_succeeded = false;
			break;
		}
	}
	if( !create_succeeded ){	
		// free memory consumed by any SMLs	
		for( seqI = 0; seqI < seq_count; seqI++ ){
			gap_list.sml_table[ seqI ]->Clear();
			delete gap_list.sml_table[ seqI ];
		}

		for( seqI = 0; seqI < seq_count; seqI++ ){
			cerr << "Creating dmSML for seqI " << seqI << endl;
			// presumably we ran out of memory and couldn't use a MemorySML.	
			// try using a FileSML with external sort
			string concat_file = CreateTempFileName("seqconcat");

			concat_file += ".raw";	// need .raw extension to tell stupid libGenome it's a raw file
			gnRAWSource::Write( *(gap_list.seq_table[ seqI ]), concat_file.c_str() );
			delete_files.push_back( concat_file );
			delete gap_list.seq_table[ seqI ];	// make sure memory gets freed!
			cerr << "Wrote raw sequence for seqI " << seqI << endl;
			gap_list.seq_table[ seqI ] = new gnSequence();
			gap_list.seq_table[ seqI ]->LoadSource( concat_file.c_str() );
			cerr << "Loaded sequence " << seqI << gap_list.seq_table[ seqI ]->length() << "b.p.\n";
			string sml_file = CreateTempFileName("dmsml");	
			DNAFileSML* sml = new DNAFileSML( sml_file.c_str() );	
			gap_list.sml_table[ seqI ] = sml;	
			sml->dmCreate( *(gap_list.seq_table[ seqI ]), default_seed );	
			delete_files.push_back( sml_file );
			delete_files.push_back( sml_file + ".coords" );
		}
	}

	//	Find all exact matches in the gap region
	nway_mh.Clear();
	nway_mh.FindMatches( gap_list );
	gap_list.MultiplicityFilter( seq_count );
//	nway_mh.GetMatchList( gap_list );

	// free memory used by SMLs!
	for( seqI = 0; seqI < seq_count; seqI++ ){
		gap_list.sml_table[ seqI ]->Clear();
		delete gap_list.sml_table[ seqI ];
	}
	
	if( debug_lcb_extension ){
		ofstream debug_extension_out( "new_extension_matches.txt" );
		WriteList( gap_list, debug_extension_out );
		debug_extension_out.close();
	}

	//	
	// If an N mask was used, transpose MUMs back into the previous	
	// sequence coordinates	
	//	
	if( !create_succeeded ){
		for( seqI = 0; seqI < seq_count; seqI++ )	
			transposeMatches( gap_list, seqI, ((FileSML*)gap_list.sml_table[ seqI ])->getUsedCoordinates() );
	}
	//
	// Transpose MUMs back into their original sequence coordinates
	//
	for( seqI = 0; seqI < seq_count; seqI++ )
		transposeMatches( gap_list, seqI, iv_regions[ seqI ] );

	EliminateOverlaps( gap_list );
	gap_list.MultiplicityFilter( seq_count );
	// filter out matches that are too short
	gap_list.LengthFilter( MIN_ANCHOR_LENGTH );

	// free memory used by sequences!
	for( seqI = 0; seqI < seq_count; seqI++ )
		delete gap_list.seq_table[ seqI ];

	for( int delI = 0; delI < delete_files.size(); delI++ )	
		boost::filesystem::remove( delete_files[delI] );

	new_matches.insert( new_matches.end(), gap_list.begin(), gap_list.end() );
}



class MatchLeftEndComparator {
public:
	MatchLeftEndComparator( unsigned seq = 0 ){
		m_seq = seq;
	}
	MatchLeftEndComparator( MatchLeftEndComparator& msc ){
		m_seq = msc.m_seq;
	}
	// TODO??  make this do a wraparound comparison if all is equal?
	boolean operator()(const AbstractMatch* a, const AbstractMatch* b) const{
		int32 start_diff = max( a->FirstStart(), m_seq ) - max( b->FirstStart(), m_seq );
		if(start_diff == 0){
			uint32 m_count = a->SeqCount();
			m_count = m_count <= b->SeqCount() ? m_count : b->SeqCount();
			for(uint32 seqI = m_seq; seqI < m_count; seqI++){
				int64 a_start = absolut( a->Start( seqI ) ), b_start = absolut( b->Start( seqI ) );
				int64 diff = a_start - b_start;
				if(a_start == (int64)NO_MATCH || b_start == (int64)NO_MATCH)
					continue;
				else if(diff == 0)
					continue;
				else
					return diff < 0;
			}
		}
		return start_diff < 0;
	}
private:
	unsigned m_seq;
};

/**
 * Transposes the coordinates of matches in mlist to correspond to the original
 * set of source sequence regions described by seq_regions, splitting matches if
 * necessary.
 */
void transposeMatches( MatchList& mlist, uint seqI, const vector< int64 >& seq_regions ){
	if( seq_regions.size() < 2 )	
		return; // no work to be done here...

	uint matchI = 0;
	MatchLeftEndComparator msc( seqI );
	sort( mlist.begin(), mlist.end(), msc );
	uint regionI = 0;
	gnSeqI region_sum = seq_regions[ 1 ] - seq_regions[ 0 ];
	gnSeqI region_start_sum = 0;
	MatchList new_matches;

	for( ; matchI < mlist.size(); matchI++ ){
		// find the translated start coordinate for this match
		int64 trans_start = mlist[ matchI ]->Start( seqI );
		int64 iv_orig_start = trans_start;
		if( trans_start == 0 )
			continue;
		while( region_sum < absolut( trans_start ) && regionI + 2 < seq_regions.size() ){
			regionI += 2;
			region_start_sum = region_sum;
			region_sum += seq_regions[ regionI + 1 ] - seq_regions[ regionI ];
		}

		if( trans_start < 0 )
			trans_start = -seq_regions[ regionI ] - ( -trans_start - region_start_sum ) + 1;
		else if( trans_start > 0 )
			trans_start = seq_regions[ regionI ] + ( trans_start - region_start_sum ) - 1;

		int64 trans_end = mlist[ matchI ]->Start( seqI );
		trans_end += trans_end > 0 ? mlist[ matchI ]->Length() - 1: -(int64)(mlist[ matchI ]->Length()) + 1;
		
		mlist[ matchI ]->SetStart( seqI, trans_start );
		
		// this bad boy may need to be split
		gnSeqI end_region_sum = region_sum;
		gnSeqI end_prev_sum = region_start_sum;
		uint end_regionI = regionI;
		Match* cur_match = mlist[ matchI ];
		while( end_region_sum < absolut( trans_end ) && end_regionI + 2 < seq_regions.size() ){
			end_regionI += 2;

			Match* left_match = new Match( *cur_match );
			// clip off the part going to the other match
			if( left_match->Start( seqI ) < 0 ){
				cur_match->CropStart( absolut( iv_orig_start ) + left_match->Length() - end_region_sum - 1);
				left_match->CropEnd( cur_match->Length() );
			}else{
				cur_match->CropEnd( absolut( iv_orig_start ) + left_match->Length() - end_region_sum - 1);
				left_match->CropStart( cur_match->Length() );
			}

			iv_orig_start += iv_orig_start > 0 ? cur_match->Length(): -(int64)cur_match->Length();

			if( trans_start < 0 )
				trans_start = -seq_regions[ end_regionI ] - ( -iv_orig_start - end_region_sum ) + 1;
			else if( trans_start > 0 )
				trans_start = seq_regions[ end_regionI ] + ( iv_orig_start - end_region_sum ) - 1;
			
			left_match->SetStart( seqI, trans_start );

			cur_match = left_match;
			new_matches.push_back( left_match );

			end_prev_sum = end_region_sum;
			end_region_sum += seq_regions[ end_regionI + 1 ] - seq_regions[ end_regionI ];

		}
//		if( end_region_sum == absolut( trans_end ) )
//			cerr << "Beware of a possible bug in transposeMatches()\n";
	}
	
	// voila... coordinates are translated
	mlist.insert( mlist.end(), new_matches.begin(), new_matches.end() );
}

void ComputeLCBs( MatchList& meml, set<uint>& breakpoints, vector<MatchList>& lcb_list, vector<int64>& weights ){

	// there must be at least one end of a block defined
	if( breakpoints.size() < 1 )
		return;
		
	lcb_list.clear();
	weights.clear();
	
	// organize the LCBs into different MatchList instances

	set<uint>::iterator break_iter = breakpoints.begin();
	uint prev_break = 0;	// prev_break is the first match in the current block
	MatchList lcb = meml;
	for( ; break_iter != breakpoints.end(); break_iter++ ){
		lcb.clear();
		lcb.insert( lcb.begin(), meml.begin() + prev_break, meml.begin() + *break_iter + 1 );
		prev_break = *break_iter + 1;
		
		// code to filter LCBs based on their coverage
		uint64 coverage;
		GetLCBCoverage( lcb, coverage );
		weights.push_back( coverage );

		// add the new MatchList to the set if it made the cut
		lcb_list.push_back( lcb );
	}
}

void Aligner::Recursion( MatchList& r_list, Match* r_begin, Match* r_end, boolean nway_only ){
	try{
	gnSeqI gap_size = 0;
	uint seqI = 0;
//	gnSeqI min_gap_size = 0;
	boolean create_ok = true;
	// create gnSequences for each intervening region
	// create a MatchList for the intervening region
	MatchList gap_list;
	
	gap_list.seq_table.reserve( seq_count );
	gap_list.sml_table.reserve( seq_count );
	vector< int64 > starts;
	uint below_cutoff_count = 0;
// 
//	Get the sequence in the intervening gaps between these two matches
//
	for( seqI = 0; seqI < seq_count; seqI++ ){
		int64 gap_end = 0;
		int64 gap_start = 0;
		getInterveningCoordinates( r_list.seq_table, r_begin, r_end, seqI, gap_start, gap_end );
		if( (r_end && r_end->Start( seqI ) == NO_MATCH) ||
			(r_begin && r_begin->Start( seqI ) == NO_MATCH )){
			below_cutoff_count++;
			cerr << "It's screwed up\n";
			gap_list.seq_table.push_back( new gnSequence() );
			gap_list.sml_table.push_back( new DNAMemorySML() );
			continue;
		}
		if( gap_end < 0 && gap_start > 0 ){
			create_ok = false;
			cerr << "It's screwed up 2\n";
			break; // bail out on directional inconsistency
		}else if( gap_end < 0 && gap_start > 0 ){
			cerr << "It's screwed up 3\n";
			create_ok = false;
			break;	// bail out on directional inconsistency
		}
		int64 diff = gap_end - gap_start;
		diff = 0 < diff ? diff : 0;
		gap_size = diff < gap_size ? gap_size : diff;

		if( gap_start == 0 )
			cerr << "scheiss\n";

		if( debug )
			cout << r_list.seq_table[ seqI ]->length() << endl;

		if( diff < min_recursive_gap_length )
			below_cutoff_count++;
		starts.push_back( gap_start );
		gnSequence* new_seq = new gnSequence( r_list.seq_table[ seqI ]->subseq( gap_start, diff ) );
		gap_list.seq_table.push_back( new_seq );
		gap_list.sml_table.push_back( new DNAMemorySML() );
	}
	
	// only perform recursive anchoring if the gapped regions are long enough
	// otherwise just let ClustalW do the work
	if( below_cutoff_count + 1 < seq_count ){
		if( nway_only )
			nway_mh.Clear();
		else
			gap_mh.get().Clear();

		multimap< uint, uint > mer_sizes;
		// calculate potential mer sizes for searches
		for( seqI = 0; seqI < seq_count; seqI++ ){
			uint search_mer_size = getDefaultSeedWeight( gap_list.seq_table[ seqI ]->length() );
			mer_sizes.insert( multimap< uint, uint >::value_type( search_mer_size, seqI ) );
		}
		multimap< uint, uint >::iterator mer_iter = mer_sizes.end();
		mer_iter--;
		vector< uint > search_seqs;
		while( mer_iter != mer_sizes.end() ){
			uint prev_mer = mer_iter->first;
			uint new_seqs = 0;
			while( true ){
				if( mer_iter->first < MIN_DNA_SEED_WEIGHT )
					break;
				if( mer_iter->first == prev_mer || search_seqs.size() < 2 ){
					search_seqs.push_back( mer_iter->second );
					new_seqs++;
					if( mer_iter == mer_sizes.begin() ){
						mer_iter = mer_sizes.end();	// signify that the scan is complete
						break;
					}
					prev_mer = mer_iter->first;
					mer_iter--;
				}else
					break;
			}

			if( search_seqs.size() < 2 )
				break;
			// look for MUMs
			
			//
			//	Create sorted mer lists for the intervening gap region
			//

			uint64 default_seed = getSeed( prev_mer );
			if( prev_mer < MIN_DNA_SEED_WEIGHT )
				break;
			for( uint seqI = 0; seqI < gap_list.seq_table.size(); seqI++ ){
				gap_list.sml_table[ seqI ]->Clear();
				gap_list.sml_table[ seqI ]->Create( *(gap_list.seq_table[ seqI ]), default_seed );
			}
			//
			//	Find all exact matches in the gap region
			//
			MatchList cur_mems = gap_list;
			cur_mems.clear();
			if( nway_only ){
				// no sense in searching for matches in subsets!!
				if( search_seqs.size() < seq_count )
					continue;
				nway_mh.ClearSequences();
				nway_mh.FindMatches( cur_mems );
			}else{
				gap_mh.get().ClearSequences();
				gap_mh.get().FindMatches( cur_mems );
			}
			for( size_t mI = 0; mI < cur_mems.size(); ++mI )
				cur_mems[mI]->Free();
			cur_mems.clear();
		}
		if( nway_only )
			nway_mh.GetMatchList( gap_list );
		else
			gap_mh.get().GetMatchList( gap_list );
		

		// delete overlaps/inclusions		
		EliminateOverlaps( gap_list );
		// mult. filter after EliminateOverlaps because e.o. may generate some subset matches
		if( nway_only )
			gap_list.MultiplicityFilter( seq_count );
		
		// for anchor accuracy, throw out any anchors that are shorter than the minimum
		// anchor length after EliminateOverlaps()
		gap_list.LengthFilter( MIN_ANCHOR_LENGTH );

	//	if( min_gap_size < search_mer_size )
	//		create_ok = false;
		if( gap_list.size() > 0 && create_ok ){

	/*		if( debug ){
				cout << "Starting mem: " << *r_begin << endl;
				cout << "Next mem: " << *r_end << endl;
				list<Match*>::iterator gappy_iter = gap_list.begin();
				while( gappy_iter != gap_list.end() ){
					cout << **gappy_iter;
					cout << endl;
					gappy_iter++;
				}
			}
	*/

			// move all the matches that were found
			vector< Match* >::iterator mum_iter = gap_list.begin();
			for( ; mum_iter != gap_list.end(); ){
				boolean add_ok = true;
				for( uint seqI = 0; seqI < (*mum_iter)->SeqCount(); ++seqI ){
					int64 gap_start;
					if( (*mum_iter)->Start( seqI ) == NO_MATCH )
						continue;
					else if( (*mum_iter)->Start( seqI ) < 0 ){
						gap_start = r_begin != NULL ? -r_begin->End( seqI ) : 0;
						if( gap_start > 0 )
	//						gap_start = -r_end->Start( seqI ) + r_end->Length() - 1;
							gap_start = r_end != NULL ? r_end->Start( seqI ) - r_end->Length() + 1 : 0;
						else if( r_begin )
							add_ok = false;
						(*mum_iter)->SetStart( seqI, (*mum_iter)->Start( seqI ) + gap_start );
					}else{
						// insert them all before mem_iter
						gap_start = r_begin != NULL ? r_begin->End( seqI ) : 0;
						if( gap_start < 0 ){
							gap_start = r_end != NULL ? r_end->Start( seqI ) - r_end->Length() + 1 : 0;
							add_ok = false;
						}
						(*mum_iter)->SetStart( seqI, (*mum_iter)->Start( seqI ) + gap_start );
					}
				}
				if( add_ok )
					r_list.push_back( *mum_iter );
				else{
					(*mum_iter)->Free();
					(*mum_iter) = NULL;
				}
				++mum_iter;
			}
	//		for( ; mum_iter != gap_list.end(); )
	//			match_allocator.Free( *mum_iter );
		}
	}
	// delete sequences and smls
	for( uint seqI = 0; seqI < gap_list.seq_table.size(); ++seqI )
		delete gap_list.seq_table[ seqI ];
	for( uint seqI = 0; seqI < gap_list.sml_table.size(); ++seqI )
		delete gap_list.sml_table[ seqI ];
		
	gap_list.seq_table.clear();
	gap_list.sml_table.clear();
	
	}catch( gnException& gne ){
		cerr << gne << endl;
	}catch( exception& e ){
		cerr << e.what() << endl;
	}catch(...){
		cerr << "When I say 'ohhh' you say 'shit'!\n";
	}
}

// compute the gapped alignments between anchors in an LCB
void AlignLCBInParallel( bool collinear_genomes, mems::GappedAligner* gal, MatchList& mlist, Interval& iv, AlnProgressTracker& apt )
{
	// check whether this function can do anything useful...
	if( !collinear_genomes && mlist.size() < 2 ){
		iv.SetMatches( mlist );
		return;
	}
	size_t galI = 0;
	vector<GappedAlignment*> gapped_alns(mlist.size()+1, NULL);
	vector<int> success(gapped_alns.size(), 0);
	gnSeqI progress_base = apt.cur_leftend;
//#pragma omp parallel for
	for( int mI = 0; mI < mlist.size()-1; mI++ )
	{
		// align the region between mI and mI+1
		GappedAlignment ga(mlist.seq_table.size(),0);
		gapped_alns[mI] = ga.Copy();

		bool align_success = gal->Align( *(gapped_alns[mI]), mlist[mI], mlist[mI+1], mlist.seq_table );
		if(align_success)
			success[mI] = 1;
		if(mI % 50 == 0 && mI > 0)
		{
			// update and print progress
			int done = 0;
			for( int i = 0; i < gapped_alns.size(); i++ )
				if(gapped_alns[i] != NULL)
					done++;
//#pragma omp critical
{
			double cur_progress = ((double)(progress_base+done) / (double)apt.total_len)*100.0;
			printProgress((uint)apt.prev_progress, (uint)cur_progress, cout);
			apt.prev_progress = cur_progress;
}
		}
	}
	apt.cur_leftend += mlist.size()-1;

	// merge the alignments and anchors back together
	vector<AbstractMatch*> merged(mlist.size()*2 + 1);
	size_t mlistI = 0;
	size_t gappedI = 0;
	bool turn = true;
	size_t mJ = 0;

	// check if genomes are collinear and get the start and end alignments if necessary
	if(collinear_genomes)
	{
		GappedAlignment ga_tmp(mlist.seq_table.size(),0);
		GappedAlignment* ga = ga_tmp.Copy();
		bool align_success = gal->Align( *ga, NULL, mlist[0], mlist.seq_table );
		if(align_success)
			merged[mJ++] = ga;
		gapped_alns[mlist.size()] = ga_tmp.Copy();
		align_success = gal->Align( *(gapped_alns[mlist.size()]), mlist.back(), NULL, mlist.seq_table );
		if(align_success)
			success[mlist.size()] = 1;
	}
	for( ; mJ < merged.size() && mlistI < mlist.size();  )
	{
		if(turn)
			merged[mJ++] = mlist[mlistI++];
		else if(success[gappedI])
			merged[mJ++] = gapped_alns[gappedI++];
		else
			gappedI++;
		turn = !turn;
	}
	// add the last alignment
	if( success[mlist.size()]==1 )
		merged[mJ++] = gapped_alns.back();
	merged.resize(mJ);

	iv.SetMatches(merged);
}

// compute the gapped alignments between anchors in an LCB
void Aligner::AlignLCB( MatchList& mlist, Interval& iv ){
	// check whether this function can do anything useful...
	if( !collinear_genomes && mlist.size() < 2 ){
		iv.SetMatches( mlist );
		return;
	}

	vector< AbstractMatch* > iv_matches;
	boolean debug_recurse = false;
	int64 config_value = 138500;
	int print_interval = 50;
	try{
	list< Match* > match_list;
	match_list.insert( match_list.end(), mlist.begin(), mlist.end() );
	mlist.clear();
	MatchList r_list = mlist;

	list< Match* >::iterator recurse_iter = match_list.begin();
	list< Match* >::iterator recurse_prev = match_list.begin();
	// scan ahead to the first n-way matches
	while( recurse_prev != match_list.end() && (*recurse_prev)->Multiplicity() != seq_count )
		++recurse_prev;

	recurse_iter = recurse_prev;
	if( !collinear_genomes ){
		if( recurse_iter != match_list.end() )
			++recurse_iter;
		while( recurse_iter != match_list.end() && (*recurse_iter)->Multiplicity() != seq_count )
			++recurse_iter;
	}else
		cout << "Assuming collinear genomes...\n";
	
	uint memI = 0;
	uint matchI = 0;
	while( true ){
		if( memI >= print_interval && memI % print_interval == 0 || debug)
			cout << "Number: " << memI << " match " << **recurse_prev << endl;
		++memI;
		if( debug_recurse ){
			cout << "Recursing on " << endl;
			if( recurse_prev != match_list.end() )
				cout << **recurse_prev << " and " << endl;
			if( recurse_iter != match_list.end() )
				cout << **recurse_iter << endl;
		}
		
		if( recurse_prev != match_list.end() && (*recurse_prev)->Start( 0 ) == config_value )
			cout << "";
		
		// recurse on a pair of matches! 
		// this function should locate all matches between the two iterators
		// and add them to r_list		
		r_list.clear();
		GappedAlignment* cr = NULL;
		boolean align_success = false;
		
		Match* r_lend = NULL;
		Match* r_rend = NULL;
		if( recurse_iter != recurse_prev )
			r_lend = *recurse_prev;
		if( recurse_iter != match_list.end() )
			r_rend = *recurse_iter;

		// attempt a clustalW alignment
		cr = new GappedAlignment();
		align_success = gal->Align( *cr, r_lend, r_rend, r_list.seq_table );

		// add the gapped alignment to the Interval
		if( r_lend != NULL )
			iv_matches.push_back( r_lend );
		if( align_success )
			iv_matches.push_back( cr );

		// scan ahead to the next pair of n-way matches
		recurse_prev = recurse_iter;
		if( recurse_iter != match_list.end() )
			++recurse_iter;
		while( recurse_iter != match_list.end() && (*recurse_iter)->Multiplicity() != seq_count )
			++recurse_iter;

		if( ( recurse_iter == match_list.end() && !collinear_genomes ) ||
				( recurse_prev == match_list.end() && collinear_genomes ) )
				break;
	}
	// get the last little bit at the end of the LCB.
	list< Match* >::iterator iter = recurse_prev;
	for( ; iter != recurse_iter; ++iter )
		iv_matches.push_back(*iter);

	mlist.insert( mlist.end(), match_list.begin(), match_list.end() );
	iv.SetMatches(iv_matches); 

	}catch( gnException& gne ){
		cerr << gne << endl;
	}catch(exception& e){
		cerr << e.what();
	}catch(...){
		cerr << "matrix exception?\n";
	}
}

// just search each intervening region once for matches, no gapped alignment...
void Aligner::SearchWithinLCB( MatchList& mlist, std::vector< search_cache_t >& new_cache, bool leftmost, bool rightmost){
	// check whether this function can do anything useful...
	if( !(leftmost || rightmost) && mlist.size() < 2 )
		return;

	boolean debug_recurse = false;
	int64 config_value = 138500;
	int print_interval = 50;

	try{
	list< Match* > match_list;
	match_list.insert( match_list.end(), mlist.begin(), mlist.end() );
	mlist.clear();
	MatchList r_list = mlist;

	list< Match* >::iterator recurse_iter = match_list.begin();
	list< Match* >::iterator recurse_prev = match_list.begin();
	if( !leftmost && recurse_iter != match_list.end() )
		++recurse_iter;
	
	uint memI = 0;
	uint matchI = 0;
	while( recurse_prev != match_list.end() ){
		if( memI >= print_interval && memI % print_interval == 0 || debug)
			cout << "Number: " << memI << " match " << **recurse_prev << endl;
		++memI;
		if( debug_recurse ){
			cout << "Recursing on " << endl;
			if( recurse_prev != match_list.end() )
				cout << **recurse_prev << " and " << endl;
			if( recurse_iter != match_list.end() )
				cout << **recurse_iter << endl;
		}
		
		
		// recurse on a pair of matches! 
		// this function should locate all matches between the two iterators
		// and add them to r_list		
		r_list.clear();
		Match* r_left = NULL;
		Match* r_right = NULL;
		if( recurse_iter == match_list.begin() && leftmost ){
			r_left = NULL;
			r_right = *recurse_iter;
		}else if( recurse_iter == match_list.end() && rightmost ){
			r_left = *recurse_prev;
			r_right = NULL;
		}else{
			r_left = *recurse_prev;
			r_right = *recurse_iter;
		}
		// check the cache to see whether this search has already been done!

		search_cache_t cacheval = make_pair( r_left, r_right );
		if( cacheval.first != NULL )
			cacheval.first = cacheval.first->Copy();
		if( cacheval.second != NULL )
			cacheval.second = cacheval.second->Copy();
		std::vector< search_cache_t >::iterator cache_entry = std::upper_bound( search_cache.begin(), search_cache.end(), cacheval, cache_comparator );
		if( cache_entry == search_cache.end() || 
			(cache_comparator( cacheval, *cache_entry ) || cache_comparator( *cache_entry, cacheval )) )
		{
			// search this region
			Recursion( r_list, r_left, r_right, true );
		}
		new_cache.push_back( cacheval );

		if( debug_recurse ){
			vector< Match* >::iterator r_iter = r_list.begin();
			cout << "Found matches " << endl;
			for(; r_iter != r_list.end(); ++r_iter )
				cout << **r_iter << endl;
		}

		// insert any n-way matches into the match list
		for( matchI = 0; matchI < r_list.size(); ++matchI ){
			if( r_list[ matchI ]->Multiplicity() == seq_count ){
				match_list.insert( recurse_iter, r_list[ matchI ] );
			}else
			{
				r_list[matchI]->Free();
				r_list[matchI] = NULL;
			}
		}

		// move ahead to the next pair of n-way matches
		recurse_prev = recurse_iter;
		if( recurse_iter != match_list.end() )
			++recurse_iter;
		
		// break early if we aren't assuming genome collinearity
		if( !rightmost && recurse_iter == match_list.end() )
			break;
			
	}

	mlist.insert( mlist.begin(), match_list.begin(), match_list.end() );

	}catch( gnException& gne ){
		cerr << gne << endl;
	}catch(exception& e){
		cerr << e.what();
	}catch(...){
		cerr << "matrix exception?\n";
	}

	// Multiplicity Filter...
	mlist.MultiplicityFilter( seq_count );
	EliminateOverlaps( mlist );
	// E.O. can create some matches of lower multiplicity
	mlist.MultiplicityFilter( seq_count );
}

void Aligner::consistencyCheck( uint lcb_count, vector< LCB >& adjacencies, vector< MatchList >& lcb_list, vector< int64 >& weights ){
	vector< LCB > tmp_adj = adjacencies;
	vector< MatchList > tmp_lcbs = lcb_list;
	vector< int64 > tmp_weights = weights;
	filterMatches( tmp_adj, tmp_lcbs, tmp_weights );
	MatchList emmlist;
	for( uint lcbI = 0; lcbI < tmp_lcbs.size(); lcbI++ )
		emmlist.insert( emmlist.end(), tmp_lcbs[ lcbI ].begin(), tmp_lcbs[ lcbI ].end() );
	set< uint > breakpoints;
	AaronsLCB( emmlist, breakpoints );
	
	// do the correct number of LCBs exist?
	if( lcb_count != tmp_lcbs.size() ){
		cerr << "lcb_count: " << lcb_count << "\ttmp_lcbs.size(): " << tmp_lcbs.size() << endl;
	}
	if( lcb_count != breakpoints.size() ){
		cerr << "lcb_count: " << lcb_count << "\tbreakpoints.size(): " << breakpoints.size() << endl;
	}
	if( tmp_lcbs.size() != breakpoints.size() ){
		cerr << "tmp_lcbs.size(): " << tmp_lcbs.size() << "\tbreakpoints.size(): " << breakpoints.size() << endl;
	}
}


/**
 * Version 2 of this algorithm:
 * each time two LCBs coalesce, repeatedly search their intervening region until
 * either a single LCB exists or all LCBs meet the current minimum_weight.
 * @returns		The weight of the minimum weight LCB that remains
 */
int64 greedyBreakpointElimination( gnSeqI minimum_weight, vector< LCB >& adjacencies, vector< int64 >& weights, ostream* status_out ){
	// repeatedly remove the low weight LCBs until the minimum weight criteria is satisfied
	uint lcbI = 0;
	vector< uint > low_weight;
	bool have_weight = false;
	gnSeqI min_weight = 0;
	gnSeqI prev_min_weight = 0;
	uint min_lcb = 0;
	uint lcb_count = adjacencies.size();
	boolean debug_bp_elimination = false;
	uint current_lcbI = 0;	/**< tracks how many of the LCBs are above the min weight */

	if( adjacencies.size() == 0 )
		return 0;	// nothing can be done
	uint seq_count = adjacencies[0].left_end.size();
	
	while( min_weight < minimum_weight ){
		if( lcb_count == 1 )
			break;	// if only a single LCB remains, don't remove it

		while(true){
			have_weight = false;
			min_weight = 0;
			current_lcbI = 0;	// always scan the entire set

			// start with current_lcbI since everything up to it has already been scanned
			for( lcbI = current_lcbI; lcbI < weights.size(); lcbI++ ){
				if( adjacencies[ lcbI ].lcb_id != lcbI ){
					// this lcb has been removed or merged with another lcb
					continue;
				}
				if( weights[ lcbI ] < min_weight || !have_weight ){
					min_weight = weights[ lcbI ];
					min_lcb = lcbI;
					have_weight = true;
					if( min_weight == prev_min_weight && current_lcbI > 0 )
						break;	// we've already found a minimum
								// weight LCB, stop here to save some searching
				}
			}
			lcbI = min_lcb;
			have_weight = false;
			// if the min weight changed then scan the entire set from the beginning
			if( prev_min_weight != min_weight ){
				if( status_out != NULL )
					*status_out << "There are " << lcb_count << " LCBs with minimum weight " << min_weight << endl;

				current_lcbI = 0;
				prev_min_weight = min_weight;
				continue;
			}

			// save time by skipping LCBs that have already been scanned
			current_lcbI = min_lcb;
			break;
		}
		
//		consistencyCheck( lcb_count, adjacencies, lcb_list, weights );
		if( min_weight >= minimum_weight )
			break;

		// actually remove the LCBs now
		// (only remove a single LCB for now -- it's easier to calculate adjacencies)

		// remove this LCB
		adjacencies[ lcbI ].lcb_id = -2;
		
		// update adjacencies
		uint seqI;
		uint left_adj;
		uint right_adj;
		for( seqI = 0; seqI < seq_count; seqI++ ){
			left_adj = adjacencies[ lcbI ].left_adjacency[ seqI ];
			right_adj = adjacencies[ lcbI ].right_adjacency[ seqI ];
			if( debug_bp_elimination ){
				if( left_adj == -2 || right_adj == -2 ){
					cerr << "improper linking\n";
				}
				// for debugging, check for consistency:
				if( left_adj != -1 && adjacencies[ left_adj ].right_adjacency[ seqI ] != lcbI )
					cerr << "Mutiny on the bounty!\n";
				// for debugging, check for consistency
				if( right_adj == adjacencies.size() )
					cerr << "Horrible Error -399a\n";
				if( right_adj != -1 && adjacencies[ right_adj ].left_adjacency[ seqI ] != lcbI )
					cerr << "Mutiny on the bounty!\n";
			}
			if( left_adj != -1 )
				adjacencies[ left_adj ].right_adjacency[ seqI ] = right_adj;
			if( right_adj != -1 && right_adj != adjacencies.size() )
				adjacencies[ right_adj ].left_adjacency[ seqI ] = left_adj;
			
		}
		// just deleted an lcb, drop the lcb count
		lcb_count--;

		// check for collapse
		for( seqI = 0; seqI < seq_count; seqI++ ){
			left_adj = adjacencies[ lcbI ].left_adjacency[ seqI ];
			right_adj = adjacencies[ lcbI ].right_adjacency[ seqI ];
			if( left_adj == -1 || right_adj == -1 )
				continue;	// can't collapse with a non-existant LCB!

			if( debug_bp_elimination ){
				if( right_adj == adjacencies.size() )
					cerr << "Horrible Error -399a\n";
				// check whether this LCB has already been merged
				if( left_adj != adjacencies[ left_adj ].lcb_id ||
					right_adj != adjacencies[ right_adj ].lcb_id ){
					// because adjacency pointers are always updated to point to the 
					// representative entry of an LCB, the lcb_id and the array index
					// should always be identical
					cerr << "improper linking\n";
					continue;
				}
				if( left_adj == -2 || right_adj == -2 ){
					cerr << "improper linking\n";
				}
			}

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

			if( seqJ != seq_count )
				continue;
			

			// these two can be collapsed
			// do it.  do it now.
			adjacencies[ right_adj ].lcb_id = left_adj;
			if( adjacencies[ right_adj ].lcb_id == -1 ||
				adjacencies[ right_adj ].lcb_id == -2 )
				cerr << "Trouble in the eleventh circle\n";
			weights[ left_adj ] += weights[ right_adj ];
			// unlink right_adj from the adjacency list and
			// update left and right ends of left_adj
			for( seqJ = 0; seqJ < seq_count; seqJ++ ){
				boolean j_orientation = adjacencies[ left_adj ].left_end[ seqJ ] > 0;
				uint rr_adj = adjacencies[ right_adj ].right_adjacency[ seqJ ];
				uint rl_adj = adjacencies[ right_adj ].left_adjacency[ seqJ ];
				if( j_orientation == orientation ){
					adjacencies[ left_adj ].right_end[ seqJ ] = adjacencies[ right_adj ].right_end[ seqJ ];
					adjacencies[ left_adj ].right_adjacency[ seqJ ] = rr_adj;
					if( rr_adj == adjacencies.size() )
						cerr << "Horrible Error -399a\n";
					if( rr_adj != -1 )
						adjacencies[ rr_adj ].left_adjacency[ seqJ ] = left_adj;
				}else{
					adjacencies[ left_adj ].left_end[ seqJ ] = adjacencies[ right_adj ].left_end[ seqJ ];
					adjacencies[ left_adj ].left_adjacency[ seqJ ] = rl_adj;
					if( rl_adj == adjacencies.size() )
						cerr << "Horrible Error -399a\n";
					if( rl_adj != -1 )
						adjacencies[ rl_adj ].right_adjacency[ seqJ ] = left_adj;
				}
				// update lcbI's adjacency links to point nowhere
				if( adjacencies[ lcbI ].left_adjacency[ seqJ ] == right_adj )
					adjacencies[ lcbI ].left_adjacency[ seqJ ] = left_adj;
				if( adjacencies[ lcbI ].right_adjacency[ seqJ ] == right_adj )
					adjacencies[ lcbI ].right_adjacency[ seqJ ] = left_adj;


			}
			// just collapsed an lcb, decrement
			lcb_count--;
		}
	}
	return min_weight;
}

class LCBLeftEndComp
{
public:
	LCBLeftEndComp() : ssc(0) {};
	bool operator()( const MatchList& a, const MatchList& b )
	{
		return ssc(a.front(), b.front());
	}
protected:
	SingleStartComparator<AbstractMatch> ssc;
};

/**
 * Takes a set of filtered LCB adjacencies and an unfiltered set of matches as input
 * returns a filtered set of matches that reflects the LCBs found
 */
void filterMatches( vector< LCB >& adjacencies, vector< MatchList >& lcb_list, vector< int64 >& weights ){
	if( lcb_list.size() < 1 )
		return;
	MatchList lcb_tmp = lcb_list[ 0 ];
	lcb_tmp.clear();
	vector< MatchList > filtered_lcbs = vector< MatchList >( lcb_list.size(), lcb_tmp );
	uint lcbI;
	for( lcbI = 0; lcbI < adjacencies.size(); lcbI++ ){
		if( adjacencies[ lcbI ].lcb_id == lcbI ){
			filtered_lcbs[ lcbI ].insert( filtered_lcbs[ lcbI ].end(), lcb_list[ lcbI ].begin(), lcb_list[ lcbI ].end() );
			continue;
		}
		if( adjacencies[ lcbI ].lcb_id == -1 ){
			cerr << "weird";
			continue; 	// this one was removed
		}
		if( adjacencies[ lcbI ].lcb_id == -2 )
			continue; 	// this one was removed

		// this one points elsewhere
		// search and update the union/find structure for the target
		stack< uint > visited_lcbs;
		visited_lcbs.push( lcbI );
		uint cur_lcb = adjacencies[ lcbI ].lcb_id;
		while( adjacencies[ cur_lcb ].lcb_id != cur_lcb ){
			visited_lcbs.push( cur_lcb );
			cur_lcb = adjacencies[ cur_lcb ].lcb_id;
			if( cur_lcb == -1 || cur_lcb == -2 ){
//				cerr << "improper hoodidge\n";
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
	}


	lcb_list.clear();
	vector< int64 > new_weights;
	for( lcbI = 0; lcbI < filtered_lcbs.size(); lcbI++ ){
		if( filtered_lcbs[ lcbI ].size() > 0 ){
			lcb_list.push_back( filtered_lcbs[ lcbI ] );
			uint64 wt = 0;
			GetLCBCoverage( filtered_lcbs[lcbI], wt );
			new_weights.push_back( wt );
//			if( new_weights[ new_weights.size() - 1 ] != weights[ lcbI ] ){
//				cerr << "Error: Have you lost weight Susan? difference: " << new_weights[ new_weights.size() - 1 ] - weights[ lcbI ] << "\n";
//			}
		}
	}

	// sort the matches inside consolidated LCBs
	MatchStartComparator<AbstractMatch> msc( 0 );
	for( lcbI = 0; lcbI < lcb_list.size(); lcbI++ ){
		sort( lcb_list[ lcbI ].begin(), lcb_list[ lcbI ].end(), msc );
	}

	// sort the LCBs themselves
	LCBLeftEndComp llec;
	std::sort( lcb_list.begin(), lcb_list.end(), llec );

	// calculate the LCB adjacencies
	weights = new_weights;
	computeLCBAdjacencies_v2( lcb_list, weights, adjacencies );

}

void Aligner::WritePermutation( vector< LCB >& adjacencies, std::string out_filename )
{
	ofstream permutation_out( out_filename.c_str() );
	if( !permutation_out.is_open() )
	{
		cerr << "Error opening " << out_filename << endl;
		return;
	}
	for( int seqI = 0; seqI < seq_count; seqI++ )
	{
		// find the left-most LCB in this genome
		int left_lcb = 0;
		for( ; left_lcb < adjacencies.size(); left_lcb++ )
		{
			uint left_adj = adjacencies[left_lcb].left_adjacency[seqI];
			if( left_adj == -1 )
				break;
		}
		// write out lcb id's in order
		for( uint lcbI = left_lcb; lcbI < adjacencies.size(); )
		{
			if( lcbI != left_lcb )
				permutation_out << '\t';
			if( adjacencies[lcbI].left_end[seqI] < 0 )
				permutation_out << "-";
			permutation_out << adjacencies[lcbI].lcb_id;
			lcbI = adjacencies[lcbI].right_adjacency[seqI];
		}
		permutation_out << endl;
	}
}

void WritePermutationCoordinates( IntervalList& perm_iv_list, std::string out_filename )
{
	ofstream perm_out( out_filename.c_str() );
	if( !perm_out.is_open() )
	{
		cerr << "Error opening \"" << out_filename << "\"\n";
		return;
	}
	perm_out << "#";
	for( size_t seqI = 0; seqI < perm_iv_list.seq_table.size(); ++seqI )
	{
		if( seqI > 0 )
			perm_out << '\t';
		perm_out << "seq" << seqI << "_leftend\tseq" << seqI << "_rightend";
	}
	perm_out << endl;
	for( size_t ivI = 0; ivI < perm_iv_list.size(); ++ivI )
	{
		for( size_t seqI = 0; seqI < perm_iv_list.seq_table.size(); ++seqI )
		{
			if( seqI > 0 )
				perm_out << '\t';
			if( perm_iv_list[ivI].Orientation(seqI) == AbstractMatch::reverse )
				perm_out << '-';
			perm_out << perm_iv_list[ivI].LeftEnd(seqI) << '\t';
			if( perm_iv_list[ivI].Orientation(seqI) == AbstractMatch::reverse )
				perm_out << '-';
			perm_out << perm_iv_list[ivI].RightEnd(seqI);
		}
		perm_out << endl;
	}
}

void Aligner::RecursiveAnchorSearch( MatchList& mlist, gnSeqI minimum_weight, vector< MatchList >& LCB_list, boolean entire_genome, ostream* status_out ){

//
// Step 4) Identify regions of collinearity (LCBs) among the remaining n-way multi-MUMs
//
	uint lcbI;
	set<uint> breakpoints;
	vector< int64 > weights;
	vector< LCB > adjacencies;
	MatchList new_matches;
	new_matches.seq_table = mlist.seq_table;
	new_matches.seq_filename = mlist.seq_filename;

	if( mlist.size() == 0 )
		return;

	AaronsLCB( mlist, breakpoints );
	if( status_out )
		*status_out << "The " << mlist.size() << " matches constitute " << breakpoints.size() << " breakpoints\n";
	// organize the LCBs into different MatchList instances (inside of LCB_list)
	ComputeLCBs( mlist, breakpoints, LCB_list, weights );
	uint weightI;
	for( weightI = 0; weightI < weights.size(); weightI++ )
		if( weights[weightI] < cur_min_coverage || cur_min_coverage == -1 )
			cur_min_coverage = weights[weightI];

	computeLCBAdjacencies_v2( LCB_list, weights, adjacencies );

	int cur_extension_round = 0;
	int64 total_weight = 0;
	int64 prev_total_weight = 0;
	weightI = 0;
	vector< vector< int64 > > prev_iv_regions;
	do {

//		for( ; weightI < weights.size(); weightI++ )
//			total_weight += weights[ weightI ];

		int64 extension_weight = total_weight;
		int64 prev_extension_weight = total_weight;

		// only search outside existing LCBs on the whole-genome scale to save time
		if( entire_genome && extend_lcbs && total_weight != 0 &&
			cur_extension_round < this->max_extension_iters )
		{
			cur_extension_round++;
			if( status_out )
				*status_out << "Performing LCB extension\n";
			vector< vector< int64 > > cur_iv_regions;
			CreateGapSearchList( adjacencies, new_matches.seq_table, cur_iv_regions, entire_genome );
			// only do the search if there's something new to search
			if( prev_iv_regions != cur_iv_regions )
			{
				int local_round = 0;
				do {
					local_round++;
					// search the gaps between the LCBs to extend the ends of LCBs
					new_matches.clear();
					vector< vector< int64 > > new_iv_regions;
					CreateGapSearchList( adjacencies, new_matches.seq_table, new_iv_regions, entire_genome );
					SearchLCBGaps( new_matches, new_iv_regions, nway_mh );
					mlist.insert( mlist.end(), new_matches.begin(), new_matches.end() );
					
					AaronsLCB( mlist, breakpoints );
					ComputeLCBs( mlist, breakpoints, LCB_list, weights );
					cur_min_coverage = *(std::min_element(weights.begin(), weights.end()));
					computeLCBAdjacencies_v2( LCB_list, weights, adjacencies );

					// calculate the new total LCB weight
					prev_extension_weight = extension_weight;
					extension_weight = 0;
					for( weightI = 0; weightI < weights.size(); weightI++ )
						extension_weight += weights[ weightI ];
					if( status_out )
						*status_out << "Previous weight: " << prev_extension_weight << " new weight: " << extension_weight << endl;
					if( prev_extension_weight > extension_weight ){
						cerr << "Error! Previous weight: " << prev_extension_weight << " new weight: " << extension_weight << endl;
					}
				}while( extension_weight > prev_extension_weight && local_round < this->max_extension_iters);
			}
			swap( prev_iv_regions, cur_iv_regions );
		}
		
		// now search within LCBs
		if( currently_recursing && total_weight != 0 ){
			vector< search_cache_t > new_cache;
			for( lcbI = 0; lcbI < LCB_list.size(); lcbI++ ){
//				if( status_out )
//					*status_out << "Searching in LCB: " << lcbI << endl;
				int prev_size = LCB_list[ lcbI ].size();
				bool leftmost = true;
				for( int i = 0; leftmost && i < adjacencies[lcbI].left_adjacency.size(); i++ )
					if(adjacencies[lcbI].left_adjacency[i] != NO_ADJACENCY)
						leftmost = false;
				bool rightmost = true;
				for( int i = 0; rightmost && i < adjacencies[lcbI].right_adjacency.size(); i++ )
					if(adjacencies[lcbI].right_adjacency[i] != NO_ADJACENCY)
						rightmost = false;
				SearchWithinLCB( LCB_list[ lcbI ], new_cache, leftmost, rightmost );
//				if( status_out )
//					*status_out << "Gained " << LCB_list[ lcbI ].size() - prev_size << " matches\n";

			}

			// delete the previous search cache
			swap( search_cache, new_cache );
			for( size_t mI = 0; mI < new_cache.size(); mI++ )
			{
				if( new_cache[mI].first != NULL )
					new_cache[mI].first->Free();
				if( new_cache[mI].second != NULL )
					new_cache[mI].second->Free();
			}
			new_cache.clear();
			std::sort( search_cache.begin(), search_cache.end(), cache_comparator );
		}
		
		mlist.clear();
		for( lcbI = 0; lcbI < LCB_list.size(); lcbI++ ){
			mlist.insert( mlist.end(), LCB_list[ lcbI ].begin(), LCB_list[ lcbI ].end() );
		}

		if( currently_recursing && total_weight != 0 ){
			// remove low weight LCBs, while searching coalesced regions
			AaronsLCB( mlist, breakpoints );
			ComputeLCBs( mlist, breakpoints, LCB_list, weights );
			computeLCBAdjacencies_v2( LCB_list, weights, adjacencies );
			cur_min_coverage = *(std::min_element(weights.begin(), weights.end()));
		}

		
		// write  alist for debugging
//		ofstream debug_match_list( "debug_match_list.txt" );
//		mlist.WriteList( debug_match_list );
//		debug_match_list.close();

//
// Step 6) Use greedy breakpoint elimination to remove low-weight LCBs
//
		int64 cur_perm_weight = permutation_weight != -1 ? permutation_weight : minimum_weight;
		do{
			vector<double> m_weights(weights.size());
			for( size_t wI = 0; wI < weights.size(); wI++ )
				m_weights[wI] = (double)weights[wI];
			SimpleBreakpointScorer sbs(adjacencies, cur_perm_weight, this->collinear_genomes);
			if( status_out )
				(*status_out) << "Performing greedy breakpoint elimination (this may take some time)\n";

			greedyBreakpointElimination_v4(adjacencies, m_weights, sbs, NULL, false);
//			cur_min_coverage = greedyBreakpointElimination( cur_perm_weight, adjacencies, weights, status_out );
//			MatchList deleted_matches;
			filterMatches( adjacencies, LCB_list, weights );
			cur_min_coverage = *(std::min_element(weights.begin(), weights.end()));
			
			mlist.clear();
			for( lcbI = 0; lcbI < LCB_list.size(); lcbI++ ){
				mlist.insert( mlist.end(), LCB_list[ lcbI ].begin(), LCB_list[ lcbI ].end() );
			}
			if( status_out )
				*status_out << "Greedy breakpoint elimination leaves " << mlist.size() << " matches constituting " << LCB_list.size() << " LCBs covering at least " << cur_min_coverage << "b.p.\n";
			
			if( permutation_weight != -1 ){
				// construct a filename
				stringstream cur_perm_filename;
				cur_perm_filename << permutation_filename << "." << cur_perm_weight / seq_count;
				// output the permutation
				WritePermutation( adjacencies, cur_perm_filename.str() );

				// also write out condensed interval data for the permutation
				cur_perm_filename << ".lcbs";
				IntervalList perm_iv_list;
				perm_iv_list.seq_filename = mlist.seq_filename;
				perm_iv_list.seq_table = mlist.seq_table;
				for( int permI = 0; permI < LCB_list.size(); permI++ ){
					vector< AbstractMatch* > perm_vector;
					perm_vector.push_back( LCB_list[permI].front() );
					if( LCB_list[permI].size() > 1 )
						perm_vector.push_back( LCB_list[permI].back() );
					Interval perm_iv(perm_vector.begin(), perm_vector.end());
					perm_iv_list.push_back(perm_iv);
				}
				WritePermutationCoordinates( perm_iv_list, cur_perm_filename.str() );

				// get the current min weight
				vector< int64 >::iterator min_w = std::min_element( weights.begin(), weights.end() );
				// increment the current weight
				cur_perm_weight = *min_w + seq_count;
			}
		}while( cur_perm_weight < minimum_weight );
		// only enable recursive anchor search once we achieve
		// the desired weight threshold once -- for speed's sake
		if( recursive && entire_genome ){
			currently_recursing = true;
		}

		// calculate the new total LCB weight
		prev_total_weight = total_weight;
		total_weight = 0;
		for( weightI = 0; weightI < weights.size(); weightI++ )
			total_weight += weights[ weightI ];
		if( status_out )
			*status_out << "Previous weight: " << prev_total_weight << " new weight: " << total_weight << endl;
	// the weight can shrink--this isn't an error condition
//		if( prev_total_weight > total_weight ){
//			cerr << "Error! Previous weight: " << prev_total_weight << " new weight: " << total_weight << endl;
			// write out the lcb lists
//		}

//
// Step 7) Repeat 4, 5 and 6 until the total weight stabilizes
//
	}while( total_weight != prev_total_weight );

	// delete the search cache
	for( size_t mI = 0; mI < search_cache.size(); mI++ )
	{
		if( search_cache[mI].first != NULL )
			search_cache[mI].first->Free();
		if( search_cache[mI].second != NULL )
			search_cache[mI].second->Free();
	}
}

/**
 * Note: this algorithm differs from the one reported in the Mauve paper
 *       The modifications should make the Mauve method more sensitive
 * Given an initial set of multi-MUMs, the alignment is an x step process:
 * 1) Eliminate overlaps among the multi-MUMs
 * 2) Compute a phylogenetic guide tree using the multi-MUMs
 * 3) Remove subset multi-MUMs
 * 4) Identify regions of collinearity (LCBs) among the remaining n-way multi-MUMs
 * 5) Perform recursive anchor search within and outside LCBs
 *    5a) search outside until weight stabilizes
 *    5b) search within LCBs
 * 6) Use greedy breakpoint elimination to remove low-weight LCBs
 *    6a) whenever two LCBs coalesce, search the intervening region for multi-MUMs
 * 7) Repeat 4, 5 and 6 until the total weight stabilizes
 * 8) Perform gapped alignment on each LCB
 * When limited area DP and POA are integrated, step 8 will become step 5c
 * 
 */

void Aligner::align( MatchList& mlist, IntervalList& interval_list, double LCB_minimum_density, double LCB_minimum_range, boolean recursive, boolean extend_lcbs, boolean gapped_alignment, string tree_filename ){
	seq_count = mlist.seq_table.size();
	this->LCB_minimum_density = LCB_minimum_density;
	this->LCB_minimum_range = LCB_minimum_range;
	this->recursive = recursive;
	this->currently_recursing = false;
	this->extend_lcbs = extend_lcbs;
	this->gapped_alignment = gapped_alignment;

	// use LCB_minimum_range == -1 to indicate that all genomes are 
	// expected to be collinear
	this->collinear_genomes = LCB_minimum_range == -1;
	if( collinear_genomes )
		cout << "\nAssuming collinear genomes...\n";

	// set the nway_mh mask
	uint64 nway_mask = 1;
	nway_mask <<= seq_count;
	nway_mask--;
	nway_mh.SetMask( nway_mask );
		
	cout << "Starting with " << mlist.size() << " MUMs\n";
	
//
// Step 1) Eliminate overlaps among the multi-MUMs
//	
	// Remove linked inclusions
	EliminateOverlaps( mlist );
	cout << "Eliminating overlaps yields " << mlist.size() << " MUMs\n";

//
// Step 2) Compute a phylogenetic guide tree using the multi-MUMs
//

	bool guide_tree_loaded = false;
	MuscleInterface& mi = MuscleInterface::getMuscleInterface();	

	if( !guide_tree_loaded && (recursive || tree_filename != "") ){
		// Make a phylogenetic tree for ClustalW
		interval_list.seq_table = mlist.seq_table;
		interval_list.seq_filename = mlist.seq_filename;
		// use the identity matrix method and convert to a distance matrix
		NumericMatrix< double > distance;
		DistanceMatrix( mlist, distance );
		if( tree_filename == "" )
			tree_filename = CreateTempFileName("guide_tree");
		mi.CreateTree( distance, tree_filename );
	}

//
// Step 3) Remove subset multi-MUMs
//
	// Multiplicity Filter...
	mlist.MultiplicityFilter( seq_count );
	cout << "Multiplicity filter gives " << mlist.size() << " MUMs\n";

	if( mlist.size() == 0 )
		return;
	
//
// Steps 4 through 7 are contained in RecursiveAnchorSearch
//
	vector< MatchList > LCB_list;
	RecursiveAnchorSearch( mlist, (gnSeqI)LCB_minimum_range, LCB_list, true, &cout );


//
// Step 8) Perform gapped alignment on each LCB using the anchors
//
	if( gapped_alignment && recursive )
		cout << "\nMaking final gapped alignment...\n";
	interval_list.clear();
	AlnProgressTracker apt;
	apt.cur_leftend = 0;
	apt.prev_progress = 0;
	apt.total_len = 0;
	for( uint lcbI = 0; lcbI < LCB_list.size(); lcbI++ )
		apt.total_len += LCB_list[lcbI].size()-1;
	for( uint lcbI = 0; lcbI < LCB_list.size(); lcbI++ ){
		Interval new_iv;
		interval_list.push_back( new_iv );
		Interval& iv = interval_list.back();
		if( !gapped_alignment || !recursive ){
			iv.SetMatches( LCB_list[lcbI] );
		}else{
//			AlignLCB( LCB_list[ lcbI ], iv );
			AlignLCBInParallel( collinear_genomes || (LCB_list.size()==1), gal, LCB_list[ lcbI ], iv, apt );
		}
	}
	
	// finally add any unaligned regions to the interval list	
	if( gapped_alignment )
		addUnalignedIntervals( interval_list );
}

}	// namespace mems

