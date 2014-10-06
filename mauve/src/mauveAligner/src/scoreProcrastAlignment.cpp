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

#include "libMems/CompactGappedAlignment.h"
#include "libMems/MatchList.h"
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include "libMems/IntervalList.h"
#include "libGenome/gnFilter.h"
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options.hpp>
#include <boost/tuple/tuple.hpp>
#include <algorithm>
namespace po = boost::program_options;

using namespace std;
using namespace genome;
using namespace mems;

// basic data structures

/** store a pair of aligned positions and the characters */
typedef struct aligned_coords_s {
	int64 pos1;
	int64 pos2;
	char char1;
	char char2;
} aligned_coords_t;


class AlignedCoordSeqIComparator {
public:
	boolean operator()(const aligned_coords_t& a, const aligned_coords_t& b) const{
		if( abs(a.pos1) == abs(b.pos1) )
			return abs(a.pos2) < abs(b.pos2);
		return abs(a.pos1) < abs(b.pos1);
	}
};

void constructCoordList( uint seqI, uint seqJ, IntervalList& iv_list, vector< aligned_coords_t >& coord_list, vector< gnSequence* >& seq_table ){

	//
	// pre-allocate the vector
	//
	gnSeqI ij_vec_size = 0;
	for( int ivI = 0; ivI < iv_list.size(); ivI++ ){
		ij_vec_size += iv_list[ivI].AlignmentLength();
	}
	coord_list = vector< aligned_coords_t >( ij_vec_size );

	//
	// fill in the vector with all aligned pairs
	//
	gnSeqI vecI = 0;	// current place in vector
	for( int ivI = 0; ivI < iv_list.size(); ivI++ ){
		GappedAlignment* aln;
		aln = dynamic_cast< GappedAlignment* >( iv_list[ ivI ].GetMatches()[0] );
		if( aln == NULL ){
			throw "Error:  expecting interval to contain a single GappedAlignment";
		}
		int64 pos1 = aln->Start( seqI );
		int64 pos2 = aln->Start( seqJ );
	
		// if rev. comp then we're starting at the other (right) side
		if( pos1 < 0 )
			pos1 -= aln->Length( seqI ) - 1;
		if( pos2 < 0 )
			pos2 -= aln->Length( seqJ ) - 1;

		
		const std::vector< std::string >& align_matrix = GetAlignment( *aln, seq_table );
		for( gnSeqI colI = 0; colI < aln->Length(); colI++ ){
			aligned_coords_t act;
			act.char1 = align_matrix[ seqI ][ colI ];
			act.char2 = align_matrix[ seqJ ][ colI ];
			act.pos1 = act.char1 == '-' ? 0 : pos1;
			act.pos2 = act.char2 == '-' ? 0 : pos2;
			
			coord_list[ vecI++ ] = act;
			
			if( act.char1 != '-' )
				pos1++;
			if( act.char2 != '-' )
				pos2++;
		}
		
	}

	//
	// sort the vector on aligned position
	//
	AlignedCoordSeqIComparator acsc;
	sort( coord_list.begin(), coord_list.end(), acsc );
}


const gnFilter* comp_filter = gnFilter::DNAComplementFilter();

template< class PairType >
class PairFirstComparator
{
public:
	bool operator()( const PairType& a, const PairType& b )
	{
		return a.first < b.first;
	}
};

void compareAlignmentsAceD( IntervalList& correct, IntervalList& procrastinated, gnSequence& concat_sequence )
{
	gnSeqI sp_truepos = 0;
	gnSeqI sp_possible = 0;

	uint seqI = 0;
	uint seqJ = 0;
	// for now, use this value to create a unique identifier for the pairwise_component_hits bitset vector
//	uint MAX_MULTIPLICITY = 1000;
	uint seq_count = concat_sequence.contigListLength();

	// create a data structure that indicates the start offset in concatenated coordinates for a given sequence
	vector< gnSeqI > concat_coords(seq_count+1, 0);
	for( size_t seqI = 0; seqI < seq_count; ++seqI )
	{
		concat_coords[seqI+1] = concat_coords[seqI] + concat_sequence.contigLength(seqI);
	}

	// tuple stores pointer to interval, the component of the interval, and the interval's index in procrastinated
	// for each position of the concatenated sequence. 
	typedef std::pair< size_t, uint > iv_tracker_t;
	typedef vector< iv_tracker_t, boost::pool_allocator< iv_tracker_t > > tracker_vector_t;
	typedef vector< tracker_vector_t, boost::pool_allocator< tracker_vector_t > > coord_iv_map_vector_t;
	// create a map from sequence position to (interval,component) for the total length of the concat sequence
	// use boost pool allocators since this never needs to get freed
	coord_iv_map_vector_t* tmp = new coord_iv_map_vector_t( concat_coords.back() + 1 );	// heap allocate to avoid destruction when the stack frame is popped
	coord_iv_map_vector_t& coord_iv_map = *tmp;
	vector< size_t > coord_iv_counts( concat_coords.back() + 1, 0 );
	// first count the number of ivs that contain each position so we know how much to allocate
	for( size_t calcI = 0; calcI < procrastinated.size(); ++calcI )
	{
		Interval& iv = procrastinated[calcI];
		for( size_t seqI = 0; seqI < iv.SeqCount(); ++seqI )
		{
			const gnSeqI lend = iv.LeftEnd(seqI);
			if( lend == NO_MATCH )
				continue;	// this shouldn't happen with procrastAligner output, but let's be safe
			const gnSeqI rend = iv.RightEnd(seqI);
			for( size_t posI = lend; posI <= rend; ++posI )
				coord_iv_counts[posI]++;
		}
	}
	// now allocate space for the map
	for( size_t mapI = 0; mapI < coord_iv_map.size(); ++mapI )
		coord_iv_map[mapI].resize( coord_iv_counts[mapI] );

	std::fill( coord_iv_counts.begin(), coord_iv_counts.end(), 0 );	// recycle this storage to count the number added thus far

	// finally, populate the map
	for( size_t calcI = 0; calcI < procrastinated.size(); ++calcI )
	{
		Interval& iv = procrastinated[calcI];
		for( size_t seqI = 0; seqI < iv.SeqCount(); ++seqI )
		{
			const gnSeqI lend = iv.LeftEnd(seqI);
			if( lend == NO_MATCH )
				continue;	// this shouldn't happen with procrastAligner output, but let's be safe
			const gnSeqI rend = iv.RightEnd(seqI);
			for( size_t posI = lend; posI <= rend; ++posI )
			{
				coord_iv_map[posI][coord_iv_counts[posI]] = make_pair( calcI, seqI );
				coord_iv_counts[posI]++;
			}
		}
	}

	size_t all_component_count = 0;
	size_t all_component_pair_count = 0;
	size_t component_pair_count = 0;
	// create a vector of bitsets for each iv to store whether their components were correctly aligned
	vector< bitset_t > component_hits( procrastinated.size() );
	// Follow Aaron's lead and store pairwise component hits in bitset_t vector
	vector< bitset_t > pairwise_component_hits( procrastinated.size() );
	for( size_t ivI = 0; ivI < component_hits.size(); ++ivI )
	{
		// make sure this value is always greater than the largest max multiplicity
//		if( MAX_MULTIPLICITY < procrastinated[ivI].SeqCount())
//			MAX_MULTIPLICITY *= 10;
		// possible pairwise component combinations for this interval
//		component_pair_count = ( ( procrastinated[ivI].SeqCount() *  (procrastinated[ivI].SeqCount() - 1) ) / 2 );
		// let this be oversized for easier indexing, but correct for it when calculating the PPV below...
		component_pair_count = procrastinated[ivI].SeqCount() * procrastinated[ivI].SeqCount();
		pairwise_component_hits[ivI].resize(component_pair_count, false);
		component_hits[ivI].resize(procrastinated[ivI].SeqCount(), false);
		all_component_pair_count += ( ( procrastinated[ivI].SeqCount() *  (procrastinated[ivI].SeqCount() - 1) ) / 2 );
		all_component_count += procrastinated[ivI].SeqCount();
	}
	
	// sort each vector of iv_tracker_t by iv memory address (first element) so we can later do set intersections
	for( size_t posI = 0; posI < coord_iv_map.size(); ++posI )
		std::sort( coord_iv_map[posI].begin(), coord_iv_map[posI].end() );

	tracker_vector_t intersect_buf1( all_component_count );	// storage for set intersections
	tracker_vector_t intersect_buf2( all_component_count );	// storage for set intersections
	PairFirstComparator< iv_tracker_t > pfc;

	// now, for each pair of aligned positions in the correct alignment, determine whether they
	// lie in a procrastAligner chain
	size_t all_pair_count = (seq_count * (seq_count - 1)) / 2;
	size_t pair_count = 0;
	
	for( seqI = 0; seqI < seq_count; seqI++ )
	{
		for( seqJ = seqI+1; seqJ < seq_count; seqJ++ )
		{
			size_t prev_count = pair_count;
			pair_count++;
			if( (pair_count * 100) / all_pair_count != (prev_count * 100) / all_pair_count )
			{
				cout << (pair_count * 100) / all_pair_count << "%..";
				cout.flush();
			}
			vector< aligned_coords_t > cor;
			
			//construct the coord list just for the correct alignment
			vector< gnSequence* > seq_table( seq_count, (gnSequence*)NULL );
			constructCoordList( seqI, seqJ, correct, cor, seq_table );

			gnSeqI corI = 0;
			// skip any gaps aligned to gaps
			while( corI < cor.size() && cor[ corI ].pos1 == 0 )
				corI++;

			for( ; corI < cor.size(); corI++ )
			{
				if( cor[ corI ].pos1 != 0 && cor[ corI ].pos2 != 0)	// don't count positions aligned to gaps
					sp_possible++;
				else
					continue;

				// which positions do the correct pair have in the concatenated alignment space?
				gnSeqI trans_pos1 = genome::absolut( cor[corI].pos1 ) + concat_coords[seqI];
				gnSeqI trans_pos2 = genome::absolut( cor[corI].pos2 ) + concat_coords[seqJ];
				
				// which chain(s) do these positions fall into?
				// are any of them the same chain?
				tracker_vector_t::iterator last_int1 = std::set_intersection( 
					coord_iv_map[trans_pos1].begin(),coord_iv_map[trans_pos1].end(),
					coord_iv_map[trans_pos2].begin(),coord_iv_map[trans_pos2].end(), 
					intersect_buf1.begin(), pfc );

				if( last_int1 == intersect_buf1.begin() )
				{
					// not contained in any chain.  false negative
				}else{
					// make a list of pairs for each position
					// set_intersection always puts elements from the first set into the output buffer,
					// since the elements in the second set may have the same iv ptr but a different
					// match component, we want a list of those as well
					tracker_vector_t::iterator last_int2 = std::set_intersection( 
						coord_iv_map[trans_pos2].begin(),coord_iv_map[trans_pos2].end(), 
						coord_iv_map[trans_pos1].begin(),coord_iv_map[trans_pos1].end(),
						intersect_buf2.begin(), pfc );
					
					size_t pcount = last_int1 - intersect_buf1.begin();
					bool found = false;	// set this to true if at least one element has different match components
					for( size_t pI = 0; pI < pcount; ++pI )
					{
						// make sure they're not in the same component (probably a very rare occurrence)
						size_t component_1 = intersect_buf1[pI].second;
						size_t component_2 = intersect_buf2[pI].second;
						size_t ivI = intersect_buf1[pI].first;
						if( component_1 == component_2 )
							continue;	// no alignment here

						// make sure the relative orientations match
						bool cor_orient = (cor[corI].pos1 > 0) == (cor[corI].pos2 > 0);
						bool calc_orient = (procrastinated[ivI].Orientation(component_1) == procrastinated[ivI].Orientation(component_2));
						if( cor_orient != calc_orient )
							continue;	// calculated alignment has the wrong strand

						// make sure they're not aligned to something else...
						CompactGappedAlignment<>* cga = dynamic_cast< CompactGappedAlignment<>* >(procrastinated[ivI].GetMatches()[0]);
						size_t col_1 = cga->SeqPosToColumn(component_1, trans_pos1);
						const vector< bitset_t >& aln_mat = cga->GetAlignment();
						// if they're not aligned, make sure they're in the same gap.
						// they might get aligned later if we were to actually align the procrastAligner chains
						// instead of just finding anchors.
						if( !aln_mat[component_2].test(col_1) )
						{
							// if we encounter any columns between col_1 and col_2 that have
							// component_1 and component_2 aligned then we wouldn't ever align
							// pos_1 and pos_2 without changing the anchoring
							size_t col_2 = cga->SeqPosToColumn(component_2, trans_pos2);
							size_t col_first = col_1;
							size_t col_last = col_2;
							if( col_first < col_last )
								swap(col_first, col_last);
							size_t colI = col_first;
							for( ; colI <= col_last; ++colI )
							{
								if( aln_mat[component_1].test(colI) && aln_mat[component_2].test(colI) )
									break;
							}
							if( colI <= col_last )
								continue;	// an anchor intervenes...  bummer.
						}

						// mark these components as good
						found = true;
						component_hits[ivI].set( component_1 );
						component_hits[ivI].set( component_2 );
						
						// Always use the smallest component first
						if( component_2 < component_1 )
							swap(component_1, component_2);

						// calculate signficand for creating double
//						double significand = (double)(component_2+1)/(double)MAX_MULTIPLICITY;
						// store merged_component
//						double merged_component = (double)(component_1+1)+significand;

						// and use as unique pairwise index for each pair to take advantage
						// of bitset_t vector
						size_t sig = component_1 * cga->SeqCount() + component_2;
						pairwise_component_hits[ivI].set( sig );
					}
					if( found )
						sp_truepos++;
				}
			}
		}
	}

	cout << "\ndone!\n";
	// yaaay! we're done.  report score.
	cout << "sp_truepos " << sp_truepos << endl;
	cout << "sp_possible " << sp_possible << endl;
	cout << "SP sensitivity: " << ((double)sp_truepos) / ((double)sp_possible) << endl;
	double components_correct = 0;
	double components_possible = 0;
	for( size_t ivI = 0; ivI < component_hits.size(); ++ivI )
	{
		components_correct += component_hits[ivI].count();
		components_possible += component_hits[ivI].size();
	}
	cout << "Match component PPV: " << components_correct / components_possible << endl;

	double pairwise_components_correct = 0;
	for( size_t ivI = 0; ivI < pairwise_component_hits.size(); ++ivI )
	{
		pairwise_components_correct  += pairwise_component_hits[ivI].count();
	}
	cout << "Pairwise match component PPV: " << pairwise_components_correct / (double)all_component_pair_count << endl;
}


/**
 * program to score alignments
 * reads in a "correct" alignment and a procrastinated alignment
 * scores the procrastinated alignment based on the correct one
 */
int main( int argc, char* argv[] )
{
	
	string correct_aln_fname;
	string procrast_aln_fname;
	string sequence_file;
	
	
	if( argc < 2 ){
		cout << "scoreProcrastAlignment <correct alignment> <procrastAligner output>\n";
		return -1;
	}
	// Declare the supported options.
	
	po::variables_map vm;
	try {

        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "get help message")
            ("correct", po::value<string>(&correct_aln_fname), "correct Alignment(XMFA)")
			("calculated", po::value<string>(&procrast_aln_fname), "procrastAligner output")
			("sequence", po::value<string>(&sequence_file), "FastA sequence file")
        ;

                
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);    

        if (vm.count("help")) {
            cout << desc << "\n";
            return 1;
        }

        
    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }
	
	

	ifstream correct_in;
	correct_in.open( correct_aln_fname.c_str() );
	if( !correct_in.is_open() ){
		cerr << "Error opening " << correct_aln_fname << endl;
		return -1;
	}
	ifstream procrast_in;
	procrast_in.open( procrast_aln_fname.c_str() );
	if( !procrast_in.is_open() ){
		cerr << "Error opening " << procrast_aln_fname << endl;
		return -1;
	}
	
try{
	IntervalList correct_ivs;
	IntervalList procrast_ivs;
	std::vector< bitset_t > align_matrix;
	vector< gnSeqI > leftend;
	cout << "Reading correct alignment into interval list...";
	correct_ivs.ReadStandardAlignment( correct_in );
	cout << " finished" << endl;
	correct_in.close();

	cout << "Reading procrastAlignment into interval list...";
    procrast_ivs.ReadStandardAlignmentCompact( procrast_in );
	cout << " finished" << endl;
	procrast_in.close();

	gnSequence concat_sequence;
	concat_sequence.LoadSource( sequence_file );	// fixme, read this filename from command line or something -- this should be unaligned sequence
	compareAlignmentsAceD( correct_ivs, procrast_ivs, concat_sequence );

}catch( gnException& gne ){
	cerr << gne << endl;
}catch( exception& e ){
	cerr << e.what() << endl;
}catch( char const* c ){
	cerr << c << endl;
}

}
