/*******************************************************************************
 * $Id: extractBackbone.cpp,v 1.2 2004/02/28 00:01:31 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/IntervalList.h"
#include "libMems/Islands.h"

using namespace std;
using namespace genome;
using namespace mems;

void print_usage( const char* pname ){
	cerr << "Usage: " << pname << " <mauve alignment> <min bb sequence length> <max bb gap size> <backbone output>\n";
}

int main( int argc, const char* argv[] ){
	if( argc <= 0 ){
		print_usage( "extractBackbone" );
		return -1;
	}
	if( argc != 5 ){
		print_usage( argv[0] );
		return -1;
	}
	
	string alignment_fname = argv[1];
	int64 min_bb_length = atol( argv[2] );
	int64 max_gap_length = atol( argv[3] );
	string output_fname = argv[4];

	ifstream alignment_in;
	alignment_in.open( alignment_fname.c_str() );
	if( !alignment_in.is_open() ){
		cerr << "Error opening " << alignment_fname << endl;
		return -1;
	}
	
	
	IntervalList aligned_ivs;
	aligned_ivs.ReadList( alignment_in );
	LoadSequences(aligned_ivs, &cout);

	vector< GappedAlignment > backbone_data;
	simpleFindBackbone( aligned_ivs, min_bb_length, max_gap_length, backbone_data );
	IntervalList backbone_ivs;
	backbone_ivs.seq_table = aligned_ivs.seq_table;
	backbone_ivs.seq_filename = aligned_ivs.seq_filename;
	// construct a new IntervalList containing only backbone regions
	for( uint bbI = 0; bbI < backbone_data.size(); bbI++ ){
		vector< AbstractMatch* > tmp( 1, &backbone_data[ bbI ] );
		backbone_ivs.push_back( Interval( tmp.begin(), tmp.end() ) );
	}
	
	ofstream output( output_fname.c_str() );
	if( !output.is_open() ){
		cerr << "Error opening " << output_fname << endl;
		return -1;
	}
	backbone_ivs.WriteList( output );

	return 0;
}
