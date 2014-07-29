/*******************************************************************************
 * $Id: calculateCoverage.cpp,v 1.5 2004/02/28 00:01:31 darling Exp $
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
	cerr << "Usage: " << pname << " <source mauve alignment>  <sequence 1>...<sequence N>\n";
}

#ifndef NELEMS
#define NELEMS(a) ( sizeof( a ) / sizeof( *a ) )
#endif

int main( int argc, const char* argv[] ){

try{
	if( argc <= 0 ){
		print_usage( "calculateCoverage" );
		return -1;
	}
	if( argc < 2 ){
		print_usage( argv[0] );
		return -1;
	}
	
	//
	// Load sequences
	//
	string alignment_fname = argv[1];
	vector< string > sequence_fname;
	vector< gnSequence* > source_seqs;
	for( uint argI = 2; argI < argc; argI++ ){
		sequence_fname.push_back( argv[ argI ] );
		cout << "Loading " << sequence_fname[ argI - 2 ];
		try{
			source_seqs.push_back( new gnSequence() );
			source_seqs[ argI - 2 ]->LoadSource( sequence_fname[ argI - 2 ] );
		}catch( gnException& gne ){
			cerr << gne << endl;
			return -1;
		}
		cout << "   " << source_seqs[ argI - 2 ]->length() << " bp\n";
	}
	
	//
	// Load IntervalList matches
	//
	ifstream alignment_in;
	alignment_in.open( alignment_fname.c_str() );
	if( !alignment_in.is_open() ){
		cerr << "Error opening " << alignment_fname << endl;
		return -1;
	}
	
	cout << "Loading alignment...\n";
	IntervalList aligned_ivs;
	aligned_ivs.ReadList( alignment_in );
	
	for( uint ivI = 0; ivI < aligned_ivs.size(); ivI++ ){
		cout << "Interval " << ivI;
		Interval& iv = aligned_ivs[ ivI ];
		for( uint seqI = 0; seqI < source_seqs.size(); seqI++ ){
			cout << '\t' << iv.Length( seqI );
		}
		cout << endl;
	}

	
}catch( gnException& gne ){
	cerr << gne << endl;
}catch( exception& e ){
	cerr << e.what() << endl;
}catch(...){

}
	return 0;
}
