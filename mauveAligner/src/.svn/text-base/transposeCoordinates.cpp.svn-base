/*******************************************************************************
 * $Id: transposeCoordinates.cpp,v 1.1 2004/02/28 00:01:31 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include "libMems/Aligner.h"

using namespace std;
using namespace genome;
using namespace mems;

void print_usage( const char* pname ){
	cout << "Usage: " << pname << " <match list> <coordinates file> <sequence ID> <match list output>\n";
}

int main( int argc, const char* argv[] ){
	if( argc != 5 ){
		print_usage("transposeCoordinates");
		return -1;
	}
	
	string match_filename = argv[1];
	ifstream match_file( match_filename.c_str() );
	if( !match_file.is_open() ){
		cerr << "Error opening \"" << match_filename << "\"" << endl;
		return -1;
	}

	string coord_filename = argv[2];
	ifstream coord_file( coord_filename.c_str() );
	if( !coord_file.is_open() ){
		cerr << "Error opening \"" << coord_filename << "\"" << endl;
		return -1;
	}
	
	int trans_seq = atoi( argv[3] );
	
	MatchList mlist;
	ReadList( mlist, match_file );
	mlist.MultiplicityFilter( mlist.seq_filename.size() );
	
	int64 coord;
	vector< int64 > coord_list;
	while( coord_file >> coord ){
		coord_list.push_back( coord );
	}
	transposeMatches( mlist, trans_seq, coord_list );
//	for( uint ivI = 0; ivI < iv_list.size(); ivI++ ){
//	}
	
	string match_outname = argv[4];
	ofstream match_out( match_outname.c_str() );
	if( !match_out.is_open() ){
		cerr << "Error opening \"" << match_outname << "\"" << endl;
		return -1;
	}
	WriteList( mlist, match_out );
	
	
	return 0;
}


