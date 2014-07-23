/*******************************************************************************
 * $Id: uniqueMerCount.cpp,v 1.1 2004/02/28 00:01:31 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/DNAFileSML.h"

using namespace std;
using namespace genome;
using namespace mems;

void print_usage( const char* pname ){
	cerr << "Usage: " << pname << " <Sorted Mer List>\n";
}

int main( int argc, const char* argv[] ){
	if( argc != 2 ){
		print_usage("uniqueMerCount");
		return -1;
	}

	string sml_filename = argv[1];
	DNAFileSML* file_sml = new DNAFileSML();
	boolean success = true;
	try{
		file_sml->LoadFile( sml_filename );
	}catch( gnException& gne ){
		success = false;
		cerr << gne << endl;
		return -1;
	}
	cout << endl << file_sml->UniqueMerCount() << endl;
}

