#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include "libGenome/gnFilter.h"
#include "libMems/IntervalList.h"
#include "libMems/MatchList.h"
#include "libMems/GappedAlignment.h"
#include "libMems/Matrix.h"
#include "libMems/MatchProjectionAdapter.h"
#include "libMems/Aligner.h"
#include "libGenome/gnFASSource.h"
#include "libMems/ProgressiveAligner.h"
using namespace std;
using namespace genome;
using namespace mems;




int main( int argc, char* argv[] )
{
	if( argc < 6 )
	{
		cerr << "Usage: alignmentProjector <input xmfa> <output xmfa> <mfa seq input> <mfa seq output> <list of seqs to include, starting at 0>\n";
		return -1;
	}
	ifstream aln_in;
	aln_in.open( argv[1] );
	if( !aln_in.is_open() ){
		cerr << "Error opening " << argv[1] << endl;
		return -1;
	}
	ofstream aln_out;
	aln_out.open( argv[2] );
	if( !aln_out.is_open() ){
		cerr << "Error writing to " << argv[2] << endl;
		return -1;
	}
	string mfa_seqs = argv[3];
	string mfa_output = argv[4];
	
	try{
		IntervalList input_ivs;
		input_ivs.ReadStandardAlignment( aln_in );
		aln_in.close();

		MatchList ml;
		ml.seq_filename = input_ivs.seq_filename;
		LoadMFASequences( ml, mfa_seqs, NULL );
		input_ivs.seq_table = ml.seq_table;

		// create a projection list
		vector< uint > projection;
		IntervalList proj_ivs;
		for( int i = 5; i < argc; ++i )
		{
			projection.push_back( atoi( argv[i] ) );
			proj_ivs.seq_filename.push_back( mfa_seqs );
			proj_ivs.seq_table.push_back( input_ivs.seq_table[projection.back()] );
		}

		vector< vector< MatchProjectionAdapter* > > LCB_list;
		vector< LCB > projected_adjs;
		projectIntervalList( input_ivs, projection, LCB_list, projected_adjs );

		cout << "projection has " << LCB_list.size() << " LCBs\n";
		proj_ivs.resize( LCB_list.size() );
		for( size_t lcbI = 0; lcbI < LCB_list.size(); ++lcbI )
			proj_ivs[lcbI].SetMatches( LCB_list[lcbI] );

		proj_ivs.WriteStandardAlignment( aln_out );

		gnSequence seq;
		seq.LoadSource( mfa_seqs );
		ofstream seq_out( mfa_output.c_str() );
		gnSequence proj_seq;
		for( size_t projI = 0; projI < projection.size(); ++projI )
			proj_seq += seq.contig(projection[projI]);
		gnFASSource::Write(proj_seq,seq_out,false,false);

	}catch( gnException& gne ){
		cerr << gne << endl;
		return -1;
	}catch( exception& e ){
		cerr << e.what() << endl;
		return -2;
	}catch( char const* c ){
		cerr << c << endl;
		return -3;
	}catch(...){
		cerr << "Unhandled exception" << endl;
		return -4;
	}
}

