#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include "libGenome/gnFilter.h"
#include "libMems/IntervalList.h"
#include "libMems/MatchList.h"
#include "libMems/GappedAlignment.h"
#include "libMems/Matrix.h"
#include "libMems/MatchProjectionAdapter.h"
#include "libMems/Aligner.h"
#include "libMems/Islands.h"
#include "libGenome/gnFASSource.h"
#include <boost/tuple/tuple.hpp>
#include "libMems/ProgressiveAligner.h"

using namespace std;
using namespace genome;
using namespace mems;

typedef boost::tuple< uint, gnSeqI, gnSeqI, vector< uint > > bbcol_t;

int main( int argc, char* argv[] )
{
	if( argc < 5 )
	{
		cerr << "Usage: projectAndStrip <input xmfa> <output xmfa> <seq1> <seq2>...<seqN>\n";
		cerr << "\nNumeric sequence identifiers start at 0.\n";
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
	vector<uint> seq_ids(argc-3);
	vector<uint> not_ids;
	for( size_t i = 3; i < argc; ++i )
		seq_ids[i - 3] = atoi(argv[i]);

	try{
		IntervalList input_ivs;
		input_ivs.ReadStandardAlignment( aln_in );
		aln_in.close();

		LoadSequences( input_ivs, NULL );

		not_ids.resize( input_ivs.seq_table.size() );
		for( size_t i = 0; i < not_ids.size(); i++ )
			not_ids[i] = i;
		for( size_t i = 0; i < seq_ids.size(); i++ )
			not_ids[seq_ids[i]] = (std::numeric_limits<size_t>::max)();
		std::sort( not_ids.begin(), not_ids.end() );
		not_ids.resize( not_ids.size() - seq_ids.size() );

		IntervalList output_ivs;
		output_ivs.seq_table = input_ivs.seq_table;
		output_ivs.seq_filename = input_ivs.seq_filename;
		
		vector< GappedAlignment* > gaga_list;
		
		for( size_t ivI = 0; ivI < input_ivs.size(); ivI++ )
		{
			Interval& iv = input_ivs[ivI];
			size_t j = 0;
			for( ; j < seq_ids.size(); j++ )
			{
				if( iv.LeftEnd( seq_ids[j] ) == NO_MATCH )
					break;
			}
			if( j == seq_ids.size() )
			{
				vector<string> aln_mat;
				GetAlignment( iv, input_ivs.seq_table, aln_mat );
				Interval new_iv;
				GappedAlignment ga(seq_ids.size(), 0);
				GappedAlignment* gaga = ga.Copy();
				vector<string> sub_mat( seq_ids.size() );
				for( size_t sI = 0; sI < seq_ids.size(); sI++ )
				{
					gaga->SetStart( sI, iv.Start(seq_ids[sI]) );
					gaga->SetLength( iv.Length(seq_ids[sI]), sI );
					swap( sub_mat[sI], aln_mat[seq_ids[sI]] );
				}
				gaga->SetAlignment(sub_mat);
				gaga_list.push_back( gaga );
			}
		}

		for( size_t gI = 0; gI < gaga_list.size(); gI++ )
			if( gaga_list[gI]->Orientation(0) == AbstractMatch::reverse )
				gaga_list[gI]->Invert();

		cout << "constructing LCBs\n";
		vector< gnSeqI > bps;
		IntervalList real_out_ivs;
		IdentifyBreakpoints(gaga_list, bps);
		vector< vector< GappedAlignment* > > coal_ivs;
		ComputeLCBs_v2(gaga_list, bps, coal_ivs);
		real_out_ivs.seq_filename.resize(seq_ids.size());
		real_out_ivs.seq_table.resize(seq_ids.size());
		for( size_t sI = 0; sI < seq_ids.size(); sI++ )
		{
			real_out_ivs.seq_filename[sI] = input_ivs.seq_filename[seq_ids[sI]];
			real_out_ivs.seq_table[sI] = input_ivs.seq_table[seq_ids[sI]];
		}
		real_out_ivs.resize( coal_ivs.size() );
		for( size_t cI = 0; cI < coal_ivs.size(); cI++ )
			real_out_ivs[cI].SetMatches(coal_ivs[cI]);
		cout << "real_out_ivs.size() " << real_out_ivs.size() << endl;



		addUnalignedIntervals( real_out_ivs );
		real_out_ivs.WriteStandardAlignment( aln_out );
		aln_out.close();
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

