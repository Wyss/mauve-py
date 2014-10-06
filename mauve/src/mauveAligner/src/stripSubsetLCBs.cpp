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
#include "libGenome/gnFASSource.h"
#include <boost/tuple/tuple.hpp>

using namespace std;
using namespace genome;
using namespace mems;

typedef boost::tuple< uint, gnSeqI, gnSeqI, vector< uint > > bbcol_t;

int main( int argc, char* argv[] )
{
	if( argc < 4 )
	{
		cerr << "Usage: stripSubsetLCBs <input xmfa> <input bbcols> <output xmfa> [min LCB size] [min genomes] [randomly subsample to X kb]\n";
		return -1;
	}
	ifstream aln_in;
	aln_in.open( argv[1] );
	if( !aln_in.is_open() ){
		cerr << "Error opening " << argv[1] << endl;
		return -1;
	}
	ifstream bbcols_in;
	bbcols_in.open( argv[2] );
	if( !bbcols_in.is_open() )
	{
		cerr << "Error opening " << argv[2] << endl;
		return -2;
	}
	ofstream aln_out;
	aln_out.open( argv[3] );
	if( !aln_out.is_open() ){
		cerr << "Error writing to " << argv[3] << endl;
		return -1;
	}

	size_t min_block_length = 0;
	if(argc>=5){
		min_block_length = atol(argv[4]);
	}
	size_t min_genome_count = -1;
	if(argc>=6){
		min_genome_count = atol(argv[5]);
	}
	size_t subsample_kb = 0;
	if(argc>=7){
		subsample_kb = atol(argv[6]);
	}
	

	try{
		IntervalList input_ivs;
		input_ivs.ReadStandardAlignment( aln_in );
		aln_in.close();

		LoadSequences( input_ivs, NULL );

		// read the bbcols file
		vector< bbcol_t > bbcols;
		string cur_line;
		while( getline( bbcols_in, cur_line ) )
		{
			stringstream line_str(cur_line);
			size_t cur_token;
			size_t tokenI = 0;
			bbcol_t bbcol;
			while( line_str >> cur_token )
			{
				switch(tokenI)
				{
					case 0:
						bbcol.get<0>() = cur_token;
						break;
					case 1:
						bbcol.get<1>() = cur_token;
						break;
					case 2:
						bbcol.get<2>() = cur_token;
						break;
					default:
						bbcol.get<3>().push_back(cur_token);
						break;
				}
				tokenI++;
			}
			bbcols.push_back(bbcol);
		}
		cout << "Read " << bbcols.size() << " backbone entries\n";

		IntervalList output_ivs;
		output_ivs.seq_table = input_ivs.seq_table;
		output_ivs.seq_filename = input_ivs.seq_filename;
/*		for( size_t i = 0; i < input_ivs.size(); ++i )
		{
			cout << "LCB " << i << " multiplicity: " << input_ivs[i].Multiplicity() << endl;
			for( size_t seqI = 0; seqI < input_ivs.seq_table.size(); ++seqI )
			{
				cout << input_ivs[i].LeftEnd(seqI) << '\t' << input_ivs[i].RightEnd(seqI) << '\t';
			}
			cout << endl;
			if( input_ivs[i].Multiplicity() == input_ivs.seq_table.size() )
				output_ivs.push_back( input_ivs[i] );
		}
*/
		cout << "seq_count is: " << input_ivs.seq_table.size() << endl;
		if(min_genome_count==-1) min_genome_count = input_ivs.seq_table.size();

		for( size_t bbI = 0; bbI < bbcols.size(); bbI++ )
		{
			if( bbcols[bbI].get<3>().size() < min_genome_count )
				continue;
			Interval* sub_iv = input_ivs[bbcols[bbI].get<0>()].Copy();
			sub_iv->CropStart( bbcols[bbI].get<1>() - 1 );
			sub_iv->CropEnd( sub_iv->Length() - bbcols[bbI].get<2>() );
			// calculate mean length
			size_t avglen = 0;
			for(size_t seqI=0; seqI < sub_iv->SeqCount(); seqI++){
				avglen += sub_iv->Length(seqI);
			}
			avglen /= sub_iv->SeqCount();
			if(avglen >= min_block_length){
				output_ivs.push_back( *sub_iv );
			}
			sub_iv->Free();
		}
		if(subsample_kb==0){
			cout << "output_ivs.size() " << output_ivs.size() << endl;
			output_ivs.WriteStandardAlignment( aln_out );
		}else{
			set<size_t> sampled;
			double cur_kb=0;
			for(; cur_kb < (double)subsample_kb && sampled.size() < output_ivs.size(); cur_kb++){
				int block = rand()%output_ivs.size();
				if(sampled.find(block)!=sampled.end()){
					continue;
				}
				sampled.insert(block);
				cur_kb += (double)(output_ivs[block].AlignmentLength()) / 1000.0;
				
			}
			IntervalList new_ivs;
			new_ivs.seq_table=output_ivs.seq_table;
			new_ivs.seq_filename=output_ivs.seq_filename;
			int i=0;
			for(set<size_t>::iterator siter = sampled.begin(); siter != sampled.end(); siter++){
				new_ivs.push_back(output_ivs[*siter]);
			}
			cout << "Writing " << cur_kb << " kb of alignment columns in " << new_ivs.size() << " blocks" << endl;
			new_ivs.WriteStandardAlignment( aln_out );
		}

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

