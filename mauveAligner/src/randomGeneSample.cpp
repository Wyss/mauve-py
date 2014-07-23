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
#include "libMems/Backbone.h"
#include "libGenome/gnFeature.h"
#include "libGenome/gnFASSource.h"

using namespace std;
using namespace genome;
using namespace mems;

typedef boost::tuple< uint, gnSeqI, gnSeqI, vector< uint > > bbcol_t;

int main( int argc, char* argv[] )
{
	if( argc < 6 )
	{
		cerr << "Usage: randomGeneSample <input xmfa> <backbone seq file> <sample genome> <number of genes> <output base name> [random seed]\n";
		return -1;
	}
	ifstream aln_in;
	aln_in.open( argv[1] );
	if( !aln_in.is_open() ){
		cerr << "Error opening " << argv[1] << endl;
		return -1;
	}
	uint gene_count = atoi( argv[4] );
	uint sgI = atoi( argv[3] );
	string output_base = argv[5];

	if( argc == 7 )
		srand(atoi(argv[6]));
	else
		srand(time(NULL));

	IntervalList input_ivs;
	input_ivs.ReadStandardAlignment( aln_in );
	aln_in.close();
	LoadSequences( input_ivs, &cout );
	
	vector< bb_seqentry_t > backbone;
	ifstream bb_in;
	bb_in.open( argv[2] );
	if( !bb_in.is_open() ){
		cerr << "Error opening \"" << argv[2] << "\"" << endl;
		return -2;
	}
	readBackboneSeqFile( bb_in, backbone );
	bb_in.close();

	gnSequence* gen0 = input_ivs.seq_table[sgI];
	vector< gnBaseFeature* > genes;
	for( size_t featI = 0; featI < gen0->getFeatureListLength(); featI++ )
	{
		gnBaseFeature* feat = gen0->getFeature(featI);
		if( feat->GetName() == "CDS" )
			genes.push_back( feat );
		else
			delete feat;
	}

	cout << genes.size() << " of the " << gen0->getFeatureListLength() << " annotated features are CDS\n";

	// pick a gene at random from the first genome, extract the alignment, and write it to a file
	for( size_t geneI = 0; geneI < gene_count; geneI++ )
	{
		cerr << "picking gene\n";
		int randy;
		do{
			randy = rand() % genes.size();
			// has this gene already been used?
			if( genes[randy] == NULL )
				continue;
			// is this gene part of N-way backbone?
			gnLocation loc = genes[randy]->GetLocation(0);
			int64 lend = loc.GetFirst();
			int64 rend = loc.GetLast();
			size_t bbI = 0;
			for( ; bbI < backbone.size(); bbI++ )
			{
				if( genome::absolut(backbone[bbI][sgI].first) <= lend && rend <= genome::absolut(backbone[bbI][sgI].second) )
					break;
			}
			size_t seqI = 0;
			for( ; bbI < backbone.size() && seqI < input_ivs.seq_table.size(); ++seqI )
			{
				if( backbone[bbI][seqI].first == 0 || backbone[bbI][seqI].second == 0 )
					break;
			}
			if( seqI == input_ivs.seq_table.size() && bbI < backbone.size() )
				break;	// found a containing segment
		}while(true);
		// print out the feature name
		for( size_t qI = 0; qI < genes[randy]->GetQualifierListLength(); qI++ )
		{
			if( genes[randy]->GetQualifierName(qI) == "gene" )
				cout << "gene:\t" << genes[randy]->GetQualifierValue(qI) << endl;
		}
		// extract the alignment
		gnLocation loc = genes[randy]->GetLocation(0);
		int64 lend = loc.GetFirst();
		int64 rend = loc.GetLast();
		cerr << "lend: " << lend << "\trend: " << rend << endl;
		size_t ivI = 0;
		for( ivI = 0; ivI < input_ivs.size(); ivI++ )
		{
			if( input_ivs[ivI].Start(sgI) != NO_MATCH )
			{
//				cerr << "iv: " << ivI << "\tstart: " << input_ivs[ivI].Start(sgI) << "\tlength: " << input_ivs[ivI].Length(sgI) << endl;
				gnSeqI iv_rend = genome::absolut(input_ivs[ivI].Start(sgI)) + input_ivs[ivI].Length(sgI);
				if(  genome::absolut(input_ivs[ivI].Start(sgI)) < lend && rend < iv_rend )
					break;
			}
		}
		if( ivI == input_ivs.size() )
			cerr << "Error: unable to assign gene to an interval!\n" << "coordinates: " << lend << '\t' << rend << endl;
		cerr << "making iv_cga\n";
		CompactGappedAlignment<> iv_cga(input_ivs[ivI]);
		CompactGappedAlignment<> col_cga;
		cerr << "getting left and right cols\n";
		gnSeqI lcol = iv_cga.SeqPosToColumn( sgI, lend );
		gnSeqI rcol = iv_cga.SeqPosToColumn( sgI, rend );
		cerr << "left col: " << lcol << "\tright_col: " << rcol << endl;
		iv_cga.copyRange(col_cga, lcol, rcol-lcol + 1);
		cerr << "getting alignment\n";
		vector< string > aln;
		GetAlignment( col_cga, input_ivs.seq_table, aln );
		gnSequence gene_aln;
		for( size_t i = 0; i < aln.size(); i++ )
		{
			gene_aln += aln[i];
			stringstream ss;
			ss << "seq" << i;
			gene_aln.setContigName(i, ss.str());
		}
		cerr << "writing fasta\n";
		stringstream of_name;
		of_name << output_base << "_" << geneI << ".fas";
		gnFASSource::Write( gene_aln, of_name.str() );

		// done with this gene
		delete genes[randy];
		genes[randy] = NULL;
	}

}

