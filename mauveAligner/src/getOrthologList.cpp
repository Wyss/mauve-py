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
#include "libMems/DistanceMatrix.h"

using namespace std;
using namespace genome;
using namespace mems;

typedef boost::tuple< uint, gnSeqI, gnSeqI, vector< uint > > bbcol_t;

void printGI( ostream& out, gnBaseFeature* f )
{
	// print out the feature GI
	size_t qI = 0;
	for( ; qI < f->GetQualifierListLength(); qI++ )
	{
		if( f->GetQualifierName(qI) == "db_xref" )
		{
			string qval = f->GetQualifierValue(qI);
			if( qval.substr(0,4) == "\"GI:" )
			{
				out << qval;
			}
		}
	}
}

double computeAvgCoverage( vector< bb_seqentry_t >& backbone, vector< size_t >& nway_bb, vector< gnBaseFeature* >& ortho_cds )
{
	vector< double > covs( ortho_cds.size() );
	double cov_sum = 0;
	for( size_t oI = 0; oI < ortho_cds.size(); oI++ )
	{
		gnLocation floc = ortho_cds[oI]->GetLocation(0);
		double intlen = 0;
		for( size_t bbI = 0; bbI < nway_bb.size(); bbI++ )
		{
			gnLocation loc;
			loc.SetStart(absolut(backbone[nway_bb[bbI]][oI].first));
			loc.SetEnd(absolut(backbone[nway_bb[bbI]][oI].second));
			gnLocation intloc = floc.GetIntersection(loc,gnLocation::determinedRegions);
			intlen += intloc.GetEnd()-intloc.GetStart();
		}
		covs[oI] = intlen / (double)(floc.GetEnd()-floc.GetStart());
		cov_sum += covs[oI];
	}
	return cov_sum / ((double)ortho_cds.size());
}


int main( int argc, char* argv[] )
{
	if( argc < 6 )
	{
		cerr << "Usage: getOrthologList <input xmfa> <backbone seq file> <reference genome> <CDS ortholog filename> <CDS alignment base name>\n";
		return -1;
	}
	ifstream aln_in;
	aln_in.open( argv[1] );
	if( !aln_in.is_open() ){
		cerr << "Error opening " << argv[1] << endl;
		return -1;
	}
	uint sgI = atoi( argv[3] );
	string ortho_fname = argv[4];
	string output_base = argv[5];

	IntervalList input_ivs;
	input_ivs.ReadStandardAlignment( aln_in );
	aln_in.close();
	LoadSequences( input_ivs, &cout );

	size_t seq_count = input_ivs.seq_table.size();
	
	vector< bb_seqentry_t > backbone;
	ifstream bb_in;
	bb_in.open( argv[2] );
	if( !bb_in.is_open() ){
		cerr << "Error opening \"" << argv[2] << "\"" << endl;
		return -2;
	}
	readBackboneSeqFile( bb_in, backbone );
	bb_in.close();

	ofstream ortho_out( ortho_fname.c_str() );
	if( !ortho_out.is_open() )
	{
		cerr << "Error opening \"" << ortho_fname << "\"\n";
		return -3;
	}

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

	size_t ortho_count = 0;
	size_t rr_count = 0;
	size_t partial_rr_annotated = 0;

	ortho_out << "OrthoID";
	for( size_t seqI = 0; seqI < seq_count; seqI++ )
		ortho_out << "\tGI_in_Genome_" << seqI;
	ortho_out << "\tCoverage\tIdentity\tRearranged\n";

	// pick a gene at random from the first genome, extract the alignment, and write it to a file
	for( size_t geneI = 0; geneI < genes.size(); geneI++ )
	{
		if( geneI == 156 )
			cerr << "watchme\n";
		// is this gene part of N-way backbone?
		gnLocation loc = genes[geneI]->GetLocation(0);
		int64 lend = loc.GetFirst();
		int64 rend = loc.GetLast();
		vector< size_t > intersecting_bb;
		size_t bbI = 0;
		for( size_t bbI = 0; bbI < backbone.size(); bbI++ )
		{
			if( (absolut(backbone[bbI][sgI].first) <= lend && lend <= absolut(backbone[bbI][sgI].second)) ||
			    (absolut(backbone[bbI][sgI].first) <= rend && rend <= absolut(backbone[bbI][sgI].second)) ||
			    (lend <= absolut(backbone[bbI][sgI].first) && absolut(backbone[bbI][sgI].first) <= rend) )
				intersecting_bb.push_back(bbI);
		}
		vector< size_t > nway_bb;
		for( size_t bbI = 0; bbI < intersecting_bb.size(); bbI++ )
		{
			size_t seqI = 0;
			for( ; seqI < input_ivs.seq_table.size(); ++seqI )
			{
				if( backbone[intersecting_bb[bbI]][seqI].first == 0 || backbone[intersecting_bb[bbI]][seqI].second == 0 )
					break;
			}
			if( seqI == input_ivs.seq_table.size() )
				nway_bb.push_back(intersecting_bb[bbI]);
		}

		// skip to the next CDS if this one wasn't part of some n-way backbone
		if( nway_bb.size() == 0 )
			continue;

		// use the alignment to find CDS that overlap in this region


		// extract the alignment
		size_t ivI = 0;
		// identify the interval that has the biggest intersection
		vector< pair< size_t, size_t > > iv_overlap;
		for( ivI = 0; ivI < input_ivs.size(); ivI++ )
		{
			if( input_ivs[ivI].Start(sgI) != NO_MATCH )
			{
				size_t inter_size = 0;
				for( size_t bbI = 0; bbI < nway_bb.size(); bbI++ )
				{
					gnLocation loc1;
					loc1.SetStart( input_ivs[ivI].LeftEnd(sgI) );
					loc1.SetEnd( input_ivs[ivI].RightEnd(sgI) );
					gnLocation loc2;
					loc2.SetStart( absolut(backbone[nway_bb[bbI]][sgI].first) );
					loc2.SetEnd( absolut(backbone[nway_bb[bbI]][sgI].second) );
					gnLocation intloc = loc1.GetIntersection( loc2, gnLocation::determinedRegions );
					gnLocation intloc2 = intloc.GetIntersection( loc, gnLocation::determinedRegions );
					inter_size += intloc2.GetEnd() - intloc2.GetStart();
				}
				if( inter_size > 0 )
					iv_overlap.push_back( make_pair( inter_size, ivI ) );
			}
		}
		bool partial_rr = false;
		std::sort( iv_overlap.begin(), iv_overlap.end() );
		if( iv_overlap.size() == 0 )
		{
			cerr << "Warning: unable to assign gene to an interval!\n" << "coordinates: " << lend << '\t' << rend << endl;
			continue;
		}else{
			ivI = iv_overlap.back().second;
			if( iv_overlap.size() > 1 )
			{
				partial_rr = true;
				rr_count++;
			}
		}
		CompactGappedAlignment<> iv_cga(input_ivs[ivI]);
		CompactGappedAlignment<> col_cga;
		gnLocation loc1;
		loc1.SetStart( input_ivs[ivI].LeftEnd(sgI) );
		loc1.SetEnd( input_ivs[ivI].RightEnd(sgI) );
		gnLocation intloc = loc1.GetIntersection( loc, gnLocation::determinedRegions );
		gnSeqI lcol = iv_cga.SeqPosToColumn( sgI, intloc.GetStart() );
		gnSeqI rcol = iv_cga.SeqPosToColumn( sgI, intloc.GetEnd() );
		if( rcol < lcol )
			swap( rcol, lcol );	// handle reverse complement
		iv_cga.copyRange(col_cga, lcol, rcol-lcol + 1);
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

		stringstream of_name;
		of_name << output_base << "_" << ortho_count << ".fas";
		gnFASSource::Write( gene_aln, of_name.str() );

		// find orthologous CDS features...
		vector< gnBaseFeature* > ortho_cds( seq_count, NULL );
		size_t ocds_count = 0;
		for( size_t seqI = 0; seqI < input_ivs.seq_table.size(); seqI++ )
		{
			gnLocation seqloc;
			seqloc.SetStart(col_cga.LeftEnd(seqI));
			seqloc.SetEnd(col_cga.RightEnd(seqI));
			vector< gnBaseFeature* > int_feats;
			vector< uint32 > indie;
			input_ivs.seq_table[seqI]->getIntersectingFeatures( seqloc, int_feats, indie );
			vector< pair< gnSeqI, size_t > > overlap_frac;
			for( size_t featI = 0; featI < int_feats.size(); featI++ )
			{
				if( int_feats[featI]->GetName() == "CDS" )
				{
					gnLocation l = seqloc.GetIntersection( int_feats[featI]->GetLocation(0), gnLocation::determinedRegions );
					size_t max_bb = 0;
					for( size_t bbI = 0; bbI < nway_bb.size(); bbI++ )
					{
						gnLocation bbloc;
						bbloc.SetBounds( absolut(backbone[nway_bb[bbI]][seqI].first), absolut(backbone[nway_bb[bbI]][seqI].second) );
						gnLocation l2 = bbloc.GetIntersection( l, gnLocation::determinedRegions );
						if( l2.GetEnd() - l2.GetStart() > max_bb )
							max_bb = l2.GetEnd() - l2.GetStart();
					}
					overlap_frac.push_back( make_pair( max_bb, featI ) );
				}else
					delete int_feats[featI];
			}
			std::sort( overlap_frac.begin(), overlap_frac.end() );
			if( overlap_frac.size() > 0 )
			{
				ortho_cds[seqI] = int_feats[ overlap_frac.back().second ];
				ocds_count++;
			}
		}

		if( ocds_count == seq_count )
		{
			if( ortho_count == 88 )
				cerr << "watchme\n";
			ortho_out << ortho_count;
			for( size_t i = 0; i < seq_count; i++ )
			{
				ortho_out << '\t';
				printGI( ortho_out, ortho_cds[i] );
			}

			double cov = computeAvgCoverage( backbone, nway_bb, ortho_cds );
			ortho_out << '\t' << cov;
			NumericMatrix<double> identity;
			vector< AbstractMatch* > amvec( 1, &col_cga );
			BackboneIdentityMatrix( amvec, input_ivs.seq_table, identity );
			double id = 0;
			for( size_t i = 0; i < seq_count; i++ )
				for( size_t j = i+1; j < seq_count; j++ )
					id += identity(i,j);
			id /= (double)(seq_count * (seq_count-1)) / 2.0;

			ortho_out << '\t' << id;
			ortho_out << '\t';
			if( partial_rr )
			{
				partial_rr_annotated++;
				ortho_out << "*";
			}
			ortho_out << endl;
			ortho_count++;
		}
		for( size_t oI = 0; oI < ortho_cds.size(); oI++ )
			if( ortho_cds[oI] != NULL )
				delete ortho_cds[oI];

	}
	cout << ortho_count << " out of " << genes.size() << " genes were at least partially conserved\n";
	cout << rr_count << " CDS appear to be broken by rearrangement, of which " << partial_rr_annotated << " are still annotated as CDS in all genomes\n";
}

