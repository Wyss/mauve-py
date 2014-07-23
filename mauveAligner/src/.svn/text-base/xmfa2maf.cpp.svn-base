#include "libMems/IntervalList.h"
#include "libMems/ProgressiveAligner.h"
#include <fstream>

using namespace mems;
using namespace std;
using namespace genome;

int main(int argc, char* argv[] ){
	if(argc != 3){
		cerr << "Usage: xmfa2maf <xmfa input> <maf output>\n";
		return -1;
	}
	ifstream ifile(argv[1]);
	if(!ifile.is_open()){
		cerr << "Error reading \"" << argv[1] << "\"\n";
		return -2;
	}
	ofstream ofile(argv[2]);
	if(!ofile.is_open()){
		cerr << "Error writing to \"" << argv[2] << "\"\n";
		return -2;
	}

	IntervalList xmfa;
	xmfa.ReadStandardAlignment(ifile);
	LoadSequences(xmfa, &cout);

	// break alignments on chromosome boundaries
	vector<AbstractMatch*> alignments;
	for(int ivI=0; ivI < xmfa.size(); ivI++){
		alignments.push_back( xmfa[ivI].Clone() );
	}

	vector< vector< gnSeqI > > chromo_bounds( xmfa.seq_table.size() );
	for(int seqI=0; seqI < xmfa.seq_table.size(); seqI++){
		for(int cI=1; cI < xmfa.seq_table[seqI]->contigListSize(); cI++){
			chromo_bounds[seqI].push_back( xmfa.seq_table[seqI]->contigStart(cI) );
		}
		SSC<AbstractMatch> msc( seqI );
		sort( alignments.begin(), alignments.end(), msc );
		AbstractMatchSeqManipulator amsm( seqI );
		applyBreakpoints( chromo_bounds[seqI], alignments, amsm );
	}

	ofile << "##maf version=1 program=progressiveMauve\n";
	for(int ivI=0; ivI < alignments.size(); ivI++ ){
		ofile << "a\n";
		vector<string> aln;
		GetAlignment( *((Interval*)(alignments[ivI])), xmfa.seq_table, aln );

		for( int seqI=0; seqI < xmfa.seq_filename.size(); seqI++ ){
			if(alignments[ivI]->LeftEnd(seqI)==0)
				continue;	// sequence not defined in this block

			// determine which contig this alignment is in
			uint32 l_contigI, r_contigI;
			gnSeqI l_baseI = alignments[ivI]->LeftEnd(seqI);
			gnSeqI r_baseI = alignments[ivI]->RightEnd(seqI)-1;
			xmfa.seq_table[seqI]->globalToLocal( l_contigI, l_baseI );
			xmfa.seq_table[seqI]->globalToLocal( r_contigI, r_baseI );
			string contig_name = xmfa.seq_table[seqI]->contigName( l_contigI );
			if(l_contigI != r_contigI){
				cerr << "interval " << ivI << " seq " << seqI << " left " << alignments[ivI]->LeftEnd(seqI) << " right " << alignments[ivI]->RightEnd(seqI) << endl;
				cerr << "l_baseI " << l_baseI << " r_baseI " << r_baseI << " l_contigI " << l_contigI << " r_contigI " << r_contigI << " name " << contig_name << endl;
				cerr << "Error, input alignment must spans multiple contigs/chromosomes. Unable to translate to MAF\n";
				return -1;
			}
			ofile << "s " << xmfa.seq_filename[seqI] << "." << contig_name;
			ofile.flush();

			int64 lend = l_baseI-1;
			if(alignments[ivI]->Orientation(seqI) == AbstractMatch::reverse){
				lend = xmfa.seq_table[seqI]->contigLength(l_contigI) - l_baseI - alignments[ivI]->Length(seqI) + 1;
			}
			ofile << " " << lend;
			ofile << " " << alignments[ivI]->Length(seqI);
			ofile << " " << (alignments[ivI]->Orientation(seqI) == AbstractMatch::reverse ? "-" : "+");
			ofile << " " << xmfa.seq_table[seqI]->contigLength(l_contigI);
			ofile << " " << aln[seqI] << endl;
		}
		ofile << endl;
	}
	ofile.close();

	return 0;
}
