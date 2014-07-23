#include "libGenome/gnFASSource.h"

using namespace std;
using namespace genome;

int main( int argc, char* argv[] )
{
	if( argc != 2 )
	{
		cerr << "Usage: gappiness <MFA file>\n";
	}
	string aln_fname = argv[1];
	gnSequence gns;
	gns.LoadSource( aln_fname );
	cout << "aln_length\t" << gns.contig(0).length() << endl;
	gnSeqI total_len = 0;
	for( uint seqI = 0; seqI < gns.contigListSize(); seqI++ )
	{
		string cur_seq = gns.contig(seqI).ToString();
		gnSeqI len = 0;
		for( size_t charI = 0; charI < cur_seq.size(); charI++ )
		{
			if( cur_seq[charI] != '-' )
				len++;
		}
		cout << "seq" << seqI << "_len\t" << len << endl;
		total_len += len;
	}
	double avg_seq_len = (double)total_len / (double)gns.contigListSize();
	cout << "avg_seq_len\t" << avg_seq_len << endl;
	cout << "gappiness\t" << (double)(gns.contig(0).length()) / avg_seq_len << endl;

	// compute average pairwise identity
	gnSeqI total_id = 0;
	gnSeqI total_possible = 0;
	for( uint seqI = 0; seqI < gns.contigListSize(); seqI++ )
		for( uint seqJ = seqI + 1; seqJ < gns.contigListSize(); seqJ++ )
	{
		string cur_seqI = gns.contig(seqI).ToString();
		string cur_seqJ = gns.contig(seqJ).ToString();
		for( size_t colI = 0; colI < cur_seqI.size(); colI++ )
		{
			if( cur_seqI[colI] == '-' || cur_seqJ[colI] == '-' )
				continue;
			total_possible++;
			if( toupper(cur_seqI[colI]) == toupper(cur_seqJ[colI]) )
				total_id++;
		}
	}
	cout << "percent_id\t" << (double)total_id / (double)total_possible << endl;
	return 0;
}

