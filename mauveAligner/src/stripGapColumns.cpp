#include "libMems/IntervalList.h"
#include "libMems/MatchList.h"
#include "libMems/GappedAlignment.h"
#include <fstream>
#include <string>
#include <vector>

using namespace std;
using namespace genome;
using namespace mems;

int main( int argc, char* argv[] )
{
	if( argc != 3 )
	{
		cerr << "Usage: stripGapColumns <input XMFA> <output XMFA>\n";
		return -1;
	}

	ifstream aln_infile( argv[1] );
	if( !aln_infile.is_open() )
	{
		cerr << "Error opening \"" << argv[1] << "\"\n";
		return -1;
	}
	IntervalList iv_list;
	iv_list.ReadStandardAlignment( aln_infile );
	LoadSequences( iv_list, &cout );
	IntervalList iv_outlist;
	iv_outlist.seq_filename = iv_list.seq_filename;
	iv_outlist.seq_table = iv_list.seq_table;
	for( uint ivI = 0; ivI < iv_list.size(); ivI++ )
	{
		Interval& cur_iv = iv_list[ivI];
		vector< string > alignment;
		GetAlignment( cur_iv, iv_list.seq_table, alignment );
		vector< string > seq_align = vector< string >( cur_iv.SeqCount() );
		for( gnSeqI colI = 0; colI < cur_iv.AlignmentLength(); colI++ )
		{
			uint seqI = 0;
			for( ; seqI < cur_iv.SeqCount(); seqI++ )
			{
				if( alignment[seqI][colI] == '-' )
					break;
			}
			if( seqI != cur_iv.SeqCount() )
				continue;
			for( seqI = 0; seqI < cur_iv.SeqCount(); seqI++ )
			{
				seq_align[seqI] += alignment[seqI][colI];
			}
		}

		GappedAlignment* new_ga = new GappedAlignment( seq_align.size(), seq_align[0].size() );
		new_ga->SetAlignment( seq_align );
		for( uint seqI = 0; seqI < cur_iv.SeqCount(); seqI++ )
		{
			new_ga->SetStart( seqI, cur_iv.Start( seqI ) );
			new_ga->SetLength( cur_iv.Length( seqI ), seqI );
		}
		vector< AbstractMatch* > am_list( 1, new_ga );
		Interval new_iv(am_list.begin(), am_list.end());
		iv_outlist.push_back( new_iv );
	}

	ofstream iv_outfile( argv[2] );
	if( !iv_outfile.is_open() )
	{
		cerr << "Error opening \"" << argv[2] << "\"\n" << endl;
		return -2;
	}
	iv_outlist.WriteStandardAlignment( iv_outfile );
	return 0;
}
