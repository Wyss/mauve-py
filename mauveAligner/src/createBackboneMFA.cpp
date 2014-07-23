#include "libMems/Interval.h"
#include "libMems/Islands.h"
#include "libGenome/gnFASSource.h"

using namespace std;
using namespace mems;
using namespace genome;

int main( int argc, char* argv[] )
{
	IntervalList iv_list;
	if( argc != 3 )
	{
		cerr << "Usage: <input interval file> <output MFA name>\n";
		return -1;
	}
	ifstream in_file( argv[1] );
	if( !in_file.is_open() )
	{
		cerr << "Error opening \"" << argv[1] << "\"\n";
		return -1;
	}
	iv_list.ReadList( in_file );
	LoadSequences(iv_list, NULL);
	string base_name = argv[2];
	cout << "Input alignment has " << iv_list.size() << " intervals\n";
	vector< string > superaln = vector<string>( iv_list.seq_table.size() );
	for( uint lcbI = 0; lcbI < iv_list.size(); lcbI++ )
	{
		// only use 1/30 LCBs
		if( lcbI % 30 != 0 )
			continue;
		gnAlignedSequences gnas;
		iv_list[lcbI].GetAlignedSequences( gnas, iv_list.seq_table );
		for( uint seqI = 0; seqI < gnas.sequences.size(); seqI++ )
			superaln[seqI] += gnas.sequences[seqI];
	}

	ofstream out_file( base_name.c_str() );
	if( !out_file.is_open() )
	{
		cerr << "Error opening \"" << base_name << "\"\n";
		return -2;
	}
	gnSequence gns;
	for( uint seqI = 0; seqI < superaln.size(); seqI++ )
	{
		stringstream seq_name;
		seq_name << seqI;
//		seq_name << "(" << iv_list[lcbI].Start(seqI) << "-" << iv_list[lcbI].Start(seqI) + iv_list[lcbI].Length(seqI) << ")";
		gns += superaln[seqI];
		gns.setContigName( gns.contigListSize()-1, seq_name.str() );
	}
	gnFASSource::Write( gns, out_file, false, false );
	return 0;
}

