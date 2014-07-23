#include "libMems/Interval.h"
#include "libMems/Islands.h"
#include "libGenome/gnFASSource.h"

using namespace std;
using namespace genome;
using namespace mems;

int main( int argc, char* argv[] )
{
	IntervalList iv_list;
	if( argc != 3 )
	{
		cerr << "Usage: <input interval file> <output base name>";
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
	for( uint lcbI = 0; lcbI < iv_list.size(); lcbI++ )
	{
		gnAlignedSequences gnas;
		iv_list[lcbI].GetAlignedSequences( gnas, iv_list.seq_table );
		stringstream lcb_filename;
		lcb_filename << base_name << ".lcb_" << lcbI;
		ofstream out_file( lcb_filename.str().c_str() );
		if( !out_file.is_open() )
		{
			cerr << "Error opening \"" << lcb_filename.str() << "\"\n";
			return -2;
		}
		// write a multi-FastA
		gnSequence gns;
		for( uint seqI = 0; seqI < gnas.sequences.size(); seqI++ )
		{
			stringstream seq_name;
			seq_name << seqI;
//			seq_name << "(" << iv_list[lcbI].Start(seqI) << "-" << iv_list[lcbI].Start(seqI) + iv_list[lcbI].Length(seqI) << ")";
			gns += gnas.sequences[seqI];
			gns.setContigName( gns.contigListSize()-1, seq_name.str() );
		}
		gnFASSource::Write( gns, out_file, false, false );
	}

	return 0;
}

