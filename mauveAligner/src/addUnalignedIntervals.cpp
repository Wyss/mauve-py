#include "libMems/Interval.h"
#include "libMems/Islands.h"

using namespace std;
using namespace genome;
using namespace mems;

int main( int argc, char* argv[] )
{
	IntervalList iv_list;
	if( argc != 3 )
	{
		cerr << "Usage: <input interval file> <output interval file>\n";
		return -1;
	}
	ifstream in_file( argv[1] );
	if( !in_file.is_open() )
	{
		cerr << "Error opening \"argv[1]\"\n";
		return -1;
	}
	iv_list.ReadStandardAlignment( in_file );
	LoadSequences(iv_list, NULL);
	addUnalignedIntervals( iv_list );
	ofstream out_file( argv[2] );
	if( !out_file.is_open() )
	{
		cerr << "Error opening \"argv[2]\"\n";
		return -2;
	}
	iv_list.WriteStandardAlignment( out_file );
	return 0;
}
