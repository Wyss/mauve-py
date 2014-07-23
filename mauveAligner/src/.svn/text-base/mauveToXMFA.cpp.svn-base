#include "libMems/IntervalList.h"
#include <fstream>
#include <string>
#include "libMems/MatchList.h"

using namespace std;
using namespace genome;
using namespace mems;

int main( int argc, char* argv[] )
{
	if( argc != 3 )
	{
		cerr << "Usage: mauveToXMFA <Mauve Alignment input> <XMFA output>\n";
		return -1;
	}
	ifstream mauve_file( argv[1] );
	if( !mauve_file.is_open() )
	{
		cerr << "Error opening \"" << argv[1] << "\"\n";
		return -2;
	}
	ofstream xmfa_file( argv[2] );
	if( !xmfa_file.is_open() )
	{
		cerr << "Error opening \"" << argv[2] << "\"\n";
		return -3;
	}

	IntervalList iv_list;
	iv_list.ReadList( mauve_file );
	LoadSequences(iv_list, &cout);
	iv_list.WriteStandardAlignment( xmfa_file );
}

