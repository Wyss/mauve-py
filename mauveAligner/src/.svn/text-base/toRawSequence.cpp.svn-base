#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnSequence.h"
#include "libGenome/gnRAWSource.h"

using namespace std;
using namespace genome;


int main( int argc, char* argv[] ){

	if( argc != 3 ){
		cout << argv[0] << " <input sequence> <output file>\n";
	}
	gnSequence seq;
	try{
		seq.LoadSource( argv[1] );
		cout << argv[1] << " is " << seq.length() << "b.p.\n";
		gnRAWSource::Write( seq, argv[2] );
	}catch( gnException& gne ){
		cerr << gne << endl;
		return -1;
	}
	return 0;
}
