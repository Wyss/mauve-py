#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnSequence.h"
#include "libGenome/gnRAWSource.h"

int main( int argc, char* argv[] ){

	if( argc != 3 ){
		cout << argv[0] << " <input sequence> <output file>\n";
	}
	gnSequence seq;
	try{
		seq.LoadSource( argv[1] );
		cout << argv[1] << " has " << seq.contigListLength() << " contigs\n";
		for( int contigI = 0; contigI < seq.contigListLength(); contigI++ ){
			gnSequence contig = seq.contig( contigI );
			string contig_name = seq.contigName( contigI );
			cout << "contig " << contig_name << " has " << contig.length() << "b.p.\n";
			gnRAWSource::Write( contig, contig_name+".raw" );
		}
	}catch( gnException& gne ){
		cerr << gne << endl;
		return -1;
	}
	return 0;
}
