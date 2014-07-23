#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnSequence.h"
#include "libGenome/gnGBKSource.h"
#include "libGenome/gnStringHeader.h"

using namespace genome;
using namespace std;

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
			// add all necessary headers
			string locus_hdr = "LOCUS       " + contig_name;
			locus_hdr += "                          DNA                  CON 27-Jan-2005";
			gnStringHeader* gnsh = new gnStringHeader( "LOCUS", locus_hdr );
			contig.addHeader( 0, gnsh, 0 );
			gnGBKSource::Write( contig, contig_name+".gbk" );
		}
//		gnRAWSource::Write( seq, argv[2] );
	}catch( gnException& gne ){
		cerr << gne << endl;
		return -1;
	}
	return 0;
}
