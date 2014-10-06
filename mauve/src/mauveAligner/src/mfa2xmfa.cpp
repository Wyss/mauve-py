#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnFASSource.h"
#include "libGenome/gnSequence.h"
#include <algorithm>
#include <fstream>

using namespace std;
using namespace genome;

void print_usage( const char* pname ){
	cerr << "Usage: " << pname << " <MFA alignment input> <XMFA alignment output> [Unaligned FastA output]\n";
}

int main( int argc, char* argv[] ) {
	if( argc < 3 ){
		if( argc == 0 )
			print_usage( "mfa2xmfa" );
		else
			print_usage( argv[0] );
		return -1;
	}
	
	gnSequence mfa_seq;
	string mfa_name = argv[1];
	try{
		mfa_seq.LoadSource( mfa_name );
	}catch( gnException& gne ){
		cerr << gne << endl;
		return -1;
	}
	ofstream xmfa_out( argv[2] );
	if( !xmfa_out.is_open() ){
		cerr << "Error opening " << argv[2] << endl;
		return -1;
	}

	// if unaligned mfa output is desired, create it now
	string mfa_output_name;
	ofstream mfa_out;
	if(argc >= 4)
	{
		mfa_output_name = argv[3];
		mfa_out.open( mfa_output_name.c_str() );
		if( !mfa_out.is_open() ){
			cerr << "Error opening " << argv[3] << endl;
			return -1;
		}

		gnSequence unaligned;
		for( size_t seqI = 0; seqI < mfa_seq.contigListSize(); seqI++ )
		{
			string cur_seq = mfa_seq.contig(seqI).ToString();
			string::iterator striter = std::remove( cur_seq.begin(), cur_seq.end(), '-' );
			cur_seq.resize( striter - cur_seq.begin() );
			unaligned += cur_seq;
			unaligned.setContigName( seqI, mfa_seq.contigName(seqI) );
		}
		gnFASSource::Write( unaligned, mfa_out, false, false );

		// create xmfa header if unaligned seq is to be written
		xmfa_out << "#FormatVersion Mauve1\n";
	}

	unsigned int seq_count = mfa_seq.contigListSize();
	// find the max length alignment entry and add gaps
	// to the ends of shorter entries for consistency
	gnSeqI max_length = 0;
	for( uint seqI = 0; seqI < seq_count; seqI++ ){
		if( mfa_seq.contigLength( seqI ) > max_length )
			max_length = mfa_seq.contigLength( seqI );
	}


	// count the number of base pairs in each sequence
	vector< gnSeqI > seq_lens( seq_count, 0 );
	for( uint seqI = 0; seqI < seq_count; seqI++ ){
		string cur_seq;
		mfa_seq.contig( seqI ).ToString(cur_seq);
		for( gnSeqI baseI = 0; baseI < cur_seq.length(); baseI++ ){
			if( cur_seq[ baseI ] != '-' )
				seq_lens[ seqI ]++;
		}
		// fill in xmfa header details if unaligned seq is to be written
		if(mfa_output_name.size() > 0)
		{
			xmfa_out << "#Sequence" << seqI + 1 << "File\t" << mfa_output_name << endl;
			xmfa_out << "#Sequence" << seqI + 1 << "Entry\t" << seqI + 1 << endl;
			xmfa_out << "#Sequence" << seqI + 1 << "Format\tFastA\n";
		}
	}

	// write xmfa body	
	for( uint seqI = 0; seqI < seq_count; seqI++ ){
		string cur_seq;
		mfa_seq.contig( seqI ).ToString(cur_seq);
		// if the alignment entry is shorter than the max length
		// then add gaps to the end for consistency
		if( cur_seq.length() < max_length )
			cur_seq += string( max_length - cur_seq.length(), '-' );

		xmfa_out << "> " << seqI + 1 << ":";
		xmfa_out << 1 << "-" << seq_lens[ seqI ] << " + " << mfa_seq.contigName( seqI );
		xmfa_out << endl;
		gnSeqI cur_pos = 0;
		for( ; cur_pos < cur_seq.length(); cur_pos += 80 ){
			gnSeqI cur_len = cur_pos + 80 < cur_seq.length() ? 80 : cur_seq.length() - cur_pos;
			xmfa_out.write( cur_seq.data() + cur_pos, cur_len );
			xmfa_out << endl;
		}
	}
		
	xmfa_out << "=" << endl;
	
}
