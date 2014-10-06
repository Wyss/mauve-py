#include "libMems/IntervalList.h"
#include "libGenome/gnFASSource.h"
#include "libMems/GappedAlignment.h"
#include <algorithm>

using namespace std;
using namespace genome;
using namespace mems;

/**
 * program to extract source sequences from an alignment
 */
int main( int argc, char* argv[] ){
	
	if( argc < 3 ){
		cout << "Sometimes you've got an alignment but you just can't seem to find the sequences that went into it." << endl;
		cout << "unalign <input alignment xmfa> <output Multi-FastA>\n";
		return -1;
	}
	
	string input_fname = argv[ 1 ];
	string output_fname = argv[ 2 ];

	ifstream alignment_in;
	alignment_in.open( input_fname.c_str() );
	if( !alignment_in.is_open() ){
		cerr << "Error opening " << input_fname << endl;
		return -1;
	}
	
	ofstream mfa_out;
	mfa_out.open( output_fname.c_str() );
	if( !mfa_out.is_open() ){
		cerr << "Error opening " << output_fname << endl;
		return -1;
	}
	
try{
	IntervalList ivs;
	cerr << "Reading " << input_fname << endl;
	ivs.ReadStandardAlignment( alignment_in );
	alignment_in.close();
	if( ivs.size() == 0 ){
		cerr << "Error! The alignment doesn't contain any intervals!\n";
		return -1;
	}
	cerr << "Successfully read " << input_fname << endl;
	cerr << "Removing gaps...\n";
	uint seq_count = ivs[ 0 ].SeqCount();
	gnSequence output_seq;
	for( uint seqI = 0; seqI < seq_count; seqI++ ){
		gnSequence cur_seq;
		AbstractMatchStartComparator<Interval> ivcomp(seqI);
		sort( ivs.begin(), ivs.end(), ivcomp );
		for( uint ivI = 0; ivI < ivs.size(); ivI++ ){
			const vector< AbstractMatch* >& matches = ivs[ivI].GetMatches();
			const vector< string >& alignment = GetAlignment(*((GappedAlignment*)matches[0]), vector<gnSequence*>(seq_count) );
			cur_seq += alignment[seqI];
			if(ivs[ivI].LeftEnd(seqI)<0)	cur_seq.setReverseComplement(true, cur_seq.contigListLength()-1);
		}
		string strseq = cur_seq.ToString();
		// strip gaps
		string gapless_seq;
		for( string::size_type charI = 0; charI < cur_seq.size(); charI++ ){
			if( strseq[ charI ] != '-' )
				gapless_seq += strseq[ charI ];
		}
		output_seq += gapless_seq;
		if(ivs.seq_filename.size()>0){
			output_seq.setContigName(seqI,ivs.seq_filename[seqI]);
			gnSequence file_seq;
			file_seq += gapless_seq;
			gnFASSource::Write( file_seq, ivs.seq_filename[seqI] );
		}
	}
	cerr << "Writing " << output_fname << endl;
	gnFASSource::Write( output_seq, mfa_out );
	
	
}catch( gnException& gne ){
	cerr << gne << endl;
	return -2;
}catch( exception& e ){
	cerr << e.what() << endl;
	return -3;
}catch( char const* c ){
	cerr << c << endl;
	return -4;
}

}
