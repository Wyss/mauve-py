#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include "libMems/IntervalList.h"

using namespace std;
using namespace genome;
using namespace mems;

class livComp {
public:
	livComp( uint seq ){ m_seq = seq; };
	bool operator()( const pair< Interval*, uint >& a, const pair< Interval*, uint >& b )
	{
		return a.first->LeftEnd(m_seq) < b.first->LeftEnd(m_seq);
	}
protected:
	uint m_seq;
};

int main( int argc, char* argv[] )
{
	if( argc < 3 )
	{
		cerr << "Usage: makeBadgerMatrix <input xmfa> <output badger file>\n";
		return -1;
	}
	ifstream aln_in;
	aln_in.open( argv[1] );
	if( !aln_in.is_open() ){
		cerr << "Error opening " << argv[1] << endl;
		return -1;
	}
	ofstream badger_out;
	badger_out.open( argv[2] );
	if( !badger_out.is_open() ){
		cerr << "Error writing to " << argv[2] << endl;
		return -1;
	}

	try{
		IntervalList input_ivs;
		input_ivs.ReadStandardAlignment( aln_in );
		aln_in.close();

		vector< pair< Interval*, uint > > labeled_ivs;
		for( size_t ivI = 0; ivI < input_ivs.size(); ivI++ )
		{
			if( input_ivs[ivI].Multiplicity() != input_ivs.seq_filename.size() )
				continue;	// not an N-way block
			labeled_ivs.push_back( make_pair( &input_ivs[ivI], ivI ) );
		}
		for( uint seqI = 0; seqI < input_ivs.seq_filename.size(); ++seqI )
		{
			badger_out << input_ivs.seq_filename[seqI];
			livComp lc(seqI);
			std::sort( labeled_ivs.begin(), labeled_ivs.end(), lc );
			if( seqI == 0 )
			{
				for( size_t ivI = 0; ivI < labeled_ivs.size(); ivI++ )
				{
					labeled_ivs[ivI].second = ivI + 1;
					if( labeled_ivs[ivI].first->Orientation(seqI) == AbstractMatch::reverse )
						labeled_ivs[ivI].first->Invert();
				}
			}
			vector< size_t > other( labeled_ivs.size() * 2 + 2 );
			for( size_t ivI = 0; ivI < labeled_ivs.size(); ivI++ )
			{
				if(labeled_ivs[ivI].first->Orientation(seqI) == AbstractMatch::forward)
				{
					other[ivI*2+1] = absolut(labeled_ivs[ivI].second)*2 - 1;
					other[ivI*2+2] = absolut(labeled_ivs[ivI].second)*2;
				}else{
					other[ivI*2+1] = absolut(labeled_ivs[ivI].second)*2;
					other[ivI*2+2] = absolut(labeled_ivs[ivI].second)*2 - 1;
				}
			}
			for( size_t ivI = 0; ivI < other.size(); ivI++ )
			{
				badger_out << "," << other[ivI];
			}
			badger_out << "\nstandard";
			for( size_t ivI = 0; ivI < labeled_ivs.size(); ivI++ )
			{
				badger_out << "," << (labeled_ivs[ivI].first->Orientation(seqI) == AbstractMatch::reverse? "-" : "") << labeled_ivs[ivI].second;
			}

			badger_out << endl;
		}

	}catch( gnException& gne ){
		cerr << gne << endl;
		return -1;
	}catch( exception& e ){
		cerr << e.what() << endl;
		return -2;
	}catch( char const* c ){
		cerr << c << endl;
		return -3;
	}catch(...){
		cerr << "Unhandled exception" << endl;
		return -4;
	}
}

