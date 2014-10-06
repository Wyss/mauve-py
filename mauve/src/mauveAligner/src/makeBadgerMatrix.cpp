#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include "libMems/IntervalList.h"
#include "libMems/MatchList.h"
#include "libMems/Aligner.h"

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
	if( argc != 4 )
	{
		cerr << "Usage: makeBadgerMatrix <input xmfa> <output badger file> <LCB coordinate file>\n";
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

	ofstream coord_out;
	coord_out.open( argv[3] );
	if( !coord_out.is_open() ){
		cerr << "Error writing to " << argv[3] << endl;
		return -2;
	}

	try{
		IntervalList input_ivs;
		input_ivs.ReadStandardAlignment( aln_in );
		aln_in.close();

		vector< pair< Interval*, uint > > labeled_ivs( input_ivs.size() );
		for( size_t ivI = 0; ivI < input_ivs.size(); ivI++ )
			labeled_ivs[ivI] = make_pair( &input_ivs[ivI], ivI );

		// write out block boundaries
		for( uint seqI = 0; seqI < input_ivs.seq_filename.size(); ++seqI )
		{
			if(seqI > 0) coord_out << '\t';
			coord_out << "seq" << seqI << "_leftend\tseq" << seqI << "_rightend";
		}
		coord_out << endl;
		for( size_t ivI = 0; ivI < input_ivs.size(); ivI++ )
		{
			if( labeled_ivs[ivI].first->Multiplicity() == 1 )
				continue;
			for( uint seqI = 0; seqI < input_ivs.seq_filename.size(); ++seqI )
			{
				if(seqI > 0) coord_out << '\t';
				string sign = labeled_ivs[ivI].first->Start(seqI) < 0 ? "-" : "";
				coord_out << sign << labeled_ivs[ivI].first->LeftEnd(seqI) << '\t' << sign << labeled_ivs[ivI].first->RightEnd(seqI);
			}
			coord_out << endl;
		}

		for( uint seqI = 0; seqI < input_ivs.seq_filename.size(); ++seqI )
		{
			badger_out << input_ivs.seq_filename[seqI];
			livComp lc(seqI);
			std::sort( labeled_ivs.begin(), labeled_ivs.end(), lc );
			for( size_t ivI = 0; ivI < labeled_ivs.size(); ivI++ )
			{
				if( labeled_ivs[ivI].first->LeftEnd(seqI) == NO_MATCH )
					continue;
				if( labeled_ivs[ivI].first->Multiplicity() == 1 )
					continue;
				int fs = labeled_ivs[ivI].first->FirstStart();
				const char* dir = labeled_ivs[ivI].first->Orientation(seqI) == labeled_ivs[ivI].first->Orientation(fs) ? "" : "-";
				badger_out << "," << dir << labeled_ivs[ivI].second + 1;
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

