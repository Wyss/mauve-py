#include "libMems/IntervalList.h"
#include "libMems/Aligner.h"
#include <fstream>
#include <string>
#include <vector>

using namespace std;
using namespace genome;
using namespace mems;


int main( int argc, char* argv[] )
{
	// 
	if( argc < 4 )
	{
		cerr << "Usage: toGrimmFormat <Mauve Alignment> <genome 1 chr lengths>...<genome N chr lengths>\n";
		return -1;
	}
	ifstream aln_file( argv[1] );
	if( !aln_file.is_open() )
	{
		cerr << "Error opening \"" << argv[1] << "\"\n";
		return -1;
	}
	vector< vector< int64 > > chr_lens;
	for( uint genomeI = 2; genomeI < argc; genomeI++ )
	{
		ifstream cur_file( argv[genomeI] );
		if( !cur_file.is_open() )
		{
			cerr << "Error opening \"" << argv[genomeI] << "\"\n";
			return -2;
		}
		int64 cur_len = 0;
		vector< int64 > len_vector;
		while( cur_file >> cur_len )
		{
			if( len_vector.size() > 0 )
				len_vector.push_back( cur_len + len_vector[ len_vector.size() - 1 ] );
			else
				len_vector.push_back( cur_len );
		}
		chr_lens.push_back( len_vector );
		cerr << "Read " << argv[genomeI] << ", " << len_vector.size() << " chromosomes covering " << len_vector[len_vector.size()-1] << " nt " << endl;
	}
try{
	IntervalList iv_list;
	iv_list.ReadList( aln_file );
	cerr << "Read " << argv[1] << endl;
	vector< int64 > weights = vector< int64 >( iv_list.size(), 1 );
	vector< LCB > adjacencies;
	cerr << "computeLCBAdjacencies\n";
	computeLCBAdjacencies_v2( iv_list, weights, adjacencies );
	uint seq_count = iv_list.seq_filename.size();
	for( uint seqI = 0; seqI < seq_count; seqI++ )
	{
		cerr << "Analyzing seq " << seqI << endl;
		cout << ">" << iv_list.seq_filename[seqI] << endl;
		uint leftmost_lcb = 0;
		for( ; leftmost_lcb < adjacencies.size(); leftmost_lcb++ )
			if( adjacencies[ leftmost_lcb ].left_adjacency[seqI] == -1 )
				break;
		uint adjI = leftmost_lcb;
		uint cur_chromosome = 0;
		while( adjI != -1 && adjI != -2 && adjI < adjacencies.size() )
		{
			if( absolut(adjacencies[ adjI ].left_end[seqI]) > chr_lens[seqI][cur_chromosome] )
			{
				cout << " $\n";
				cur_chromosome++;
			}else if( adjI != leftmost_lcb )
				cout << " ";
			if( adjacencies[ adjI ].left_end[seqI] < 0 )
				cout << "-";
			cout << adjacencies[ adjI ].lcb_id + 1;
			adjI = adjacencies[ adjI ].right_adjacency[seqI];
		}
		cout << " $" << endl;
	}
}catch(gnException& gne){
	cerr << gne << endl;
}
}
