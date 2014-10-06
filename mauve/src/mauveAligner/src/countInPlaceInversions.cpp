#include "libMems/IntervalList.h"
#include "libMems/Aligner.h"
#include <fstream>
#include <string>
#include <vector>
#include <utility>

using namespace std;
using namespace genome;
using namespace mems;

int main( int argc, char* argv[] )
{
	if( argc != 2 )
	{
		cerr << "Usage: countInPlaceInversions <Mauve Alignment>\n";
		return -1;
	}
	ifstream aln_file( argv[1] );
	if( !aln_file.is_open() )
	{
		cerr << "Error opening \"" << argv[1] << "\"\n";
		return -1;
	}

	IntervalList iv_list;
	iv_list.ReadList( aln_file );
	vector< int64 > weights = vector< int64 >( iv_list.size(), 1 );
	vector< LCB > adjacencies;
	computeLCBAdjacencies_v2( iv_list, weights, adjacencies );
	uint seq_count = iv_list.seq_filename.size();
	vector< pair< uint, uint > > inv_seqs;

	for( uint adjI = 0; adjI < adjacencies.size(); adjI++ )
	{
		// find in place inversions
		uint seqI = 1;
		for( ; seqI < seq_count; seqI++ )
		{
			if( adjacencies[adjI].left_adjacency[0] != adjacencies[adjI].left_adjacency[seqI] ||
				adjacencies[adjI].right_adjacency[0] != adjacencies[adjI].right_adjacency[seqI] )
				break;
		}
		if( seqI == seq_count )
		{
			// in place inversion
			// count forward
			uint forward_count = 0;
			for( seqI = 0; seqI < seq_count; seqI++ )
			{
				if( adjacencies[adjI].left_end[seqI] > 0 )
					forward_count++;
			}
			for( seqI = 0; seqI < seq_count; seqI++ )
			{
				if( forward_count * 2 > seqI && adjacencies[adjI].left_end[seqI] < 0 )
					inv_seqs.push_back( make_pair( adjI, seqI ) );
				if( forward_count * 2 < seqI && adjacencies[adjI].left_end[seqI] > 0 )
					inv_seqs.push_back( make_pair( adjI, seqI ) );
			}
		}
	}
	for( uint invI = 0; invI < inv_seqs.size(); invI++ )
	{
		cout << "In-place inversion in seq " << inv_seqs[invI].second;
		cout << "\tlend: " << adjacencies[inv_seqs[invI].first].left_end[inv_seqs[invI].second];
		cout << "\trend: " << adjacencies[inv_seqs[invI].first].right_end[inv_seqs[invI].second] << endl;
	}
}
