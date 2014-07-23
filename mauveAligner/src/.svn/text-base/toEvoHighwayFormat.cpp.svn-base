#include "libMems/IntervalList.h"
#include "libMems/Aligner.h"
#include <fstream>
#include <string>
#include <vector>

using namespace std;
using namespace genome;
using namespace mems;

// find the chromosome that a given coordinate belongs to
int getChromosome( vector< int64 >& chr_lens, int64 pos )
{
	int chrI = 0;
	for( ; chrI < chr_lens.size(); chrI++ )
		if( chr_lens[chrI] > pos )
			break;
	return chrI;
}

// convert a number to a four letter base 26 number
string getAlphabetID( uint chromo_counter )
{
	string rval = "aaaa";
	int charI = 3;
	while( charI > 0 && chromo_counter > 0 )
	{
		int rem1 = chromo_counter % 26;
		chromo_counter /= 26;
		rval[charI--] = (char)(rem1 + 97);
	}
	return rval;
}

int main( int argc, char* argv[] )
{
	// 
	if( argc < 4 )
	{
		cerr << "Usage: toEvoHighwayFormat <Mauve Alignment> <reference genome id> <genome 1 chr lengths>...<genome N chr lengths>\n";
		return -1;
	}
	ifstream aln_file( argv[1] );
	if( !aln_file.is_open() )
	{
		cerr << "Error opening \"" << argv[1] << "\"\n";
		return -1;
	}
	uint ref_id = atoi( argv[2] );

	vector< vector< int64 > > chr_lens;
	for( uint genomeI = 3; genomeI < argc; genomeI++ )
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
		if( seqI == ref_id )
			continue;
		uint leftmost_lcb = 0;
		for( ; leftmost_lcb < adjacencies.size(); leftmost_lcb++ )
			if( adjacencies[ leftmost_lcb ].left_adjacency[seqI] == -1 )
				break;
		uint adjI = leftmost_lcb;
		uint cur_chromosome = 0;
		uint chromo_counter = 0;

		while( adjI != -1 && adjI != -2 && adjI < adjacencies.size() )
		{
			if( absolut(adjacencies[adjI].left_end[seqI]) > chr_lens[seqI][cur_chromosome] )
			{
				cur_chromosome++;
				chromo_counter = 0;
			}

			// write out a row for an evo highway synteny block
			// write ref name
			cout << iv_list.seq_filename[ref_id];
			// write ref chromosome
			int ref_chr = getChromosome( chr_lens[ref_id], absolut(adjacencies[adjI].left_end[ref_id]) );
			cout << '\t' << ref_chr + 1;

			// write ref interval
			if( ref_chr > 0 )
			{
				cout << '\t' << absolut(adjacencies[ adjI ].left_end[ref_id]) - chr_lens[ref_id][ref_chr - 1];
				cout << '\t' << absolut(adjacencies[ adjI ].right_end[ref_id]) - chr_lens[ref_id][ref_chr - 1];
			}else{
				cout << '\t' << absolut(adjacencies[ adjI ].left_end[ref_id]);
				cout << '\t' << absolut(adjacencies[ adjI ].right_end[ref_id]);
			}

			// write species chromosome
			cout << '\t' << cur_chromosome + 1;
			cout << getAlphabetID( chromo_counter );
			// write species interval
			if( cur_chromosome > 0 )
			{
				cout << '\t' << absolut(adjacencies[ adjI ].left_end[seqI]) - chr_lens[seqI][cur_chromosome - 1];
				cout << '\t' << absolut(adjacencies[ adjI ].right_end[seqI]) - chr_lens[seqI][cur_chromosome - 1];
			}else{
				cout << '\t' << absolut(adjacencies[ adjI ].left_end[seqI]);
				cout << '\t' << absolut(adjacencies[ adjI ].right_end[seqI]);
			}
			// write strand
			cout << '\t';
			if( adjacencies[ adjI ].left_end[ref_id] > 0 && adjacencies[ adjI ].left_end[seqI] < 0 ||
				adjacencies[ adjI ].left_end[ref_id] < 0 && adjacencies[ adjI ].left_end[seqI] > 0 )
				cout << "-";
			cout << 1;
			// write target name
			cout << '\t' << iv_list.seq_filename[seqI];
			// write lcb id
			cout << '\t' << adjacencies[adjI].lcb_id + 1 << endl;
			adjI = adjacencies[adjI].right_adjacency[seqI];
			chromo_counter++;
		}
	}
}catch(gnException& gne){
	cerr << gne << endl;
}
}
