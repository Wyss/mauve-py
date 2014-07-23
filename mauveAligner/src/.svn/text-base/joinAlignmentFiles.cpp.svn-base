#include "libMems/IntervalList.h"
#include <fstream>
#include <vector>
#include <sstream>
#include "libMems/SlotAllocator.h"
#include "libMems/Match.h"
#include "libMems/GappedAlignment.h"

using namespace std;
using namespace genome;
using namespace mems;

int main( int argc, char* argv[] )
{
	if( argc != 4 )
	{
		cerr << "joinAlignments <mauve .mln base name> <number of files> <mauve output file>\n";
		return -1;
	}
	string base_name = argv[1];
	stringstream aln_count_str(argv[2]);
	string out_fname = argv[3];
	uint aln_count;
	aln_count_str >> aln_count;
	cerr << "aln_count is: " << aln_count << endl;
	cerr << "fix this trash code\n";
	throw "shit";
/*
try{
	SlotAllocator< Match >& sa = SlotAllocator< Match >::GetSlotAllocator();
	IntervalList all_iv_list;
	for( uint alnI = 1; alnI <= aln_count; alnI++ )
	{
		IntervalList cur_iv_list;
		try{
			stringstream aln_fname;
			aln_fname << base_name << alnI << ".mln";
			ifstream cur_aln_file( aln_fname.str().c_str() );
			if( !cur_aln_file.is_open() )
			{
				cerr << "Couldn't open: \"" << aln_fname.str() << "\"\n";
				return -1;
			}
			cur_iv_list.ReadList( cur_aln_file );
			// hack: trim out all gapped alignments
			for( uint ivI = 0; ivI < cur_iv_list.size(); ivI++ )
			{
				Interval& cur_iv = cur_iv_list[ivI];
				vector<AbstractMatch*> new_matches;
				for( uint mI = 0; mI < cur_iv.matches.size(); mI++ )
				{
					GappedAlignment* ga = dynamic_cast<GappedAlignment*>(cur_iv.matches[mI]);
					if( ga == NULL )
					{
						if( mI < 5 || mI > cur_iv.matches.size() - 5 )
							new_matches.push_back( cur_iv.matches[mI] );
						else
							sa.Free(static_cast<Match*>(cur_iv.matches[mI]));
						continue;
					}
					delete ga;
				}
				cur_iv.matches = new_matches;
				cur_iv.CalculateOffset();
			}
		}catch(gnException& gne){
			// try reading the .alignment file instead of the .mln
			stringstream aln_fname;
			aln_fname << base_name << alnI << ".alignment";
			ifstream cur_aln_file( aln_fname.str().c_str() );
			if( !cur_aln_file.is_open() )
			{
				cerr << "Couldn't open: \"" << aln_fname.str() << "\"\n";
				return -1;
			}
			cur_iv_list.ReadStandardAlignment( cur_aln_file );
			for( uint ivI = 0; ivI < cur_iv_list.size(); ivI++ )
			{
				cout << ((GappedAlignment*)cur_iv_list[ivI].matches[0])->Start(0) << endl;
			}
		}
		if( alnI == 0 )
		{
			all_iv_list = cur_iv_list;
		}else{
			all_iv_list.insert( all_iv_list.end(), cur_iv_list.begin(), cur_iv_list.end() );
		}
		// progress update
		if( (alnI*100)/aln_count != ((alnI*100)-1)/aln_count ){
			cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bRead " << (alnI*100)/aln_count << "% of data";
			cout.flush();
		}
	}
	cout << endl << "Writing output\n";
	ofstream out_file( out_fname.c_str() );
	if( !out_file.is_open() )
	{
		cerr << "Error opening \"" << out_fname << "\"\n";
		return -2;
	}
	all_iv_list.WriteList( out_file );
}catch( gnException& gne )
{
        cerr << gne << endl;
}
*/
	return 0;
}
