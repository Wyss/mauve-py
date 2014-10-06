// coordinateTranslate
// (c) Aaron Darling 2011
// Licensed under the GPL

#include <libMems/IntervalList.h>
#include <fstream>
#include <libMems/CompactGappedAlignment.h>
#include <libMems/MatchList.h>

using namespace mems;
using namespace std;
using namespace genome;

int main( int argc, char* argv[] ){
	if(argc != 3){
		cerr << "Usage: coordinateTranslate <XMFA alignment> <alignment coordinate file>\n";
		cerr << "Alignment coordinate file should be structured into 2 tab-delimited columns: <block ID> <column>\n";
		cerr << "Output will be the nearest aligned position for each genome in the block, with 0 entries for genomes undefined in the block\n";
		return -1;
	}
	ifstream in_aln( argv[1] );
	if(!in_aln.is_open() ){
		cerr << "Error opening alignment file \"" << argv[1] << "\"\n";
		return -2;
	}
	IntervalList iv_list;
	iv_list.ReadStandardAlignment(in_aln);
	LoadSequences( iv_list, NULL );

	ifstream in_coords( argv[2] );
	if(!in_coords.is_open() ){
		cerr << "Error opening coordinate file \"" << argv[2] << "\"\n";
		return -2;
	}
	int block_id;
	while( in_coords >> block_id ){
		int block_col;
		in_coords >> block_col;
		std::vector<gnSeqI> pos;
		std::vector<bool> column;
		iv_list[block_id].GetColumn( block_col, pos, column );
		for(int i=0; i<pos.size(); i++){
			if(i>0) cout << "\t";
			cout << (column[i] ? pos[i] : 0); 
		}
		cout << "\n";
	}
	return 0;
}


