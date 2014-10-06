#include "libMems/IntervalList.h"
#include "libGenome/gnFASSource.h"
#include <sstream>

using namespace std;
using namespace genome;
using namespace mems;

void extractSubAlignment( IntervalList& iv_list, IntervalList& sub_list, uint seqI, gnSeqI lend, gnSeqI length )
{
	// find the relevant interval
	uint ivI = 0;
	for( ; ivI < iv_list.size(); ivI++ )
		if( iv_list[ivI].LeftEnd(seqI) <= lend &&
			lend < iv_list[ivI].LeftEnd(seqI) + iv_list[ivI].Length(seqI) )
			break;

	// we've now got the starting interval, crop appropriately
	gnSeqI crop_left_amt = lend - iv_list[ivI].LeftEnd(seqI);
	Interval iv(iv_list[ivI]);
	iv.CropLeft(crop_left_amt, seqI);
	gnSeqI crop_right_amt = length < iv.Length(seqI) ?  iv.Length(seqI) - length : 0;
	iv.CropRight(crop_right_amt, seqI);
	iv.CalculateOffset();
	sub_list.push_back(iv);
}

int main( int argc, char* argv[] )
{
	if( argc != 4 )
	{
		cerr << "Usage: extractSubAlignment <XMFA alignment input> <Multi-FastA base name> <sub-alignment spec file>\n";
		cerr << "where subalignment spec file is tab delimited text of the form:\n";
		cerr << "<genome id>\t<left end>\t<length>\n";
	}

	string alignment_infilename = argv[1];
	string alignment_outfilename = argv[2];
	string spec_filename = argv[3];

	ifstream alignment_infile( alignment_infilename.c_str() );
	if( !alignment_infile.is_open() )
	{
		cerr << "Error opening \"" << alignment_infilename << "\"\n";
		return -1;
	}
	
	IntervalList iv_list, iv_sublist;
	iv_list.ReadStandardAlignment( alignment_infile );

	ifstream spec_infile( spec_filename.c_str() );
	if( !spec_infile.is_open() )
	{
		cerr << "Error opening \"" << spec_filename << "\"\n";
		return -1;
	}
	vector< gnSequence* > seq_table( iv_list.seq_filename.size(), new gnSequence() );
	size_t ivI = 0;
	string cur_line;
	while( getline( spec_infile, cur_line ) )
	{
		stringstream line_str( cur_line );
		uint seqI;
		int64 lend;
		gnSeqI length;
		if( !(line_str >> seqI) )
			break;
		if( !(line_str >> lend) )
			break;
		if( !(line_str >> length) )
			break;
		extractSubAlignment( iv_list, iv_sublist, seqI, lend, length );

		gnAlignedSequences gnas;
		iv_sublist[0].GetAlignedSequences( gnas, seq_table );
		stringstream ss;
		ss << alignment_outfilename << ".interval_" << ivI;
		ofstream out_file( ss.str().c_str() );
		if( !out_file.is_open() )
		{
			cerr << "Error opening \"" << ss.str() << "\"\n";
			return -1;
		}
		gnSequence mfa;
		for( uint seqI = 0; seqI < seq_table.size(); seqI++ )
		{
			mfa += gnas.sequences[seqI];
			stringstream cname;
			cname << seqI << "(" << iv_sublist[0].Start(seqI) << ":" << iv_sublist[0].Start(seqI) + iv_sublist[0].Length(seqI) << ")";
			mfa.setContigName( mfa.contigListLength() - 1, cname.str() );
		}
		gnFASSource::Write( mfa, out_file, false, false );
		iv_sublist.clear();
		ivI++;
	}
}