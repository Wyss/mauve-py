/*******************************************************************************
 * $Id: getAlignmentWindows.cpp,v 1.2 2004/02/28 00:01:31 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/IntervalList.h"
#include "libMems/Islands.h"
#include "libGenome/gnFASSource.h"
#include "libMems/GappedAlignment.h"
#include "libMems/Interval.h"
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string/erase.hpp>

using namespace std;
using namespace genome;
using namespace mems;

void print_usage( const char* pname ){
	cerr << "Usage: " << pname << " <XMFA alignment> <window length> <window shift amount> <base output filename>\n";
}

int main( int argc, const char* argv[] ){
	if( argc <= 0 ){
		print_usage( "getAlignmentWindows" );
		return -1;
	}
	if( argc != 5 ){
		print_usage( argv[0] );
		return -1;
	}
	
	string alignment_fname = argv[1];
	int64 window_length = atol( argv[2] );
	int64 shift_length = atol( argv[3] );
	string output_basename = argv[4];
	
	ifstream alignment_in;
	alignment_in.open( alignment_fname.c_str() );
	if( !alignment_in.is_open() ){
		cerr << "Error opening " << alignment_fname << endl;
		return -1;
	}
	
	IntervalList aligned_ivs;
	aligned_ivs.ReadStandardAlignment( alignment_in );
	cout << "Read " << aligned_ivs[0].SeqCount() << " sequences with " << aligned_ivs.size() << " aligned intervals from " << alignment_fname << endl;
	cout.flush();
	MatchList mlist;
	mlist.seq_filename = aligned_ivs.seq_filename;
	if( mlist.seq_filename.size() > 0 )
		LoadSequences(mlist, &cout);
	else if( aligned_ivs.size() == 1 )
	{
		mlist.seq_filename.resize( aligned_ivs[0].SeqCount() );	
		mlist.seq_table.resize( aligned_ivs[0].SeqCount() );
		std::vector< mems::AbstractMatch* > matches;
		aligned_ivs[0].StealMatches(matches);
		std::vector< string > seqs = mems::GetAlignment( *((mems::GappedAlignment*)matches[0]), mlist.seq_table );
		for( size_t seqI = 0; seqI < mlist.seq_table.size(); ++seqI )
		{
			boost::algorithm::erase_all( seqs[seqI], std::string("-") );
			mlist.seq_table[seqI] = new gnSequence( seqs[seqI] );
		}
		aligned_ivs[0].SetMatches( matches );
	}else{
		cerr << "Error, source sequence file references not given\n";
	}
	// for each interval, extract sliding windows and write them to Multi-FastA files
	for( uint ivI = 0; ivI < aligned_ivs.size(); ivI++ )
	{
		vector< string > alignment;
		GetAlignment( aligned_ivs[ivI], mlist.seq_table, alignment );
		Interval& iv = aligned_ivs[ivI];
		stringstream ivnum;
		ivnum << ivI;
		boost::filesystem::path base_path = output_basename;
		boost::filesystem::create_directory( base_path );
		boost::filesystem::path iv_path = output_basename;
		iv_path /= "interval_" + ivnum.str();
		boost::filesystem::create_directory( iv_path );
		for( gnSeqI window_leftend = 0; window_leftend < iv.AlignmentLength(); window_leftend += shift_length )
		{
			gnSeqI cur_window_size = window_leftend + window_length < iv.AlignmentLength() ? window_length : iv.AlignmentLength() - window_leftend;

			stringstream window_filename;
			window_filename << "window_" << window_leftend << "_to_" << window_leftend + cur_window_size - 1 << ".mfa";
			boost::filesystem::path window_path = iv_path;
			window_path /= window_filename.str();
			ofstream out_file( window_path.string().c_str() );
			if( !out_file.is_open() )
			{
				cerr << "Error opening \"" << window_filename.str() << "\"\n";
				return -2;
			}
			// write a multi-FastA
			gnSequence gns;
			for( uint seqI = 0; seqI < iv.SeqCount(); seqI++ )
			{
				stringstream seq_name;
				seq_name << seqI;
				gns += alignment[seqI].substr(window_leftend, cur_window_size);
				gns.setContigName( gns.contigListSize()-1, seq_name.str() );
			}
			gnFASSource::Write( gns, out_file, false, false );
			if( cur_window_size < window_length )
				break;
		}
		// now write the whole interval as a single MFA
		boost::filesystem::path lcb_path = iv_path;
		lcb_path /= "lcb.mfa";
		ofstream lcb_out( lcb_path.string().c_str() );
		if( !lcb_out.is_open() )
		{
			cerr << "Error opening " << lcb_path.string() << endl;
			return -3;
		}
		gnSequence fns;
		for( uint seqI = 0; seqI < iv.SeqCount(); seqI++ )
		{
			stringstream seq_name;
			seq_name << seqI;
			fns += alignment[seqI];
			fns.setContigName( fns.contigListSize()-1, seq_name.str() );
		}
		gnFASSource::Write( fns, lcb_out, false, false );

	}
	return 0;
}

