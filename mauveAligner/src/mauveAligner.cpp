/*******************************************************************************
 * $Id: memsApp.cpp,v 1.49 2004/04/23 00:18:45 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "mauveAligner.h"
#include "getopt.h"
#include <sstream>
#include <stdexcept>
#include "libGenome/gnSequence.h"
#include "libMems/Matrix.h"
#include "libMems/NumericMatrix.h"
#include "libMems/DNAFileSML.h"
#include "libMems/MemHash.h"
#include "libMems/MaskedMemHash.h"
#include "libMems/Aligner.h"
#include "libMems/MatchList.h"
#include "libMems/RepeatHash.h"
#include "libMems/Interval.h"
#include "libMems/IntervalList.h"
#include "libMems/gnAlignedSequences.h"
#include "libMems/Islands.h"
#include "libMems/MuscleInterface.h"
#include "libMems/DistanceMatrix.h"

#include "boost/filesystem/operations.hpp"

using namespace std;
using namespace genome;
using namespace mems;

class MLDeleter {
public:
	MLDeleter( MatchList& ml ) : mlist( ml ) {}
	~MLDeleter(){ mlist.Clear(); }
private:
	MatchList& mlist;
};

#define NELEMS(a) ( sizeof( a ) / sizeof( *a ) )

int main( int argc, char* argv[] ){
#if	WIN32
// Multi-tasking does not work well in CPU-bound
// console apps running under Win32.
// Reducing the process priority allows GUI apps
// to run responsively in parallel.
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	return doAlignment(argc, argv);
}

/**
 * This application uses libMems to produce full scale multiple
 * genomic alignments.  First the command line is parsed to get the names of data files
 * and user specified options.  Next each sequence and its corresponding sorted mer list
 * are loaded.  If the sorted mer list fails to load a new one is created.  
 * If it is necessary to find matches in the sequences instead of loading them, each 
 * sequence and SML are added to a MemHash which searches for exact matches.  
 * Then LCBs are found if the user requested it.  Finally, either the MatchList or the
 * LCB list is written to disk.
 */
int doAlignment( int argc, char* argv[] ){
try{
	if( argc <= 0 ){
		print_usage( "mauveAligner" );
		return -1;
	}
	if( argc == 1 ){
		print_usage( argv[0] );
		return -1;
	}

	// set the Muscle path
	MuscleInterface& mi = MuscleInterface::getMuscleInterface();
	mi.ParseMusclePath( argv[0] );

	//
	// definitions of the variables that can be set by the user on the command line:
	//
	vector<string> seq_files;
	vector<string> sml_files;
	vector<gnSequence*> seq_table;
	vector<DNAFileSML*> sml_table;
	uint seed_size = 0;	// Use default settings
	int seed_rank = 0;
	boolean recursive = true;
	boolean lcb_extension = true;
	boolean gapped_alignment = true;
	boolean create_LCBs = true;
	boolean calculate_coverage = false;
	int64 LCB_size = -1;
	string output_file = "";
	boolean read_matches = false;
	boolean read_lcbs = false;
	boolean find_repeats = false;
	boolean print_stats = false;
	boolean eliminate_overlaps = false;
	boolean nway_filter = false;
	boolean collinear_genomes = false;
	string match_input_file = "";
	string lcb_stats_file = "";
	string island_file = "";
	string lcb_file = "";
	string tree_filename = "";
	string coverage_list_file = "";
	boolean output_alignment = false;
	string alignment_output_dir = "";
	string alignment_output_format = "";
	string alignment_output_file = "";
	string match_log = "";
	string offset_log = "";
	string merge_log = "";
	// island related
	uint island_size = 0;
	uint island_break_min = 0;
	// backbone related
	uint backbone_size = 0;
	uint max_backbone_gap = 0;
	int64 min_r_gap_length = -1;
	string backbone_file = "";
	boolean output_backbone = false;
	// for parallelization of LCB alignment
	vector< int > realign_lcbs;
	string muscle_args = "";
	string gapped_aligner;

	string permutation_filename;
	int64 permutation_weight = -1;

	boolean lcb_match_input_format = false;
	int opt_max_extension_iters = -1;

	uint seqI;
	boolean print_version = false;
	int max_gapped_alignment_length = -1;
	
	ostream* detail_list_out = NULL;	/**< output stream for detail list */

	//
	// parse command line with gnu getopt
	//
	int opt;
	int config_opt;
	int ac = argc;
	char** av = argv;
	// 'm' mer size
	// 'r' recursive
	const char* short_args= "";
	enum opt_names{
		opt_mums,
		opt_no_recursion,
		opt_no_lcb_extension,
		opt_no_gapped_alignment,
		opt_seed_size,
		opt_seed_type,
		opt_weight,
		opt_output,
		opt_eliminate_overlaps,
		opt_n_way_filter,
		opt_match_input,
		opt_lcb_input,
		opt_output_alignment,
		opt_id_matrix,
		opt_island_size,
		opt_island_output,
		opt_island_break_min,
		opt_backbone_size,
		opt_max_backbone_gap,
		opt_backbone_output,
		opt_coverage_output,
		opt_repeats,
		opt_gapped_aligner,
		opt_max_gapped_aligner_length,
		opt_min_recursive_gap_length,
		opt_output_guide_tree,
		opt_alignment_output_dir,
		opt_alignment_output_format,
		opt_match_log,
		opt_offset_log,
		opt_merge_match_log,
		opt_version,
		opt_scratch_path,
		opt_realign_lcb,
		opt_id_matrix_input,
		opt_collinear,
		opt_muscle_args,
		opt_permutation_matrix_output,
		opt_permutation_matrix_min_weight,
		opt_lcb_match_input,
		opt_max_extension_iterations,
	};
	struct option long_opts[] = {
		{"mums", no_argument, &config_opt, opt_mums},
		{"no-recursion", no_argument, &config_opt, opt_no_recursion},
		{"no-lcb-extension", no_argument, &config_opt, opt_no_lcb_extension},
		{"no-gapped-alignment", no_argument, &config_opt, opt_no_gapped_alignment},
		{"seed-size", required_argument, &config_opt, opt_seed_size},
		{"seed-type", required_argument, &config_opt, opt_seed_type},
		{"weight", required_argument, &config_opt, opt_weight},
		{"output", required_argument, &config_opt, opt_output},
		{"eliminate-overlaps", no_argument, &config_opt, opt_eliminate_overlaps},
		{"n-way-filter", no_argument, &config_opt, opt_n_way_filter},
		{"match-input", required_argument, &config_opt, opt_match_input},
		{"lcb-input", required_argument, &config_opt, opt_lcb_input},
		{"output-alignment", optional_argument, &config_opt, opt_output_alignment},
		{"id-matrix", optional_argument, &config_opt, opt_id_matrix},
		{"island-size", required_argument, &config_opt, opt_island_size},
		{"island-output", required_argument, &config_opt, opt_island_output},
		{"island-break-min", required_argument, &config_opt, opt_island_break_min},
		{"backbone-size", required_argument, &config_opt, opt_backbone_size},
		{"max-backbone-gap", required_argument, &config_opt, opt_max_backbone_gap},
		{"backbone-output", optional_argument, &config_opt, opt_backbone_output},
		{"coverage-output", optional_argument, &config_opt, opt_coverage_output},
		{"repeats", no_argument, &config_opt, opt_repeats},
		{"max-gapped-aligner-length", required_argument, &config_opt, opt_max_gapped_aligner_length},
		{"min-recursive-gap-length", required_argument, &config_opt, opt_min_recursive_gap_length},
		{"output-guide-tree", required_argument, &config_opt, opt_output_guide_tree},
		{"alignment-output-dir", required_argument, &config_opt, opt_alignment_output_dir},
		{"alignment-output-format", required_argument, &config_opt, opt_alignment_output_format},
		{"match-log", required_argument, &config_opt, opt_match_log},
		{"offset-log", required_argument, &config_opt, opt_offset_log},
		{"merge-match-log", required_argument, &config_opt, opt_merge_match_log},
		{"version", no_argument, &config_opt, opt_version},
		{"scratch-path", required_argument, &config_opt, opt_scratch_path},
		{"realign-lcb", required_argument, &config_opt, opt_realign_lcb},
		{"id-matrix-input", required_argument, &config_opt, opt_id_matrix_input},
		{"collinear", no_argument, &config_opt, opt_collinear},
		{"muscle-args", required_argument, &config_opt, opt_muscle_args},
		{"permutation-matrix-output", required_argument, &config_opt, opt_permutation_matrix_output},
		{"permutation-matrix-min-weight", required_argument, &config_opt, opt_permutation_matrix_min_weight},
		{"lcb-match-input", no_argument, &config_opt, opt_lcb_match_input},
		{"max-extension-iterations", required_argument, &config_opt, opt_max_extension_iterations},

		{0, 0, 0, 0}	// for correct termination of option list
						// getopt_long can segfault without this
	};

	int indexptr;
	while( (opt = getopt_long( ac, av, short_args, long_opts, &indexptr )) != EOF ){
		switch( opt ){
			case 0:
				switch(config_opt){
					case opt_mums:
						create_LCBs = false;
						break;
					case opt_no_recursion:
						recursive = false;
						break;
					case opt_no_lcb_extension:
						lcb_extension = false;
						break;
					case opt_no_gapped_alignment:
						gapped_alignment = false;
						break;
					case opt_seed_size:
						seed_size = atoi( optarg );
						break;
					case opt_seed_type:
						if( strcmp( "solid", optarg ) == 0 )
							seed_rank = SOLID_SEED;
						else if( strcmp( "coding", optarg ) == 0 )
							seed_rank = CODING_SEED;
						else if( strcmp( "spaced", optarg ) == 0 )
							seed_rank = 0;
						else if( strcmp( "spaced1", optarg ) == 0 )
							seed_rank = 1;
						else if( strcmp( "spaced2", optarg ) == 0 )
							seed_rank = 2;
						else
							cerr << "Warning: --seed-type parameter not understood.  Using default spaced seeds\n";
						break;
					case opt_weight:
						LCB_size = atol( optarg );
						break;
					case opt_output:
						output_file = optarg;
						break;
					case opt_eliminate_overlaps:
						eliminate_overlaps = true;
						break;
					case opt_n_way_filter:
						nway_filter = true;
						break;
					case opt_match_input:
						read_matches = true;
						match_input_file = optarg;
						break;
					case opt_lcb_input:
						lcb_file = optarg;
						read_lcbs = true;
						break;
					case opt_output_alignment:
						output_alignment = true;
						if( optarg != NULL )
							alignment_output_file = optarg;
						break;
					case opt_id_matrix:
						break;
					case opt_island_size:
						island_size = atoi( optarg );
						break;
					case opt_island_output:
						island_file = optarg;
						break;
					case opt_island_break_min:
						island_break_min = atoi( optarg );
						break;
					case opt_backbone_size:
						backbone_size = atoi( optarg );
						break;
					case opt_max_backbone_gap:
						max_backbone_gap = atoi( optarg );
						break;
					case opt_backbone_output:
						backbone_file = optarg;
						output_backbone = true;
						break;
					case opt_coverage_output:
						if( optarg != NULL )
							coverage_list_file = optarg;
						calculate_coverage = true;
						break;
					case opt_repeats:
						find_repeats = true;
						break;
					case opt_gapped_aligner:
						gapped_aligner = optarg;
						break;
					case opt_max_gapped_aligner_length:
						max_gapped_alignment_length = atoi( optarg );
						break;
					case opt_min_recursive_gap_length:
						min_r_gap_length = atol( optarg );
						break;
					case opt_output_guide_tree:
						tree_filename = optarg;
						break;
					case opt_alignment_output_dir:
						alignment_output_dir = optarg;
						break;
					case opt_alignment_output_format:
						alignment_output_format = optarg;
						break;
					case opt_match_log:
						match_log = optarg;
						break;
					case opt_offset_log:
						offset_log = optarg;
						break;
					case opt_merge_match_log:
						merge_log = optarg;
						break;
					case opt_version:
						print_version = true;
						break;
					case opt_scratch_path:
						FileSML::registerTempPath( optarg );
						break;
					case opt_realign_lcb:
						realign_lcbs.push_back( atoi( optarg ) );
						break;
					case opt_id_matrix_input:
					case opt_collinear:
						collinear_genomes = true;
						break;
					case opt_muscle_args:
						muscle_args = optarg;
						mi.SetExtraMuscleArguments( muscle_args );
						break;
					case opt_permutation_matrix_output:
						permutation_filename = optarg;
						break;
					case opt_permutation_matrix_min_weight:
						permutation_weight = atol(optarg);
						break;
					case opt_lcb_match_input:
						lcb_match_input_format = true;
						break;
					case opt_max_extension_iterations:
						opt_max_extension_iters = atoi(optarg);
						break;
					default:
						print_usage( argv[0] );
						return -1;
				}
				break;
			default:
				print_usage( argv[0] );
				return -1;
		}
	}
	// now read in the seq and sml file names from av
	boolean seq_name_arg = true;
	for( int optI = optind; optI < argc; optI++ ){
		if( seq_name_arg )
			seq_files.push_back( av[ optI ] );
		else
			sml_files.push_back( av[ optI ] );
		seq_name_arg = !seq_name_arg;
	}
	
	// print the version if the user requested it
	if( print_version ){
		cerr << "mauveAligner " << " build date " << __DATE__ << " at " << __TIME__ << endl;
	}


	//
	// check validity of command line option combinations
	//
	if( ( island_size != 0 && island_file == "" ) || ( island_size == 0 && island_file != "" ) ){
		cerr << "Error: Both --island-output and --island-size must be specified to generate islands\n";
		return -1;
	}

	if( (alignment_output_dir == "" && alignment_output_format != "") || 
		(alignment_output_dir != "" && alignment_output_format == "") ){
		cerr << "Error: Both --alignment-output-dir and --alignment-output-format must be specified in order to generate alignment output in a custom format\n";
		return -1;
	}
	
	if( alignment_output_format != "" ){
		if( !gnAlignedSequences::isSupportedFormat( alignment_output_format ) ){
			cerr << "Error:  " << alignment_output_format << " is not a supported alignment format.\n";
			return -1;
		}
	}

	if( find_repeats ){
		if( create_LCBs || read_matches || read_lcbs || calculate_coverage || 
		    island_file != "" || island_size != 0 || recursive || lcb_stats_file != "" ){
			cerr << "A paramater has been specified that is incompatible with repeat list generation\n";
			return -1;
		}
	}

	//
	// done parsing and checking command line options
	// Start doing the work
	//

	MatchList match_list;
	MLDeleter deleter( match_list );
	
	if( seq_files.size() == 1 && sml_files.size() == 0 ){
		LoadMFASequences( match_list, seq_files[0], &cout);
		if( find_repeats || ( !read_lcbs && !read_matches ) )
			match_list.CreateMemorySMLs(seed_size, &cout, seed_rank);
	}else if( seq_files.size() != sml_files.size() ){
		cerr << "Error: Each sequence file must have a corresponding SML file specified.\n";
		return -1;
	}else{
		match_list.seq_filename = seq_files;
		match_list.sml_filename = sml_files;
		LoadSequences( match_list, &cout );
		if( find_repeats || !read_matches || ( !read_lcbs && !read_matches ) )
			match_list.LoadSMLs( seed_size, &cout, seed_rank );
	}

	ostream* match_out;
	if( output_file != "" ){
		ofstream* match_out_file = new ofstream( output_file.c_str() );
		if( !match_out_file->is_open() ){
			cerr << "Error opening " << output_file << endl;
			return -2;
		}
		match_out = match_out_file;
	}else
		match_out = &cout;
	
	// search for repetitive regions
	if( find_repeats ){
		RepeatHash repeat_finder;
		repeat_finder.LogProgress( &cout );
		repeat_finder.FindMatches( match_list );
		WriteList( match_list, *match_out );
		match_out->flush();
		return 0;
	}
	
	// read matches if the user requested it
	if( read_matches ){
		ifstream match_in( match_input_file.c_str() );
		if( !match_in.is_open() ){
			cerr << "Error opening " << match_input_file << endl;
			return -2;
		}
		if( !lcb_match_input_format )
		{
			try{
				ReadList( match_list, match_in );
			}catch( gnException& gne ){
				cerr << "Error reading " << match_input_file << "\nPossibly corrupt file or invalid file format\n";
				return -2;
			}
		}else{
			IntervalList m_iv_list;
			m_iv_list.ReadList( match_in );
			for( int ivI = 0; ivI < m_iv_list.size(); ivI++ ){
				for( int mI = 0; mI < m_iv_list[ivI].GetMatches().size(); mI++ ){
					Match* m = dynamic_cast< Match* >(m_iv_list[ivI].GetMatches()[mI]);
					if( m != NULL && m->Multiplicity() > 1)
						match_list.push_back(m->Copy());
				}
			}
		}
		if( seq_files.size() > 1 )
			match_list.seq_filename = seq_files;
		else if( match_list.seq_table.size() == 0 )
			// fill seq_table with empty sequences
			for( seqI = 0; seqI < match_list.seq_filename.size(); seqI++ )
				match_list.seq_table.push_back( new gnSequence() );
	}else if ( !read_lcbs ){
		// get full subset matches
		MaskedMemHash match_finder;

		if( nway_filter ){
			// only find the n-way matches
			uint64 nway_mask = 1;
			nway_mask <<= match_list.seq_table.size();
			nway_mask--;
			match_finder.SetMask( nway_mask );
		}
		match_finder.LogProgress( &cout );
		fstream match_log_out;
		if( match_log != "" ){
			match_log_out.open( match_log.c_str(), ios::in | ios::out );
			if( !match_log_out.is_open() ){
				cerr << "Error opening " << match_log << endl;
				return -1;
			}
			match_finder.SetMatchLog( &match_log_out );
			// append to whatever's already in the file
			match_log_out.seekg( 0, ios::end );
		}
		fstream offset_log_out;
		vector< gnSeqI > offset_start;
		for( seqI = 0; seqI < match_list.seq_table.size(); seqI++ )
			offset_start.push_back( 0 );
		
		if( offset_log != "" ){
			offset_log_out.open( offset_log.c_str(), ios::in | ios::out );
			if( !offset_log_out.is_open() ){
				cerr << "Error opening " << offset_log << endl;
				return -1;
			}
			match_finder.SetOffsetLog( &offset_log_out );
			string last_line;
			string cur_line;
			while( getline( offset_log_out, cur_line ) ){
				last_line = cur_line;
			}
			if( last_line != "" ){
				stringstream cur_off_stream( last_line );
				for( seqI = 0; seqI < match_list.seq_table.size(); seqI++ )
					cur_off_stream >> offset_start[ seqI ];
			}
			offset_log_out.clear();
		}
		ifstream merge_log_in;
		if( merge_log != "" ){
			merge_log_in.open( merge_log.c_str() );
			if( !merge_log_in.is_open() ){
				cerr << "Error opening " << merge_log << endl;
				return -1;
			}

			for( seqI = 0; seqI < match_list.seq_table.size(); seqI++ ){
				if( !match_finder.AddSequence( match_list.sml_table[ seqI ], match_list.seq_table[ seqI ] ) ){
					ErrorMsg( "Error adding " + match_list.seq_filename[seqI] + "\n");
					return -1;
				}
			}
			match_finder.LoadFile( merge_log_in );
			match_finder.GetMatchList( match_list );
		}else{
			match_finder.FindMatchesFromPosition( match_list, offset_start );
		}
		match_log_out.close();
		offset_log_out.close();
		match_finder.Clear();
	}
		

	// write out a match list if the user doesn't want LCBs
	if( !create_LCBs && !read_lcbs){
		if( eliminate_overlaps ){
			EliminateOverlaps( match_list );
		}

		if( nway_filter ){
			match_list.MultiplicityFilter( match_list.seq_table.size() );
		}
		
		WriteList( match_list, *match_out );
		match_out->flush();
		
		// output a guide tree or a coverage list if necessary
		// beware that selecting the nway filter above will cause the guide tree
		// and coverage lists to be incorrect
		vector< pair< uint64, uint64 > > coverage_list;
		if( tree_filename != "" || calculate_coverage ){
			// only count each base pair once!
			if( !eliminate_overlaps )
				EliminateOverlaps( match_list );
		}

		if( tree_filename != "" ){
			NumericMatrix< double > distance;
			DistanceMatrix( match_list.seq_table.size(), coverage_list, distance );
			MuscleInterface& mi = MuscleInterface::getMuscleInterface();
			if( tree_filename == "" )
				tree_filename = CreateTempFileName("guide_tree");
			mi.CreateTree( distance, tree_filename );
		}

		return 0;
	}
	
	// check whether the input sequences were masked to eliminate excess NNNNNs
	for( seqI = 0; seqI < match_list.sml_table.size(); seqI++ ){
		FileSML* cur_sml = dynamic_cast< FileSML* >(match_list.sml_table[ seqI ]);
		if( cur_sml != NULL ){
			const vector< int64 >& seq_coords = cur_sml->getUsedCoordinates();
			if( seq_coords.size() > 0 ){
				transposeMatches( match_list, seqI, seq_coords );
			}
		}
	}
	
	// at this point any SortedMerLists used to identify the initial set of MUMs
	// are no longer necessary.  Free them
	for( uint smlI = 0; smlI < match_list.sml_table.size(); smlI++ ){
		match_list.sml_table[ smlI ]->Clear();
		delete match_list.sml_table[ smlI ];
	}
	match_list.sml_table.clear();
	
	// Align the sequences if necessary
	if( LCB_size < 0 ){
		// calculate a default LCB weight, 3 times the mer size times the seq. count
		if( seed_size <= 0 )	
			seed_size = MatchList::GetDefaultMerSize( match_list.seq_table );
		LCB_size = seed_size * 3 * match_list.seq_table.size();
	}else{
		// adjust the LCB weight for the number of sequences being aligned
		LCB_size *= match_list.seq_table.size();
	}

	// check that LCB_size can be set appropriately
	if( create_LCBs && LCB_size < 0) {
		cerr << "A minimum LCB size greater than 0 must be specified in order to create LCBs.\n";
		return -1;
	}

	// hack to communicate that the genomes are collinear
	if( collinear_genomes )
		LCB_size = -1;
	
	Aligner aligner( match_list.seq_table.size() );

	if( min_r_gap_length >= 0 ){
		aligner.SetMinRecursionGapLength( min_r_gap_length );
	}

	aligner.SetGappedAligner( MuscleInterface::getMuscleInterface() );
	if( max_gapped_alignment_length != -1 )
		aligner.SetMaxGappedAlignmentLength( max_gapped_alignment_length );
	
	if( permutation_weight != -1 && permutation_filename == "" )
		cerr << "A permutation output file must be specified to generate signed permutations\n";
	if( permutation_weight == -1 && permutation_filename != "" )
		permutation_weight = LCB_size;
	if( permutation_weight != -1 )
	{
		permutation_weight *= match_list.seq_table.size();
		aligner.SetPermutationOutput( permutation_filename, permutation_weight );
	}
	if( opt_max_extension_iters != -1 )
	{
		aligner.SetMaxExtensionIterations(opt_max_extension_iters);
	}

	IntervalList interval_list;
	interval_list.seq_table = match_list.seq_table;
	interval_list.seq_filename = match_list.seq_filename;
	if( lcb_file == "" ){

		try{
			aligner.align( match_list, interval_list, 0, LCB_size, recursive, lcb_extension, gapped_alignment, tree_filename );
		}catch( gnException& gne ){
			cerr << gne << endl;
		}
		interval_list.WriteList( *match_out );
		match_out->flush();

	}else if( read_lcbs ){
		ifstream lcb_input( lcb_file.c_str() );
		if( !lcb_input.is_open() ){
			cerr << "Error opening " << lcb_file << endl;
			return -2;
		}
		try{

			interval_list.seq_table = match_list.seq_table;
			interval_list.seq_filename = match_list.seq_filename;
			interval_list.ReadList( lcb_input );
//			addUnalignedIntervals( interval_list );
		}catch( gnException& gne ){
			cerr << gne << endl;
			cerr << "Error reading " << lcb_file << "\nPossibly corrupt file or invalid file format\n";
			return -2;
		}
	}
	if( realign_lcbs.size() > 0 ){
		// set up a new IntervalList
		IntervalList realigned_intervals;
		realigned_intervals.seq_table = interval_list.seq_table;
		realigned_intervals.seq_filename = interval_list.seq_filename;
		for( int realignI = 0; realignI < realign_lcbs.size(); realignI++ ){
			// extract a match list from the interval list for this LCB
			Interval& iv = interval_list[ realignI ];
			// clear any matches from the current match_list
			match_list.clear();
			for( int matchI = 0; matchI < iv.GetMatches().size(); matchI++ ){
				AbstractMatch* m = iv.GetMatches()[ matchI ];
				Match* match = dynamic_cast< Match* >( m );
				if( match != NULL && m->Multiplicity() > 1)
					match_list.push_back( match->Copy() );
			}
			aligner.align( match_list, realigned_intervals, 0, LCB_size, recursive, false, gapped_alignment, tree_filename );
		}
		
		// once all intervals have been realigned reset the interval_list
		interval_list = realigned_intervals;
	}
	
	if( output_alignment ){
		if( !gapped_alignment )
			addUnalignedIntervals( interval_list );
		if( alignment_output_file == "" || alignment_output_file == "-" ){
			interval_list.WriteStandardAlignment( cout );
		}else{
			ofstream align_out( alignment_output_file.c_str() );
			if( !align_out.is_open() ){
				cerr << "Error opening " << alignment_output_file << endl;
				return -1;
			}
			interval_list.WriteStandardAlignment( align_out );
			align_out.close();
		}
	}
	uint lcbI;
	
	// output alignments in another format if the user asked for it
	if( alignment_output_dir != "" ){
		boost::filesystem::path output_dir = alignment_output_dir;
		boost::filesystem::create_directory( output_dir );

		for( lcbI = 0; lcbI < interval_list.size(); lcbI++ ){
			gnAlignedSequences gnas;
			interval_list[ lcbI ].GetAlignedSequences( gnas, match_list.seq_table );
			ostringstream oss;
			oss << "lcb_" << lcbI << ".txt";
			boost::filesystem::path outtie = output_dir / oss.str();
			ofstream alignment_lcb_out( outtie.string().c_str(), ios::trunc );
			if( !alignment_lcb_out.is_open() ){
				cerr << "Error opening " << oss.str() << endl;
				return -1;
			}
			gnas.output( alignment_output_format, alignment_lcb_out );
		}
	}

	//
	// output an identity matrix if requested
	//
	if( print_stats ){
		ostream* stats_out;
		if( lcb_stats_file == "" || lcb_stats_file == "-" ){
			stats_out = &cout;
		}else{
			ofstream* stats_out_file = new ofstream( lcb_stats_file.c_str() );
			if( !stats_out_file->is_open() ){
				cerr << "Error opening " << lcb_stats_file << endl;
				return -1;
			}
			stats_out = stats_out_file;
		}
		NumericMatrix< double > identity;
		IdentityMatrix( interval_list, identity );
		identity.print( *stats_out );
		if( lcb_stats_file == "" || lcb_stats_file == "-" ){
			delete stats_out;
		}
	}

	//
	// output backbone if it was requested
	//
	if( output_backbone ){
		ostream* backbone_out;
		if( backbone_file != "" ){
			ofstream* backbone_out_file = new ofstream( backbone_file.c_str() );
			if( !backbone_out_file->is_open() ){
				cerr << "Error opening " << backbone_file << endl;
				return -1;
			}
			backbone_out = backbone_out_file;
		}else
			backbone_out = &cout;

		vector< GappedAlignment > backbone_data;
		simpleFindBackbone( interval_list, backbone_size, max_backbone_gap, backbone_data );
		outputBackbone( backbone_data, *backbone_out );
		if( backbone_file != "" ){
			delete backbone_out;
		}
	}

	//
	// output islands if they were requested
	//
	if( island_file != "" ){
		ostream* island_out;
		if( island_file == "-" )
			island_out = &cout;
		else{
			ofstream* island_out_file = new ofstream( island_file.c_str() );
			if( !island_out_file->is_open() ){
				cerr << "Error opening " << island_file << endl;
				return -1;
			}
			island_out = island_out_file;
		}
		simpleFindIslands( interval_list, island_size, *island_out );
		findIslandsBetweenLCBs( interval_list, island_size, *island_out );

		if( island_file != "-" ){
			delete island_out;
		}
	}
	match_list.clear();	// bad.  leaks memory.
}catch( gnException& gne ) {
	cerr << "Unhandled gnException: " << gne << endl;
	return -10;
}catch( exception& e ) {
	cerr << "Unhandled exception: " << e.what() << endl;
	return -11;
}catch( char* message ){
	cerr << "Unhandled exception: " << message << endl;
	return -12;
}catch(...){
	cerr << "Unknown exception occurred.\n";
	return -13;
}

	return 0;
}

void print_usage( const char* pname ){
	cerr << "Usage:" << endl;
	cerr << pname << " [options] <seq1 filename> <sml1 filename> ... "
		<< " <seqN filename> <smlN filename>" << endl;
	cerr << "Options:" << endl;
	cerr << "\t    --output=<file> Output file name.  Prints to screen by default" << endl;
	cerr << "\t    --mums Find MUMs only, do not attempt to determine locally collinear blocks (LCBs)\n";
	cerr << "\t    --no-recursion Don't perform recursive anchor identification (implies --no-gapped-alignment)" << endl;
	cerr << "\t    --no-lcb-extension If determining LCBs, don't attempt to extend the LCBs\n";
	cerr << "\t    --seed-size=<number> Initial seed match size, default is log_2( average seq. length )" << endl;
	cerr << "\t    --max-extension-iterations=<number> Limit LCB extensions to this number of attempts, default is 4\n";
	cerr << "\t    --eliminate-inclusions Eliminate linked inclusions in subset matches.\n";
	cerr << "\t    --weight=<number> Minimum LCB weight in base pairs per sequence" << endl;
	cerr << "\t    --match-input=<file> Use specified match file instead of searching for matches\n";
	cerr << "\t    --lcb-match-input  Indicates that the match input file contains matches that have been clustered into LCBs\n";
	cerr << "\t    --lcb-input=<file> Use specified lcb file instead of constructing LCBs (skips LCB generation)\n";
	cerr << "\t    --scratch-path=<path>  For large genomes, use a directory for storage of temporary data.  Should be given two or more times to with different paths.\n";
	cerr << "\t    --id-matrix=<file> Generate LCB stats and write them to the specified file\n";
	cerr << "\t    --island-size=<number> Find islands larger than the given number\n";
	cerr << "\t    --island-output=<file> Output islands the given file (requires --island-size)\n";
	cerr << "\t    --backbone-size=<number> Find stretches of backbone longer than the given number of b.p.\n";
	cerr << "\t    --max-backbone-gap=<number> Allow backbone to be interrupted by gaps up to this length in b.p.\n";
	cerr << "\t    --backbone-output=<file> Output islands the given file (requires --island-size)\n";
	cerr << "\t    --coverage-output=<file> Output a coverage list to the specified file (- for stdout)\n";
	cerr << "\t    --repeats Generates a repeat map.  Only one sequence can be specified\n";
	cerr << "\t    --output-guide-tree=<file> Write out a guide tree to the designated file\n";
	cerr << "\t    --collinear Assume that input sequences are collinear--they have no rearrangements\n";
	cerr << "\nGapped alignment controls:\n";
	cerr << "\t    --no-gapped-alignment Don't perform a gapped alignment\n";
	cerr << "\t    --max-gapped-aligner-length=<number> Maximum number of base pairs to attempt aligning with the gapped aligner\n";
	cerr << "\t    --min-recursive-gap-length=<number> Minimum size of gaps that Mauve will perform recursive MUM anchoring on (Default is 200)\n";
	cerr << "\nSigned permutation matrix options:\n";
	cerr << "\t    --permutation-matrix-output=<file> Write out the LCBs as a signed permutation matrix to the given file\n";
	cerr << "\t    --permutation-matrix-min-weight=<number> A permutation matrix will be written for every set of LCBs with weight between this value and the value of --weight\n";
	cerr << "\nAlignment output options:\n";
	cerr << "\t    --alignment-output-dir=<directory> Outputs a set of alignment files (one per LCB) to a given directory\n";
	cerr << "\t    --alignment-output-format=<directory> Selects the output format for --alignment-output-dir\n";
	cerr << "\t    --output-alignment=<file> Write out an XMFA format alignment to the designated file\n";
	cerr << endl;
	
	const vector< string >& formats = gnAlignedSequences::getSupportedFormats();
	cerr << "Supported alignment output formats are: ";
	for( int formatI = 0; formatI < formats.size(); formatI++ ){
		if( formatI > 0 )
			cerr << ", ";
		cerr << formats[ formatI ];
	}
	cerr << endl;
	cerr << endl;
}

