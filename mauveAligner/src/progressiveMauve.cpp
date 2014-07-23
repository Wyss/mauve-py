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
#include "libMems/Matrix.h"
#include "libMems/NumericMatrix.h"
#include "libGenome/gnSequence.h"
#include "libMems/DNAFileSML.h"
#include "libMems/MemHash.h"
#include "libMems/MatchList.h"
#include "libMems/Interval.h"
#include "libMems/IntervalList.h"
#include "libMems/gnAlignedSequences.h"
#include "libMems/Islands.h"
#include "libMems/MuscleInterface.h"
#include "libMems/Backbone.h"
//#include "libMems/twister.h"

#include "libMems/ProgressiveAligner.h"
#include "libMems/PairwiseMatchFinder.h"
#include "libMems/HomologyHMM/parameters.h"
#include "UniqueMatchFinder.h"

#include <boost/filesystem.hpp>

#include "libMems/Memory.h"

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

class OptionList;

class MauveOption : public option
{
public:
	MauveOption( OptionList& ol, const char* name, int has_arg, const std::string& usage_info);

	boolean set;
	std::string arg_value;
	std::string usage_info;
};


class OptionList : public vector< MauveOption* >
{
public:
	OptionList() : opt_list(NULL){};
	~OptionList()
	{
		if( opt_list != NULL )
			delete[] opt_list;
	}
	struct option* getOptions()
	{
		if( opt_list == NULL )
		{
			opt_list = new option[ this->size() + 1 ];
			int i = 0;
			for( ; i < this->size(); i++ ){
				opt_list[i] = *(*this)[i];
			}
			struct option empty = {0,0,0,0};
			opt_list[i] = empty;
		}
		return opt_list;
	}
	int config_opt;
protected:
	struct option* opt_list;
};

MauveOption::MauveOption( OptionList& ol, const char* name, int has_arg, const std::string& usage_info ) :
	set( false ),
	usage_info( usage_info )
{
	this->name = name;
	this->has_arg = has_arg;
	this->flag = &ol.config_opt;
	this->val = ol.size();
	ol.push_back(this);
}

void print_usage( const char* pname, OptionList& option_list )
{
	cerr << "progressiveMauve usage:\n\n";
	cerr << "When each genome resides in a separate file:" << endl;
	cerr << pname << " [options] <seq1 filename> ... <seqN filename>" << endl << endl;
	cerr << "When all genomes are in a single file:" << endl;
	cerr << pname << " [options] <seq filename>" << endl << endl;
	cerr << "Options:" << endl;
	for( size_t optionI = 0; optionI < option_list.size(); optionI++ )
	{
		cerr << "\t" << "--" << option_list[optionI]->name;
		cerr << (option_list[optionI]->has_arg == no_argument ? " " : "=");
		cerr << option_list[optionI]->usage_info << endl;
	}
	cerr << endl << endl;
	cerr << "Examples:\n";
	cerr << pname << " --output=my_seqs.xmfa my_genome1.gbk my_genome2.gbk my_genome3.fasta\n";
	cerr << "\nIf genomes are in a single file and have no rearrangement:\n";
	cerr << pname << " --collinear --output=my_seqs.xmfa my_genomes.fasta\n";
}

void printMatchSizes()
{
	UngappedLocalAlignment< HybridAbstractMatch<> > ula;
	UngappedLocalAlignment< SparseAbstractMatch<> > sula;
	CompactGappedAlignment<> cga;
	MatchHashEntry	mhe;
	bitset_t bitset;
	Match m;
	cerr << "sizeof(UngappedLocalAlignment< HybridAbstractMatch<> >) " << sizeof(ula) << endl;
	cerr << "sizeof(UngappedLocalAlignment< SparseAbstractMatch<> >) " << sizeof(sula) << endl;
	cerr << "sizeof(m) " << sizeof(m) << endl;
	cerr << "sizeof(CompactGappedAlignment<>) " << sizeof(cga) << endl;
	cerr << "sizeof(boost::dynamic_bitset) " << sizeof(bitset) << endl;
	cerr << "sizeof(MatchHashEntry) " << sizeof(mhe) << endl;
}

#ifndef WIN32
#include <signal.h>
#endif

/**
 * Aborts the running progressiveMauve program
 */
void terminateProgram( int sig )
{
	std::cerr << "Caught signal " << sig << std::endl;
	std::cerr << "Cleaning up and exiting!\n";
	deleteRegisteredFiles();
	std::cerr << "Temporary files deleted.\n";
	exit(sig);	
}

#ifdef WIN32
BOOL WINAPI handler(DWORD dwCtrlType)
{
	switch(dwCtrlType)
	{
	case CTRL_C_EVENT:
	case CTRL_BREAK_EVENT:
	case CTRL_CLOSE_EVENT:
	case CTRL_LOGOFF_EVENT:
	case CTRL_SHUTDOWN_EVENT:
		terminateProgram(dwCtrlType);
	default:
		break;
	}
	return true;
}
#endif

int main( int argc, char* argv[] )
{
#if	WIN32
// Multi-tasking does not work well in CPU-bound
// console apps running under Win32.
// Reducing the process priority allows GUI apps
// to run responsively in parallel. (thanks Bob Edgar!)
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
// also register a handler to clean up during abnormal shutdown
	BOOL status = SetConsoleCtrlHandler(handler, TRUE);
#else
// register a signal handler to catch errors and control-c and clean up...
	signal( SIGINT, terminateProgram );
	signal( SIGTERM, terminateProgram );
	signal( SIGSEGV, terminateProgram );
#endif
	// delete temp files at program exit!
	atexit( &deleteRegisteredFiles );

	return doAlignment(argc, argv);
}

void getPatternText( int64 seed_pattern, char pattern[65] )
{
	char pat[65] = {
		'0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0',
		'0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0',
		'0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0',
		'0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0',
		'\0'};
	int lastone = 64;
	for( int i = 63; i >= 0; i-- )
	{
		pat[i] = seed_pattern & 0x1 ? '1' : '0';
		lastone = pat[i] == '1' ? i : lastone;
		seed_pattern >>= 1;
	}
	memcpy( pattern, pat + lastone, 65 - lastone );
}

void getDefaultSmlFileNames( const vector< string >& seq_files, vector< string >& sml_files, int seed_weight, int seed_rank )
{
	int64 seed_pattern = getSeed(seed_weight, seed_rank);
	// convert seed pattern to text;
	char pattern[65];
	getPatternText(seed_pattern, pattern);
	sml_files.resize(seq_files.size());
	for( int seqI = 0; seqI < seq_files.size(); seqI++ )
		sml_files[seqI] = seq_files[seqI] + "." + pattern + ".sslist";
}

void applyBackbone( IntervalList& iv_list, string& bbcols_fname, string& bb_fname, size_t island_gap_size, double hmm_identity, double pgh, double pgu )
{
	ofstream bb_out( bb_fname.c_str() );
	backbone_list_t bb_list;
	// adapt to the GC of the sequences
	double gc_content = computeGC( iv_list.seq_table );
	std::cout << "Organisms have " << std::setprecision(3) << gc_content*100 << "% GC\n";

	Params hmm_params = getAdaptedHoxdMatrixParameters( gc_content );
	hmm_params.iGoHomologous = pgh;
	hmm_params.iGoUnrelated = pgu;
	adaptToPercentIdentity( hmm_params, hmm_identity );

	detectAndApplyBackbone(iv_list, bb_list, hmm_params);
	bb_list.clear();

	BigGapsDetector bgd( island_gap_size );
	detectBackbone( iv_list, bb_list, &bgd );

	writeBackboneSeqCoordinates( bb_list, iv_list, bb_out );
	std::vector< bb_seqentry_t > bb_seq_list;
	bb_out.close();
	std::ifstream bbseq_input( bb_fname.c_str() );
	readBackboneSeqFile( bbseq_input, bb_seq_list );

	mergeAdjacentSegments( bb_seq_list );
	addUniqueSegments( bb_seq_list );
	bbseq_input.close();
	bb_out.open(bb_fname.c_str());
	writeBackboneSeqFile( bb_out, bb_seq_list );

	ofstream bbcols_out( bbcols_fname.c_str() );
	writeBackboneColumns( bbcols_out, bb_list );
	iv_list.backbone_filename = bbcols_fname;
}

/**
 * progressive alignment.  wheee.
 */
int doAlignment( int argc, char* argv[] ){
//try{
	OptionList mauve_options;
	MauveOption opt_island_gap_size( mauve_options, "island-gap-size", required_argument, "<number> Alignment gaps above this size in nucleotides are considered to be islands [20]" );
	MauveOption opt_profile( mauve_options, "profile", required_argument, "<file> (Not yet implemented) Read an existing sequence alignment in XMFA format and align it to other sequences or alignments" );
	MauveOption opt_apply_backbone( mauve_options, "apply-backbone", required_argument, "<file> Read an existing sequence alignment in XMFA format and apply backbone statistics to it" );
	MauveOption opt_disable_backbone( mauve_options, "disable-backbone", no_argument, "Disable backbone detection" );
	MauveOption opt_mums( mauve_options, "mums", no_argument, "Find MUMs only, do not attempt to determine locally collinear blocks (LCBs)" );
	MauveOption opt_seed_weight( mauve_options, "seed-weight", required_argument, "<number> Use the specified seed weight for calculating initial anchors" );
	MauveOption opt_output( mauve_options, "output", required_argument, "<file> Output file name.  Prints to screen by default" );
	MauveOption opt_backbone_output( mauve_options, "backbone-output", required_argument, "<file> Backbone output file name (optional)." );
	MauveOption opt_match_input( mauve_options, "match-input", required_argument, "<file> Use specified match file instead of searching for matches" );
	MauveOption opt_input_id_matrix( mauve_options, "input-id-matrix", required_argument, "<file> An identity matrix describing similarity among all pairs of input sequences/alignments" );
	MauveOption opt_max_gapped_aligner_length( mauve_options, "max-gapped-aligner-length", required_argument, "<number> Maximum number of base pairs to attempt aligning with the gapped aligner" );
	MauveOption opt_input_guide_tree( mauve_options, "input-guide-tree", required_argument, "<file> A phylogenetic guide tree in NEWICK format that describes the order in which sequences will be aligned" );
	MauveOption opt_output_guide_tree( mauve_options, "output-guide-tree", required_argument, "<file> Write out the guide tree used for alignment to a file" );
	MauveOption opt_version( mauve_options, "version", no_argument, "Display software version information" );
	MauveOption opt_debug( mauve_options, "debug", no_argument, "Run in debug mode (perform internal consistency checks--very slow)" );
	MauveOption opt_scratch_path_1( mauve_options, "scratch-path-1", required_argument, "<path> Designate a path that can be used for temporary data storage.  Two or more paths should be specified." );
	MauveOption opt_scratch_path_2( mauve_options, "scratch-path-2", required_argument, "<path> Designate a path that can be used for temporary data storage.  Two or more paths should be specified." );
	MauveOption opt_collinear( mauve_options, "collinear", no_argument, "Assume that input sequences are collinear--they have no rearrangements" );
	MauveOption opt_scoring_scheme( mauve_options, "scoring-scheme", required_argument, "<ancestral|sp_ancestral|sp> Selects the anchoring score function.  Default is extant sum-of-pairs (sp)." );
	MauveOption opt_no_weight_scaling( mauve_options, "no-weight-scaling", no_argument, "Don't scale LCB weights by conservation distance and breakpoint distance" );
	MauveOption opt_max_breakpoint_distance_scale( mauve_options, "max-breakpoint-distance-scale", required_argument, "<number [0,1]> Set the maximum weight scaling by breakpoint distance.  Defaults to 0.5" );
	MauveOption opt_conservation_distance_scale( mauve_options, "conservation-distance-scale", required_argument, "<number [0,1]> Scale conservation distances by this amount.  Defaults to 0.5" );
	MauveOption opt_muscle_args( mauve_options, "muscle-args", required_argument, "<arguments in quotes> Additional command-line options for MUSCLE.  Any quotes should be escaped with a backslash" );
	MauveOption opt_skip_refinement( mauve_options, "skip-refinement", no_argument, "Do not perform iterative refinement" );
	MauveOption opt_skip_gapped_alignment( mauve_options, "skip-gapped-alignment", no_argument, "Do not perform gapped alignment" );
	MauveOption opt_bp_dist_estimate_min_score( mauve_options, "bp-dist-estimate-min-score", required_argument, "<number> Minimum LCB score for estimating pairwise breakpoint distance" );
	MauveOption opt_mem_clean( mauve_options, "mem-clean", no_argument, "Set this to true when debugging memory allocations" );
	MauveOption opt_gap_open( mauve_options, "gap-open", required_argument, "<number> Gap open penalty" );
	MauveOption opt_penalize_repeats( mauve_options, "repeat-penalty", required_argument, "<negative|zero> Sets whether the repeat scores go negative or go to zero for highly repetitive sequences.  Default is negative." );
	MauveOption opt_gap_extend( mauve_options, "gap-extend", required_argument, "<number> Gap extend penalty" );
	MauveOption opt_substitution_matrix( mauve_options, "substitution-matrix", required_argument, "<file> Nucleotide substitution matrix in NCBI format" );
	MauveOption opt_weight( mauve_options, "weight", required_argument, "<number> Minimum pairwise LCB score" );
	MauveOption opt_min_scaled_penalty( mauve_options, "min-scaled-penalty", required_argument, "<number> Minimum breakpoint penalty after scaling the penalty by expected divergence" );
	MauveOption opt_go_homologous( mauve_options, "hmm-p-go-homologous", required_argument, "<number> Probability of transitioning from the unrelated to the homologous state [0.00001]" );
	MauveOption opt_go_unrelated( mauve_options, "hmm-p-go-unrelated", required_argument, "<number> Probability of transitioning from the homologous to the unrelated state [0.000000001]" );
	MauveOption opt_hmm_identity( mauve_options, "hmm-identity", required_argument, "<number> Expected level of sequence identity among pairs of sequences, ranging between 0 and 1 [0.7]" );
	MauveOption opt_seed_family( mauve_options, "seed-family", no_argument, "Use a family of spaced seeds to improve sensitivity" );
	MauveOption opt_solid_seeds( mauve_options, "solid-seeds", no_argument, "Use solid seeds. Do not permit substitutions in anchor matches." );
	MauveOption opt_coding_seeds( mauve_options, "coding-seeds", no_argument, "Use coding pattern seeds. Useful to generate matches coding regions with 3rd codon position degeneracy." );
	MauveOption opt_disable_cache( mauve_options, "disable-cache", no_argument, "Disable recursive anchor search cacheing to workaround a crash bug" );
	MauveOption opt_recursive( mauve_options, "no-recursion", no_argument, "Disable recursive anchor search" );

	if( argc <= 0 ){
		print_usage( "mauveAligner", mauve_options );
		return -1;
	}
	if( argc == 1 ){
		print_usage( argv[0], mauve_options );
		return -1;
	}

	// default values for homology HMM transitions
	double pgh = 0.00001;
	double pgu = 0.000000001;
	double hmm_identity = 0.7;	// percent identity modeled by the HMM homologous state
	size_t island_gap_size = 20;

	// set the Muscle path
	MuscleInterface& mi = MuscleInterface::getMuscleInterface();
	mi.ParseMusclePath( argv[0] );

	// parse the options
	//
	// parse command line with gnu getopt
	//
	int opt;
	int ac = argc;
	char** av = argv;
	int indexptr;
	while( (opt = getopt_long( ac, av, "", mauve_options.getOptions(), &indexptr )) != EOF ){
		if( opt == 0 )
		{
			mauve_options[mauve_options.config_opt]->set = true;
			if( optarg != NULL )
				mauve_options[mauve_options.config_opt]->arg_value = optarg;
		}else{
			print_usage( argv[0], mauve_options );
			return -1;
		}
	}
	
	if( opt_scratch_path_1.set )
		FileSML::registerTempPath( opt_scratch_path_1.arg_value.c_str() );
	if( opt_scratch_path_2.set )
		FileSML::registerTempPath( opt_scratch_path_2.arg_value.c_str() );

	// set the random number generator to a fixed seed for repeatability
	// this should be changed if the algorithm ever depends on true pseudo-randomness
	SetTwisterSeed(37);

	if( opt_go_homologous.set )
		pgh = strtod( opt_go_homologous.arg_value.c_str(), NULL );
	if( opt_go_unrelated.set )
		pgu = strtod( opt_go_unrelated.arg_value.c_str(), NULL );
	if( opt_hmm_identity.set )
		hmm_identity = strtod( opt_hmm_identity.arg_value.c_str(), NULL );
	if( opt_island_gap_size.set )
		island_gap_size = atoi( opt_island_gap_size.arg_value.c_str() );

	// for debugging only:
	if( opt_apply_backbone.set )
	{
		IntervalList iv_list;
		ifstream in_file( opt_apply_backbone.arg_value.c_str() );
		ofstream out_file( opt_output.arg_value.c_str() );
		iv_list.ReadStandardAlignment(in_file);
		MatchList ml;
		ml.seq_filename = iv_list.seq_filename;
		if( ml.seq_filename[0] != ml.seq_filename[1] )
			LoadSequences(ml, &cout);
		else
			LoadMFASequences(ml, ml.seq_filename[0], &cout);
		iv_list.seq_table = ml.seq_table;
		string bb_fname = opt_output.arg_value + ".backbone";
		string bbcols_fname = opt_output.arg_value + ".bbcols";
		applyBackbone( iv_list, bbcols_fname, bb_fname, island_gap_size, hmm_identity, pgh, pgu );
		iv_list.WriteStandardAlignment(out_file);
		return 0;
	}

	//
	// definitions of the variables that can be set by the user on the command line:
	//
	vector<string> seq_files;
	vector<string> sml_files;
	vector<gnSequence*> seq_table;
	vector<DNAFileSML*> sml_table;
	uint mer_size = 0;	// Use default settings
	boolean create_LCBs = true;
	string output_file = "";
	string tree_filename = "";

	boolean lcb_match_input_format = false;

	uint seqI;
	
	ostream* detail_list_out = NULL;	/**< output stream for detail list */

	// now read in the seq file names from av
	boolean seq_name_arg = true;
	for( int optI = optind; optI < argc; optI++ )
		seq_files.push_back( av[ optI ] );

	// set sml_names
	for( size_t seq_fileI = 0; seq_fileI < seq_files.size(); seq_fileI++ )
		sml_files.push_back( seq_files[seq_fileI] + ".sslist" );
	
	// print the version if the user requested it
	if( opt_version.set ){
		cerr << "progressiveMauve " << " build date " << __DATE__ << " at " << __TIME__ << endl;
	}

	if( seq_files.size() == 0 )
	{
		if( !opt_version.set )
			print_usage( argv[0], mauve_options );
		return 0;
	}

	//
	// done parsing and checking command line options
	// Start doing the work
	//

	MatchList pairwise_match_list;
	if( opt_seed_weight.set )
	{
		mer_size = atoi( opt_seed_weight.arg_value.c_str() );
	}
	
	if( seq_files.size() == 1 ){
		LoadMFASequences( pairwise_match_list, seq_files[0], &cout );
		pairwise_match_list.CreateMemorySMLs( mer_size, &cout );
	}else{
		pairwise_match_list.seq_filename = seq_files;
		pairwise_match_list.sml_filename = sml_files;
		// testing: rewrite seq files in RAW format
		LoadAndCreateRawSequences( pairwise_match_list, &cout );
//		LoadSequences( pairwise_match_list, &cout );
		if(opt_solid_seeds.set)
			pairwise_match_list.LoadSMLs( mer_size, &cout, SOLID_SEED, true );
		else if(opt_coding_seeds.set)
			pairwise_match_list.LoadSMLs( mer_size, &cout, CODING_SEED );
		else
			pairwise_match_list.LoadSMLs( mer_size, &cout, CODING_SEED );
	}

	ostream* match_out;
	if( opt_output.set ){
		ofstream* match_out_file = new ofstream( opt_output.arg_value.c_str() );
		if( !match_out_file->is_open() ){
			cerr << "Error opening \"" << opt_output.val << "\"\n";
			return -2;
		}
		match_out = match_out_file;
	}else
		match_out = &cout;
	
	if(opt_mem_clean.set)
		debugging_memory = true;

	// read matches if the user requested it
	if( opt_match_input.set ){
		ifstream match_in( opt_match_input.arg_value.c_str() );
		if( !match_in.is_open() ){
			cerr << "Error opening " << opt_match_input.arg_value << endl;
			return -2;
		}
		try{
			ReadList( pairwise_match_list, match_in );
		}catch( gnException& gne ){
			cerr << gne << endl;
			cerr << "Error reading " << opt_match_input.arg_value << "\nPossibly corrupt file or invalid file format\n";
			return -2;
		}

		if( seq_files.size() > 1 )
			pairwise_match_list.seq_filename = seq_files;
		else if( pairwise_match_list.seq_table.size() == 0 )
			// fill seq_table with empty sequences
			for( seqI = 0; seqI < pairwise_match_list.seq_filename.size(); seqI++ )
				pairwise_match_list.seq_table.push_back( new gnSequence() );
	}else if( !opt_seed_family.set ){
		if( pairwise_match_list.seq_table.size() > 4 )
		{
			UniqueMatchFinder umf;
			umf.LogProgress( &cout );
			umf.FindMatches( pairwise_match_list );
			umf.Clear();
		}else{
			PairwiseMatchFinder pmf;
			pmf.LogProgress( &cout );
			pmf.FindMatches( pairwise_match_list );
			pmf.Clear();
		}
		cout << "done.\n";
	}else{
		// use an entire seed family to do the search
		if( mer_size == 0 )
		{
			size_t avg = 0;
			for( int seqI = 0; seqI < pairwise_match_list.seq_table.size(); seqI++ )
				avg += pairwise_match_list.seq_table[seqI]->length();
			avg /= pairwise_match_list.seq_table.size();
			mer_size = getDefaultSeedWeight( avg );
		}
		// search with the longest seeds first so that overlapping matches tend to get contained
		vector< pair< int, int > > length_ranks(3);
		length_ranks[0] = make_pair( getSeedLength( getSeed(mer_size, 0) ), 0 );
		length_ranks[1] = make_pair( getSeedLength( getSeed(mer_size, 1) ), 1 );
		length_ranks[2] = make_pair( getSeedLength( getSeed(mer_size, 2) ), 2 );
		std::sort( length_ranks.begin(), length_ranks.end() );

		UniqueMatchFinder umf;
		for( int seedI = 2; seedI >= 0; seedI-- )
		{
			umf.LogProgress( &cout );
			int64 seed_pattern = getSeed(mer_size, length_ranks[seedI].second );
			char pattern[65];
			getPatternText( seed_pattern, pattern );
			cout << "\nSearching with seed pattern " << pattern << "\n";
			MatchList cur_list;
			cur_list.seq_filename = pairwise_match_list.seq_filename;
			cur_list.seq_table = pairwise_match_list.seq_table;
			if( seq_files.size() == 1 )
				cur_list.CreateMemorySMLs( mer_size, &cout, length_ranks[seedI].second );
			else
			{
				getDefaultSmlFileNames( cur_list.seq_filename, cur_list.sml_filename, mer_size, length_ranks[seedI].second );
				cur_list.LoadSMLs(mer_size, &cout, length_ranks[seedI].second);
			}
			umf.FindMatches( cur_list );
			umf.ClearSequences();
			for( size_t smlI = 0; smlI < cur_list.sml_table.size(); smlI++ )
				delete cur_list.sml_table[smlI];	// free memory
			for( size_t curI = 0; curI < cur_list.size(); curI++ )
				cur_list[curI]->Free();	// free more memory!
		}
		umf.GetMatchList(pairwise_match_list);
		cout << "done\n";
		umf.Clear();
	}
	
	if( opt_mums.set )
	{
		WriteList(pairwise_match_list, *match_out);
		for( size_t seqI = 0; seqI < pairwise_match_list.seq_table.size(); seqI++ )
			delete pairwise_match_list.seq_table[seqI];	// an auto_ptr or shared_ptr could be great for this
		for( size_t seqI = 0; seqI < pairwise_match_list.sml_table.size(); seqI++ )
			delete pairwise_match_list.sml_table[seqI];
		return 0;
	}

	// check whether the input sequences were masked to eliminate excess NNNNNs
	for( seqI = 0; seqI < pairwise_match_list.sml_table.size(); seqI++ ){
		FileSML* cur_sml = dynamic_cast< FileSML* >(pairwise_match_list.sml_table[ seqI ]);
		if( cur_sml != NULL ){
			const vector< int64 >& seq_coords = cur_sml->getUsedCoordinates();
			if( seq_coords.size() > 0 ){
				transposeMatches( pairwise_match_list, seqI, seq_coords );
			}
		}
	}
	
	// free any match search memory
	SlotAllocator<MatchHashEntry>& allocator = SlotAllocator<MatchHashEntry>::GetSlotAllocator();
	allocator.Purge();
	
	ProgressiveAligner aligner( pairwise_match_list.seq_table.size() );
	if( opt_skip_gapped_alignment.set )
		aligner.setGappedAlignment(false);
	if( opt_skip_refinement.set )
		aligner.setRefinement(false);
	if( opt_debug.set )
		debug_aligner = true;

	// check that LCB_size can be set appropriately
	if( opt_weight.set )
	{
		double lcb_weight = strtod( opt_weight.arg_value.c_str(), NULL );
		if( lcb_weight < 0 )
		{
			cerr << "A minimum LCB size greater than 0 must be specified in order to create LCBs.\n";
			return -1;
		}else
			aligner.setBreakpointPenalty( lcb_weight );
	}

	if( opt_collinear.set )
		aligner.setCollinear(true);

	if( opt_max_gapped_aligner_length.set )
	{
		int64 mgal = atol( opt_max_gapped_aligner_length.arg_value.c_str() );
		aligner.SetMaxGappedAlignmentLength( mgal );
	}

	if( opt_seed_family.set )
		aligner.setUseSeedFamilies(true);

	penalize_repeats = true;
	if(opt_penalize_repeats.set && opt_penalize_repeats.arg_value == "zero")
		penalize_repeats = false;

	if( opt_scoring_scheme.set )
	{
		if( opt_scoring_scheme.arg_value == "ancestral" )
			aligner.setLcbScoringScheme(ProgressiveAligner::AncestralScoring);
		else if( opt_scoring_scheme.arg_value == "ancestral_sp" )
			aligner.setLcbScoringScheme(ProgressiveAligner::AncestralSumOfPairsScoring);
		else if( opt_scoring_scheme.arg_value == "sp" )
			aligner.setLcbScoringScheme(ProgressiveAligner::ExtantSumOfPairsScoring);
		else
		{
			cerr << "Unrecognized scoring scheme: " << opt_scoring_scheme.arg_value << endl;
			return -2;
		}
	}else	// default to extant sp
		aligner.setLcbScoringScheme(ProgressiveAligner::ExtantSumOfPairsScoring);
	if( opt_no_weight_scaling.set )
		aligner.setUseLcbWeightScaling(false);
	if( opt_max_breakpoint_distance_scale.set )
	{
		double d = strtod( opt_max_breakpoint_distance_scale.arg_value.c_str(), NULL );
		aligner.setBreakpointDistanceScale(d);
	}
	if( opt_conservation_distance_scale.set )
	{
		double d = strtod( opt_conservation_distance_scale.arg_value.c_str(), NULL );
		aligner.setConservationDistanceScale(d);
	}
	if( opt_bp_dist_estimate_min_score.set )
	{
		double d = strtod( opt_bp_dist_estimate_min_score.arg_value.c_str(), NULL );
		aligner.setBpDistEstimateMinScore(d);
	}
	if( opt_disable_cache.set )
	{
		aligner.SetUseCacheDb(false);
	}

	if( opt_min_scaled_penalty.set )
	{
		aligner.setMinimumBreakpointPenalty(strtod( opt_min_scaled_penalty.arg_value.c_str(), NULL ) );
	}
	if( pairwise_match_list.seq_table.size() != 0 )
	{
		aligner.setPairwiseMatches( pairwise_match_list );
	}
	if( opt_muscle_args.set )
	{
		MuscleInterface& mi = MuscleInterface::getMuscleInterface();
		mi.SetExtraMuscleArguments(opt_muscle_args.arg_value);
	}
	if( opt_recursive.set )
		aligner.SetRecursive(false);
	else
		aligner.SetRecursive(true);

	PairwiseScoringScheme pss;
	if( opt_gap_open.set )
	{
		pss.gap_open = atoi(opt_gap_open.arg_value.c_str());
	}
	if( opt_gap_extend.set )
	{
		pss.gap_extend = atoi(opt_gap_open.arg_value.c_str());
	}
	if( opt_substitution_matrix.set )
	{
		ifstream sub_in( opt_substitution_matrix.arg_value.c_str() );
		if( !sub_in.is_open() )
		{
			cerr << "Error opening substitution matrix file: \"" << opt_substitution_matrix.arg_value << "\"\n";
			return -1;
		}
		score_t matrix[4][4];
		readSubstitutionMatrix( sub_in, matrix );
		pss = PairwiseScoringScheme(matrix, pss.gap_open, pss.gap_extend);
	}
	aligner.setPairwiseScoringScheme(pss);

	if( opt_input_guide_tree.set )
		aligner.setInputGuideTreeFileName( opt_input_guide_tree.arg_value );
	if( opt_output_guide_tree.set )
		aligner.setOutputGuideTreeFileName( opt_output_guide_tree.arg_value );

	// if we will be doing a profile-profile or profile-sequence alignment
	// then read in the profile
	IntervalList profile_1;
	IntervalList profile_2;
	if( opt_profile.set ){
		cerr << "Profile-profile alignment not yet implemented\n";
		return -3;
	}

	IntervalList interval_list;
	interval_list.seq_table = pairwise_match_list.seq_table;
	interval_list.seq_filename = pairwise_match_list.seq_filename;

	if( opt_profile.set )
		; //aligner.alignPP(profile_1, profile_2, interval_list );
	else
		aligner.align( interval_list.seq_table, interval_list );

	if( !opt_disable_backbone.set )
	{

		string bbcols_fname = opt_output.arg_value + ".bbcols";
		string bb_seq_fname = opt_backbone_output.arg_value;
		if( !opt_backbone_output.set )
			bb_seq_fname = opt_output.arg_value + ".backbone";
		applyBackbone( interval_list, bbcols_fname, bb_seq_fname, island_gap_size, hmm_identity, pgh, pgu );
	}

	interval_list.WriteStandardAlignment(*match_out);
	match_out->flush();

	for( size_t seqI = 0; seqI < pairwise_match_list.seq_table.size(); seqI++ )
		delete pairwise_match_list.seq_table[seqI];	// an auto_ptr or shared_ptr could be great for this
	for( size_t seqI = 0; seqI < pairwise_match_list.sml_table.size(); seqI++ )
		delete pairwise_match_list.sml_table[seqI];

// only explicitly free memory if absolutely necessary
// since free() is very slow and the OS will reclaim it at program exit anyways
	if(opt_mem_clean.set)
	{
		// free memory used by pairwise matches
		for( size_t mI = 0; mI < pairwise_match_list.size(); mI++ )
			pairwise_match_list[mI]->Free();

		if( opt_output.set )
			delete match_out;
	}

/*
}catch( gnException& gne ) {
	cerr << "Unhandled gnException: " << gne << endl;
	throw gne;
	return -10;
}catch( exception& e ) {
	cerr << "Unhandled exception: " << e.what() << endl;
	throw e;
	return -11;
}catch( char* message ){
	cerr << "Unhandled exception: " << message << endl;
	throw message;
	return -12;
}catch( const char* message ){
	cerr << "Unhandled exception: " << message << " (const)\n";
	throw message;
	return -14;
}catch(...){
	cerr << "Unknown exception occurred.\n";
	throw;
	return -13;
}
*/
	return 0;
}

