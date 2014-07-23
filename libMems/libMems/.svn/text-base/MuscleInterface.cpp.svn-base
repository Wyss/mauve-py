/*******************************************************************************
 * $Id: MuscleInterface.cpp,v 1.27 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing, copying, and modification
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/MuscleInterface.h"

#include "libGenome/gnFilter.h"
#include "libGenome/gnFASSource.h"
#include "libGenome/gnStringTools.h"
#include "libMUSCLE/muscle.h"
#include "libMUSCLE/params.h"
#include "libMUSCLE/msa.h"
#include "libMUSCLE/seq.h"
#include "libMUSCLE/seqvect.h"
#include "libMUSCLE/tree.h"
#include "libMUSCLE/clust.h"
#include "libMUSCLE/profile.h"
#include "libMUSCLE/distfunc.h"
#include "libMUSCLE/clustsetdf.h"
#include "libMUSCLE/textfile.h"
#include "libMUSCLE/types.h"
#include "boost/algorithm/string/erase.hpp"
#include "boost/algorithm/string/case_conv.hpp"

#include <sstream>
#include <fstream>

using namespace std;
using namespace genome;

// this gets defined in muscle.cpp, but not declared in any headers
namespace muscle {
extern void MUSCLE(SeqVect &v, MSA &msaOut);
extern void RefineW(const MSA &msaIn, MSA &msaOut);
}

using namespace muscle;

namespace mems {

bool debug_muscle = false;

bool pipeExec( char** cmd_argv, const string& command, const string& input, string& output, string& error );

char** parseCommand( const string& cmd );
char** parseCommand( const string& cmd ){
	// first count tokens

	// tokenize on "
	stringstream qs( cmd );
	string cur_str;
	boolean in_quote = true;
	int token_count = 0;
	vector< string > cmd_tokens;
	while( getline( qs, cur_str, '"' ) ){
		// never start out in a quote
		in_quote = !in_quote;
		if( cur_str.length() == 0 )
			continue;
		if( in_quote ){
			cmd_tokens.push_back( cur_str );
		}else{
			stringstream ss( cur_str );
			string asdf;
			while( ss >> asdf )
				cmd_tokens.push_back( asdf );
		}
	}
	char ** cmd_array = new char*[ cmd_tokens.size() + 1 ];
	for( int tokI = 0; tokI < cmd_tokens.size(); tokI++ ){
		cmd_array[ tokI ] = new char[ cmd_tokens[ tokI ].length() + 1 ];
		strcpy( cmd_array[ tokI ], cmd_tokens[ tokI ].c_str() );
	}
	cmd_array[ cmd_tokens.size() ] = NULL;
	return cmd_array;
}

#if !defined(WIN32)
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>

// unix pipelined execution code
bool pipeExec( char** cmd_argv, const string& command, const string& input, string& output, string& error ){
	int stdin_pipe[2], stdout_pipe[2], stderr_pipe[2];
	boolean success = false;
	pid_t sid;
	pid_t pid1;
	const char* fail;
	char buf[1024];
	ssize_t bread = 0;
	int rval = 0;

	if((sid = setsid()) < 0)sid = getpgrp();
	
	if((sid < 0 && (fail = "sid"))
	 || (pipe(stdin_pipe) < 0 && (fail = "stdin"))
	 || (pipe(stdout_pipe) < 0 && (fail = "stdout"))
//	 || (pipe(stderr_pipe) < 0 && (fail = "stderr"))
	)
    {
		fprintf(stderr, "Ouch, the world just collapsed (%s).\n", fail);
		perror("muscle:");
		goto cleanup;
	}
	
	fcntl(stdin_pipe[0], F_SETFL, fcntl(stdin_pipe[0], F_GETFL) & ~O_NONBLOCK);
	fcntl(stdin_pipe[1], F_SETFL, fcntl(stdin_pipe[1], F_GETFL) & ~O_NONBLOCK);
	fcntl(stdout_pipe[0], F_SETFL, fcntl(stdout_pipe[0], F_GETFL) & ~O_NONBLOCK);
	fcntl(stdout_pipe[1], F_SETFL, fcntl(stdout_pipe[1], F_GETFL) & ~O_NONBLOCK);
/*	fcntl(stderr_pipe[0], F_SETFL, fcntl(stderr_pipe[0], F_GETFL) & ~O_NONBLOCK);
	fcntl(stderr_pipe[1], F_SETFL, fcntl(stderr_pipe[1], F_GETFL) & ~O_NONBLOCK);
*/	
	if((pid1 = fork()) < 0)goto cleanup;	
	if(pid1)
		setpgid(pid1, sid);
	else
	{
		dup2(stdin_pipe[0], 0);
		dup2(stdout_pipe[1], 1);
//		dup2(stderr_pipe[1], 2);
		close( stdin_pipe[0] );
		close( stdin_pipe[1] );
		close( stdout_pipe[0] );
		close( stdout_pipe[1] );
//		close( stderr_pipe[0] );
//		close( stderr_pipe[1] );
		execvp(cmd_argv[0], cmd_argv);
		_exit(errno);
	}
	rval = write( stdin_pipe[1], input.c_str(), input.size() );
	if( rval == -1 )
		perror( "write: " );
	if( close( stdin_pipe[1] ) )
		perror( "close stdin_w: " );
	if( close( stdin_pipe[0] ) )
		perror( "close stdin_r: " );

	close( stdout_pipe[1] );
	// read the alignment
	while(true){
		bzero( buf, sizeof(buf) );
		bread = read( stdout_pipe[0], buf, 1023 );
		if( bread == 0 )
			break;	// reached EOF
		if( bread == -1 ){
			perror("muscle read: " );
		}
		output += buf;
	}
	wait( NULL );
	success = true;
	
cleanup:
	close(stdin_pipe[0]);
	close(stdin_pipe[1]);
	close(stdout_pipe[0]);
	close(stdout_pipe[1]);
//	close(stderr_pipe[0]);
//	close(stderr_pipe[1]);
	return success;
};


#else

//windows piping code
#include <windows.h>
#define bzero(a) memset(a,0,sizeof(a)) //easier -- shortcut

bool IsWinNT()  //check if we're running NT
{
  OSVERSIONINFO osv;
  osv.dwOSVersionInfoSize = sizeof(osv);
  GetVersionEx(&osv);
  return (osv.dwPlatformId == VER_PLATFORM_WIN32_NT);
}

void ErrorMessage(char *str)  //display detailed error info
{
  LPVOID msg;
  FormatMessage(
    FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM,
    NULL,
    GetLastError(),
    MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), // Default language
    (LPTSTR) &msg,
    0,
    NULL
  );
  printf("%s: %s\n",str,msg);
  LocalFree(msg);
}

bool pipeExec( char** cmd_argv, const string& command, const string& input, string& output, string& error ){

	char buf[1024];           //i/o buffer

	STARTUPINFO si;
	SECURITY_ATTRIBUTES sa;
	SECURITY_DESCRIPTOR sd;               //security information for pipes
	PROCESS_INFORMATION pi;
	HANDLE newstdin_w,newstdout_w,newstderr_w,newstdin_r,newstdout_r,newstderr_r;
	HANDLE read_stdout,read_stderr,write_stdin;  //pipe handles
	boolean success = false;

	if (IsWinNT())        //initialize security descriptor (Windows NT)
	{
		InitializeSecurityDescriptor(&sd,SECURITY_DESCRIPTOR_REVISION);
		SetSecurityDescriptorDacl(&sd, true, NULL, false);
		sa.lpSecurityDescriptor = &sd;
	}
	else sa.lpSecurityDescriptor = NULL;

	sa.nLength = sizeof(SECURITY_ATTRIBUTES);
	sa.bInheritHandle = true;         //allow inheritable handles

	if (!CreatePipe(&newstdin_r,&newstdin_w,&sa,0))   //create stdin pipe
	{
		ErrorMessage("CreatePipe");
		goto finito;
	}
	if (!CreatePipe(&newstdout_r,&newstdout_w,&sa,0))  //create stdout pipe
	{
		ErrorMessage("CreatePipe");
		goto finito;
	}
	if (!CreatePipe(&newstderr_r,&newstderr_w,&sa,0))  //create stdout pipe
	{
		ErrorMessage("CreatePipe");
		goto finito;
	}
	// Duplicate the write handle to the pipe so it is not inherited. 
	boolean fSuccess = DuplicateHandle(GetCurrentProcess(), newstdin_w, 
		GetCurrentProcess(), &write_stdin, 0, 
		FALSE,                  // not inherited 
		DUPLICATE_SAME_ACCESS); 
	if (! fSuccess){
		ErrorMessage("DuplicateHandle failed"); 
		goto finito;
	}
	CloseHandle(newstdin_w); 
	newstdin_w = INVALID_HANDLE_VALUE;

	// Duplicate the read handle to the pipe so it is not inherited. 
	fSuccess = DuplicateHandle(GetCurrentProcess(), newstdout_r, 
		GetCurrentProcess(), &read_stdout, 0, 
		FALSE,                  // not inherited 
		DUPLICATE_SAME_ACCESS); 
	if (! fSuccess){
		ErrorMessage("DuplicateHandle failed"); 
		goto finito;
	}
	CloseHandle(newstdout_r); 
	newstdout_r = INVALID_HANDLE_VALUE;

	// Duplicate the read handle to the pipe so it is not inherited. 
	fSuccess = DuplicateHandle(GetCurrentProcess(), newstderr_r, 
		GetCurrentProcess(), &read_stderr, 0, 
		FALSE,                  // not inherited 
		DUPLICATE_SAME_ACCESS); 
	if (! fSuccess){
		ErrorMessage("DuplicateHandle failed"); 
		goto finito;
	}
	CloseHandle(newstderr_r); 
	newstderr_r = INVALID_HANDLE_VALUE;

	GetStartupInfo(&si);      //set startupinfo for the spawned process
	/*
		The dwFlags member tells CreateProcess how to make the process.
		STARTF_USESTDHANDLES validates the hStd* members. STARTF_USESHOWWINDOW
		validates the wShowWindow member.
	*/
	si.dwFlags = STARTF_USESTDHANDLES|STARTF_USESHOWWINDOW;
	si.wShowWindow = SW_HIDE;
	si.hStdOutput = newstdout_w;
	si.hStdError = newstderr_w;     //set the new handles for the child process
	si.hStdInput = newstdin_r;

	//spawn the child process
	char* cmd = new char[ command.length() + 1 ];
	strcpy( cmd, command.c_str() );
	if (!CreateProcess(NULL,cmd,NULL,NULL,TRUE,0,
						NULL,NULL,&si,&pi))
	{
		delete cmd;
		ErrorMessage("CreateProcess");
		goto finito;
	}
	delete cmd;

	unsigned long exit=0;  //process exit code
	unsigned long bread;   //bytes read
	unsigned long avail;   //bytes available

	WriteFile(write_stdin, input.c_str(), input.size(), &bread, NULL); //send data to stdin
	CloseHandle(write_stdin);
	write_stdin = INVALID_HANDLE_VALUE;

	// Wait until child process exits.
	while( true ){
		GetExitCodeProcess( pi.hProcess, &exit );
		if( exit != STILL_ACTIVE )
			WaitForSingleObject( pi.hProcess, INFINITE );
			
		// read anything that came to stdout
		PeekNamedPipe(read_stdout,buf,1023,&bread,&avail,NULL);
		if( avail == 0 )
			Sleep(5);	// didn't get anything, so take a break to avoid hogging the CPU...
		while( avail > 0 ){
			bzero(buf);
			int read_size = 1023 < avail ? 1023 : avail;
			ReadFile(read_stdout,buf,read_size,&bread,NULL);  //read the stdout pipe
			avail -= bread;
			output += buf;
		}

		// read anything that came to stderr
		PeekNamedPipe(read_stderr,buf,1023,&bread,&avail,NULL);
		while( avail > 0 ){
			bzero(buf);
			int read_size = 1023 < avail ? 1023 : avail;
			ReadFile(read_stderr,buf,read_size,&bread,NULL);  //read the stdout pipe
			avail -= bread;
			error += buf;
		}

		if( exit != STILL_ACTIVE )
			break;
	}
	// Wait until child process exits.
    WaitForSingleObject( pi.hProcess, INFINITE );
	success = true;

	//clean up and exit
finito:
    if( pi.hThread != INVALID_HANDLE_VALUE )
		CloseHandle(pi.hThread);
    if( pi.hProcess != INVALID_HANDLE_VALUE )
		CloseHandle(pi.hProcess);
    if( newstdin_r != INVALID_HANDLE_VALUE )
		CloseHandle(newstdin_r);
    if( newstdout_w != INVALID_HANDLE_VALUE )
		CloseHandle(newstdout_w);
    if( newstderr_w != INVALID_HANDLE_VALUE )
		CloseHandle(newstderr_w);
    if( read_stdout != INVALID_HANDLE_VALUE )
		CloseHandle(read_stdout);
    if( read_stderr != INVALID_HANDLE_VALUE )
		CloseHandle(read_stderr);
    if( write_stdin != INVALID_HANDLE_VALUE )
		CloseHandle(write_stdin);
	return success;
}

#endif



MuscleInterface& MuscleInterface::getMuscleInterface()
{
        static MuscleInterface m_ci;

        return m_ci;
}

MuscleInterface::MuscleInterface() : GappedAligner() {
	muscle_path = "muscle_aed";
	muscle_arguments = "-stable -quiet -seqtype DNA";
	muscle_cmdline = parseCommand( muscle_path + " " + muscle_arguments );
	max_alignment_length = 12500;
}

void MuscleInterface::ParseMusclePath( const char* argv0 ){
	// get the execution path
	string path_str = argv0;
	// trim quotes
	if( path_str[0] == '"' )
		path_str = path_str.substr( 1, path_str.size() - 2 );
	standardizePathString( path_str );
	string::size_type i = path_str.rfind('/');
	if( i != string::npos )
		path_str.erase(i+1, path_str.length() - (i+1));
	else
		path_str.clear();
	SetMusclePath( '"' + path_str + "muscle_aed\"");
}

void MuscleInterface::SetMusclePath( const string& path ){
	muscle_path = path;
	ClearCommandLine();
	muscle_cmdline = parseCommand( muscle_path + " " + muscle_arguments );
}

void MuscleInterface::SetExtraMuscleArguments( const string& args )
{
	extra_muscle_arguments = args;
}

void MuscleInterface::SetMuscleArguments( const string& args )
{
	ClearCommandLine();
	muscle_arguments = args + " " + extra_muscle_arguments;
	muscle_cmdline = parseCommand( muscle_path + " " + args + " " + extra_muscle_arguments );
}

MuscleInterface& MuscleInterface::operator=( const MuscleInterface& ci ){
	GappedAligner::operator =( ci );
	return *this;
}

//tjt: not the best way of doing this, should have just one Align function that takes an AbstractMatch*,
	//     not both Match* & AbstractMatch* in separate, nearly identical functions..
	//     Such a change would involve changes to GappedAligner, and would require some additional care taken
	//     with SeqCount & Multiplicity, as well as seq_table[ seqI ]->length()/seq_table[ 0 ]->length(i),
	//     for now, leave like this. hopefully sooner than later, make pretty!
boolean MuscleInterface::Align( GappedAlignment& cr, Match* r_begin, Match* r_end, vector< gnSequence* >& seq_table ){
	gnSeqI gap_size = 0;
	boolean create_ok = true;
	uint seq_count = seq_table.size();
	//seq_count = r_begin->Multiplicity();
	uint seqI;
	uint align_seqs = 0;
	vector< string > tmp_mat = vector< string >( seq_count );
try{

// 
//	Get the sequence in the intervening gaps between these two matches
//
	vector< string > seq_data;
	vector< int64 > starts;
	vector< uint > seqs;
	const gnFilter* rc_filter = gnFilter::DNAComplementFilter();
	
	for( seqI = 0; seqI < seq_count; seqI++ ){

		// skip this sequence if it's undefined
		if( (r_end != NULL && r_end->Start( seqI ) == NO_MATCH ) ||
			(r_begin != NULL && r_begin->Start( seqI ) == NO_MATCH) ){
			starts.push_back( NO_MATCH );
			continue;
		}

		// determine the size of the gap
		int64 gap_start = 0;
		int64 gap_end = 0;
		getInterveningCoordinates( seq_table, r_begin, r_end, seqI, gap_start, gap_end );

		int64 diff = gap_end - gap_start;
		if( diff <= 0 || diff > max_alignment_length ){
			starts.push_back( NO_MATCH );
			continue;	// skip this sequence if it's either too big or too small
		}
		seqs.push_back( seqI );
// the gnSequence pointers are shared across threads and have a common ifstream
		// extract sequence data
		if( r_end == NULL || r_end->Start( seqI ) > 0 ){
			starts.push_back( gap_start );
			seq_data.push_back( seq_table[ seqI ]->ToString( diff , gap_start ) );
		}else{
			// reverse complement the sequence data.
			starts.push_back( -gap_start );
			string cur_seq_data = seq_table[ seqI ]->ToString( diff , gap_start );
			rc_filter->ReverseFilter( cur_seq_data );
			seq_data.push_back( cur_seq_data );
		}
	}

	if( seqs.size() <= 1 )
		create_ok = false;

	if( create_ok ){
//		SetMuscleArguments( " -quiet -stable -seqtype DNA " );
		vector< string > aln_matrix;
		if( !CallMuscleFast( aln_matrix, seq_data, 0, 0 ) ){
			cout << "Muscle was unable to align:\n";
			if( r_begin )
				cout << "Left match: " << *r_begin << endl;
			if( r_end )
				cout << "Right match: " << *r_end << endl;
			return false;
		}

		gnSeqI aln_length = aln_matrix.size() == 0 ? 0 : aln_matrix[0].length();
		cr = GappedAlignment( seq_count, aln_length );
		vector< string > aln_mat = vector< string >( seq_count );

		// set sequence starts
		for( uint seqI = 0; seqI < seqs.size(); seqI++ ){
			cr.SetLength( seq_data[ seqI ].size(), seqs[ seqI ] );
			aln_mat[ seqs[ seqI ] ] = aln_matrix[ seqI ];
		}
		for( uint seqI = 0; seqI < seq_count; seqI++ ){
			cr.SetStart( seqI, starts[ seqI ] );
			if( aln_mat[ seqI ].length() != aln_length )
				aln_mat[ seqI ] = string( aln_length, '-' );
		}

		cr.SetAlignment( aln_mat );

		return true;
	}
}catch(exception& e){
	cerr << "At: " << __FILE__ << ":" << __LINE__ << endl;
	cerr << e.what();
}
	return false;
}

static int failure_count = 0;

boolean MuscleInterface::Align( GappedAlignment& cr, AbstractMatch* r_begin, AbstractMatch* r_end, vector< gnSequence* >& seq_table){
	gnSeqI gap_size = 0;
	boolean create_ok = true;
	//tjt: set the seq_count to a match m's multiplicity
	//     even though all components n of match m could be 
	//     less than the k sequences
	//     if n == k, then perhaps there is 1 match component per sequence
	//     if k = 1, n == repeat match multiplicity, where n >= 2
	//     
	uint seq_count = r_begin->Multiplicity();
	uint seqI;
	uint align_seqs = 0;
	vector< string > tmp_mat = vector< string >( seq_count );
try{

// 
//	Get the sequence in the intervening gaps between these two matches
//
	vector< string > seq_data;
	vector< int64 > starts;
	vector< uint > seqs;
	const gnFilter* rc_filter = gnFilter::DNAComplementFilter();
	
	//std::cout << "getting regions between match components to align" << std::endl;
	for( seqI = 0; seqI < seq_count; seqI++ ){

		// skip this sequence if it's undefined
		if( (r_end != NULL && r_end->Start( seqI ) == NO_MATCH ) ||
			(r_begin != NULL && r_begin->Start( seqI ) == NO_MATCH) ){
			starts.push_back( NO_MATCH );
			continue;
		}

		// determine the size of the gap
		int64 gap_start = 0;
		int64 gap_end = 0;

		// determine the size of the gap
		gap_end = r_end != NULL ? r_end->Start( seqI ) : seq_table[ seqI ]->length() + 1;
		gap_start = r_begin != NULL ? r_begin->End( seqI ) + 1 : 1;
		if( gap_end < 0 || gap_start < 0 ){
			gap_end   = r_begin != NULL ? -r_begin->Start( seqI ) : seq_table[ 0 ]->length() + 1;
			gap_start = r_end != NULL ? -r_end->Start( seqI ) + r_end->Length( seqI ) : 1;
		}
		if( gap_end <= 0 || gap_start <= 0 ){
			// if either is still < 0 then there's a problem...
			genome::ErrorMsg( "Error constructing intervening coordinates" );
		}

		int64 diff = gap_end - gap_start;
		
		//diff <= 0 ||
		if( diff <= 0 || diff > max_alignment_length ){
			starts.push_back( NO_MATCH );
			continue;	// skip this sequence if it's either too big or too small
		}

		seqs.push_back( seqI );

		// extract sequence data
		if (0 )
		{
			starts.push_back( gap_start );
			seq_data.push_back( "A" );
			std::cout << "A" << std::endl;
			diff = 1;
		}
// the gnSequence pointers are shared across threads and have a common ifstream
		if( r_end == NULL || r_end->Start( seqI ) > 0 ){
			starts.push_back( gap_start );
			//std::cout << seq_table[ 0 ]->ToString( diff , gap_start ) << std::endl;
			//tjt: all sequences are concatenated together into 1 seq_table entry
			//
			seq_data.push_back( seq_table[ 0 ]->ToString( diff , gap_start ) );
		}else{
			// reverse complement the sequence data.
			starts.push_back( -gap_start );
			//tjt: all sequences are concatenated together into 1 seq_table entry
			//     
			string cur_seq_data = seq_table[ 0 ]->ToString( diff , gap_start );
			rc_filter->ReverseFilter( cur_seq_data );
			seq_data.push_back( cur_seq_data );
			//std::cout << cur_seq_data << std::endl;
		}
	}

    //no seqs able to be aligned..
    if( seqs.size() == 0)
        create_ok = false;


	if( create_ok ){
//		SetMuscleArguments( " -quiet -stable -seqtype DNA " );
		vector< string > aln_matrix;
		if( !CallMuscleFast( aln_matrix, seq_data, 0, 0 ) ){
			cout << "Muscle was unable to align:\n";
			return false;
		}
        
        //fill in regions between adjacent seeds with gaps
        //if aln_matrix is smaller than multiplicity, then we know 
        //that there are some regions between seeds that have len == 0
        if (aln_matrix.size() != r_begin->Multiplicity() && 0)
        {
            for( uint seqI = 0; seqI < starts.size(); seqI++ )
            {
                //if this a position between two adjacent matches..
                if (starts.at(seqI) == NO_MATCH)
                {
                    //calculate the number of gaps to fill in
                    int64 gap_end = r_end != NULL ? r_end->Start( seqI ) : seq_table[ seqI ]->length() + 1;
		            int64 gap_start = r_begin != NULL ? r_begin->End( seqI ) + 1 : 1;
                    if( r_end == NULL || r_end->Start( seqI ) > 0 ){
			            starts[seqI] = 0;//gap_start;
			            seq_data.insert(seq_data.begin()+(seqI),"");
		            }else{
			            starts[seqI] = 0;//-gap_start;
			            seq_data.insert(seq_data.begin()+(seqI),"");
		            }
                    string tmp(aln_matrix[0].length(), '-');
                    aln_matrix.insert(aln_matrix.begin()+(seqI), tmp);
                    seqs.insert(seqs.begin()+(seqI),seqI);
                }
            }
        }
		gnSeqI aln_length = aln_matrix.size() == 0 ? 0 : aln_matrix[0].length();
		cr = GappedAlignment( seq_count, aln_length );
		vector< string > aln_mat = vector< string >( seq_count );

		// set sequence starts
		for( uint seqI = 0; seqI < seqs.size(); seqI++ ){
			cr.SetLength( seq_data[ seqI ].size(), seqs[ seqI ] );
			aln_mat[ seqs[ seqI ] ] = aln_matrix[ seqI ];
		}
		for( uint seqI = 0; seqI < seq_count; seqI++ ){
			cr.SetStart( seqI, starts[ seqI ] );
			if( aln_mat[ seqI ].length() != aln_length )
				aln_mat[ seqI ] = string( aln_length, '-' );
		}

		cr.SetAlignment( aln_mat );

		return true;
	}
}catch(exception& e){
	cerr << "At: " << __FILE__ << ":" << __LINE__ << endl;
	cerr << e.what();
}
	return false;
}

boolean MuscleInterface::CallMuscle( vector< string >& aln_matrix, const vector< string >& seq_table )
{
	gnSequence seq;

	try{
		ostringstream input_seq_stream;
		//istringstream muscle_input_seq_stream;
		for( uint seqI = 0; seqI < seq_table.size(); seqI++ ){
			seq += seq_table[ seqI ];
			seq.setContigName( seqI, "seq" );
		}
		gnFASSource::Write( seq, input_seq_stream, false, true );		
		// now open a pipe to Muscle
		string muscle_cmd = muscle_path + " " + muscle_arguments;
		string output;
		string error;
		boolean success = pipeExec( muscle_cmdline, muscle_cmd, input_seq_stream.str(), output, error );
		if( !success || output.size() == 0 )
		{
			throw "b0rk3d";
		}

		istringstream output_aln_stream( output );
		string cur_line;

		// parse the fasta output
		while( getline( output_aln_stream, cur_line ) ){
			if( cur_line[0] == '>' ){
				aln_matrix.push_back( "" );
				continue;
			}
			gnSeqI len = cur_line.size();
			len = cur_line[ len - 1 ] == '\r' ? len - 1 : len;
			uint seqI = aln_matrix.size() - 1;
			aln_matrix[ seqI ] += cur_line.substr( 0, len );
		}

		return true;
	}catch( gnException& gne ){
	}catch( exception& e ){
	}catch(...){
	}
	cerr << "muscle failed!  saving failed input data to muscle_failure_" << failure_count << ".txt\n";
	cerr << "Please contact the Mauve developers about this problem\n";
	stringstream debug_fname;
	debug_fname << "muscle_failure_" << failure_count++ << ".txt";
	ofstream debug_file( debug_fname.str().c_str() );
	gnFASSource::Write(seq, debug_file, false);
	debug_file.close();
	return false;
}

// version 2 of this code: attempt to call muscle without performing costly disk I/O!!
boolean MuscleInterface::CallMuscleFast( vector< string >& aln_matrix, const vector< string >& seq_table, int gap_open, int gap_extend )
{
	if (gap_open != 0)
		g_scoreGapOpen.get() = gap_open;
	if (gap_extend != 0)
		g_scoreGapExtend.get() = gap_extend;
	g_SeqType.get() = SEQTYPE_DNA;	// we're operating on DNA
	g_uMaxIters.get() = 1;			// and we don't want to refine the alignment...yet
	g_bStable.get() = true;			// we want output seqs in the same order as input
	g_bQuiet.get() = true;			// and don't print anything to the console
	g_SeqWeight1.get() = SEQWEIGHT_ClustalW;	// not sure what weighting scheme works best for DNA

	SetMaxIters(g_uMaxIters.get());
	SetSeqWeightMethod(g_SeqWeight1.get());

	// now construct a SeqVect containing input sequences
	SeqVect sv;
	const char* seqname = "seq00000";
	for( size_t seqI = 0; seqI < seq_table.size(); seqI++ )
	{
		Seq curseq;
		curseq.SetId(seqI);
		curseq.SetName(seqname);
		curseq.resize(seq_table[seqI].size());
		std::copy(seq_table[seqI].begin(), seq_table[seqI].end(), curseq.begin());
		sv.AppendSeq(curseq);
	}

	MSA msaTmp;
	MUSCLE(sv,msaTmp);

	// now extract the alignment
	aln_matrix.clear();
	aln_matrix.resize(msaTmp.GetSeqCount());
	for( size_t seqI = 0; seqI < msaTmp.GetSeqCount(); seqI++ )
	{
		unsigned indie = msaTmp.GetSeqIndex(seqI);
		const char* buf = msaTmp.GetSeqBuffer(indie);
		string curseq(buf, msaTmp.GetColCount());
		swap(aln_matrix[seqI],curseq);
	}
	return true;	// how can it possibly fail? :)
}

bool MuscleInterface::Refine( GappedAlignment& ga, size_t windowsize )
{
	const vector< string >& seq_table = GetAlignment( ga, vector< gnSequence* >() );
	vector< string > aln_table;
	for( uint seqI = 0; seqI < ga.SeqCount(); seqI++ )
	{
		if( ga.LeftEnd(seqI) != NO_MATCH )
		{
			aln_table.push_back( seq_table[seqI] );
		}
	}
	vector< string > aln_matrix;
	if( windowsize == 0 )
		SetMuscleArguments( " -quiet -refine -seqtype DNA " );
	else
	{
		stringstream sstr;
		sstr << " -quiet -seqtype DNA -refinew -refinewindow " << windowsize << " ";
		SetMuscleArguments( sstr.str() );
	}
	bool success = CallMuscle( aln_matrix, aln_table );
	if( success )
	{
		aln_table.clear();
		uint alnI = 0;
		for( uint seqI = 0; seqI < ga.SeqCount(); seqI++ )
		{
			if( ga.LeftEnd(seqI) != NO_MATCH )
				aln_table.push_back( aln_matrix[alnI++] );
			else
				aln_table.push_back( string( aln_matrix[0].size(), '-' ) );
		}
		ga.SetAlignment( aln_table );
	}
	return success;
}

void msaFromSeqTable(MSA& msa, const vector< string >& seq_table, unsigned id_base = 0)
{
	msa.SetSize(seq_table.size(), seq_table[0].size());
	for( uint seqI = 0; seqI < seq_table.size(); seqI++ )
	{
		stringstream ss;
		ss << "seq" << seqI;
		msa.SetSeqName(seqI, ss.str().c_str());
		msa.SetSeqId(seqI,seqI+id_base);
		for(size_t i = 0; i < seq_table[seqI].size(); i++)
			msa.SetChar(seqI, i, seq_table[seqI][i]);
	}
}


bool MuscleInterface::RefineFast( GappedAlignment& ga, size_t windowsize )
{
	const vector< string >& seq_table = GetAlignment( ga, vector< gnSequence* >() );
	vector< string > aln_table;
	for( uint seqI = 0; seqI < ga.SeqCount(); seqI++ )
	{
		if( ga.LeftEnd(seqI) != NO_MATCH )
		{
			aln_table.push_back( seq_table[seqI] );
		}
	}

	g_SeqType.get() = SEQTYPE_DNA;	// we're operating on DNA
	g_uMaxIters.get() = 1;			// and we don't want to refine the alignment...yet
	g_bStable.get() = true;			// we want output seqs in the same order as input
	g_bQuiet.get() = true;			// and don't print anything to the console
	g_SeqWeight1.get() = SEQWEIGHT_ClustalW;	// not sure what weighting scheme works best for DNA

	g_uRefineWindow.get() = windowsize;
	g_uWindowTo.get() = 0;

	SetMaxIters(g_uMaxIters.get());
	SetSeqWeightMethod(g_SeqWeight1.get());

	MSA::SetIdCount(seq_table.size());

	// create an MSA
	MSA msa;
	msaFromSeqTable(msa, seq_table);

	SetAlpha(ALPHA_DNA);
	msa.FixAlpha();
	SetPPScore(PPSCORE_SPN);
	SetMuscleInputMSA(msa);

	Tree GuideTree;
	TreeFromMSA(msa, GuideTree, g_Cluster2.get(), g_Distance2.get(), g_Root2.get());
	SetMuscleTree(GuideTree);

	MSA msaOut;
	MSA* finalMsa;

	if(windowsize == 0)
	{
		if (g_bAnchors.get())
			RefineVert(msa, GuideTree, g_uMaxIters.get());
		else
			RefineHoriz(msa, GuideTree, g_uMaxIters.get(), false, false);
		finalMsa = &msa;
	}else{
		RefineW(msa, msaOut);
		finalMsa = &msaOut;
	}


	ValidateMuscleIds(*finalMsa);
	ValidateMuscleIds(GuideTree);

	// now extract the alignment
	vector< string > aln_matrix;
	aln_matrix.resize(finalMsa->GetSeqCount());
	for( size_t seqI = 0; seqI < finalMsa->GetSeqCount(); seqI++ )
	{
		unsigned indie = finalMsa->GetSeqIndex(seqI);
		const char* buf = finalMsa->GetSeqBuffer(indie);
		string curseq(buf, finalMsa->GetColCount());
		swap(aln_matrix[seqI],curseq);
	}

	ga.SetAlignment( aln_matrix );
	return true;
}


void stripGapColumns( std::vector< std::string >& aln )
{
	size_t cur_col = 0;
	size_t gap_seq = 0;
	for( size_t colI = 0; colI < aln[0].size(); colI++ )
	{
		gap_seq = 0;
		for( ; gap_seq < aln.size(); gap_seq++ )
			if( aln[gap_seq][colI] != '-' )
				break;
		if( gap_seq != aln.size() )
		{
			for( gap_seq = 0; gap_seq < aln.size(); gap_seq++ )
				aln[gap_seq][cur_col] = aln[gap_seq][colI];
			cur_col++;
		}
	}
	for( gap_seq = 0; gap_seq < aln.size(); gap_seq++ )
		aln[gap_seq].resize(cur_col);
}

void stripGaps( std::string& str )
{
	std::string::iterator striter = std::remove(str.begin(), str.end(), '-');
	str.resize(striter - str.begin());
}

bool MuscleInterface::ProfileAlign( const GappedAlignment& ga1, const GappedAlignment& ga2, GappedAlignment& aln, bool anchored )
{
	try{
		const vector< string >& aln1 = GetAlignment( ga1, vector< gnSequence* >() );
		const vector< string >& aln2 = GetAlignment( ga2, vector< gnSequence* >() );
		vector< uint > order;
		ostringstream input_seq_stream;
		gnSequence seq;
		vector< string > aln11( ga1.Multiplicity() );
		vector< string > aln22( ga2.Multiplicity() );
		size_t curI = 0;
		for( uint seqI = 0; seqI < aln1.size(); seqI++ )
		{
			if( ga1.LeftEnd(seqI) != NO_MATCH )
			{
				aln11[curI++] = aln1[seqI];
				order.push_back(seqI);
			}
		}
		curI = 0;
		for( uint seqI = 0; seqI < aln2.size(); seqI++ )
		{
			if( ga2.LeftEnd(seqI) != NO_MATCH )
			{
				aln22[curI++] = aln2[seqI];
				order.push_back(seqI);
			}
		}
// strip the gap columns only if we're doing unanchored PP alignment
		if( !anchored )
		{
			stripGapColumns(aln11);
			stripGapColumns(aln22);
		}
		for( uint seqI = 0; seqI < aln11.size(); seqI++ )
		{
			seq += aln11[ seqI ];
			seq.setContigName( seq.contigListLength()-1, "seq" );
		}

		gnFASSource::Write( seq, input_seq_stream, false, true );
		input_seq_stream << "=\n";

		gnSequence seq2;
		for( uint seqI = 0; seqI < aln22.size(); seqI++ )
		{
			seq2 += aln22[ seqI ];
			seq2.setContigName( seq2.contigListLength()-1, "seq" );
		}

		gnFASSource::Write( seq2, input_seq_stream, false, true );
		input_seq_stream << "=\n";

		if( debug_muscle )
		{
			// for debugging: write the anchored profiles to a file
			stringstream debug_fname;
			debug_fname << "muscle_debug_" << failure_count++ << ".txt";
			ofstream debug_file( debug_fname.str().c_str() );
			debug_file << input_seq_stream.str();
			debug_file.close();
		}

		// now open a pipe to Muscle
		string musc_args = "-quiet -seqtype DNA -profile -ProfileOnStdIn ";
		if( anchored )
			musc_args += "-AnchoredPP ";
		SetMuscleArguments( musc_args );
		string output;
		string error;
		string muscle_cmd = muscle_path + " " + muscle_arguments;
		if( debug_muscle )
		{
			cerr << "Running " << muscle_cmd << endl;
		}
		boolean success = pipeExec( muscle_cmdline, muscle_cmd, input_seq_stream.str(), output, error );
		if( !success || output.size() == 0 )
		{
			if( output.size() == 0 )
				cerr << "\nmuscle nothing\n";
			else
				cerr << "\nunsuccessful muscle\n";
			return false;
		}

		istringstream output_aln_stream( output );
		string cur_line;

		// parse the fasta output
		vector< string > aln_matrix( ga1.SeqCount() );
		int ordI = -1;
		while( getline( output_aln_stream, cur_line ) ){
			if( cur_line[0] == '>' ){
				ordI++;
				continue;
			}
			gnSeqI len = cur_line.size();
			len = cur_line[ len - 1 ] == '\r' ? len - 1 : len;
			uint seqI = aln_matrix.size() - 1;
			aln_matrix[ order[ordI] ] += cur_line.substr( 0, len );
		}
		for( size_t i = 0; i < aln_matrix.size(); i++ )
		{
			if( aln_matrix[i].size() == 0 )
				aln_matrix[i].resize( aln_matrix[order[0]].size(), '-' );
		}

		aln.SetAlignment( aln_matrix );
		for( uint seqI = 0; seqI < ga1.SeqCount(); seqI++ )
			if( ga1.LeftEnd(seqI) != NO_MATCH )
			{
				aln.SetLeftEnd(seqI, ga1.LeftEnd(seqI));
				aln.SetLength(ga1.Length(seqI), seqI);
			}
		for( uint seqI = 0; seqI < ga2.SeqCount(); seqI++ )
			if( ga2.LeftEnd(seqI) != NO_MATCH )
			{
				aln.SetLeftEnd(seqI, ga2.LeftEnd(seqI));
				aln.SetLength(ga2.Length(seqI), seqI);
			}
		return true;
	}catch( gnException& gne ){
	}catch( exception& e ){
	}catch(...){
	}
	return false;
}


bool MuscleInterface::ProfileAlignFast( const GappedAlignment& ga1, const GappedAlignment& ga2, GappedAlignment& aln, bool anchored )
{
	try{
		const vector< string >& aln1 = GetAlignment( ga1, vector< gnSequence* >() );
		const vector< string >& aln2 = GetAlignment( ga2, vector< gnSequence* >() );
		vector< uint > order;
		vector< string > aln11( ga1.Multiplicity() );
		vector< string > aln22( ga2.Multiplicity() );
		size_t curI = 0;
		for( uint seqI = 0; seqI < aln1.size(); seqI++ )
		{
			if( ga1.LeftEnd(seqI) != NO_MATCH )
			{
				aln11[curI++] = aln1[seqI];
				order.push_back(seqI);
			}
		}
		curI = 0;
		for( uint seqI = 0; seqI < aln2.size(); seqI++ )
		{
			if( ga2.LeftEnd(seqI) != NO_MATCH )
			{
				aln22[curI++] = aln2[seqI];
				order.push_back(seqI);
			}
		}
// strip the gap columns only if we're doing unanchored PP alignment
		if( !anchored )
		{
			stripGapColumns(aln11);
			stripGapColumns(aln22);
		}

		g_SeqType.get() = SEQTYPE_DNA;	// we're operating on DNA
		g_uMaxIters.get() = 1;			// and we don't want to refine the alignment...yet
		g_bStable.get() = true;			// we want output seqs in the same order as input
		g_bQuiet.get() = true;			// and don't print anything to the console
		g_SeqWeight1.get() = SEQWEIGHT_ClustalW;	// not sure what weighting scheme works best for DNA

		SetMaxIters(g_uMaxIters.get());
		SetSeqWeightMethod(g_SeqWeight1.get());

		MSA::SetIdCount(order.size());

		MSA msa1;
		MSA msa2;
		MSA msaOut;
		msaFromSeqTable(msa1, aln11);
		msaFromSeqTable(msa2, aln22, msa1.GetSeqCount());

		SetAlpha(ALPHA_DNA);
		msa1.FixAlpha();
		msa2.FixAlpha();
		SetPPScore(PPSCORE_SPN);

		if(anchored)
		{
			AnchoredProfileProfile(msa1, msa2, msaOut);
		}else{
			ProfileProfile(msa1, msa2, msaOut);
		}

		// get the output
		vector< string > aln_matrix( aln1.size() );
		for( size_t seqI = 0; seqI < msaOut.GetSeqCount(); seqI++ )
		{
			unsigned indie = msaOut.GetSeqIndex(seqI);
			const char* buf = msaOut.GetSeqBuffer(indie);
			string curseq(buf, msaOut.GetColCount());
			swap(aln_matrix[order[indie]],curseq);

			// debugging, check that sequences came out in the same order they went in!
/*			string inseq = aln1[order[indie]];
			string outseq = aln_matrix[order[indie]];
			stripGaps(inseq);
			stripGaps(outseq);
			if(inseq != outseq)
			{
				unsigned indie = msaOut.GetSeqIndex(seqI);
				cerr << "bad indie " << indie << endl;
				genome::breakHere();
			}
*/
		}
		// fill empty seqs with gaps
		for( size_t seqI = 0; seqI < aln_matrix.size(); seqI++ )
			if(aln_matrix[seqI].size() == 0)
				aln_matrix[seqI].resize(msaOut.GetColCount(), '-');

		aln.SetAlignment( aln_matrix );
		for( uint seqI = 0; seqI < ga1.SeqCount(); seqI++ )
			if( ga1.LeftEnd(seqI) != NO_MATCH )
			{
				aln.SetLeftEnd(seqI, ga1.LeftEnd(seqI));
				aln.SetLength(ga1.Length(seqI), seqI);
			}
		for( uint seqI = 0; seqI < ga2.SeqCount(); seqI++ )
			if( ga2.LeftEnd(seqI) != NO_MATCH )
			{
				aln.SetLeftEnd(seqI, ga2.LeftEnd(seqI));
				aln.SetLength(ga2.Length(seqI), seqI);
			}
		return true;

	}catch( gnException& gne ){
	}catch( exception& e ){
	}catch(...){
	}
	return false;
}


void MuscleInterface::CreateTree( const NumericMatrix<double>& distances, const std::string& tree_filename  )
{
	g_bQuiet.get() = true;			// don't print anything to the console!
	DistFunc df;
	df.SetCount( distances.rows() );
	for( size_t i = 0; i < distances.rows(); i++ )
		for( size_t j = 0; j < distances.rows(); j++ )
			df.SetDist( i, j, distances(i,j) );

	for( size_t i = 0; i < distances.rows(); i++ )
	{
		stringstream ss;
		ss << "seq";
		ss << i + 1;
		df.SetName( i, ss.str().c_str() );
		df.SetId( i, i );
	}
	ClustSetDF csdf( df );
	Clust crusty;
	crusty.Create( csdf, CLUSTER_NeighborJoining );
	Tree tt;
	tt.FromClust( crusty );
	TextFile tf( tree_filename.c_str(), true );
	tt.ToFile( tf );
}


}
