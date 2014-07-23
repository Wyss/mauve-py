/*******************************************************************************
 * $Id: MatchList.cpp,v 1.22 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing, copying, and modification
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/RepeatMatchList.h"
#include "libMems/DNAFileSML.h"
#include "libMems/DNAMemorySML.h"
#include "libMems/MemHash.h"
#include <map>
#include <sstream>
#include <ctime>

using namespace std;
using namespace genome;
namespace mems {

typedef void* MatchID_t;


RepeatMatchList::RepeatMatchList() :MatchList()
{
}


void RepeatMatchList::LoadSequences( ostream* log_stream ){
	
	if( seq_filename.size() == 0 )
		return;

	gnSeqI total_len = 0;
	for( uint seqI = 0; seqI < seq_filename.size(); seqI++ ){
		gnSequence* file_sequence = new gnSequence();
		// Load the sequence and tell the user if it loaded successfully
		try{
			file_sequence->LoadSource( seq_filename[ seqI ] );
		}catch( gnException& gne ){
			delete file_sequence;
			if( gne.GetCode() == FileNotOpened() )
				cerr << "Error loading " << seq_filename[ seqI ] << endl;
			else
				cerr << gne;
			return;
		}catch( exception& e ){
			delete file_sequence;
			cerr << "Unhandled exception loading " << seq_filename[ seqI ] << endl;
			cerr << "At: " << __FILE__ << ":" << __LINE__ << endl;
			cerr << e.what();
			return;
		}catch( ... ){
			delete file_sequence;
			cerr << "Unknown exception when loading " << seq_filename[ seqI ] << endl;
			return;
		}
		
		total_len += file_sequence->length();
		seq_table.push_back( file_sequence );
		if( log_stream != NULL ){
			(*log_stream) << "Sequence loaded successfully.\n";
			(*log_stream) << seq_filename[ seqI ] << " " << file_sequence->length() << " base pairs.\n";
		}
	}

}

void RepeatMatchList::LoadSMLs( uint mer_size, ostream* log_stream ){

	// if the mer_size parameter is 0 then calculate a default mer size for these sequences
	if( mer_size == 0 ){
		mer_size = GetDefaultMerSize( seq_table );
		if( log_stream != NULL ){
			(*log_stream) << "Using weight " << mer_size << " mers for initial seeds\n";
		}
	}

	// load and creates SMLs as necessary
	//punt: tjt
	//uint64 default_seed = getSeed( mer_size );
	uint64 default_seed = getSolidSeed( mer_size );
	vector< uint > create_list;
	uint seqI = 0;
	for( seqI = 0; seqI < seq_table.size(); seqI++ ){
		// define a DNAFileSML to store a sorted mer list
		DNAFileSML* file_sml = new DNAFileSML();
		sml_table.push_back( file_sml );

		boolean success = true;
		try{
			file_sml->LoadFile( sml_filename[ seqI ] );
		}catch( gnException& gne ){
			success = false;
			create_list.push_back( seqI );
		}
		boolean recreate = false;
		if(success && (file_sml->Seed() != default_seed )){
			if( log_stream != NULL )
				(*log_stream) << "Default seed mismatch.  A new sorted mer list will be created.\n";
			recreate = true;
			create_list.push_back( seqI );
		}

		if( success && !recreate && log_stream != NULL )
			(*log_stream) << "Sorted mer list loaded successfully\n";
	}

	// free up memory before creating any SMLs
	if( create_list.size() > 0 )
		for( seqI = 0; seqI < sml_table.size(); seqI++ ){
			sml_table[ seqI ]->Clear();
			delete sml_table[ seqI ];
			sml_table[ seqI ] = NULL;
		}
	
	// create any SMLs that need to be created
	for( uint createI = 0; createI < create_list.size(); createI++ ){
		if( log_stream != NULL )
			(*log_stream) << "Creating sorted mer list\n";
		try{

		time_t start_time = time(NULL);
		sml_table[ create_list[ createI ] ] = new DNAFileSML( sml_filename[ create_list[ createI ] ] );
		sml_table[ create_list[ createI ] ]->Create( *seq_table[ create_list[ createI ] ], default_seed );
		time_t end_time = time(NULL);
	 	if( log_stream != NULL )
			(*log_stream) << "Create time was: " << end_time - start_time << " seconds.\n";
		
		}catch(...){
			cerr << "Error creating sorted mer list\n";
			throw;
		}
	}
	
	// reload the other SMLs now that creation has completed
	if( create_list.size() > 0 ){
		for( seqI = 0; seqI < sml_filename.size(); seqI++ ){
			if( sml_table[ seqI ] != NULL )
				continue;
			sml_table[ seqI ] = new DNAFileSML( sml_filename[ seqI ] );
			try{
				((DNAFileSML*)sml_table[ seqI ])->LoadFile( sml_filename[ seqI ] );
			}catch( gnException& gne ){
				cerr << "Error loading sorted mer list\n";
				throw;
			}
		}
	}
}
void RepeatMatchList::ReadList(istream& match_file){
	string tag;
	gnSeqI len;
	int64 start;
	unsigned int seq_count;
	
	match_file >> tag;	//format version tag
	if( tag != "FormatVersion" ){
		Throw_gnEx(InvalidFileFormat());
	}
	match_file >> tag;	//format version
	if( tag != "3" ){
		Throw_gnEx(InvalidFileFormat());
	}
	match_file >> tag;	//sequence count tag
	if( tag != "SequenceCount" ){
		Throw_gnEx(InvalidFileFormat());
	}
	match_file >> seq_count;	//sequence count
	if(seq_count < 2){
		Throw_gnEx(InvalidFileFormat());
	}
	
	// read the sequence file names and lengths
	for( unsigned int seqI = 0; seqI < seq_count; seqI++ ){		
		match_file >> tag;	// name tag
		getline( match_file, tag );
		// skip the tab character
		tag = tag.substr( 1 );
		seq_filename.push_back(tag);
//		try{
//			gnSequence *new_seq = new gnSequence();
//			new_seq->LoadSource(tag);
//			seq_table.push_back( new_seq );
//		}catch( gnException& gne );
		match_file >> tag;	// length tag
		gnSeqI seq_len;
		match_file >> seq_len;	// length
		if( seqI < seq_table.size() )
			if( seq_table[ seqI ]->length() != seq_len ){
				cerr << "Warning: Genome sizes in the match list differ.\n";
				cerr << "seq_table[ " << seqI << " ]->length() " << seq_table[ seqI ]->length() << " seq_len: " << seq_len << endl;
			}
	}

	// read the number of matches
	unsigned int match_count;
	match_file >> tag;	// match count tag
	match_file >> match_count;	// match count
		
	// read the matches
	map< MatchID_t, Match* > match_map;
	string cur_line;
	getline( match_file, cur_line );
	while( getline( match_file, cur_line ) ){
		MatchHashEntry mhe( seq_count, 0 );
		stringstream line_stream( cur_line );
		
		line_stream >> len;
		mhe.SetLength(len);

		for(uint32 seqI = 0; seqI < seq_count; seqI++){
			line_stream >> start;
			mhe.SetStart(seqI, start);
		}
		
		mhe.CalculateOffset();
		
		MatchID_t match_id;
		line_stream >> match_id;
		
		uint sub_count;
		boolean bad_stream = false;
		line_stream >> sub_count;
		if(sub_count > 0)
			throw "Unable to read file, invalid format, cannot read subset information\n";

		if( bad_stream )
			break;

		uint sup_count;
		line_stream >> sup_count;
		if(sup_count > 0)
			throw "Unable to read file, invalid format, cannot read superset information\n";
		if( bad_stream )
			break;
		
		Match* new_match = mhe.Copy();
		push_back( new_match );
		match_map.insert( map< MatchID_t, Match* >::value_type( match_id, new_match ));
	}
	if( match_count != size() ){
		Throw_gnEx(InvalidFileFormat());
	}
}

void RepeatMatchList::WriteList(ostream& match_file) const{
	if( size() == 0 )
		return;
	Match* first_mem = *(begin());
	unsigned int seq_count = first_mem->SeqCount();

	match_file << "FormatVersion" << '\t' << 3 << "\n";
	match_file << "SequenceCount" << '\t' << seq_count << "\n";
	for(unsigned int seqI = 0; seqI < seq_count; seqI++){
		match_file << "Sequence" << seqI << "File" << '\t';
		if( seq_filename.size() > seqI )
			match_file << seq_filename[seqI];
		else
			match_file << "null";
		match_file << "\n";
		match_file << "Sequence" << seqI << "Length" << '\t';
		if( seq_table.size() > seqI )
			match_file << seq_table[seqI]->length();
		else
			match_file << "0";
		match_file << "\n";
	}

	match_file << "MatchCount" << '\t' << size() << endl;

	//get all the mems out of the hash table and write them out
    vector<Match*>::const_iterator match_iter;
	match_iter = begin();
	set<Match*> cur_set;
	set<Match*>::iterator set_iter;
	for(; match_iter != end(); match_iter++){
		// print the match
		match_file << **match_iter << '\t';

		// print the Multiplicity
		match_file << (*match_iter)->Multiplicity() << '\t';

		// print the match address
		match_file << (MatchID_t)(*match_iter) << '\t';
				
		// print subset id's
		match_file << 0;

		// print superset id's
		match_file << '\t' << 0;
		match_file << endl;
	}
}

}	// namespace mems
