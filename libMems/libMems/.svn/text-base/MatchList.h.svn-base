/*******************************************************************************
 * $Id: MatchList.h,v 1.10 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _MatchList_h_
#define _MatchList_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <list>
#include "libMems/SortedMerList.h"
#include "libMems/DNAFileSML.h"
#include "libMems/DNAMemorySML.h"
#include "libGenome/gnSequence.h"
#include "libMems/Match.h"
#include "libMems/gnRAWSequence.h"
#include "libGenome/gnRAWSource.h"
#include "libMems/Files.h"
#include <sstream>
#include <map>
#include <ctime>

namespace mems {

template< typename MatchPtrType >
class GenericMatchList : public std::vector< MatchPtrType > 
{
public:
	GenericMatchList(){};
	GenericMatchList( const GenericMatchList& ml );
	GenericMatchList& operator=( const GenericMatchList& ml );


	/**
	 * Attempts to load SMLs designated by the
	 * elements of the sml_filename vector.  This
	 * method will create the sorted mer lists if they do not exist.
	 * The DNAFileSML objects are created on the heap
	 * and are not deallocated when this class is destroyed.  They should
	 * be manually destroyed when no longer in use.
	 * @param	seed_rank	The rank of the seed to use, 0-2 are ranked spaced seeds, 
	 *						other options include CODING_SEED and SOLID_SEED
	 */
	void LoadSMLs( uint mer_size, std::ostream* log_stream, int seed_rank = 0, bool solid = false, bool force_recreate = false );

	/**
	 * Loads sequences to align from a Multi-FastA file and constructs a SML
	 * for each sequence entry in the file.
	 * The genome::gnSequence and SortedMerList objects are created on the heap
	 * and are not deallocated when this class is destroyed.  They should
	 * be manually destroyed when no longer in use.
	 *
	 * @param mer_size		The seed size to use when constructing the sorted mer lists
	 * @param log_stream	An output stream to log messages to.  If NULL no logging is done
	 * @param load_smls		Specifies whether sorted mer lists should be created 
	 * 						for each sequence entry
	 */
	void CreateMemorySMLs( uint mer_size, std::ostream* log_stream, int seed_rank = 0 );

	/**
	 * Calculates a default search mer size for the given set of sequences
	 * @param seq_table		The vector of sequences to calculate a default mer size for
	 */
	static uint GetDefaultMerSize( const std::vector< genome::gnSequence* >& seq_table );
	
	/**
	 * Deletes the genome::gnSequence, SortedMerList, and Match objects associated
	 * with this GenericMatchList.
	 */
	void Clear();
	
	/**
	 * Removes all matches that have a multiplicity lower than the specified level
	 * @param mult	The multiplicity filter threshold
	 */
	void MultiplicityFilter( unsigned mult );

	/**
	 * Removes all matches that shorter than the specified length
	 * @param length	The minimum length
	 */
	void LengthFilter( gnSeqI length );

	/**
	 * Removes matches that do not match in exactly the sequences specified in filter_spec
	 * @param filter_spec 	The specification of the exact filter, true designates that the
	 *						match must exist in that sequence.  filter_spec must contain
	 *						one boolean entry for every sequence.
	 */
//	void ExactFilter( valarray< bool >& filter_spec );
	/**
	 * Removes matches that do not intersect with the sequences specified in filter_spec
	 * @param filter_spec 	The specification of the intersection filter, true designates
	 *						match must exist in that sequence.  filter_spec must contain
	 *						one boolean entry for every sequence.
	 */
//	void IntersectFilter( valarray< bool >& filter_spec );

	
	std::vector<std::string> sml_filename;		/**< The file names of the sorted mer list for each sequence, may be empty or null */
	std::vector<std::string> seq_filename;		/**< The file names of the sequence data, may be empty or null */
	std::vector<SortedMerList*> sml_table;	/**< The sorted mer list associated with each sequence, may be empty or null */
	std::vector<genome::gnSequence*> seq_table;		/**< The actual sequences associated with the matches stored in this list.  Should not be empty or null. */

protected:

};

typedef GenericMatchList< Match* > MatchList;

CREATE_EXCEPTION( InvalidArgument );

/**
 * Thrown when a file being read is invalid
 */
CREATE_EXCEPTION(InvalidFileFormat)


/**
 * Reads a GenericMatchList from an input stream
 * Sequence and SML file names are read into the seq_filename
 * and sml_filename vectors, but the actual files are not
 * opened.  The calling function should load them after
 * using this method.
 * @param match_stream The input stream to read from
 */
void ReadList( MatchList& mlist, std::istream& match_stream );

/**
 *  Writes a GenericMatchList to the designated output stream
 * @param match_stream The output stream to write to
 */
void WriteList( const MatchList& mlist, std::ostream& match_stream );

typedef void* MatchID_t;

template< typename MatchPtrType >
GenericMatchList< MatchPtrType >::GenericMatchList( const GenericMatchList< MatchPtrType >& ml ){
	*this = ml;
}

template< typename MatchPtrType >
GenericMatchList< MatchPtrType >& GenericMatchList< MatchPtrType >::operator=( const GenericMatchList< MatchPtrType >& ml ){
	std::vector< MatchPtrType >::operator=( ml );
	sml_filename = ml.sml_filename;
	seq_filename = ml.seq_filename;
	sml_table = ml.sml_table;
	seq_table = ml.seq_table;
	return *this;
}

/**
 * Attempts to load the sequences designated by the
 * elements of the seq_filename vector.
 * The genome::gnSequence objects are created on the heap
 * and are not deallocated when this class is destroyed.  They should
 * be manually destroyed when no longer in use.
 */
template< typename MatchListType >
void LoadSequences( MatchListType& mlist, std::ostream* log_stream ){
	
	if( mlist.seq_filename.size() == 0 )
		return;

	for( uint seqI = 0; seqI < mlist.seq_filename.size(); seqI++ ){
		genome::gnSequence* file_sequence = new genome::gnSequence();
		// Load the sequence and tell the user if it loaded successfully
		try{
			file_sequence->LoadSource( mlist.seq_filename[ seqI ] );
		}catch( genome::gnException& gne ){
			delete file_sequence;
			if( gne.GetCode() == genome::FileNotOpened() )
				std::cerr << "Error loading " << mlist.seq_filename[ seqI ] << std::endl;
			else
				std::cerr << gne;
			return;
		}catch( std::exception& e ){
			delete file_sequence;
			std::cerr << "Unhandled exception loading " << mlist.seq_filename[ seqI ] << std::endl;
			std::cerr << "At: " << __FILE__ << ":" << __LINE__ << std::endl;
			std::cerr << e.what();
			return;
		}catch( ... ){
			delete file_sequence;
			std::cerr << "Unknown exception when loading " << mlist.seq_filename[ seqI ] << std::endl;
			return;
		}
		
		mlist.seq_table.push_back( file_sequence );
		if( log_stream != NULL ){
			(*log_stream) << "Sequence loaded successfully.\n";
			(*log_stream) << mlist.seq_filename[ seqI ] << " " << file_sequence->length() << " base pairs.\n";
		}
	}

}

/**
 * Loads the sequences designated by the elements of the seq_filename vector and
 * creates temporary RAW sequence files.  The resulting gnSequences are gnRAWSequences.
 * The genome::gnRAWSequence objects are created on the heap
 * and are not deallocated when this class is destroyed.  They should
 * be manually destroyed when no longer in use.
 */
template< typename MatchListType >
void LoadAndCreateRawSequences( MatchListType& mlist, std::ostream* log_stream ){
	
	if( mlist.seq_filename.size() == 0 )
		return;

	for( uint seqI = 0; seqI < mlist.seq_filename.size(); seqI++ ){
		genome::gnSequence* file_sequence = new genome::gnSequence();
		// Load the sequence and tell the user if it loaded successfully
		try{
			file_sequence->LoadSource( mlist.seq_filename[ seqI ] );
		}catch( genome::gnException& gne ){
			delete file_sequence;
			if( gne.GetCode() == genome::FileNotOpened() )
				std::cerr << "Error loading " << mlist.seq_filename[ seqI ] << std::endl;
			else
				std::cerr << gne;
			return;
		}catch( std::exception& e ){
			delete file_sequence;
			std::cerr << "Unhandled exception loading " << mlist.seq_filename[ seqI ] << std::endl;
			std::cerr << "At: " << __FILE__ << ":" << __LINE__ << std::endl;
			std::cerr << e.what();
			return;
		}catch( ... ){
			delete file_sequence;
			std::cerr << "Unknown exception when loading " << mlist.seq_filename[ seqI ] << std::endl;
			return;
		}

		// now create a temporary raw sequence
		std::string tmpfilename = "rawseq";
		tmpfilename = CreateTempFileName("rawseq");
		genome::gnRAWSource::Write( *file_sequence, tmpfilename );
		delete file_sequence;
		registerFileToDelete( tmpfilename );

		if( log_stream != NULL )
			(*log_stream) << "Storing raw sequence at " << tmpfilename << std::endl;	
		genome::gnRAWSequence* raw_seq = new genome::gnRAWSequence( tmpfilename );
		mlist.seq_table.push_back( raw_seq );
		if( log_stream != NULL ){
			(*log_stream) << "Sequence loaded successfully.\n";
			(*log_stream) << mlist.seq_filename[ seqI ] << " " << raw_seq->length() << " base pairs.\n";
		}
	}
}


template< typename MatchPtrType >
void GenericMatchList< MatchPtrType >::LoadSMLs( uint mer_size, std::ostream* log_stream, int seed_rank, bool solid, bool force_create ){

	// if the mer_size parameter is 0 then calculate a default mer size for these sequences
	if( mer_size == 0 ){
		mer_size = GetDefaultMerSize( seq_table );
		if( log_stream != NULL ){
			(*log_stream) << "Using weight " << mer_size << " mers for initial seeds\n";
		}
	}

	// load and creates SMLs as necessary
	uint64 default_seed = getSeed( mer_size, seed_rank );
	if (solid)
		uint64 default_seed = getSolidSeed( mer_size );
	std::vector< uint > create_list;
	uint seqI = 0;
	for( seqI = 0; seqI < seq_table.size(); seqI++ ){
		// define a DNAFileSML to store a sorted mer list
		DNAFileSML* file_sml = new DNAFileSML();
		sml_table.push_back( file_sml );

		boolean success = true;
		try{
			file_sml->LoadFile( sml_filename[ seqI ] );
		}catch( genome::gnException& gne ){
			success = false;
			create_list.push_back( seqI );
		}
		boolean recreate = false;
		if(success && force_create){
			if( log_stream != NULL )
				(*log_stream) << "SML exists, but forcefully recreating.  A new sorted mer list will be created.\n";
			recreate = true;
			create_list.push_back( seqI );
		}
		else if(success && (file_sml->Seed() != default_seed )){
			if( log_stream != NULL )
				(*log_stream) << "Default seed mismatch.  A new sorted mer list will be created.\n";
			recreate = true;
			create_list.push_back( seqI );
		}

		if( success && !recreate && log_stream != NULL && !force_create )
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
			std::cerr << "Error creating sorted mer list\n";
			throw;
		}
	}
	
	// reload the other SMLs now that creation has completed
	if( create_list.size() > 0 ){
		for( seqI = 0; seqI < seq_filename.size(); seqI++ ){
			if( sml_table[ seqI ] != NULL )
				continue;
			sml_table[ seqI ] = new DNAFileSML( sml_filename[ seqI ] );
			try{
				((DNAFileSML*)sml_table[ seqI ])->LoadFile( sml_filename[ seqI ] );
			}catch( genome::gnException& gne ){
				std::cerr << "Error loading sorted mer list\n";
				throw;
			}
		}
	}
}

template< typename MatchPtrType >
uint GenericMatchList< MatchPtrType >::GetDefaultMerSize( const std::vector< genome::gnSequence* >& seq_table ){
	gnSeqI total_len = 0;
	for( uint seqI = 0; seqI < seq_table.size(); seqI++ )
		total_len += seq_table[ seqI ]->length();
	return getDefaultSeedWeight( total_len / seq_table.size() );
}


/**
 * Loads sequences to align from a Multi-FastA file 
 * The genome::gnSequence and SortedMerList objects are created on the heap
 * and are not deallocated when this class is destroyed.  They should
 * be manually destroyed when no longer in use.
 *
 * @param mfa_filename	The name of the Multi-FastA file to read in.  Each 
 *						sequence entry will be treated as a separate sequence to 
 *						be aligned.
 * @param log_stream	An output stream to log messages to.  If NULL no logging is done
 */
template< typename MatchListType >
void LoadMFASequences( MatchListType& mlist, const std::string& mfa_filename, std::ostream* log_stream ) {
	genome::gnSequence file_sequence;
	// Load the sequence and tell the user if it loaded successfully
	try{
		file_sequence.LoadSource( mfa_filename );
	}catch( genome::gnException& gne ){
		if( gne.GetCode() == genome::FileNotOpened() )
			std::cerr << "Error loading " << mfa_filename << std::endl;
		else
			std::cerr << gne;
		return;
	}catch( std::exception& e ){
		std::cerr << "Unhandled exception loading " << mfa_filename << std::endl;
		std::cerr << "At: " << __FILE__ << ":" << __LINE__ << std::endl;
		std::cerr << e.what();
		return;
	}catch( ... ){
		std::cerr << "Unknown exception when loading " << mfa_filename << std::endl;
		return;
	}

	mlist.seq_filename.clear();
	gnSeqI total_len = 0;
	for( uint contigI = 0; contigI < file_sequence.contigListSize(); contigI++ ){
		genome::gnSequence* contig_seq = new genome::gnSequence( file_sequence.contig( contigI ) );
		mlist.seq_filename.push_back( mfa_filename );
//		mlist.seq_filename.push_back( file_sequence.contigName( contigI ) );
		if( log_stream != NULL ){
			(*log_stream) << "Sequence loaded successfully.\n";
			(*log_stream) << mlist.seq_filename[ contigI ] << " " << contig_seq->length() << " base pairs.\n";
		}
		mlist.seq_table.push_back( contig_seq );
	}
}

template< typename MatchPtrType >
void GenericMatchList< MatchPtrType >::CreateMemorySMLs( uint mer_size, std::ostream* log_stream, int seed_rank ) 
{
	// if the mer_size parameter is 0 then calculate a default mer size for these sequences
	if( mer_size == 0 ){
		mer_size = GetDefaultMerSize( seq_table );
		if( log_stream != NULL ){
			(*log_stream) << "Using " << mer_size << "-mers for initial seeds\n";
		}
	}

	uint64 default_seed = getSeed( mer_size, seed_rank );

	// define a DNAMemorySML to store a sorted mer list
	for( uint contigI = 0; contigI < seq_table.size(); contigI++ )
	{
		DNAMemorySML* contig_sml = new DNAMemorySML();
		boolean success = true;
		if( log_stream != NULL )
			(*log_stream) << "Creating sorted mer list\n";
		time_t start_time = time(NULL);
		contig_sml->Create( *seq_table[contigI], default_seed );
		time_t end_time = time(NULL);
	 	if( log_stream != NULL )
			(*log_stream) << "Create time was: " << end_time - start_time << " seconds.\n";
		
		sml_table.push_back( contig_sml );
	}
}

template< typename MatchPtrType >
void GenericMatchList< MatchPtrType >::Clear() {
	for( uint seqI = 0; seqI < seq_table.size(); seqI++ ){
		if( seq_table[ seqI ] != NULL )
			delete seq_table[ seqI ];
	}
	for( uint seqI = 0; seqI < sml_table.size(); seqI++ ){
		if( sml_table[ seqI ] != NULL )
			delete sml_table[ seqI ];
	}
	typename std::vector<MatchPtrType>::iterator match_iter = this->begin();
	for(; match_iter != this->end(); match_iter++ ){
		(*match_iter)->Free();
		(*match_iter) = NULL;
	}
	seq_table.clear();
	sml_table.clear();
	this->clear();
	seq_filename.clear();
	sml_filename.clear();
}

/**
 * Use this to update linkage pointers after copying an entire set of Matches
 */
template< class FromType, class ToType, class MatchListType >
void RemapSubsetMatchAddresses( std::map<FromType, ToType>& old_to_new_map, MatchListType& match_list );


template< class FromType, class ToType, class MatchListType >
void RemapSubsetMatchAddresses( std::map<FromType, ToType>& old_to_new_map, MatchListType& match_list )
{
	// now remap the subset and superset links
	typename MatchListType::iterator match_iter = match_list.begin();
	//typedef typename MatchListType::value_type MatchType;
	//typedef typename Match MatchType;
	typename std::map<FromType, ToType>::iterator map_iter;
	for(; match_iter != match_list.end(); ++match_iter ){
		// remap all subsets
		std::set< Match* >& subsets = (*match_iter)->Subsets();
		std::set< Match* > new_subsets;
		std::set< Match* >::iterator sub_iter = subsets.begin();
		for(; sub_iter != subsets.end(); ++sub_iter ){
			map_iter = old_to_new_map.find( (FromType)*sub_iter );
			new_subsets.insert( map_iter->second );
		}
		subsets = new_subsets;

		// remap all supersets
		std::set< Match* >& supersets = (*match_iter)->Supersets();
		std::set< Match* > new_supersets;
		std::set< Match* >::iterator super_iter = supersets.begin();
		for(; super_iter != supersets.end(); ++super_iter ){
			map_iter = old_to_new_map.find( (FromType)*super_iter );
			new_supersets.insert( map_iter->second );
		}
		supersets = new_supersets;
	}
}

inline
void ReadList(MatchList& mlist, std::istream& match_file)
{
	std::string tag;
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
		std::getline( match_file, tag );
		// skip the tab character
		tag = tag.substr( 1 );
		mlist.seq_filename.push_back(tag);
		match_file >> tag;	// length tag
		gnSeqI seq_len;
		match_file >> seq_len;	// length
		if( seqI < mlist.seq_table.size() )
			if( mlist.seq_table[ seqI ]->length() != seq_len ){
				std::cerr << "Warning: Genome sizes in the match list differ.\n";
				std::cerr << "seq_table[ " << seqI << " ]->length() " << mlist.seq_table[ seqI ]->length() << " seq_len: " << seq_len << std::endl;
			}
	}

	// read the number of matches
	unsigned int match_count;
	match_file >> tag;	// match count tag
	match_file >> match_count;	// match count
		
	// read the matches
	std::map< MatchID_t, Match* > match_map;
	std::string cur_line;
	std::getline( match_file, cur_line );
	while( getline( match_file, cur_line ) ){
		Match mhe( seq_count );
		std::stringstream line_stream( cur_line );
		
		line_stream >> len;
		mhe.SetLength(len);

		for(uint32 seqI = 0; seqI < seq_count; seqI++){
			line_stream >> start;
			mhe.SetStart(seqI, start);
		}
		
		MatchID_t match_id;
		line_stream >> match_id;
		
		uint sub_count;
		boolean bad_stream = false;
		line_stream >> sub_count;
		if(sub_count > 0)
			throw "Unable to read file, invalid format, cannot read subset data\n";

		if( bad_stream )
			break;

		uint sup_count;
		line_stream >> sup_count;
		if(sub_count > 0)
			throw "Unable to read file, invalid format, cannot read superset data\n";
		if( bad_stream )
			break;
		
		Match* new_match = mhe.Copy();
		mlist.push_back( new_match );
		match_map.insert( std::map< MatchID_t, Match* >::value_type( match_id, new_match ));
	}
	if( match_count != mlist.size() ){
		Throw_gnEx(InvalidFileFormat());
	}
}

inline
void WriteList( const MatchList& mlist, std::ostream& match_file)
{
	if( mlist.size() == 0 )
		return;
	Match* first_mem = *(mlist.begin());
	unsigned int seq_count = first_mem->SeqCount();

	match_file << "FormatVersion" << '\t' << 3 << "\n";
	match_file << "SequenceCount" << '\t' << seq_count << "\n";
	for(unsigned int seqI = 0; seqI < seq_count; seqI++){
		match_file << "Sequence" << seqI << "File" << '\t';
		if( mlist.seq_filename.size() > seqI )
			match_file << mlist.seq_filename[seqI];
		else
			match_file << "null";
		match_file << "\n";
		match_file << "Sequence" << seqI << "Length" << '\t';
		if( mlist.seq_table.size() > seqI )
			match_file << mlist.seq_table[seqI]->length();
		else
			match_file << "0";
		match_file << "\n";
	}

	match_file << "MatchCount" << '\t' << mlist.size() << std::endl;

	//get all the mems out of the hash table and write them out
	std::vector<Match*>::const_iterator match_iter;
	match_iter = mlist.begin();
	std::set<Match*> cur_set;
	std::set<Match*>::iterator set_iter;
	for(; match_iter != mlist.end(); match_iter++){
		// print the match
		match_file << **match_iter << '\t';

		// print the match address
		match_file << (MatchID_t)(*match_iter) << '\t';
		
		// print subset id's
		match_file << 0;

		// print superset id's
		match_file << '\t' << 0;
		match_file << std::endl;
	}
}

template< typename MatchPtrType >
void GenericMatchList< MatchPtrType >::MultiplicityFilter( unsigned mult ){

	size_t cur = 0;
	for( uint memI = 0; memI < this->size(); memI++ ){
		if( (*this)[ memI ]->Multiplicity() == mult )
			(*this)[cur++] = (*this)[memI];
		else{
			(*this)[ memI ]->Free();
			(*this)[ memI ] = NULL;
		}
	}
	this->resize(cur);
}

template< typename MatchPtrType >
void GenericMatchList< MatchPtrType >::LengthFilter( gnSeqI length ){

	size_t cur = 0;
	for( size_t memI = 0; memI < this->size(); memI++ ){
		if( (*this)[ memI ]->Length() >= length )
			(*this)[cur++] = (*this)[memI];
		else{
			(*this)[ memI ]->Free();
			(*this)[ memI ] = NULL;
		}
	}
	this->resize(cur);
}

}

#endif	//_MatchList_h_
