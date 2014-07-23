/*******************************************************************************
 * $Id: FileSML.cpp,v 1.22 2004/04/26 21:13:58 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing, copying, and modification
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/


#include "libMems/FileSML.h"
// for CreateTempFileName():
#include "libMems/Aligner.h"
#include "libGenome/gnFilter.h"
#include "libGenome/gnRAWSource.h"
#include <algorithm>
#include <cmath>
#include "boost/filesystem/operations.hpp"

using namespace std;
using namespace genome;
namespace mems {

FileSML& FileSML::operator=(const FileSML& sa)
{
 	SortedMerList::operator=( sa );
 	filename = sa.filename;
	sarray_start_offset = sa.sarray_start_offset;
	seq_coords = sa.seq_coords;
	sarfile.open(filename.c_str(), ios::binary | ios::in );
	if(!sarfile.is_open()){
		DebugMsg("FileSML::=: Unable to open suffix array file.\n");
		sarfile.clear();
		return *this;
	}
	return *this;
}

void FileSML::Clear() {
	SortedMerList::Clear();
	filename = "";
	sarfile.close();
	sarray_start_offset = 0;
	seq_coords.clear();
}

void FileSML::LoadFile(const string& fname){
	filename = fname;
	sarfile.open(fname.c_str(), ios::binary | ios::in );
	if(!sarfile.is_open()){
		sarfile.clear();
		Throw_gnExMsg( FileNotOpened(), "Unable to open file.\n");
	}
	// read the header into a temporary header struct just
	// in case it's bogus
	SMLHeader tmp_header;
	sarfile.read((char*)&tmp_header, sizeof(struct SMLHeader));
	if(sarfile.gcount() < (int)sizeof(struct SMLHeader)){
		sarfile.clear();
		Throw_gnExMsg( FileUnreadable(), "Unable to read file.");
	}
	if(tmp_header.version != FormatVersion()){
		Throw_gnExMsg( FileUnreadable(), "Unsupported file format.");
	}
	header = tmp_header;
	
	SetMerMaskSize( header.seed_weight );
	seed_mask = mer_mask;
	SetMerMaskSize( header.seed_length );

	//header is ok.  read the sequence.
	gnSeqI seq_len = header.length;
	if(header.circular)
		seq_len += header.seed_length - 1;
	binary_seq_len = ((uint64)seq_len * (uint64)header.alphabet_bits) / 32;
	if(((uint64)seq_len * (uint64)header.alphabet_bits) % 32 != 0)
		binary_seq_len++;
	binary_seq_len+=2;	//fix for access violations.

	if(sequence != NULL)
		delete[] sequence;
	sequence = new uint32[binary_seq_len];

	sarfile.read((char*)sequence, binary_seq_len*sizeof(uint32));
	if(sarfile.gcount() < (int64)(binary_seq_len*sizeof(uint32))){
		sarfile.clear();
		Throw_gnExMsg( FileUnreadable(), "Error reading sequence data.");
	}

	sarray_start_offset = sarfile.tellg();
	sarfile.seekg(sarray_start_offset + sizeof(gnSeqI) * header.length);
	if(!sarfile.good()){
		sarfile.clear();
		Throw_gnExMsg( FileUnreadable(), "Premature end of file.");
	}
	filename = fname;

	// create a memory-map to the data of interest
	sardata.open(fname);
	
	// check whether there is a .coords mask file to read
	string coordfile = filename + ".coords";
	ifstream coord_in( coordfile.c_str() );
	if( coord_in.is_open() ){
		seq_coords.clear();
		int64 cur_coord;
		while( coord_in >> cur_coord ){
			seq_coords.push_back( cur_coord );
		}
	}
}

void FileSML::OpenForWriting( boolean truncate ){
	// Open smlfile for writing
	boolean was_open = sarfile.is_open();
	if(was_open)
		sarfile.close();
	if( truncate )
		sarfile.open(filename.c_str(), ios::binary | ios::in | ios::out | ios::trunc );
	else
		sarfile.open(filename.c_str(), ios::binary | ios::in | ios::out );
	if(!sarfile.is_open() || !sarfile.good()){
		sarfile.clear();
		if(was_open)
			sarfile.open(filename.c_str(), ios::binary | ios::in );
		Throw_gnExMsg(FileNotOpened(), "Unable to open file for writing.");
	}
}

boolean FileSML::WriteHeader(){
	if(!sarfile.is_open()){
		Throw_gnExMsg(IOStreamFailed(), "File is not valid.");
	}
	boolean success = true;
	const char* errormsg = "";
	// Open sarfile for writing and write new header.
	OpenForWriting( false );
	sarfile.write((char*)&header, sizeof(struct SMLHeader));
	if(!sarfile.good()){
		errormsg = "Error writing header to disk.";
		success = false;
	}

	// reopen the sorted mer list file read-only
	sarfile.close();
	sarfile.open(filename.c_str(), ios::binary | ios::in );
	if(!sarfile.is_open()){
		errormsg = "Error opening sorted mer list file.";
		success = false;
	}
	if(!success)
		Throw_gnExMsg(IOStreamFailed(), errormsg);
	return success;
}

gnSeqI FileSML::UniqueMerCount(){
	gnSeqI tmp_count = header.unique_mers;
	SortedMerList::UniqueMerCount();
	if(tmp_count != header.unique_mers)
		WriteHeader();
	return header.unique_mers;
}

//change the description in memory and on disk
void FileSML::SetDescription(const string& d){
	strncpy(header.description, d.c_str(), DESCRIPTION_SIZE-1);
	WriteHeader();
} 

void FileSML::SetID(const sarID_t d){
	header.id = d;
	WriteHeader();
}


extern "C" {
#include "libMems/dmSML/dmsort.h"
}

char** FileSML::tmp_paths = NULL;

void FileSML::registerTempPath( const string& path ) {
	string tmp_path = path;
	// add trailing path separator if necessary
#ifdef WIN32
	if( tmp_path[ tmp_path.size() - 1 ] != '\\' )
		tmp_path += "\\";
#else
	if( tmp_path[ tmp_path.size() - 1 ] != '/' )
		tmp_path += "/";
#endif
		
	if( tmp_paths == NULL ){
		tmp_paths = new char*[1];
		tmp_paths[ 0 ] = NULL;
	}
	
	int path_count = 0;
	while( tmp_paths[ path_count ] != NULL )
		path_count++;

	// create a new array with room for another element
	char** tmp_tmp_paths = new char*[ path_count+1 ];
	// copy old elements
	for( int pathI = 0; pathI < path_count; pathI++ )
		tmp_tmp_paths[ pathI ] = tmp_paths[ pathI ];
	// add new element
	tmp_tmp_paths[ path_count ] = new char[ tmp_path.size() + 1 ];
	strncpy( tmp_tmp_paths[ path_count ], tmp_path.c_str(), tmp_path.size() + 1 );
	// set null terminator element
	tmp_tmp_paths[ path_count + 1 ] = NULL;
	
	// set new paths
	char** old_paths = tmp_paths;
	tmp_paths = tmp_tmp_paths;
	
	// free old array
	delete[] old_paths;
}

const char* FileSML::getTempPath( int pathI ){
	return tmp_paths[ pathI ];
}

int FileSML::getTempPathCount(){
	int path_count = 0;
	while( tmp_paths && tmp_paths[ path_count ] != NULL )
		path_count++;
	return path_count;
}


void maskNNNNN( const gnSequence& in_seq, gnSequence& out_seq, vector< int64 >& seq_coords, int mask_n_length ) {
	
	gnSeqI seqI = 1;
	gnSeqI read_length = 1024*1024;
	string cur_seq;
	gnSeqI n_count = 0;
	gnSeqI n_stretch_start = 0;
	gnSeqI n_stretch_end = 1;

	while( seqI <= in_seq.length() ){
		read_length = seqI + read_length < in_seq.length() ? read_length : in_seq.length() - seqI + 1;
		in_seq.ToString( cur_seq, read_length, seqI );
		
		uint charI = 0;
		for( ; charI < cur_seq.size(); charI++ ){
			if( cur_seq[ charI ] == 'N' || cur_seq[ charI ] == 'n' ){
				if( n_count == 0 ){
					n_stretch_start = seqI + charI;
				}
				n_count++;
			}else{
				if( n_count > mask_n_length ){
					if( n_stretch_start - n_stretch_end != 0 ){
						// Add the sequence region to the output sequence
						out_seq += in_seq.subseq( n_stretch_end, n_stretch_start - n_stretch_end );
						// add the masked coordinates
						seq_coords.push_back( n_stretch_end );
						seq_coords.push_back( n_stretch_start - 1 );
					}
					// update n_stretch_end to the first non N character
					n_stretch_end = seqI + charI;
				}
				n_count = 0;
			}
		}
		seqI += read_length;
	}
	out_seq += in_seq.subseq( n_stretch_end, seqI - n_stretch_end );

	// add the masked coordinates
	seq_coords.push_back( n_stretch_end );
	seq_coords.push_back( seqI - 1 );
}

	// use dmSML to construct the SML
	// then read it in using LoadFile()
void FileSML::dmCreate(const gnSequence& seq, const uint64 seed){
	// Filter NNNNNs
	gnSequence masked_seq;
	seq_coords.clear();
	maskNNNNN( seq, masked_seq, seq_coords, 0 );
	
	// write a raw sequence to a tmp file stored in the first scratch path
	string rawfile = CreateTempFileName("dm_rawseq");
	gnRAWSource::Write( masked_seq, rawfile.c_str() );
	
	// write a sequence coordinate file
	if( seq_coords.size() > 0 ){
		string coordfile = filename + ".coords";
		ofstream coord_out( coordfile.c_str() );
		if( !coord_out.is_open() ){
			cerr << "Could not open " << coordfile << endl;
			throw "";
		}
		
		for( int coordI = 0; coordI < seq_coords.size(); coordI+=2 ){
			coord_out << seq_coords[ coordI ] << '\t' << seq_coords[ coordI + 1 ] << endl;
		}
		coord_out.close();
	}
	
	
	// run dmSML
	const char* const* scratch_paths = (const char* const*)tmp_paths;
	sarfile.close();
	int rval = dmSML( rawfile.c_str(), filename.c_str(), scratch_paths, seed );
	if( rval != 0 )
		cerr << "Crap.  It's broke, return value " << rval << endl;
	
	boost::filesystem::remove( rawfile );
	// load the sorted mer list
	LoadFile( filename );
}

void FileSML::Create(const gnSequence& seq, const uint64 seed){

	vector<bmer> sml_array;
	bool is_spaced_seed = getSeedWeight(seed) != getSeedLength(seed);
	OpenForWriting( true );

	try{
		SortedMerList::Create( seq, seed );
		
		if( is_spaced_seed )
			FillDnaSeedSML(seq, sml_array);
		else
			FillSML(seq, sml_array);

	}catch(...){	
		// if there was a memory allocation error then
		// try using dmSML to do an external sort
		sarfile.clear();
		sarfile.close();
		sarfile.clear();
		if( sequence != NULL )
			delete[] sequence;
		binary_seq_len = 0;

		dmCreate( seq, seed );
	}

//	RadixSort(s_array);
	sort(sml_array.begin(), sml_array.end(), &bmer_lessthan);
	
	/* now write out the file header */
	sarfile.write((char*)&header, sizeof(struct SMLHeader));

	if(!sarfile.good()){
		sarfile.clear();
		Throw_gnExMsg( IOStreamFailed(), "Error writing sorted mer list header to disk.\n");
	}

	/* write out the actual sequence */
	sarfile.write((char*)sequence, binary_seq_len*sizeof(uint32));
	sarray_start_offset = sarfile.tellg();

	/* write out the sorted mer list */
	for(gnSeqI suffixI=0; suffixI < sml_array.size(); suffixI++)
		sarfile.write((char*)&(sml_array[suffixI].position), sizeof(smlSeqI_t));
	
	sarfile.flush();
	if(!sarfile.good()){
		sarfile.clear();
		Throw_gnExMsg( IOStreamFailed(), "Error writing sorted mer list to disk.\n");
	}
	// reopen the sorted mer list file read-only
	sarfile.close();
	sarfile.open(filename.c_str(), ios::binary | ios::in );
	if(!sarfile.is_open())
		Throw_gnExMsg( FileNotOpened(), "FileSML::Create: Error opening sorted mer list file.\n");

	sardata.open(filename);
}

bmer FileSML::operator[](gnSeqI index)
{
	bmer tmp_mer;
	tmp_mer.position = base()[index];
	tmp_mer.mer = GetSeedMer(tmp_mer.position);
	return tmp_mer;
}


boolean FileSML::Read(vector<bmer>& readVector, gnSeqI size, const gnSeqI offset)
{
	if(!sarfile.is_open()){
		DebugMsg("FileSML::Read: Error sar file not open.\n");
		return false;
	}

	gnSeqI total_len = SMLLength();
	if(offset >= total_len){
		readVector.clear();
		return false;
	}
	gnSeqI readlen = offset + size < total_len ? size : total_len - offset;
	
	readVector.resize( readlen );

	//copy data to the vector
	for(gnSeqI j=0; j < readlen; j++){
		bmer tmp_mer;
		tmp_mer.position = base()[offset+j];
		if( tmp_mer.position > header.length ){
			string errmsg = "Corrupted SML, position ";
			errmsg += tmp_mer.position + " is out of range\n";
			ErrorMsg( errmsg );
			cerr << errmsg;
		}else
			tmp_mer.mer = GetSeedMer(tmp_mer.position);
		readVector[ j ] = tmp_mer;
	}
	return true;
}

void FileSML::BigCreate(const gnSequence& seq, const uint32 split_levels, const uint32 mersize){
//	unsigned long freemem = wxGetFreeMemory();	//get the amount of free memory.
//	unsigned long neededmem = GetNeededMemory(seq.length());
//	if(neededmem >= freemem && neededmem > MEMORY_MINIMUM){ // divide and conquer
	if(split_levels > 0){	// split_levels defines the number of times to divide and conquer
		uint64 midpoint = seq.length() / 2;
		midpoint = (midpoint * header.alphabet_bits) / 32;
		midpoint = (midpoint / header.alphabet_bits) * 32;
		gnSequence seqA = seq.subseq(1, midpoint);
		gnSequence seqB = seq.subseq(1 + midpoint, seq.length() - midpoint);
		seqA.setCircular(false);
		seqB.setCircular(false);
		cout << "Splitting " << seq.length() << " to " << seqA.length() << " and " << seqB.length() << "\n";

		//create the first sar
		string temp_sarfile = CreateTempFileName("bdsa_split");
		FileSML* temp_sar = this->Clone();
		temp_sar->filename = temp_sarfile.c_str();
		temp_sar->BigCreate(seqA, split_levels - 1, mersize);

		//create the second sar
		string temp_sarfile2 = CreateTempFileName("bdsa_split");
		FileSML* temp_sar2 = this->Clone();
		temp_sar2->filename = temp_sarfile2.c_str();
		temp_sar2->BigCreate(seqB, split_levels - 1, mersize);

		//merge them to this file
		cout << "Merging " << seqA.length() << " and " << seqB.length() << "\n";
		Merge(*temp_sar, *temp_sar2);
		//free up RAM
		delete temp_sar;
		delete temp_sar2;
		//erase the temp files.
		boost::filesystem::remove(temp_sarfile);
		boost::filesystem::remove(temp_sarfile2);
	}else{
		Create(seq, mersize);
	}
}

void FileSML::RadixSort(vector<bmer>& s_array){
	vector<bmer> *source_buckets;
	vector<bmer> *tmp_buckets;
	vector<bmer> *buckets;
	uint32 radix_size = 11;
	uint64 radix_mask = 0xFFFFFFFF;
	radix_mask <<= 32;
	radix_mask |= 0xFFFFFFFF;
	radix_mask >>= 64 - radix_size;
	
	uint32 bucket_count = (uint32) pow((double)2, (double)radix_size);
	uint32 cur_shift_bits = 0;
	buckets = new vector<bmer>[bucket_count];
	source_buckets = new vector<bmer>[bucket_count];
	uint64 cur_bucket;
	for(uint32 merI = 0; merI < s_array.size(); merI++){
		cur_bucket = s_array[merI].mer & radix_mask;
		source_buckets[cur_bucket].push_back(s_array[merI]);
	}
	s_array.clear();
	cur_shift_bits += radix_size;
	radix_mask <<= radix_size;
	while(cur_shift_bits < 64){
		for(uint32 bucketI = 0; bucketI < bucket_count; bucketI++){
			for(uint32 merI = 0; merI < source_buckets[bucketI].size(); merI++){
				cur_bucket = source_buckets[bucketI][merI].mer & radix_mask;
				cur_bucket >>= cur_shift_bits;
				buckets[cur_bucket].push_back(source_buckets[bucketI][merI]);
			}
			source_buckets[bucketI].clear();
		}
		cur_shift_bits += radix_size;
		radix_mask <<= radix_size;
		tmp_buckets = source_buckets;
		source_buckets = buckets;
		buckets = tmp_buckets;
	}
	s_array.clear();
	for(uint32 bucketI = 0; bucketI < bucket_count; bucketI++){
		for(uint32 merI = 0; merI < source_buckets[bucketI].size(); merI++){
			s_array.push_back(source_buckets[bucketI][merI]);
		}
		source_buckets[bucketI].clear();
	}
	delete[] source_buckets;
	delete[] buckets;
}

//Merges the supplied sorted mer lists into this one, overwriting the existing sml.
//KNOWN BUG:  The first sorted mer list must have (length * alphabet_bits) / word_bits == 0
//for Merge to work properly.
void FileSML::Merge(SortedMerList& sa, SortedMerList& sa2){
STACK_TRACE_START
	SMLHeader sa_head = sa.GetHeader();
	SMLHeader sa_head2 = sa2.GetHeader();
	
	//basic copying
	header = sa_head;
	//take the smaller mer_size
	if(sa_head.seed_length < sa_head2.seed_length){
		header.seed_length = sa_head.seed_length;
		mer_mask = sa.GetMerMask();
	}else{
		header.seed_length = sa_head2.seed_length;
		mer_mask = sa2.GetMerMask();
	}
	header.unique_mers = NO_UNIQUE_COUNT;
	header.length += sa_head2.length;

	//allocate some memory
	const uint32 SEQ_BUFFER_SIZE = 200000;
	Array<uint32> seq_buf ( SEQ_BUFFER_SIZE + header.seed_length );

	//do some sanity checks on the sars we're merging.
	if(sa_head.alphabet_bits != sa_head2.alphabet_bits ||
	  sa_head.version != sa_head2.version ||
	  memcmp(sa_head.translation_table, sa_head2.translation_table, UINT8_MAX)){
		Throw_gnExMsg(SMLMergeError(), "Incompatible sorted mer lists.");
	}
	
	OpenForWriting( true );

	//write the header
	sarfile.write((char*)&header, sizeof(struct SMLHeader));
	if(!sarfile.good()){
		sarfile.clear();
		sarfile.close();
		sarfile.open(filename.c_str(), ios::binary | ios::in );
		Throw_gnExMsg(IOStreamFailed(), "Error writing sorted mer list header to disk.");
	}

	//copy sequence data into memory.
	uint32 binary_seq_len = (header.length * header.alphabet_bits) / 32;
	if((header.length * header.alphabet_bits) % 32 > 0)
		binary_seq_len++;

	//The +1 is to avoid access violations when copying in the
	//binary sequence before shifting.
	if( sequence != NULL )
		delete[] sequence;
	sequence = new uint32[binary_seq_len+1];
	sa.GetBSequence(sequence, sa_head.length, 0);

	uint32 bseq_len1 = (sa_head.length * sa_head.alphabet_bits) / 32;
	uint32 bseq_remainder = (sa_head.length * sa_head.alphabet_bits) % 32;
	if(bseq_remainder > 0){
		sa2.GetBSequence(&(sequence[bseq_len1]), sa_head2.length, 0);
		//mask off the end of the first sequence
		uint32 end_mask = 0xFFFFFFFF;
		end_mask <<= bseq_remainder;
		sequence[bseq_len1] &= end_mask;

		//shift the second sequence over.
		for(uint32 i=bseq_len1; i < binary_seq_len; i++){
			uint32 tmp = sequence[i+1];
			tmp >>= 32 - bseq_remainder;
			sequence[i] |= tmp;
			sequence[i+1] <<= bseq_remainder;
		}
	}else
		sa2.GetBSequence(&(sequence[bseq_len1]), sa_head2.length, 0);
	
	//write the sequence
	sarfile.write((char*)sequence, binary_seq_len * sizeof(uint32));
	sarray_start_offset = sarfile.tellg();

	//get new mers in the middle
	vector<bmer> middle_mers;
	bmer mid_mer;
	for(uint32 midI = sa_head.length - header.seed_length + 1; midI < sa_head.length; midI++){
		mid_mer.position = midI;
		mid_mer.mer = GetMer(midI);
		middle_mers.push_back(mid_mer);
	}
	sort(middle_mers.begin(), middle_mers.end(), &bmer_lessthan);
	//put a special mer at the end which will never go into the sorted mer list
	//since every possible mer is less than it.
	mid_mer.mer = 0xFFFFFFFF;
	mid_mer.mer <<= 32;
	mid_mer.mer |= 0xFFFFFFFF;
	mid_mer.position = GNSEQI_END;
	middle_mers.push_back(mid_mer);
	//merge and write the sorted mer lists
	vector<bmer> array1, array2;
	uint32 SAR_BUFFER_SIZE = SEQ_BUFFER_SIZE/2;  //actual size is this number * 13 bytes
	uint32 k=0, l=0, midI=0;
	uint32 m = 0, n = 0;
	gnSeqI bufferI=0;
	do{
		//mergesort them
		while(m < array1.size() && n < array2.size()){
			if(array1[m].mer <= array2[n].mer){
				if(array1[m].mer <= middle_mers[midI].mer){
					seq_buf.data[bufferI] = array1[m].position;
					m++;
					bufferI++;
				}else{
					seq_buf.data[bufferI] = middle_mers[midI].position;
					midI++;
					bufferI++;
				}
			}else if(array2[n].mer <= middle_mers[midI].mer){
				seq_buf.data[bufferI] = array2[n].position + sa_head.length;
				n++;
				bufferI++;
			}else{
				seq_buf.data[bufferI] = middle_mers[midI].position;
				midI++;
				bufferI++;
			}
			if(bufferI == SEQ_BUFFER_SIZE){
				sarfile.write((char*)seq_buf.data, bufferI * sizeof(uint32));
				bufferI = 0;
			}
		}
		if(m == array1.size()){
			sa.Read(array1, SAR_BUFFER_SIZE, k);
			k += array1.size();
			m = 0;
		}
		if(n == array2.size()){
			sa2.Read(array2, SAR_BUFFER_SIZE, l);
			l += array2.size();
			n = 0;
		}
	}while(array1.size() != 0 && array2.size() != 0);
	if(bufferI > 0)
		sarfile.write((char*)seq_buf.data, (bufferI)*sizeof(uint32));
	//consolidate the remaining mers to a known vector
	vector<bmer> remaining_mers;
	for(;m < array1.size(); m++)
		remaining_mers.push_back(array1[m]);
	for(;n < array2.size(); n++){
		remaining_mers.push_back(array2[n]);
		remaining_mers[remaining_mers.size()-1].position += sa_head.length;
	}
	for(;midI < middle_mers.size() - 1; midI++)
		remaining_mers.push_back(middle_mers[midI]);
	//merge them with the remaining middle_mers
	sort(remaining_mers.begin(), remaining_mers.end(), &bmer_lessthan);
	uint32 remI = 0;
	for(;remI < remaining_mers.size(); remI++)
		seq_buf.data[remI] = remaining_mers[remI].position;
	if(remI > 0)
		sarfile.write((char*)seq_buf.data, (remI)*sizeof(uint32));

	if(!sarfile.good()){
		sarfile.clear();
		sarfile.close();
		sarfile.open(filename.c_str(), ios::binary | ios::in );
		Throw_gnExMsg(IOStreamFailed(), "Error writing position array.");
	}
	// reopen the sorted mer list file read-only
	sarfile.close();
	sarfile.open(filename.c_str(), ios::binary | ios::in );
	if(!sarfile.is_open()){
		sarfile.clear();
		Throw_gnExMsg(FileNotOpened(), "Error opening sorted mer list file.");
	}
STACK_TRACE_END
}

}
