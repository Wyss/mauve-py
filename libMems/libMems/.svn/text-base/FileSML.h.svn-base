/*******************************************************************************
 * $Id: FileSML.h,v 1.11 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _FileSML_h_
#define _FileSML_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#pragma warning(push)
#pragma warning(disable : 4996)
#pragma warning(pop)

#include "libGenome/gnSequence.h"
#include "libMems/SortedMerList.h"
#include <boost/iostreams/device/mapped_file.hpp>
#include <fstream>
#include <vector>
#include <string>

namespace mems {

//sequence database size will be
//base_count / 4 + base_count * 12 bytes

#define DEFAULT_MEMORY_MINIMUM 20971520  //~20 Megabytes

class FileSML : public SortedMerList
{
public:
	FileSML() : SortedMerList() {
//		file_mutex = new wxMutex();
	};
	FileSML& operator=(const FileSML& sa);
	virtual FileSML* Clone() const = 0;
	
	virtual void Clear();
	
	/**
	 * Loads an existing sorted mer list from a file on disk.
	 * @param fname The name of the file to load
	 * @throws FileNotOpened thrown if the file could not be opened
	 * @throws FileUnreadable thrown if the file was corrupt or not a sorted mer list
	 */
	virtual void LoadFile(const std::string& fname);
	/**
	 * Creates large sorted mer lists which do not fit entirely in memory.
	 * BigCreate uses an external mergesort to create large sorted mer lists.
	 * It will divide the data a number of times specified by the split_levels
	 * parameter.  Each split is written to temp files on disk and merged.
	 * @param seq The sequence to create an SML for.
	 * @param split_levels The number of times to divide the sequence in half.
	 * @param mersize The size of the mers to sort on.
	 * @see FileSML::Create
	 */
	virtual void BigCreate(const genome::gnSequence& seq, const uint32 split_levels, const uint32 mersize = DNA_MER_SIZE);
	virtual void Create(const genome::gnSequence& seq, const uint64 seed );
	virtual boolean Read(std::vector<bmer>& readVector, gnSeqI size, gnSeqI offset = 0);
	virtual void Merge(SortedMerList& sa, SortedMerList& sa2);

	virtual bmer operator[]( gnSeqI index );

	virtual gnSeqI UniqueMerCount();
	virtual void SetDescription(const std::string& d);
	virtual void SetID(const sarID_t d);
	
	virtual uint32 FormatVersion();
	static uint64 MemoryMinimum();
	virtual void RadixSort(std::vector<bmer>& s_array);

	void dmCreate(const genome::gnSequence& seq, const uint64 seed);
	static void registerTempPath( const std::string& tmp_path );

	static const char* getTempPath( int pathI );

	static int getTempPathCount();
	
	const std::vector< int64 >& getUsedCoordinates() const { return seq_coords; };

protected:
	/**
	 * Reopens the sarfile fstream in read/write mode
	 * @throws FileNotOpened thrown if the file could not be opened for writing
	 */
	virtual void OpenForWriting( boolean truncate = false );
	/**
	 * Writes the SML header to disk
	 * @throws FileNotOpened thrown if the file could not be opened for writing
	 * @throws IOStreamFailed thrown if an error occurred writing the data
	 */
	virtual boolean WriteHeader();
	/**
	 * Calculates and returns the amount of memory needed to create a sorted
	 * mer list for a sequence of the specified length.
	 * @param len The length of the sequence
	 * @return The amount of memory needed in bytes.
	 */
	virtual uint64 GetNeededMemory(gnSeqI len) = 0;

	std::string filename;
	std::fstream sarfile;
	uint64 sarray_start_offset;

	boost::iostreams::mapped_file_source sardata;
	smlSeqI_t* base(){ return (smlSeqI_t*)(sardata.data()+sarray_start_offset); }
	
	static char** tmp_paths;	/**< paths to scratch disk space that can be used for an external sort */
	std::vector< int64 > seq_coords;	/**< If Ns are masked, contains coordinates of regions without Ns */
};

// versions 2 and 5 were previous
// jump to 100 to avoid confusion with DNAFileSML
inline
uint32 FileSML::FormatVersion(){
	static uint32 f_version = 100;
	return f_version;
}

inline
uint64 FileSML::MemoryMinimum(){
	static uint32 m_minimum = DEFAULT_MEMORY_MINIMUM;
	return m_minimum;
}

void maskNNNNN( const genome::gnSequence& in_seq, genome::gnSequence& out_seq, std::vector< int64 >& seq_coords, int mask_n_length );

}

#endif   //_FileSML_h_
