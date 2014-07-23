/////////////////////////////////////////////////////////////////////////////
// File:            gnRAWSource.h
// Purpose:         Implements gnBaseSource for raw data files
// Description:     
// Changes:        
// Version:         libGenome 0.5.1 
// Author:          Aaron Darling 
// Modified by:     
// Copyright:       (c) Aaron Darling 
// Licenses:        See COPYING file for details 
/////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef _gnRAWSource_h_
#define _gnRAWSource_h_

#include "libGenome/gnDefs.h"

#include <string>
#include <fstream>
#include <vector>
#include "libGenome/gnFileSource.h"
#include "libGenome/gnFileContig.h"
#include "libGenome/gnSourceSpec.h"
#include "libGenome/gnSequence.h"

namespace genome {


/**
 * gnRAWSource reads raw sequence data from a text file.
 * This class reads and writes raw sequence to and from files.
 * A raw sequence does not contain any newlines, fragment delimiters, or other
 * type of annotation.
 * gnRAWSource is used by gnSourceFactory to read files and should only be used 
 * directly.when writing out raw files by calling gnRAWSource::Write( mySpec, "C:\\myFile.txt");
 */

class GNDLLEXPORT gnRAWSource : public gnFileSource
{
public:
	/**
	 * Empty Constructor, does nothing.
	 */
	gnRAWSource();	
	/**
	 * Clone Constructor copies the specified gnRAWSource.
	 * @param s The gnRAWSource to copy.
	 */
	gnRAWSource( const gnRAWSource& s );
	/**
	 * Destructor, frees memory.
	 */
	~gnRAWSource();
	/**
	 * Returns an exact copy of this class.
	 */
	gnRAWSource* Clone() const;
// Contig Access methods	
	uint32 GetContigListLength() const;
	boolean HasContig( const std::string& name ) const;
	uint32 GetContigID( const std::string& name ) const;
	std::string GetContigName( const uint32 i ) const;
	gnSeqI GetContigSeqLength( const uint32 i ) const;

	boolean SeqRead( const gnSeqI start, char* buf, gnSeqI& bufLen, const uint32 contigI=ALL_CONTIGS );

	/**
	 * Writes the specified gnSequence to a raw file named "filename".
	 * @param sequence The gnSequence to write out.
	 * @param filename The name of the file to write.
	 * @return True if successful, false otherwise.
	 */
	static boolean Write(gnSequence& sequence, const std::string& filename);
	/**
	 * Writes the specified source to a raw file named "filename".
	 * @param source The source to write out.
	 * @param filename The name of the file to write.
	 * @return True if successful, false otherwise.
	 */
	static boolean Write(gnBaseSource *source, const std::string& filename);
	gnGenomeSpec *GetSpec() const;
	gnFileContig* GetFileContig( const uint32 contigI ) const;

	static boolean CheckRawData( boolean set = false, boolean check = false );
private:
	boolean SeqSeek( const gnSeqI start, const uint32& contigI, uint64& startPos, uint64& readableBytes );
	boolean SeqStartPos( const gnSeqI start, gnFileContig& contig, uint64& startPos, uint64& readableBytes );
	boolean ParseStream( std::istream& fin );
	
	gnFileContig* m_contig;	
	gnGenomeSpec* m_spec;
};// class gnRAWSource
// Clone	
inline
gnRAWSource* gnRAWSource::Clone() const
{
	return new gnRAWSource( *this );
}
// Contig Access methods	
inline
uint32 gnRAWSource::GetContigListLength() const
{
	return m_contig == NULL? 0 : 1;
}
inline
boolean gnRAWSource::Write(gnBaseSource *source, const std::string& filename){
	gnSequence gns(*source->GetSpec());
	return Write(gns, filename);
}

inline
boolean gnRAWSource::CheckRawData( boolean set, boolean check ){
	static boolean check_raw_data = false;
	if( set ){
		check_raw_data = check;
	}
	return check_raw_data;
}


}	// end namespace genome

#endif
	// _gnRAWSource_h_
