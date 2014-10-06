/////////////////////////////////////////////////////////////////////////////
// File:            gnFASSource.h
// Purpose:         Implements gnBaseSource for .FAS files
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

#ifndef _gnFASSource_h_
#define _gnFASSource_h_

#include "libGenome/gnDefs.h"

#include <string>
#include <fstream>
#include <vector>
#include "libGenome/gnFileSource.h"
#include "libGenome/gnSequence.h"


namespace genome {

#define FAS_LINE_WIDTH 80

/**
 *	gnFASSource reads and writes FastA files.
 * gnFASSource is used by gnSourceFactory to read files. 
 * Files can be written in the FastA file format by calling
 * gnFASSource::Write( mySpec, "C:\\myFasFile.fas");
 */

class GNDLLEXPORT gnFASSource : public gnFileSource
{
public:
	/**
	 * Empty Constructor, does nothing.
	 */
	gnFASSource();
	/**
	 * Clone Constructor copies the specified gnFASSource.
	 * @param s The gnFASSource to copy.
	 */
	gnFASSource( const gnFASSource& s );
	/**
	 * Destructor, frees memory.
	 */
	~gnFASSource();
	/**
	 * Returns an exact copy of this class.
	 */
	gnFASSource* Clone() const;

	uint32 GetContigListLength() const;
	boolean HasContig( const std::string& name ) const;
	uint32 GetContigID( const std::string& name ) const;
	std::string GetContigName( const uint32 i ) const;
	gnSeqI GetContigSeqLength( const uint32 i ) const;
	gnFileContig* GetContig( const uint32 i ) const;

	boolean SeqRead( const gnSeqI start, char* buf, gnSeqI& bufLen, const uint32 contigI=ALL_CONTIGS ) ;
	
	/**
	 * Write the given gnSequence to a FastA file.
	 * @param sequence The gnSequence to write out.
	 * @param filename The name of the file to write.
	 * @param write_coords If true each entry's name will be followed by the coordinates of the entry
	 *                     in the context of the entrire file.
	 * @param enforce_unique_names If true each entry's name will be recorded as they are written.  Each
	 *                              successive duplicate name that is found will have an underscore and a
	 *   						   number appended to it, indicating the number of entries by the same
	 *							   name which have already been written.
	 *                             Turning this off will yield a slight performance improvement when writing
	 *                             files with a large number of entries.  (More than 1000)
	 * @throws A FileNotOpened() exception may be thrown.
	 */
	static void Write(gnSequence& sequence, const std::string& filename, boolean write_coords = true, boolean enforce_unique_names = true);

	/**
	 * Write the given gnSequence to an ostream.
	 * @param sequence The gnSequence to write out.
	 * @param m_ostream The output stream to write to.
	 * @param write_coords If true each entry's name will be followed by the coordinates of the entry
	 *                     in the context of the entrire file.
	 * @param enforce_unique_names If true each entry's name will be recorded as they are written.  Each
	 *                              successive duplicate name that is found will have an underscore and a
	 *   						   number appended to it, indicating the number of entries by the same
	 *							   name which have already been written.
	 *                             Turning this off will yield a slight performance improvement when writing
	 *                             files with a large number of entries.  (More than 1000)
	 */
	static void Write(gnSequence& sequence, std::ostream& m_ostream, boolean write_coords = true, boolean enforce_unique_names = true);

	/**
	 * Deprecated - do not use.
	 * Write the given source to a FastA file.
	 * @param source The spec to write out.
	 * @param filename The name of the file to write.
	 */
	static boolean Write(gnBaseSource *source, const std::string& filename);
	
	gnGenomeSpec *GetSpec() const;
	
	gnFileContig* GetFileContig( const uint32 contigI ) const;
private:
	boolean SeqReadImpl( const gnSeqI start, char* buf, gnSeqI& bufLen, const uint32 contigI=ALL_CONTIGS ) ;
	boolean SeqSeek( const gnSeqI start, const uint32 contigI, uint64& startPos, uint64& readableBytes );
	boolean SeqStartPos( const gnSeqI start, gnFileContig& contig, uint64& startPos, uint64& readableBytes );
	boolean ParseStream( std::istream& fin );
	
	std::vector< gnFileContig* > m_contigList;
};// class gnFASSource

inline
gnFASSource* gnFASSource::Clone() const
{
	return new gnFASSource( *this );
}

inline
uint32 gnFASSource::GetContigListLength() const
{
	return m_contigList.size();
}


}	// end namespace genome

#endif
	// _gnFASSource_h_
