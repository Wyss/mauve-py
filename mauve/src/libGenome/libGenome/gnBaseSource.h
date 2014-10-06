/////////////////////////////////////////////////////////////////////////////
// File:            gnBaseSource.h
// Purpose:         Abstract source class
// Description:     Basic interface to all source objects
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

#ifndef _gnBaseSource_h_
#define _gnBaseSource_h_

#include "libGenome/gnDefs.h"

#include <string>

#include "libGenome/gnClone.h"

namespace genome {

class gnGenomeSpec;
class gnFilter;

/**
 * gnBaseSource defines a standard interface for derived classes to
 * provide access to file, database, and network sources of genetic data.
 * First the data source is opened and parsed, using the open() function.
 * The source class then creates and a gnGenomeSpec containing references
 * to its sequence data.  The gnGenomeSpec can then be used by gnSequence
 * to give programmer friendly access to the sequence data.
 */
class GNDLLEXPORT gnBaseSource : public gnClone
{
public:
	gnBaseSource(){}
	/**
	 * Destructor, frees memory used by this source class.
	 */
	virtual ~gnBaseSource(){}
	virtual gnBaseSource* Clone() const = 0;

	/**
	 * Opens the source given in "openString" for reading.
	 * @param openString The name of the source (file, network URL, or 
	 * database) to open.
	 * @throws Will throw a FileNotOpened exception if the file was not found
	 * or was not accessible.  Will propagate a FileUnreadable exception if the
	 * file format was invalid.
	 */
	virtual void Open( std::string openString ) = 0;
	/**
	 * Opens this source for reading.
	 * @throws Will throw a FileNotOpened exception if the file was not found
	 * or was not accessible.
	 */
	virtual void Open() = 0;
	/**
	 * Closes the file or connection this source is reading from.
	 * @throws IOStreamError if an error occurs closing the file.
	 */
	virtual void Close() = 0;
	/**
	 * Get the location of the source that is being used.
	 * @return The location std::string describing this source, usually a file
	 * name.
	 */
	virtual std::string GetOpenString() const = 0;

	/**
	 * Get the number of sequence contigs in this source.
	 * @return The number of contigs in this source.
	 */
	virtual uint32 GetContigListLength() const = 0;
	/**
	 * Looks for a contig by name.
	 * Returns true if it finds the contig, otherwise false.
	 * @param name The name of the contig to look for.
	 * @return True if the named contig exists, false otherwise.
	 */
	virtual boolean HasContig( const std::string& name ) const = 0;
	/**
	 * Get a contig index by name.
	 * If the source does not contain a contig by the specified name
	 * GetContigID returns UINT32_MAX.
	 * @param name The name of the contig to look for.
	 * @return The index of the named contig or UINT32_MAX.
	 */
	virtual uint32 GetContigID( const std::string& name ) const = 0;
	/**
	 * Get the name of the specified contig.
	 * Returns an empty std::string if the specified contig is out of range.
	 * @param i The index of the contig or ALL_CONTIGS.
	 * @return The name of the contig or an empty std::string.
	 */
	virtual std::string GetContigName( const uint32 i ) const = 0;
	/**
	 * Get the total number of base pairs in the specified contig.
	 * @param i The index of the contig or ALL_CONTIGS.
	 * @return The length in base pairs of the specified contig.
	 */
	virtual gnSeqI GetContigSeqLength( const uint32 i ) const = 0;

	/**
	 * Get the filter currently being used to filter unwanted characters out of read sequences.
	 * @return A pointer to the gnFilter currently in use.
	 */
	virtual const gnFilter* GetFilter() const = 0;
	/**
	 * Set the filter that will be used to filter unwanted characters out of the sequence data.
	 * @param filter The filter to remove unwanted characters from the sequence.
	 * @throws NullPointer is thrown if the specified filter pointer is null.
	 */
	virtual void SetFilter( gnFilter* filter ) = 0;

	/**
	 * Gets raw input from this source.
	 * Read will attempt to read "bufLen" bytes starting at "pos" directly from the source.
	 * It stores the data in "buf", and returns the actual number of bytes read in bufLen.
	 * Read will return false if a serious error occurs.
	 * @param pos The position in the file to start reading.
	 * @param buf The character array to store data into.
	 * @param bufLen The number of bytes to read.
	 * @return True if the operation was successful.
	 */
	virtual boolean Read( const uint64 pos, char* buf, gnSeqI& bufLen) = 0;
	/**
	 * Gets sequence data from this source.
	 * SeqRead will attempt to read "bufLen" base pairs starting at "start", an offset into the sequence.
	 * Reading inside a specific contig can be accomplished by supplying the "contigI" parameter with
	 * a valid contig index.
	 * SeqRead stores the sequence data in "buf" and returns the actual number of bases read in "bufLen".
	 * SeqRead will return false if a serious error occurs.
	 * @param start The base pair to start reading at.
	 * @param buf The character array to store base pairs into.
	 * @param bufLen The number of base pairs to read.
	 * @param contigI The index of the contig to read or ALL_CONTIGS by default.
	 * @return True if the operation was successful.
	 */
	virtual boolean SeqRead( const gnSeqI start, char* buf, gnSeqI& bufLen, const uint32 contigI=ALL_CONTIGS ) = 0;
	/**
	 * Get the annotated sequence data as a gnGenomeSpec.
	 * GetSpec returns a gnGenomeSpec which contains the sequence, header,
	 * and feature data contained by this source.
	 * @return The annotated sequence data.
	 */
	virtual gnGenomeSpec *GetSpec() const = 0;
private:
};// class gnBaseSource


}	// end namespace genome

#endif
	// _gnBaseSource_h_
