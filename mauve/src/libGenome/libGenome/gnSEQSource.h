/////////////////////////////////////////////////////////////////////////////
// File:            gnSEQSource.h
// Purpose:         Implements gnBaseSource for .SEQ files
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

#ifndef _gnSEQSource_h_
#define _gnSEQSource_h_

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
 *	gnSEQSource is a SEQ file reader.
 * This class reads and writes the DNAStar SEQ file format.
 * gnSEQSource is used by gnSourceFactory to read files and should only be used 
 * directly.when writing out files in SEQ file format by calling
 * gnSEQSource::Write( mySpec, "C:\\mySeqFile.seq");
 */

class GNDLLEXPORT gnSEQSource : public gnFileSource
{
public:
	/**
	 * Empty Constructor, does nothing.
	 */
	gnSEQSource();	
	/**
	 * Clone Constructor copies the specified gnSEQSource.
	 * @param s The gnSEQSource to copy.
	 */
	gnSEQSource( const gnSEQSource& s );
	/**
	 * Destructor, frees memory.
	 */
	~gnSEQSource();
	/**
	 * Returns an exact copy of this class.
	 */
	gnSEQSource* Clone() const;
// Contig Access methods	
	uint32 GetContigListLength() const;
	boolean HasContig( const std::string& name ) const;
	uint32 GetContigID( const std::string& name ) const;
	std::string GetContigName( const uint32 i ) const;
	gnSeqI GetContigSeqLength( const uint32 i ) const;

	boolean SeqRead( const gnSeqI start, char* buf, gnSeqI& bufLen, const uint32 contigI=ALL_CONTIGS );

	/**
	 * Writes the specified gnSequence to a .SEQ file named "filename".
	 * @param sequence The gnSequence to write out.
	 * @param filename The name of the file to write.
	 * @return True if successful, false otherwise.
	 */
	static boolean Write(gnSequence& sequence, const std::string& filename);
	/**
	 * Writes the specified source to a .SEQ file named "filename".
	 * @param source The source to write out.
	 * @param filename The name of the file to write.
	 * @return True if successful, false otherwise.
	 */
	static boolean Write(gnBaseSource *source, const std::string& filename);
	/**
	 * Writes the given spec to a .SEQ file named "filename".
	 * @param spec The spec to write out.
	 * @param filename The name of the file to write.
	 * @return True if successful, false otherwise.
	 */
	static boolean Write(gnGenomeSpec *spec, const std::string& filename);
	gnGenomeSpec *GetSpec() const;
	gnFileContig* GetFileContig( const uint32 contigI ) const;
private:
	boolean SeqSeek( const gnSeqI start, const uint32& contigI, uint64& startPos, uint64& readableBytes );
	boolean SeqStartPos( const gnSeqI start, gnFileContig& contig, uint64& startPos, uint64& readableBytes );
	boolean ParseStream( std::istream& fin );

	static std::string& Filler(uint32 length);
	static void FormatString(std::string& data, uint32 offset, uint32 width);
	static void BaseCount(const std::string& bases, gnSeqI& a_count, gnSeqI& c_count, gnSeqI& g_count, gnSeqI& t_count, gnSeqI& other_count);
//	gnSeqI m_seqLength;
	
	gnGenomeSpec *m_spec;
	std::vector< gnFileContig* > m_contigList;	
};// class gnSEQSource
// Clone	
inline
gnSEQSource* gnSEQSource::Clone() const
{
	return new gnSEQSource( *this );
}
// Contig Access methods	
inline
uint32 gnSEQSource::GetContigListLength() const
{
	return m_contigList.size();
}
inline
boolean gnSEQSource::Write(gnSequence& sequence, const std::string& filename){
	return Write(sequence.GetSpec(), filename);
}
inline
boolean gnSEQSource::Write(gnBaseSource *source, const std::string& filename){
	return Write(source->GetSpec(), filename);
}
inline
gnGenomeSpec *gnSEQSource::GetSpec() const{
	return m_spec->Clone();
}


}	// end namespace genome

#endif
	// _gnSEQSource_h_
	
