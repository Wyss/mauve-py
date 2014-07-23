/////////////////////////////////////////////////////////////////////////////
// File:            gnABISource.h
// Purpose:         Implements gnBaseSource for ABI files
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

#ifndef _gnABISource_h_
#define _gnABISource_h_

#include "libGenome/gnDefs.h"

#include <string>
#include "libGenome/gnSequence.h"
#include "libGenome/gnFileSource.h"
#include "libGenome/gnFileContig.h"
#include "libGenome/gnSourceSpec.h"
#include "libGenome/gnGenomeSpec.h"
#include "libGenome/gnFilter.h"

namespace genome {

/**
 * gnABISource is not implemented.
 */
class GNDLLEXPORT gnABISource : public gnFileSource
{
public:
	gnABISource();	
	/**
	 * Clone Constructor copies the specified gnSEQSource.
	 * @param s The gnABISource to copy.
	 */
	gnABISource( const gnABISource& s );
	~gnABISource();
	gnABISource* Clone() const;

	uint32 GetContigListLength() const;
	boolean HasContig( const std::string& name ) const;
	uint32 GetContigID( const std::string& name ) const;
	std::string GetContigName( const uint32 i ) const;
	gnSeqI GetContigSeqLength( const uint32 i ) const;
	gnFileContig* GetContig( const uint32 i ) const;

	boolean SeqRead( const gnSeqI start, char* buf, gnSeqI& bufLen, const uint32 contigI=ALL_CONTIGS );

	/**
	 * Writes the specified gnSequence to an ABI file named "filename".
	 * @param sequence The gnSequence to write out.
	 * @param filename The name of the file to write.
	 * @return True if successful, false otherwise.
	 */
	static boolean Write(gnSequence& sequence, const std::string& filename);
	gnGenomeSpec *GetSpec() const;

	gnFileContig* GetFileContig( const uint32 contigI ) const;
private:
	boolean SeqSeek( const gnSeqI start, const uint32& contigI, uint64& startPos, uint64& readableBytes );
	boolean SeqStartPos( const gnSeqI start, gnFileContig& contig, uint64& startPos, uint64& readableBytes );
	boolean ParseStream( std::istream& fin );

	gnGenomeSpec *m_spec;
	std::vector< gnFileContig* > m_contigList;	

};// class gnABISource

inline
gnABISource* gnABISource::Clone() const
{
	return new gnABISource( *this );
}
inline
uint32 gnABISource::GetContigListLength() const
{
	return m_contigList.size();
}
inline
gnGenomeSpec *gnABISource::GetSpec() const
{
	return m_spec->Clone();
}


}	// end namespace genome

#endif
	// _gnABISource_h_
