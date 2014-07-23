/////////////////////////////////////////////////////////////////////////////
// File:            gnGBKSource.h
// Purpose:         Implements gnBaseSource for GenBank files
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

#ifndef _gnGBKSource_h_
#define _gnGBKSource_h_

#include "libGenome/gnDefs.h"

#include <string>
#include <fstream>
#include <vector>
#include "libGenome/gnFileSource.h"
#include "libGenome/gnFileContig.h"
#include "libGenome/gnSourceSpec.h"
#include "libGenome/gnSequence.h"


namespace genome {

const uint32 SEQ_COLUMN_WIDTH = 80;
const uint32 SEQ_HEADER_NAME_LENGTH = 11;
const uint32 SEQ_SUBTAG_COLUMN = 5;
const uint32 SEQ_LOCUS_CIRCULAR_COLUMN = 43;
const uint32 SEQ_LOCUS_NAME_COLUMN = 13;
const uint32 SEQ_LOCUS_NAME_LENGTH = 10;
const uint32 SEQ_LOCUS_SIZE_LENGTH = 10;
const uint32 SEQ_LOCUS_DNATYPE_OFFSET = 33;
const uint32 SEQ_LOCUS_DNATYPE_LENGTH = 7;
const uint32 SEQ_LOCUS_DIVCODE_OFFSET = 52;
const uint32 SEQ_LOCUS_DIVCODE_LENGTH = 3;
const uint32 SEQ_LOCUS_DATE_OFFSET = 62;
const uint32 SEQ_LOCUS_DATE_LENGTH = 11;
const uint32 SEQ_FEATURE_LOC_OFFSET = 21;
const uint32 SEQ_BASES_INDEX_END = 9;

/**
 *	gnGBKSource is a GenBank sequence file reader.
 * This class reads and writes the GenBank file format.
 * gnGBKSource is used by gnSourceFactory to read files and should only be used
 * directly when writing out files in GBK file format by calling
 * gnGBKSource::Write( mySpec, "C:\\mySeqFile.gbk");
 */

class GNDLLEXPORT gnGBKSource : public gnFileSource
{
public:
	/**
	 * Empty Constructor, does nothing.
	 */
	gnGBKSource();	
	/**
	 * Clone Constructor copies the specified gnGBKSource.
	 * @param s The gnGBKSource to copy.
	 */
	gnGBKSource( const gnGBKSource& s );
	/**
	 * Destructor, frees memory.
	 */
	~gnGBKSource();
	/**
	 * Returns an exact copy of this class.
	 */
	gnGBKSource* Clone() const;
// Contig Access methods	
	uint32 GetContigListLength() const;
	boolean HasContig( const std::string& name ) const;
	uint32 GetContigID( const std::string& name ) const;
	std::string GetContigName( const uint32 i ) const;
	gnSeqI GetContigSeqLength( const uint32 i ) const;

	boolean SeqRead( const gnSeqI start, char* buf, gnSeqI& bufLen, const uint32 contigI=ALL_CONTIGS );

	/**
	 * Writes the specified gnSequence to a .GBK file named "filename".
	 * @param seq The gnSequence to write out.
	 * @param filename The name of the file to write.
	 * @return True if successful, false otherwise.
	 */
	static boolean Write(gnSequence& seq, const std::string& filename);
	/**
	 * Writes the specified source to a .GBK file named "filename".
	 * @param source The source to write out.
	 * @param filename The name of the file to write.
	 * @return True if successful, false otherwise.
	 */
	static boolean Write(gnBaseSource *source, const std::string& filename);
	gnGenomeSpec *GetSpec() const;
	gnFileContig* GetFileContig( const uint32 contigI ) const;
private:
	boolean SeqReadImpl( const gnSeqI start, char* buf, gnSeqI& bufLen, const uint32 contigI=ALL_CONTIGS );
	boolean SeqSeek( const gnSeqI start, const uint32& contigI, uint64& startPos, uint64& readableBytes );
	boolean SeqStartPos( const gnSeqI start, gnFileContig& contig, uint64& startPos, uint64& readableBytes );
	boolean ParseStream( std::istream& fin );

	static std::string& Filler(uint32 length);
	static void FormatString(std::string& data, uint32 offset, uint32 width);
//	gnSeqI m_seqLength;
	
	gnGenomeSpec *m_spec;
	std::vector< gnFileContig* > m_contigList;	
};// class gnGBKSource

template< class SubSpec >
void WriteHeader(gnMultiSpec< SubSpec >* spec, const std::string& hdr, std::ofstream& m_ofstream);

// Clone	
inline
gnGBKSource* gnGBKSource::Clone() const
{
	return new gnGBKSource( *this );
}
// Contig Access methods	
inline
uint32 gnGBKSource::GetContigListLength() const
{
	return m_contigList.size();
}
inline
boolean gnGBKSource::Write(gnBaseSource *source, const std::string& filename){
	gnSequence gns(*source->GetSpec());
	return Write(gns, filename);
}
inline
gnGenomeSpec *gnGBKSource::GetSpec() const{
	return m_spec->Clone();
}


}	// end namespace genome

#endif
	// _gnGBKSource_h_
