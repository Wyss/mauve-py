/////////////////////////////////////////////////////////////////////////////
// File:            gnContigSpec.h
// Purpose:         Abstract Contig Spec class
// Description:     Defines an interface for contig specs
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

#ifndef _gnContigSpec_h_
#define _gnContigSpec_h_

#include "libGenome/gnDefs.h"

#include <vector>
#include <string>

#include "libGenome/gnBaseSpec.h"
#include "libGenome/gnBaseFeature.h"
#include "libGenome/gnDebug.h"

namespace genome {

/**
 * gnContigSpec is an interface for classes which store contigs, or reads,
 * of DNA or protein sequence
 */
class GNDLLEXPORT gnContigSpec : public gnBaseSpec
{
public:
	gnContigSpec(){}
	/**
	 * Destructor, frees memory.
	 */
	virtual ~gnContigSpec(){}
	virtual gnContigSpec* Clone() const = 0;
	virtual gnContigSpec* CloneRange( const gnSeqI startI, const gnSeqI len ) const = 0;
	/**
	 * Get the name of the source associated with this spec.
	 * @return The source name or an empty std::string if none exists.
	 */
	virtual std::string GetSourceName() const;
	/**
	 * Get the base pair index where this contig starts inside of the sequence data.
	 * @return The starting base pair.
	 */
	virtual gnSeqI GetStart() const;
	/**
	 * Get the length of this contig inside of the sequence data.
	 * @return This contigs length.
	 */
	virtual gnSeqI GetLength() const;
	/**
	 * Get the length of the source for this spec.
	 * @return The source length.
	 */
	virtual gnSeqI GetSourceLength() const = 0;
	/**
	 * Returns this contig's index in its source sequence.
	 * @return This contig's index in the source sequence.
	 */
	virtual uint32 GetSourceContigIndex() const;
	/**
	 * Sets the name of the source associated with this contig.
	 * @param sourceName The new name.
	 * @return True if successful.
	 */
	virtual void SetSourceName( const std::string& sourceName );
	/**
	 * Sets the starting base pair to read from in the contig's sequence.
	 * This does not affect the actual sequence data but only where to begin
	 * using it in this contig.
	 * @param start The new starting base pair.
	 * @return True if successful.
	 */
	virtual void SetStart( const gnSeqI start );
	/**
	 * Sets the length of reads into this sequence.
	 * This does not affect the actual sequence data but only how much of it
	 * is used in this contig.
	 * @param len The new sequence length.
	 * @return True if successful.
	 */
	virtual void SetLength( const gnSeqI len );
	/**
	 * Sets this contig's index in its source sequence.
	 * @param contigI This contig's index in the source sequence.
	 */
	virtual void SetSourceContigIndex( const uint32 contigI );
	/**
	 * Sets the reverse complement bit for this contig.
	 * This routine will translate the start index to the reverse base pair.
	 * @param value True for reverse complement, false otherwise.
	 * @return True if successful.
	 */
	virtual void SetReverseComplement( const boolean value );

	virtual void CropStart( gnSeqI cropLen );
	virtual void CropEnd( gnSeqI cropLen );

	virtual boolean SeqRead(const gnSeqI start, gnSeqC* buf, gnSeqI& bufLen, const uint32 contigI ) const;
	virtual void Clear();
protected:
	gnSeqI m_start;  //start within the genome.
	gnSeqI m_length;
	uint32 m_SourceContigIndex;

	/**
	 * all derived classes must implement this!
	 * it simply reads the specified bases into buf, disregarding circularity and reverse complement.
	 */
	virtual boolean Read(const gnSeqI start, gnSeqC* buf, gnSeqI& bufLen ) const = 0;

private:

}; // class gnContigSpec


inline
std::string gnContigSpec::GetSourceName() const{
	return m_sourceName;
}
inline
void gnContigSpec::SetSourceName(const std::string& sourceName){
	m_sourceName = sourceName;
}

inline
gnSeqI gnContigSpec::GetStart() const
{
	return m_start;
}
inline
gnSeqI gnContigSpec::GetLength() const
{
	return m_length;
}
// SET
inline
void gnContigSpec::SetStart( const gnSeqI start )
{
	m_start = start;
}
inline
void gnContigSpec::SetLength( const gnSeqI len )
{
	m_length = len;
}
inline
uint32 gnContigSpec::GetSourceContigIndex() const{
	return m_SourceContigIndex;
}
inline
void gnContigSpec::SetSourceContigIndex( const uint32 contigI ){
	m_SourceContigIndex = contigI;
}


}	// end namespace genome

#endif
	// _gnContigSpec_h_
