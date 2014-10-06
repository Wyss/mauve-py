/////////////////////////////////////////////////////////////////////////////
// File:            gnSourceSpec.h
// Purpose:         implements gnBaseSpec for source specs
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

#ifndef _gnSourceSpec_h_
#define _gnSourceSpec_h_

#include "libGenome/gnDefs.h"

#include <string>
#include "libGenome/gnContigSpec.h"
#include "libGenome/gnBaseSource.h"

namespace genome {


/**
 * gnSourceSpec stores sequence and annotation data from another source.
 * @see gnContigSpec
 * @see gnBaseSpec
 */
class GNDLLEXPORT gnSourceSpec : public gnContigSpec
{
public:
	/**
	 * Empty constructor.
	 */
	gnSourceSpec();
	/**
	 * Constructor, creates a gnSourceSpec using sequence data in the given source.
	 * A circular spec will be created if the end index is greater than the start.
	 * @param m_pSource The source to read base pairs from.
	 * @param m_ContigIndex The index of the contig in the source, ALL_CONTIGS by default.
	 * @param startI The index to start reading base pairs from the source.
	 * @param endI The index to stop reading base pairs from the source.
	 * @param revComp True if the sequence is read reverse complement.
	 */
	gnSourceSpec( gnBaseSource* m_pSource, const uint32 m_ContigIndex=ALL_CONTIGS, const gnSeqI startI=0, const gnSeqI endI=GNSEQI_END, const boolean revComp = false);
	/**
	 * Copy constructor.
	 * @param s the gnSourceSpec to copy.
	 */
	gnSourceSpec( const gnSourceSpec& s );
	~gnSourceSpec();

	gnSourceSpec* Clone() const;
// Source Spec Specific functions
	virtual void Clear();
	virtual gnSeqI GetSourceLength() const;	
	virtual gnBaseSource *GetSource() const;

	/**
	 * Copies a specified range of bases and returns a pointer to
	 * the resulting gnSourceSpec.  You must delete the copy when you
	 * are finished with it.
	 * @param startI The first base pair to copy
	 * @param len The length of the piece to copy
	 * @return A copy of the gnSourceSpec containing only the specified bases
	 */
	virtual gnSourceSpec* CloneRange( const gnSeqI startI, const gnSeqI len ) const;

protected:
	virtual boolean Read(const gnSeqI start, gnSeqC* buf, gnSeqI& bufLen ) const;
	gnBaseSource* m_pSource;
	
}; // class gnSourceSpec

inline
gnSourceSpec* gnSourceSpec::Clone() const
{
	return new gnSourceSpec( *this );
}
inline
gnSeqI gnSourceSpec::GetSourceLength() const{
	return m_pSource->GetContigSeqLength(m_SourceContigIndex);
}
inline
gnBaseSource* gnSourceSpec::GetSource() const
{
	return m_pSource;
}

inline
boolean gnSourceSpec::Read(const gnSeqI start, gnSeqC* buf, gnSeqI& bufLen) const{
	return m_pSource->SeqRead(start, buf, bufLen, m_SourceContigIndex);
}


}	// end namespace genome

#endif
	// _gnSourceSpec_h_
