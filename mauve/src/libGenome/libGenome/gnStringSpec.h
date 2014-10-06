/////////////////////////////////////////////////////////////////////////////
// File:            gnStringSpec.h
// Purpose:         implements gnContigSpec for strings
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

#ifndef _gnStringSpec_h_
#define _gnStringSpec_h_

#include "libGenome/gnDefs.h"

#include <cstdlib>
#include <cstring>
#include <string>
#include "libGenome/gnContigSpec.h"
#include "libGenome/gnBaseSource.h"

namespace genome {


/**
 * gnStringSpec stores a sequence and annotation data in memory.
 * For a more complete description see the gnBaseSpec documentation.
 */
class GNDLLEXPORT gnStringSpec : public gnContigSpec
{
public:
	/**
	 * Empty constructor.
	 */
	gnStringSpec();
	/**
	 * Constructor, creates a gnStringSpec using sequence data in the given std::string.
	 * A circular spec will be created if the end index is greater than the start.
	 * @param m_string The std::string to read base pairs from.
	 * @param startI The index to start reading base pairs from the std::string.
	 * @param endI The index to stop reading base pairs from the std::string.
	 * @param revComp True if the sequence is read reverse complement.
	 */
	gnStringSpec( const std::string& m_string, const gnSeqI startI=0, const gnSeqI endI=GNSEQI_END, const boolean revComp = false);
	/**
	 * Copy constructor.
	 * @param s the gnStringSpec to copy.
	 */
	gnStringSpec( const gnStringSpec& s );
	~gnStringSpec();
// Clone
	gnStringSpec* Clone() const;
	virtual void Clear();
// Value Access methods

	virtual gnSeqI GetSourceLength() const;

	
// Source Spec Specific functions
	virtual gnBaseSource *GetSource() const;

	/**
	 * Copies a specified range of bases and returns a pointer to
	 * the resulting gnStringSpec.  You must delete the copy when you
	 * are finished with it.
	 * @param startI The first base pair to copy
	 * @param len The length of the piece to copy
	 * @return A copy of the gnStringSpec containing only the specified bases
	 */
	virtual gnStringSpec* CloneRange( const gnSeqI startI, const gnSeqI len ) const;

protected:
	virtual boolean Read(const gnSeqI start, gnSeqC* buf, gnSeqI& bufLen ) const;
	std::string m_seqString;
	
}; // class gnStringSpec

inline
gnStringSpec* gnStringSpec::Clone() const
{
	return new gnStringSpec( *this );
}
inline
gnSeqI gnStringSpec::GetSourceLength() const{
	return m_seqString.length();
}
inline
gnBaseSource* gnStringSpec::GetSource() const
{
	return NULL;
}

inline
boolean gnStringSpec::Read(const gnSeqI start, gnSeqC* buf, gnSeqI& bufLen) const{
	memcpy(buf, m_seqString.data() + start, bufLen);
	return true;
}


}	// end namespace genome

#endif
	// _gnStringSpec_h_
