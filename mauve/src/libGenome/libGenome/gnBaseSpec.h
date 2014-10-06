/////////////////////////////////////////////////////////////////////////////
// File:            gnBaseSpec.h
// Purpose:         abstract Spec class
// Description:     Defines a basic interface for all specs
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

#ifndef _gnBaseSpec_h_
#define _gnBaseSpec_h_

#include "libGenome/gnDefs.h"

#include <vector>
#include <string>

#include "libGenome/gnClone.h"
#include "libGenome/gnBaseFeature.h"
#include "libGenome/gnBaseHeader.h"


namespace genome {

/**
 * gnBaseSpec is the class which stores genetic information and is best accessed using gnSequence.
 */
class GNDLLEXPORT gnBaseSpec : public gnClone
{
public:
	gnBaseSpec(){}
	/**
	 * Destructor, frees memory.
	 */
	virtual ~gnBaseSpec(){}
	virtual gnBaseSpec* Clone() const = 0;
	virtual gnBaseSpec* CloneRange( const gnSeqI startI, const gnSeqI len ) const = 0;
	/**
	 * Get the name of the contig associated with this spec.
	 * @return The contig name or an empty std::string if none exists.
	 */
	virtual std::string GetName() const;
	/**
	 * Sets the name for this contig
	 * @param name The new name.
	 * @return True if successful.  Honestly, I don't know how this could be unsuccessful...
	 */
	virtual void SetName( const std::string& name );
	/**
	 * Get the length of all the sequence data covered by this spec.
	 * @return This spec's length in base pairs.
	 */
	virtual gnSeqI GetLength() const = 0;
	/**
	 * Returns true if this spec is read reverse complement.
	 * @return True if this spec is read reverse complement.
	 */
	virtual boolean IsReverseComplement() const;
	/**
	 * Returns true if this spec's sequence is circular.
	 * @return True if this spec's sequence is circular.
	 */
	virtual boolean IsCircular() const;
	/**
	 * Sets the reverse complement bit for this spec.
	 * @param value True for reverse complement, false otherwise.
	 */
	virtual void SetReverseComplement( const boolean value ) = 0;
	/**
	 * Sets whether this spec should be read circular.
	 * If circular is set, reads beyond the end of this spec will pick up
	 * at the beginning and read up to the start index.
	 * @param value True for circular, false otherwise.
	 */
	virtual void SetCircular( const boolean value );
	/**
	 * Crop the first cropLen bases from the sequence. CropStart will
	 * delete features and headers associated with the cropped bases.
	 * @param cropLen The number of base pairs to delete from the beginning.
	 */
	virtual void CropStart( gnSeqI cropLen ) = 0;
	/**
	 * Crop the last cropLen bases from the sequence. CropEnd will
	 * delete features and headers associated with the cropped bases.
	 * @param cropLen The number of base pairs to delete from the end.
	 */
	virtual void CropEnd( gnSeqI cropLen ) = 0;

	/**
	 * Reads sequence data from this spec.
	 * SeqRead will attempt to read "bufLen" base pairs starting at "start", an offset into the sequence.
	 * Reading inside a specific contig can be accomplished by supplying the "contigI" parameter with
	 * a valid contig index.
	 * SeqRead stores the sequence data in "buf" and returns the actual number of bases read in "bufLen".
	 * SeqRead will return false if a serious error occurs.
	 * @param start The base pair to start reading at.
	 * @param buf The character array to store base pairs into.
	 * @param bufLen The number of base pairs to read.  This will be modified to reflect the actual number of bases read.
	 * @param contigI The index of the subspec to read or ALL_CONTIGS by default.
	 * @return True if the operation was successful.
	 */
	virtual boolean SeqRead(const gnSeqI start, gnSeqC* buf, gnSeqI& bufLen, const uint32 contigI ) const = 0;
	/**
	 * Clears all data from this spec.
	 */
	virtual void Clear();
protected:
	boolean m_reverseComplement;
	boolean m_circular;
	
	std::string m_name;
	std::string m_sourceName;
	
}; // class gnBaseSpec

inline
std::string gnBaseSpec::GetName() const{
	return m_name;
}
inline
void gnBaseSpec::SetName( const std::string& name ){
	m_name = name;
}
inline
boolean gnBaseSpec::IsReverseComplement() const{
	return m_reverseComplement;
}
inline
boolean gnBaseSpec::IsCircular() const{
	return m_circular;
}
inline
void gnBaseSpec::SetCircular( const boolean value ){
	m_circular = value;
}
inline
void gnBaseSpec::Clear(){
	m_sourceName = "";
	m_name = "";
	m_reverseComplement = false;
	m_circular = false;
}


}	// end namespace genome

#endif
	// _gnBaseSpec_h_
