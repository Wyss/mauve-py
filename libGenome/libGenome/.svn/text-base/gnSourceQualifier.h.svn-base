/////////////////////////////////////////////////////////////////////////////
// File:            gnSourceQualifier.h
// Purpose:         Source Qualifier class
// Description:     Provides an interface for gnBaseQualifier in files.
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

#ifndef _gnSourceQualifier_h_
#define _gnSourceQualifier_h_

#include "libGenome/gnDefs.h"

#include <string>
#include "libGenome/gnBaseQualifier.h"
#include "libGenome/gnBaseSource.h"

#include <utility>

namespace genome {


/**
 * gnSourceQualifier is used to store a sequence qualifier which resides 
 * in another source.  They are created by gnGBKSource and other sources 
 * and placed into a gnFeature.
 */
class GNDLLEXPORT gnSourceQualifier : public gnBaseQualifier
{
public:
	/**
	 * Empty constructor.
	 */
	gnSourceQualifier();
	/**
	 * Constructor, records the location and name of the qualifier in the source.
	 * @param source The source which contains the qualifier.
	 * @param name The name of the qualifier.
	 * @param begin The offset into the source where the qualifier starts.
	 * @param length The length of the qualifier.
	 */
	gnSourceQualifier( gnBaseSource* source, std::string& name, uint32 begin, uint32 length );
	/**
	 * Copy constructor.
	 * @param s The gnSourceQualifier to copy.
	 */
	gnSourceQualifier( const gnSourceQualifier& s );
	/**
	 * Destructor, frees memory.
	 */
	~gnSourceQualifier();

	gnSourceQualifier* Clone() const;

	std::string GetName() const;
	std::string GetValue() const;
	
	uint32 GetNameLength() const;
	uint32 GetValueLength() const;

	/**
	 * Get the qualifier's start position within the source.
	 * @return The qualifier start position.
	 */
	uint32 GetValueStart() const;
	/**
	 * Get the qualifier's start position and length within the source.
	 * @return The qualifier start position and length as a pair.
	 */
	std::pair<uint32, uint32> GetValueLoc() const;

	/**
	 * Set the qualifier's start position within the source.
	 * @param start The qualifier start position.
	 */
	void SetValueStart(const uint32 start);
	/**
	 * Get the qualifier's length within the source.
	 * @param length The qualifier length.
	 */
	void SetValueLength(const uint32 length);
	/**
	 * Set the qualifier's start position and length within the source.
	 * @param startLen The qualifier start position and length as a pair.
	 */
	void SetValueLoc(const std::pair<uint32, uint32> startLen);

private:
	std::string m_name;
	uint32 m_start, m_length;
	gnBaseSource *m_source;
}; //class gnSourceQualifier

inline 
gnSourceQualifier* gnSourceQualifier::Clone() const{
	return new gnSourceQualifier(*this);
}
inline
std::string gnSourceQualifier::GetName() const{
	return m_name;
}
inline
uint32 gnSourceQualifier::GetNameLength() const{
	return m_name.length();
}
inline
uint32 gnSourceQualifier::GetValueStart() const{
	return m_start;
}
inline
uint32 gnSourceQualifier::GetValueLength() const{
	return m_length;
}
inline
std::pair<uint32, uint32> gnSourceQualifier::GetValueLoc() const{
	std::pair<uint32, uint32> p;
	p.first = m_start;
	p.second = m_length;
	return p;
}
inline
void gnSourceQualifier::SetValueStart(const uint32 start){
	m_start = start;
}
inline
void gnSourceQualifier::SetValueLength(const uint32 length){
	m_length = length;
}
inline
void gnSourceQualifier::SetValueLoc(const std::pair<uint32, uint32> startLen){
	m_start = startLen.first;
	m_length = startLen.second;
}


}	// end namespace genome

#endif
	// _gnSourceQualifier_h_
