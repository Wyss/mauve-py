/////////////////////////////////////////////////////////////////////////////
// File:            gnStringQualifier.h
// Purpose:         abstract Qualifier class
// Description:     Provides an interface for Qualifiers in memory and on disk.
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

#ifndef _gnStringQualifier_h_
#define _gnStringQualifier_h_

#include "libGenome/gnDefs.h"

#include <string>
#include "libGenome/gnClone.h"
#include "libGenome/gnBaseQualifier.h"

namespace genome {


/**
 * gnStringQualifier stores a sequence qualifier in memory.
 * Use gnStringQualifier for a general purpose qualifier class.
 */
class GNDLLEXPORT gnStringQualifier : public gnBaseQualifier
{
public:
	/**
	 * Empty constructor.
	 */
	gnStringQualifier();
	/**
	 * Create a gnStringQualifier.
	 * @param name The qualifier name.
	 * @param value The qualifier.
	 */
	gnStringQualifier(const std::string& name, const std::string& value);
	/**
	 * Copy constructor.
	 * @param s The gnStringQualifier to copy.
	 */
	gnStringQualifier(const gnStringQualifier& s);
	/**
	 * Destructor, frees memory.
	 */
	~gnStringQualifier(){}

	gnStringQualifier* Clone() const;

	std::string GetName() const;
	std::string GetValue() const;
	/**
	 * Set the name of qualifier stored in this class.
	 * @param name The qualifier name as a std::string.
	 */
	void SetName(const std::string& name);
	/**
	 * Set the qualifier stored in this class.
	 * @param value The header as a std::string.
	 */
	void SetValue(const std::string& value);
	
	uint32 GetNameLength() const;
	uint32 GetValueLength() const;
private:
	std::string m_name;
	std::string m_value;
}; //class gnStringQualifier

inline
gnStringQualifier::gnStringQualifier(){
	m_name = "";
	m_value = "";
}
inline
gnStringQualifier::gnStringQualifier(const std::string& name, const std::string& value){
	m_name = name;
	m_value = value;
}
inline
gnStringQualifier::gnStringQualifier(const gnStringQualifier& s){
	m_name = std::string(s.m_name);
	m_value = std::string(s.m_value);
}
inline
gnStringQualifier* gnStringQualifier::Clone() const{
	return new gnStringQualifier(*this);
}
inline
std::string gnStringQualifier::GetName() const{
	return m_name;
}
inline
std::string gnStringQualifier::GetValue() const{
	return m_value;
}
inline
void gnStringQualifier::SetName(const std::string& name){
	m_name = name;
}
inline
void gnStringQualifier::SetValue(const std::string& value){
	m_value = value;
}
inline
uint32 gnStringQualifier::GetNameLength() const{
	return m_name.length();
}
inline
uint32 gnStringQualifier::GetValueLength() const{
	return m_value.length();
}



}	// end namespace genome

#endif
	// _gnStringQualifier_h_
