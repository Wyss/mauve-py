/////////////////////////////////////////////////////////////////////////////
// File:            gnStringHeader.h
// Purpose:         abstract Header class
// Description:     Provides an interface for Headers in memory.
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

#ifndef _gnStringHeader_h_
#define _gnStringHeader_h_

#include "libGenome/gnDefs.h"

#include <string>
#include "libGenome/gnClone.h"
#include "libGenome/gnLocation.h"
#include "libGenome/gnBaseHeader.h"

namespace genome {


/**
 * gnStringHeader stores sequence related header information in memory.
 * Use gnStringHeader for a general purpose headers.
 * @see gnBaseHeader
 */
class GNDLLEXPORT gnStringHeader : public gnBaseHeader
{
public:
	/**
	 * Empty constructor.
	 */
	gnStringHeader();
	/**
	 * Create a gnStringHeader.
	 * @param name The header name.
	 * @param header The header.
	 */
	gnStringHeader(const std::string& name, const std::string& header);
	/**
	 * Copy constructor.
	 * @param s The gnStringHeader to copy.
	 */
	gnStringHeader(const gnStringHeader& s);
	/**
	 * Destructor, frees memory.
	 */
	~gnStringHeader() {};

	gnStringHeader* Clone() const;

	std::string GetHeader() const;
	std::string GetHeaderName() const;

	/**
	 * Set the header stored in this class.
	 * @param header The header as a std::string.
	 */
	void SetHeader(const std::string& header);
	/**
	 * Set the header's name stored in this class.
	 * @param name The header name as a std::string.
	 */
	void SetHeaderName(const std::string& name);
	
	uint32 GetLength() const;
private:
	std::string m_name;
	std::string m_header;
}; //class gnStringHeader

inline
gnStringHeader::gnStringHeader(){
	m_header = std::string();
}
inline
gnStringHeader::gnStringHeader(const std::string& name, const std::string& header){
	m_name = name;
	m_header = header;
}
inline
gnStringHeader::gnStringHeader(const gnStringHeader& s){
	m_header = std::string(s.m_header);
}
inline
gnStringHeader* gnStringHeader::Clone() const{
	return new gnStringHeader(*this);
}
inline
std::string gnStringHeader::GetHeader() const{
	return m_header;
}
inline
std::string gnStringHeader::GetHeaderName() const{
	return m_name;
}
inline
void gnStringHeader::SetHeader(const std::string& header){
	m_header = header;
}
inline
void gnStringHeader::SetHeaderName(const std::string& name){
	m_name = name;
}
inline
uint32 gnStringHeader::GetLength() const{
	return m_header.length();
}


}	// end namespace genome

#endif
	// _gnStringHeader_h_
