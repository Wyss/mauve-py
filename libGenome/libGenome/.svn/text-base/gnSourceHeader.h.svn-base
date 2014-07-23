/////////////////////////////////////////////////////////////////////////////
// File:            gnSourceHeader.h
// Purpose:         Source Header class
// Description:     Provides an interface for Headers on disk.
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

#ifndef _gnSourceHeader_h_
#define _gnSourceHeader_h_

#include "libGenome/gnDefs.h"

#include <string>
#include "libGenome/gnBaseHeader.h"
#include "libGenome/gnBaseSource.h"

#include <utility>

namespace genome {


/**
 * gnSourceHeader is used to store a sequence header which resides in another source.
 * gnSourceHeaders are created by gnSEQSource and other sources and placed into a spec.
 */
class GNDLLEXPORT gnSourceHeader : public gnBaseHeader
{
public:
	/**
	 * Empty constructor.
	 */
	gnSourceHeader();
	/**
	 * Constructor, records the header's name and byte offset in the source.
	 * @param source The source which contains the header.
	 * @param name The name of the header.
	 * @param begin The offset into the source where the header starts.
	 * @param length The length of the header.
	 */
	gnSourceHeader( gnBaseSource* source, const std::string& name, const uint32 begin, const uint32 length );
	/**
	 * Copy constructor.
	 * @param s The gnSourceHeader to copy.
	 */
	gnSourceHeader( const gnSourceHeader& s );
	/**
	 * Destructor, frees memory.
	 */
	~gnSourceHeader();

	gnSourceHeader* Clone() const;

	std::string GetHeader() const;
	std::string GetHeaderName() const;
	/// Not implemented, use gnStringHeader. 
//	void SetHeader(const std::string& header){};
	/// Not implemented, use gnStringHeader. 
//	void SetHeaderName(const std::string& name){};
	
	uint32 GetLength() const;
	/**
	 * Get the header's start position within the source.
	 * @return The byte offset of the header's start position.
	 */
	uint32 GetHeaderStart() const;
	/**
	 * Get the header's start position and length within the source.
	 * @return The header start position and length as a pair.
	 */
	std::pair<uint32, uint32> GetHeaderLoc() const;

	/**
	 * Set the header's start position within the source.
	 * @param start The header start position.
	 */
	void SetHeaderStart(const uint32 start);
	/**
	 * Get the header's length within the source.
	 * @param length The header length.
	 */
	void SetHeaderLength(const uint32 length);
	/**
	 * Set the header's start position and length within the source.
	 * @param startLen The header start position and length as a pair.
	 */
	void SetHeaderLoc(const std::pair<uint32, uint32> startLen);

private:
	std::string m_name;
	uint32 m_start, m_length;
	gnBaseSource *m_source;
}; //class gnSourceHeader

inline 
gnSourceHeader* gnSourceHeader::Clone() const{
	return new gnSourceHeader(*this);
}
inline
std::string gnSourceHeader::GetHeaderName() const{
	return m_name;
}
inline
uint32 gnSourceHeader::GetHeaderStart() const{
	return m_start;
}
inline
uint32 gnSourceHeader::GetLength() const{
	return m_length;
}
inline
std::pair<uint32, uint32> gnSourceHeader::GetHeaderLoc() const{
	std::pair<uint32, uint32> p;
	p.first = m_start;
	p.second = m_length;
	return p;
}
inline
void gnSourceHeader::SetHeaderStart(const uint32 start){
	m_start = start;
}
inline
void gnSourceHeader::SetHeaderLength(const uint32 length){
	m_length = length;
}
inline
void gnSourceHeader::SetHeaderLoc(const std::pair<uint32, uint32> startLen){
	m_start = startLen.first;
	m_length = startLen.second;
}



}	// end namespace genome

#endif
	// _gnSourceHeader_h_
