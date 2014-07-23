/////////////////////////////////////////////////////////////////////////////
// File:            gnFileSource.h
// Purpose:         Implements gnBaseSource for .File files
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

#ifndef _gnFileSource_h_
#define _gnFileSource_h_

#include "libGenome/gnDefs.h"

#include <string>
#include <fstream>
#include "libGenome/gnBaseSource.h"
#include "libGenome/gnFileContig.h"
#include "libGenome/gnException.h"

namespace genome {

/**
 * gnFileSource is a standard interface to all file based sources of genetic information.
 * All file source classes are derived from this class.
 */
class GNDLLEXPORT gnFileSource : public gnBaseSource
{
public:
	gnFileSource();
	gnFileSource(const gnFileSource& gnfs);
	virtual ~gnFileSource(){}
	virtual gnFileSource* Clone() const = 0;
	// Open, Close	
	virtual void Open( std::string openString );
	virtual void Open();
	virtual void Close();
	virtual std::string GetOpenString() const ;
	  // Filter
	virtual const gnFilter* GetFilter() const ;
	virtual void SetFilter( gnFilter* filter );

	virtual boolean Read( const uint64 pos, char* buf, gnSeqI& bufLen );
	/**
	 * Returns a pointer to the file contig corresponding to contigI or
	 * null if none exists.
	 */
	virtual gnFileContig* GetFileContig( const uint32 contigI ) const = 0;
protected:
	void DetermineNewlineType();

	std::string m_openString;
	std::ifstream m_ifstream;
	const gnFilter* m_pFilter;
	gnNewlineType m_newlineType;
	uint32 m_newlineSize;

private:
	virtual boolean ParseStream( std::istream& fin ) = 0;
};// class gnFileSource

inline
std::string gnFileSource::GetOpenString( ) const
{
	return m_openString;
}
// Filter
inline
const gnFilter* gnFileSource::GetFilter() const
{
	return m_pFilter;
}

inline
void gnFileSource::SetFilter( gnFilter* filter )
{
	if(filter == NULL){
		Throw_gnEx(NullPointer());
	}
	m_pFilter = filter;
}


}	// end namespace genome

#endif
	// _gnFileSource_h_
