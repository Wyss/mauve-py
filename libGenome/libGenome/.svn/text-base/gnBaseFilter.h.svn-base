/////////////////////////////////////////////////////////////////////////////
// File:            gnBaseFilter.h
// Purpose:         Generic filter interface
// Description:     Filters sequences, translates, reverse complement, converts
//                   additions, etc.
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

#ifndef _gnBaseFilter_h_
#define _gnBaseFilter_h_

#include "libGenome/gnDefs.h"

#include <string>
#include "libGenome/gnClone.h"
#include "libGenome/gnDefs.h"


namespace genome {

class GNDLLEXPORT gnBaseFilter : public gnClone
{
public:		
	virtual gnBaseFilter* Clone() const = 0;
	
	/**
	 * Gets the name of this filter
	 * @return the filter name 
	 */
	virtual std::string GetName() const;
	/**
	 * Sets the name of this filter
	 * @param name the new filter name
	 */
	virtual void SetName( std::string name );
	
	/**
	 * Filter the given character
	 * @param ch The character to filter
	 * @return The filtered character
	 */
	virtual gnSeqC Filter( const gnSeqC ch ) const = 0;

	/**
	 * Filter the given character array
	 * @param seq A pointer to the character array
	 * @param len the length of the character array to filter
	 * @return The filtered character
	 */
	virtual void Filter( gnSeqC** seq, gnSeqI& len ) const = 0;

	/**
	 * Filters the given std::string
	 * @param seq The std::string to filter
	 */
	virtual void Filter( std::string &seq ) const = 0;

protected:
	std::string m_name;

};//class gnBaseFilter

inline
std::string gnBaseFilter::GetName() const
{
	return m_name;
}
inline
void gnBaseFilter::SetName( std::string name )
{
	m_name = name;
}


}	// end namespace genome

#endif
	// _gnBaseFilter_h_
