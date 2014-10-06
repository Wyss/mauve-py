/////////////////////////////////////////////////////////////////////////////
// File:            gnBaseQualifier.h
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

#ifndef _gnBaseQualifier_h_
#define _gnBaseQualifier_h_

#include "libGenome/gnDefs.h"

#include <string>
#include "libGenome/gnClone.h"


namespace genome {

/**
 * gnBaseQualifier is a general interface to sequence qualifiers.
 * gnBaseFeature uses qualifiers to store annotated sequences.
 * Use gnStringQualifier for a general purpose qualifier class.
 */
class GNDLLEXPORT gnBaseQualifier : public gnClone
{
public:
	gnBaseQualifier(){}
	virtual ~gnBaseQualifier(){}
	virtual gnBaseQualifier* Clone() const = 0;
	/**
	 * Get the name of qualifier stored in this class.
	 * @return The qualifier name as a std::string.
	 */
	virtual std::string GetName() const = 0;
	/**
	 * Get the qualifier stored in this class.
	 * @return The qualifier as a std::string.
	 */
	virtual std::string GetValue() const = 0;
	
	/**
	 * Get the length of the qualifier name stored in this class.
	 * @return The length of the qualifier name.
	 */
	virtual uint32 GetNameLength() const = 0;
	/**
	 * Get the length of the qualifier stored in this class.
	 * @return The length of the qualifier.
	 */
	virtual uint32 GetValueLength() const = 0;
private:
}; //class gnBaseQualifier


}	// end namespace genome

#endif
	// _gnBaseQualifier_h_
	
