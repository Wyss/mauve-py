/////////////////////////////////////////////////////////////////////////////
// File:            gnBaseHeader.h
// Purpose:         abstract Header class
// Description:     Provides an interface for Headers in memory and on disk.
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

#ifndef _gnBaseHeader_h_
#define _gnBaseHeader_h_

#include "libGenome/gnDefs.h"

#include <string>
#include "libGenome/gnClone.h"
#include "libGenome/gnLocation.h"

namespace genome {


/**
 * This class provides a general interface to sequence related headers.
 * Headers commonly precede sequence data in several file formats.  In
 * FastA files, the header is on the > line before a contig.
 * In GenBank files, each contig has headers describing authors and other
 * information.  In GenBank files, the header name corresponds to the
 * name of the header field.  A few genBank header names are: DEFINITION
 * ACCESSION, VERSION, KEYWORDS, SEGMENT, SOURCE, REFERENCE, and COMMENT.
 * Eventually, individual classes may be implemented for each header type.
 */
class GNDLLEXPORT gnBaseHeader : public gnClone
{
public:
	gnBaseHeader(){}
	virtual ~gnBaseHeader(){}
	virtual gnBaseHeader* Clone() const = 0;
	/**
	 * Get the header.
	 * @return The header as a std::string.
	 */
	virtual std::string GetHeader() const = 0;
	/**
	 * Get the header's name, if any.
	 * @return The header name as a std::string.
	 */
	virtual std::string GetHeaderName() const = 0;	
	/**
	 * Get the header's length in bytes.
	 * @return The length of the header.
	 */
	virtual uint32 GetLength() const = 0;
private:
}; //class gnBaseHeader


}	// end namespace genome

#endif
	// _gnBaseHeader_h_
	
