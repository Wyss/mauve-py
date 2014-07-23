/////////////////////////////////////////////////////////////////////////////
// File:            gnClone.h
// Purpose:         Abstract Clone class
// Description:     Define the a Clonable interface for other objects
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

#ifndef _gnClone_h_
#define _gnClone_h_

#include "libGenome/gnDefs.h"


namespace genome {

class GNDLLEXPORT gnClone
{
	public:
		virtual ~gnClone(){}
		virtual gnClone* Clone() const = 0;
		gnClone(){}
	private:
};// class gnClone


}	// end namespace genome

#endif
	// _gnClone_h_
