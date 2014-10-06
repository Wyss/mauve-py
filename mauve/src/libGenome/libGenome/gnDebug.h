/////////////////////////////////////////////////////////////////////////////
// File:            libGenome/gnDebug.h
// Purpose:         Debug header used for libGenome
// Description:     Debug defines all debuging tools used for libGenome
// Rev:             A
// Author:          Aaron Darling 
// Modified by:     
// Copyright:       (c) Aaron Darling 
// Licenses:        See COPYING file for details 
/////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef _gnDebug_h_
#define _gnDebug_h_

#include "libGenome/gnDefs.h"
#include <string>
#include <iostream>

#if(defined(WIN32))
#include "intrin.h"
#endif

namespace genome {

/** this little function traps into the MSVC debugger */
inline
void breakHere()
{
#if(defined(WIN32))
__debugbreak();
#endif
}

#if defined(COMMAND_LINE) || defined(_CONSOLE)

const boolean USE_COMMAND_LINE = true;
const boolean USE_GUI = false;

#elif defined(GN_GUI)

const boolean USE_COMMAND_LINE = false;
const boolean USE_GUI = true;

#else

const boolean USE_COMMAND_LINE = false;
const boolean USE_GUI = false;

#endif
GNDLLEXPORT void DebugMsg(std::string a);
inline
void DebugMsg(std::string a){
	if(USE_COMMAND_LINE){
		std::cout << a;
	}else if(USE_GUI){
	}
}

GNDLLEXPORT void ErrorMsg(std::string a);
inline
void ErrorMsg(std::string a){
	if(USE_COMMAND_LINE){
		std::cout << a;
	}else if(USE_GUI){
	}
}


}	// end namespace genome

#endif
	//_gnDebug_h_
