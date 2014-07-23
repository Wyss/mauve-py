/////////////////////////////////////////////////////////////////////////////
// File:            libGenome/gnSetup.h
// Purpose:         libGenome setup
// Description:     Defines os/compiler specific constants, etc.
//                  Included in libGenome/gnDefs.h.
// Rev:             A
// Author:          Aaron Darling 
// Modified by:     
// Copyright:       (c) Aaron Darling 
// Licenses:        
/////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef _gnSetup_h_
#define _gnSetup_h_

#include <stdlib.h>

#ifdef GNMAKINGDLL   // build the libgenome dll
#define GNDLLEXPORT __declspec(dllexport)
#define GNDLLEXPORT_DATA(type) __declspec(dllexport) type
#elif defined(GNUSINGDLL)  
// the project uses a libGenome as a dll
#define GNDLLEXPORT __declspec(dllimport)
#define GNDLLEXPORT_DATA(type) __declspec(dllimport) type
#else	// static linking
#define GNDLLEXPORT
#define GNDLLEXPORT_DATA
#endif

// attempt to auto-link the genome library on windows
#if defined(WIN64)&&defined(NDEBUG)&&!defined(FASTDEBUG)&&defined(_OPENMP)
#pragma comment(lib, "genome64omp.lib")
#endif
#if defined(WIN64)&&defined(FASTDEBUG)&&defined(_OPENMP)
#pragma comment(lib, "genome64fdomp.lib")
#endif
#if defined(WIN32)&&!defined(WIN64)&&defined(NDEBUG)&&!defined(FASTDEBUG)&&defined(_OPENMP)
#pragma comment(lib, "genomeomp.lib")
#endif
#if defined(WIN32)&&!defined(WIN64)&&defined(FASTDEBUG)&&defined(_OPENMP)
#pragma comment(lib, "genomefdomp.lib")
#endif
#if defined(WIN64)&&defined(NDEBUG)&&!defined(FASTDEBUG)&&!defined(_OPENMP)
#pragma comment(lib, "genome64.lib")
#endif
#if defined(WIN64)&&defined(FASTDEBUG)&&!defined(_OPENMP)
#pragma comment(lib, "genome64fd.lib")
#endif
#if defined(WIN32)&&!defined(WIN64)&&defined(NDEBUG)&&!defined(FASTDEBUG)&&!defined(_OPENMP)
#pragma comment(lib, "genome.lib")
#endif
#if defined(WIN32)&&!defined(WIN64)&&defined(FASTDEBUG)&&!defined(_OPENMP)
#pragma comment(lib, "genomefd.lib")
#endif


#endif
	//_gnSetup_h_
