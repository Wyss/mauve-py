/////////////////////////////////////////////////////////////////////////////
// File:            libGenome/gnDefs.h
// Purpose:         Defines common constants in libGenome.
// Description:     Defines, consts, typedef etc for libGenome
// Rev:             A
// Author:          Aaron Darling 
// Modified by:     
// Copyright:       (c) Aaron Darling 
// Licenses:        See COPYING file for details 
/////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef _gnDefs_h_
#define _gnDefs_h_

#ifdef __GNUG__
    #pragma interface "libGenome/gnDefs.h"
#endif

#include <limits.h>
#include <float.h>
#include "libGenome/gnSetup.h"

// bool	
typedef unsigned char        boolean;
typedef unsigned char        int1;
// signed
typedef signed char          int8;
typedef signed short int     int16;
typedef signed int           int32;
typedef signed long long     int64;
typedef signed char          sint8;
typedef signed short int     sint16;
typedef signed int           sint32;
typedef signed long int      sint64;
// unsigned
typedef unsigned char        uint8;
typedef unsigned short int   uint16;
typedef unsigned int         uint32;
typedef unsigned long long   uint64;
typedef unsigned int         uint;
// float
typedef float                float32;
typedef double               float64;

// defines
#ifndef TRUE
#define TRUE                 1
#endif
#ifndef FALSE
#define FALSE                0
#endif
#ifndef PI
#define PI                   3.1415926535897932384626433832795028
#endif

#ifndef BOOLEAN_MIN
#define BOOLEAN_MIN	         0
#endif
#ifndef BOOLEAN_MAX
#define BOOLEAN_MAX          1
#endif
#ifndef INT1_MIN
#define INT1_MIN	         0
#endif
#ifndef INT1_MAX
#define INT1_MAX             1
#endif

#define INT8_BYTE_SIZE       1
#define INT16_BYTE_SIZE      2
#define INT32_BYTE_SIZE      4
#define INT64_BYTE_SIZE      8

#define INT8_BIT_SIZE        8
#define INT16_BIT_SIZE       16
#define INT32_BIT_SIZE       32
#define INT64_BIT_SIZE       64

#define UINT8_BYTE_SIZE      1
#define UINT16_BYTE_SIZE     2
#define UINT32_BYTE_SIZE     4
#define UINT64_BYTE_SIZE     8

#define UINT8_BIT_SIZE       8
#define UINT16_BIT_SIZE      16
#define UINT32_BIT_SIZE      32
#define UINT64_BIT_SIZE      64

// limits.h
#ifndef INT8_MIN
#define INT8_MIN             SCHAR_MIN   //0x00
#endif
#ifndef INT8_MAX
#define INT8_MAX             SCHAR_MAX   //0x7f
#endif
#ifndef INT16_MIN
#define INT16_MIN            SHRT_MIN    //0x0000
#endif
#ifndef INT16_MAX
#define INT16_MAX            SHRT_MAX    //0x7fff
#endif
#ifndef INT32_MIN
#define INT32_MIN            INT_MIN     //0x0000000
#endif
#ifndef INT32_MAX
#define INT32_MAX            INT_MAX     //0x7fffffff
#endif
#ifndef INT64_MIN
#define INT64_MIN            LONG_MIN    //0x0000000000000000
#endif
#ifndef INT64_MAX
#define INT64_MAX            LONG_MAX    //0x7fffffffffffffff
#endif

#ifndef UINT8_MIN
#define UINT8_MIN            0           //0x00
#endif
#ifndef UINT8_MAX
#define UINT8_MAX            UCHAR_MAX   //0xff
#endif
#ifndef UINT16_MIN
#define UINT16_MIN           0           //0x0000
#endif
#ifndef UINT16_MAX
#define UINT16_MAX           USHRT_MAX   //0xffff
#endif
#ifndef UINT32_MIN
#define UINT32_MIN           0           //0x00000000
#endif
#ifndef UINT32_MAX
#define UINT32_MAX           UINT_MAX    //0xffffffff
#endif
#ifndef UINT64_MIN
#define UINT64_MIN           0           //0x0000000000000000
#endif
#ifndef UINT64_MAX
#define UINT64_MAX           ULONG_MAX   //0xffffffffffffffff
#endif

// float.h
#define FLOAT32_MIN          FLT_MIN
#define FLOAT32_MAX          FLT_MAX
#define FLOAT32_MIN_EXP      FLT_MIN_EXP
#define FLOAT32_MAX_EXP      FLT_MAX_EXP
#define FLOAT32_MIN_10_EXP   FLT_MIN_10_EXP
#define FLOAT32_MAX_10_EXP   FLT_MAX_10_EXP
#define FLOAT32_DIGIT        FLT_DIG
#define FLOAT32_RADIX        FLT_RADIX
#define FLOAT32_MIN_FRACTION FLT_EPSILON

#define FLOAT64_MIN          DBL_MIN
#define FLOAT64_MAX          DBL_MAX
#define FLOAT64_MIN_EXP      DBL_MIN_EXP
#define FLOAT64_MAX_EXP      DBL_MAX_EXP
#define FLOAT64_MIN_10_EXP   DBL_MIN_10_EXP
#define FLOAT64_MAX_10_EXP   DBL_MAX_10_EXP
#define FLOAT64_DIGIT        DBL_DIG
#define FLOAT64_RADIX        DBL_RADIX
#define FLOAT64_MIN_FRACTION DBL_EPSILON

// Sequence Types
typedef char				gnSeqC;		// Sequence Character
typedef uint64				gnSeqI;		// Sequence Index

#define GNSEQI_ERROR		UINT32_MAX	// return value
#define GNSEQI_END			UINT32_MAX	// argument value
#define GNSEQI_BEGIN		UINT64_MIN	// argument value
#define GNSEQC_NULL			0
#define GNSEQC_MIN			INT8_MIN
#define GNSEQC_MAX			INT8_MAX


#define CONTIG_SECTION_SIZE 3

// some compilers don't have abs() for 64 bit ints
#if (defined(__GNUG__) && ( __GNUC__ <= 2 )) || defined(__INTEL_COMPILER) || (defined _MSC_VER && defined __cplusplus)

int64 abs( int64 a );

#endif

#ifdef __cplusplus
namespace genome {
#endif

enum gnContigSection{
	gnContigHeader = 0,
	gnContigAnnotation = 1,
	gnContigSequence = 2
};

enum gnNewlineType{
	gnNewlineUnix = 0,
	gnNewlineWindows = 1,
	gnNewlineMac = 2
};

static const uint32 ALL_CONTIGS = UINT32_MAX;
static const uint32 BUFFER_SIZE = 100000;


#ifdef __cplusplus

template< typename T >
T absolut( const T& t )
{ 
	return t < 0 ? -t : t; 
};

/**
 * Class used to ensure correct deallocation of arrays.
 */
template<class T>
class Array{
public:
	/** Allocate a new array of length "bufsize".  The array can be
	 * accessed through the data member
	 */
	Array( uint64 bufsize ){
		data = new T[bufsize];
	}
	~Array() { 
		if (data != NULL)
			delete[] data;
	}
	/** The actual array */
	T* data;
private:
	Array( const Array& a );
	Array& operator=( const Array& a );
	Array(){};
};

#endif //__cplusplus

#ifndef __PRETTY_FUNCTION__
#define __PRETTY_FUNCTION__ "Unknown( ) "
#endif

#ifdef __cplusplus
}	// end namespace genome
#endif

#endif
	//_gnDefs_h_
