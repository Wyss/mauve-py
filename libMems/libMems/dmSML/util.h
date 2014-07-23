#ifndef _util_h_
#define _util_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdarg.h>

// these just let you get a temporary string -- don't hang onto it for very
// long though -- it will get overwritten at some point.  This is useful for
// passing as parms and such.
const char * Fmt( const char * fmt, ... );
const char * VFmt( const char * fmt, va_list args );


/// shifts a 64-bit value (in two 32 bit parts) either right or left.
/// amt negative -> left, positive -> right
void Shift64( int amt, int * hi, int * lo );


void AddTo64( unsigned int amt, unsigned int *hi, unsigned int *lo );

/** cross-platform file deletion */
int removeFile( const char* filename, int verbose );


#endif /* _util_h_ */
