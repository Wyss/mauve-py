#ifndef _aPOSIXaio_h_
#define _aPOSIXaio_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/dmSML/asyncio.h"

int OpenPAIO( aFILE * file, const char *path, int mode );
int ClosePAIO( aFILE * file );

int WritePAIO( aFILE * file, aIORec * rec );
int ReadPAIO( aFILE * file, aIORec * rec );

int QueryLastCompletePAIO( aFILE * file );

#endif /* _aPOSIXaio_h_ */
