#ifndef _awin32_h_
#define _awin32_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/dmSML/asyncio.h"

int OpenWIN32( aFILE * file, const char *path, int mode );
int CloseWIN32( aFILE * file );

int WriteWIN32( aFILE * file, aIORec * rec );
int ReadWIN32( aFILE * file, aIORec * rec );

int QueryLastCompleteWIN32( aFILE * file );

#endif /* _awin32_h_ */
