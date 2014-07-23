#ifndef _alibc_h_
#define _alibc_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

int OpenLibC( aFILE * file, const char *path, int mode );
int CloseLibC( aFILE * file );
int WriteLibC( aFILE * file, aIORec * rec );
int ReadLibC( aFILE * file, aIORec * rec );

/* Line ending test modification... */

#endif /* _alibc_h_ */
