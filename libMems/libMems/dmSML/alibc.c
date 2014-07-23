#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/dmSML/asyncio.h"
#include "libMems/dmSML/alibc.h"

#if defined USE_LIBC

int OpenLibC( aFILE * file, const char *path, int mode ) {
    FILE * result = fopen( path, mode == A_READ ? "rb" : "wb" );
    file->libchandle = result;
    if( result == NULL ) {
        return( 0 );
    }
    return( 1 );
}


int CloseLibC( aFILE * file ) {
    fclose( file->libchandle );
    return( 1 );
}


int WriteLibC( aFILE * file, aIORec * rec ) {
    fwrite( rec->buf, rec->size, rec->count, file->libchandle );
    return( 1 );
}

int ReadLibC( aFILE * file, aIORec * rec ) {
    fread( rec->buf, rec->size, rec->count, file->libchandle );
    return( 1 );
}


int OperationCompleteLibC( aFILE * file ) {
    // libc operations are atomic
    return( 1 );
}

int FileBusyLibC( aFILE * file ) {
    // libc operations are atomic
    return( 1 );
}

#endif /* USE_LIBC */
