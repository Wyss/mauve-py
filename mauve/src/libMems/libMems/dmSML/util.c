#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include "libMems/dmSML/util.h"

#define FMT_BUFFER_SIZE     (32)
#define FMT_MAX_STRING      (1024)

static char FmtBuffer[FMT_BUFFER_SIZE][FMT_MAX_STRING];
static int FmtIdx;

const char * Fmt( const char * fmt, ... ) {
    const char * ret;
    va_list args;
    va_start( args, fmt );
    ret = VFmt( fmt, args );
    va_end( args );
    return( ret );
}


const char * VFmt( const char * fmt, va_list args ) {
    if( ++FmtIdx >= FMT_BUFFER_SIZE ) {
        FmtIdx = 0;
    }
    // silly windows....
#ifdef WIN32
    _vsnprintf( FmtBuffer[FmtIdx], sizeof( FmtBuffer[FmtIdx] ), fmt, args );
#else
    vsnprintf( FmtBuffer[FmtIdx], sizeof( FmtBuffer[FmtIdx] ), fmt, args );
#endif
    FmtBuffer[FmtIdx][FMT_MAX_STRING-1] = '\0';
    return( FmtBuffer[FmtIdx] );
}


/// shifts a 64-bit value (in two 32 bit parts) either right or left.
/// amt negative -> left, positive -> right
void Shift64( int amt, int * hi, int * lo ) {
    if( amt == 0 ) {
        return;
    }
    if( amt > 0 ) {
        *lo >>= amt;
        *lo |= *hi << ((sizeof( *hi ) * 8) - amt);
        *hi >>= amt;
    } else {
        amt = -amt;
        *hi <<= amt;
        *hi |= *lo >> ((sizeof( *lo ) * 8) - amt);
        *lo <<= amt;
    }
}




void AddTo64( unsigned int amt, unsigned int *hi, unsigned int *lo ) {

    int i;
    // holds each byte value.
    int in[8], out[8], tmp[8];
    int carry;

    for( i = 0; i < 8; i++ ) {
        in[i] = out[i] = tmp[i] = 0;
    }

    in[0] = amt & 0xFF;
    in[1] = (amt >> 8) & 0xFF;
    in[2] = (amt >> 16) & 0xFF;
    in[3] = (amt >> 24) & 0xFF;

    tmp[0] = *lo & 0xFF;
    tmp[1] = (*lo >> 8) & 0xFF;
    tmp[2] = (*lo >> 16) & 0xFF;
    tmp[3] = (*lo >> 24) & 0xFF;
    tmp[4] = *hi & 0xFF;
    tmp[5] = (*hi >> 8) & 0xFF;
    tmp[6] = (*hi >> 16) & 0xFF;
    tmp[7] = (*hi >> 24) & 0xFF;
    

    /*
    out[0] = (tmp[0] + in[0]);
    carry = out[0] >> 8;
    out[0] &= 0xFF;

    out[1] = (tmp[1] + in[1] + carry);
    carry = out[1] >> 8;
    out[1] &= 0xFF;

    out[2] = (tmp[2] + in[2] + carry);
    carry = out[2] >> 8;
    out[2] &= 0xFF;

    out[3] = (tmp[3] + in[3] + carry);
    carry = out[3] >> 8;
    out[3] &= 0xFF;
    */

    carry = 0;
    for( i = 0; i < 8; i++ ) {
        out[i] = in[i] + tmp[i] + carry;
        carry = out[i] >> 8;
        out[i] &= 0xFF;
    }

    // convert back to 2-int form.
    *lo = out[0] + (out[1] << 8) + (out[2] << 16) + (out[3] << 24);
    *hi = out[4] + (out[5] << 8) + (out[6] << 16) + (out[7] << 24);

}

/** Utility function to delete a file */
int removeFile( const char* filename, int verbose )
{
#ifdef WIN32
		return remove( filename );
//        return !DeleteFile( filename );
#else
        char* rm_cmd;
        if( verbose )
                rm_cmd = Fmt( "/bin/rm -fv %s", filename );
        else
                rm_cmd = Fmt( "/bin/rm -f %s", filename );
        return system( rm_cmd );
#endif
}

