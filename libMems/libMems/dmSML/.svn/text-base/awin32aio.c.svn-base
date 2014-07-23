#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/dmSML/awin32aio.h"
#include "libMems/dmSML/util.h"
#ifdef USE_WIN32


#define WIN32_LEAN_AND_MEAN
#include <windows.h>


static VOID CALLBACK DummyCompletionRoutine( DWORD err, DWORD nbytes, LPOVERLAPPED lpo ) {
    // we poll for completion, so this is just a dummy to make windows happy.
    printf( "completion routine!\n" );
}


int OpenWIN32( aFILE * file, const char *path, int mode ) {
    HANDLE result;
    DWORD access = mode == A_READ ? GENERIC_READ : GENERIC_WRITE;
    DWORD disposition = mode == A_READ ? OPEN_EXISTING : CREATE_ALWAYS;
    result = CreateFile( 
        path, 
        access, 
        FILE_SHARE_DELETE | FILE_SHARE_READ | FILE_SHARE_WRITE,
        NULL,
        disposition,
        FILE_FLAG_OVERLAPPED,
        NULL );
    if( result == INVALID_HANDLE_VALUE ) {
    	access = GetLastError();
    	printf( "Error opening %s, code %d\n", path, access );
        return( 0 );
    }
    file->w32handle = result;
    return( 1 );
}


int CloseWIN32( aFILE * file ) {
    return( CloseHandle( file->w32handle ) );
}


int WriteWIN32( aFILE * file, aIORec * rec ) {

    static offset_t total_bytes = 0;
    DWORD err;
    if( file->mode != A_WRITE ) {
        return( 0 );
    }

    rec->w32overlapped = malloc( sizeof( *(rec->w32overlapped) ) );
    memset( rec->w32overlapped, 0, sizeof( *(rec->w32overlapped) ) );
	
	if( rec->pos != CURRENT_POS ){
		offset_t tmppos = rec->pos;
		tmppos >>= 32;
		file->filep_high = tmppos;
		// clear high bits.  Is this really necessary?
		tmppos = rec->pos;
		tmppos <<= 32;
		tmppos >>= 32;
		file->filep_low = tmppos;
	}

    rec->w32overlapped->OffsetHigh = file->filep_high;
    rec->w32overlapped->Offset = file->filep_low;

    //printf( "issuing write -- first few bytes of buffer are\n" );
    //for( i = 0; i < 20; i++ ) {
    //    printf( "%c", rec->buf[i] );
    //}
    //printf( "\n" );
    total_bytes += rec->size * rec->count;
    //printf( "total bytes: %d\n", total_bytes );
    if( WriteFileEx( 
        file->w32handle, 
        rec->buf, 
        rec->size*rec->count, 
        rec->w32overlapped,
        DummyCompletionRoutine ) == 0 ) {
        err = GetLastError();
        printf( "error with WriteFileEx: %d\n", err );
        return( 0 );
    }
    return( 1 );
}


int ReadWIN32( aFILE * file, aIORec * rec ) {
    DWORD err;
    if( file->mode != A_READ ) {
        return( 0 );
    }
    rec->w32overlapped = malloc( sizeof( *(rec->w32overlapped) ) );
    memset( rec->w32overlapped, 0, sizeof( *(rec->w32overlapped) ) );

	if( rec->pos != CURRENT_POS ){
		offset_t tmppos = rec->pos;
		tmppos >>= 32;
		file->filep_high = tmppos;
		// clear high bits.  Is this really necessary?
		tmppos = rec->pos;
		tmppos <<= 32;
		tmppos >>= 32;
		file->filep_low = tmppos;
	}

    rec->w32overlapped->OffsetHigh = file->filep_high;
    rec->w32overlapped->Offset = file->filep_low;
    if( ReadFileEx( 
        file->w32handle, 
        rec->buf, 
        rec->size*rec->count, 
        rec->w32overlapped,
        DummyCompletionRoutine ) == 0 ) {
        err = GetLastError();
        switch( err ) {
        case ERROR_HANDLE_EOF:
            printf( "readfileex says EOF -- we'll pretend it worked\n" );
            return( 1 );
        default:
            printf( "error with ReadFileEx -- Last Error: %d\n", GetLastError() );
            printf( "called:  ReadFileEx( %d, %d, %d, %d, %d )\n", 
                file->w32handle, 
                rec->buf, 
                rec->size*rec->count, 
                rec->w32overlapped,
                DummyCompletionRoutine );
            return( 0 );
        }
    }
    return( 1 );
}


int QueryLastCompleteWIN32( aFILE * file ) {
    DWORD result;
    // this operation may not have ever been executed yet (the case
    // where w32overlapped is NULL) so we must detect this.
    if( file->queuetail && file->queuetail->w32overlapped ) {
        // this is a simple poll, because we're waiting for 0 msec.
        result = WaitForSingleObject( file->w32handle, 0 );
        if( result != WAIT_TIMEOUT ) {
            return( 1 );
        } else {
            return( 0 );
        }
        //return( HasOverlappedIoCompleted( file->queuetail->w32overlapped ) );
    } else {
        return( 0 );
    }
}



#endif /* USE_WIN32 */
