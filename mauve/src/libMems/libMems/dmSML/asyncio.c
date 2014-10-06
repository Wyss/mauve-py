#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/dmSML/asyncio.h"

#include "libMems/dmSML/alibc.h"
#include "libMems/dmSML/awin32aio.h"
#include "libMems/dmSML/aPOSIXaio.h"

#include "libMems/dmSML/util.h"

#include "libMems/dmSML/buffer.h"
#include <string.h>

#if defined(USE_POSIX_AIO)||defined(USE_LINUX_AIO)
#include <unistd.h>
#include <sys/stat.h>
#endif

static int OperationNumber = 0;

int QueueEmpty( aFILE * file );
void RemoveOperation( aFILE * file );
void FreeQueue( aFILE * file );
int ExecuteWrite( aFILE * file, aIORec * rec );
int ExecuteRead( aFILE * file, aIORec * rec );
int EnqueueOperation( char * buffer, offset_t size, offset_t count, aFILE * file, offset_t pos );
void ExecuteOperation( aFILE * file );
int QueryOpComplete( aFILE * file );
int aAct( void * buffer, offset_t size, offset_t count, aFILE * file, offset_t pos );

int QueueEmpty( aFILE * file ) {
    return( (file->queuehead == file->queuetail) && (file->queuehead == NULL) );
}

void RemoveOperation( aFILE * file ) {
    aIORec * tofree;
    if( !QueueEmpty( file ) ) {
        tofree = file->queuetail;
        if( file->queuetail == file->queuehead ) {
            file->queuehead = file->queuetail = NULL;
        } else {
            file->queuetail->next->last = NULL;
            file->queuetail = file->queuetail->next;
        }
        // FIXME: ack hack from my poor design for win32
#if defined USE_WIN32
        free( tofree->w32overlapped );
#elif defined(USE_POSIX_AIO)||defined(USE_LINUX_AIO)
		free( tofree->aio_cb );
#endif
        free( tofree );
    }        
}


void FreeQueue( aFILE * file ) {
    while( !QueueEmpty( file ) ) {
        RemoveOperation( file );
    }
}



// opens a file
aFILE * aOpen( const char * path, int mode ) {
    int err = 0;
    aFILE *ret = malloc( sizeof( *ret ) );
    
    memset( ret, 0, sizeof( *ret ) );
    ret->mode = mode;
    ret->busy = 0;
#if defined USE_LINUX_AIO
    err = !OpenLinux( ret, path, mode );
#elif defined USE_POSIX_AIO
    err = !OpenPAIO( ret, path, mode );
#elif defined USE_LIBC
    err = !OpenLibC( ret, path, mode );
#elif defined USE_WIN32
    err = !OpenWIN32( ret, path, mode );
#endif
    if( err ) {
        free( ret );
        ret = NULL;
    }
    return( ret );
}


// helper to close a file
int aClose( aFILE * file ) {
    int err = 0;
    // block until the file is no longer busy.
    aWaitNotBusy( file );
#if defined USE_LINUX_AIO
    err = CloseLinux( file );
#elif defined USE_POSIX_AIO
    err = ClosePAIO( file );
#elif defined USE_LIBC
    err = CloseLibC( file );
#elif defined USE_WIN32
    err = CloseWIN32( file );
#endif
    FreeQueue( file );
    free( file );
    return( err );
}



int ExecuteWrite( aFILE * file, aIORec * rec ) {
    int err = 0;
#if defined USE_LINUX_AIO
    err = !WriteLinux( file, rec );
#elif defined USE_POSIX_AIO
    err = !WritePAIO( file, rec );
#elif defined USE_LIBC
    err = !WriteLibC( file, rec );
#elif defined USE_WIN32
    err = !WriteWIN32( file, rec );
#endif
    if( err ) {
        //printf( "error in ExecuteWrite\n" );
    } else {
        file->busy = 1;
    }
    return( err );
}


int ExecuteRead( aFILE * file, aIORec * rec ) {
    int err = 0;
#if defined USE_LINUX_AIO
    err = !ReadLinux( file, rec );
#elif defined USE_POSIX_AIO
    err = !ReadPAIO( file, rec );
#elif defined USE_LIBC
    err = !ReadLibC( file, rec );
#elif defined USE_WIN32
    err = !ReadWIN32( file, rec );
#endif
    if( err ) {
        //printf( "error in ExecuteRead\n" );
    } else {
        file->busy = 1;
    }
    return( err );
}



int EnqueueOperation( char * buffer, offset_t size, offset_t count, aFILE * file, offset_t pos ) {
    if( QueueEmpty( file ) ) {
        
		file->queuehead = file->queuetail = malloc( sizeof( *file->queuehead ) );
		memset( file->queuehead, 0, sizeof( *(file->queuehead) ) );
		file->queuehead->last = NULL;
    } else {
		file->queuehead->next = malloc( sizeof( *file->queuehead->next ) );
		memset( file->queuehead->next, 0, sizeof( *(file->queuehead->next) ) );
		file->queuehead->next->last = file->queuehead;
		file->queuehead = file->queuehead->next;
    }
	file->queuehead->buf = buffer;
	file->queuehead->size = size;
	file->queuehead->count = count;
	file->queuehead->pos = pos;
	file->queuehead->operation = ++OperationNumber;        
	file->queuehead->next = NULL;
    return( file->queuehead->operation );
}



void ExecuteOperation( aFILE * file ) {
    // if file is busy or there are no pending ops, we can't do much.
    if( !QueueEmpty( file ) && !file->busy ) {
        int err = 0;
        if( file->mode == A_WRITE ) {
            err = ExecuteWrite( file, file->queuetail );
        } else {
            err = ExecuteRead( file, file->queuetail );
        }
        if( !err ) {
            // advance file pointer
            AddTo64( file->queuetail->size * file->queuetail->count, &(file->filep_high), &(file->filep_low) );
        }
        
    }
}


int QueryOpComplete( aFILE * file ) {
    // check to see if the last operation was completed.
    if( file->queuetail != NULL ) {
#if defined USE_LINUX_AIO
	return( QueryLastCompleteLinux( file ) );
#elif defined USE_POSIX_AIO
	return( QueryLastCompletePAIO( file ) );
#elif defined USE_LIBC
        return( 1 );
#elif defined USE_WIN32
        return( QueryLastCompleteWIN32( file ) );
#endif
    }
    return( 1 );
}




// for files open for writing, ensures that all data is
// safely on disk (flushes buffer cache).
void aFlush( aFILE *file ) {
#if defined(USE_POSIX_AIO)||defined(USE_LINUX_AIO)
	if( fsync( file->file_descriptor ) )
		perror("fsync");
#elif defined USE_LIBC
    if( fflush( file->libchandle ) ) {
        printf( "error flushing stdio libc file\n" );
    }
#elif defined USE_WIN32
    if( !FlushFileBuffers( file->w32handle ) ) {
        printf( "error flushing win32 file\n" );
    }
#endif
}

// get the size in bytes of a particular file
unsigned long long aStatFileSize( const char * path ) {
#if defined(USE_POSIX_AIO)||defined(USE_LINUX_AIO)
	struct stat stat_data;
	if( stat( path , &stat_data) ){
		perror(path);
		return 0;
	}
	return stat_data.st_size;
#elif defined USE_LIBC
#error "libc aStatSize not implemented"
#elif defined USE_WIN32
	WIN32_FILE_ATTRIBUTE_DATA file_data;
	unsigned long long f_size;
	GetFileAttributesEx( path, GetFileExInfoStandard, (void*)&file_data );
	f_size = file_data.nFileSizeHigh;
	f_size <<= 32;
	f_size += file_data.nFileSizeLow;
	return f_size;
//#error "Implement me!  WIN32 aStatSize"
#endif
}


// get the size in records of a particular file
// used when skipping the binning phase
unsigned long aStatSize( const char * path ) {
#if defined(USE_POSIX_AIO)||defined(USE_LINUX_AIO)
	struct stat stat_data;
	if( stat( path , &stat_data) ){
		perror(path);
		return 0;
	}
	return stat_data.st_size / sizeof(record_t);
#elif defined USE_LIBC
#error "libc aStatSize not implemented"
#elif defined USE_WIN32
	return aStatFileSize( path ) / sizeof(record_t);
	printf("Implement me!  WIN32 aStatSize");
//#error "Implement me!  WIN32 aStatSize"
#endif
}


void aUpdateOperations( aFILE * file ) {
    int op_complete;
    // if we are busy, see if the last thing has completed.
    op_complete = QueryOpComplete( file );
    if( !op_complete ) {
        //printf( "op not yet complete on file 0x%X\n", file );
    }
    if( !QueueEmpty( file ) && file->busy && op_complete ) {
        RemoveOperation( file );
        file->busy = 0;
    }
    // if the queue is still not empty, start the next one up.
    if( !QueueEmpty( file ) ) {
        ExecuteOperation( file );
    }
        
}




int aAct( void * buffer, offset_t size, offset_t count, aFILE * file, offset_t pos ) {
    int operation = 0;
    // enter the operation in the queue, and then
    // try to execute what we can.
    // enqueue the op.
    operation = EnqueueOperation( buffer, size, count, file, pos );
    // execute operations
    ExecuteOperation( file );
    return( operation );
}


// these allow you to queue reads and writes.
// these return 0 for a failure, or an operation
// code that can be checked for completion with
// a_OperationComplete
int aWrite( void * buffer, offset_t size, offset_t count, aFILE * file, offset_t pos ) {
    return( aAct( buffer, size, count, file, pos ) );
}
int aRead( void * buffer, offset_t size, offset_t count, aFILE * file, offset_t pos ) {
    return( aAct( buffer, size, count, file, pos ) );
}


// returns 1 if the operation was completed, 0 otherwise.
int aOperationComplete( aFILE * file, int operation ) {
    aIORec *qp;
    // scan through the queue until we find the op
    // or we get to the end.  if we get to the end
    // and don't find it, it must have completed,
    // otherwise it hasn't.
    for( qp = file->queuetail; qp != NULL; qp = qp->next ) {
        if( qp->operation == operation ) {
            return( 0 );
        }
    }
    return( 1 );
}


// returns 1 if the file is doing IO, 0 otherwise.
int aFileBusy( aFILE * file ) {
    return( file->busy );
}


// blocks and waits for the specified operation to
// complete.
void aWaitComplete( aFILE * file, int operation ) {
    while( !aOperationComplete( file, operation ) ) {
        aUpdateOperations( file );
    }
}


// blocks and waits for the file to not be busy
// and for all IO operations to complete.
void aWaitNotBusy( aFILE * file ) {
    while( file->busy ) {
        aUpdateOperations( file );
    }
}


