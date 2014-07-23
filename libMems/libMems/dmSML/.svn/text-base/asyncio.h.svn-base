#ifndef _asyncio_h_
#define _asyncio_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

//#define USE_LINUX_AIO
//#define USE_LIBC_AIO	// don't use kaio

#ifdef WIN32
#   define WIN32_LEAN_AND_MEAN
#   include <windows.h>
#   define USE_WIN32
#else 
#   ifndef _LARGEFILE64_SOURCE
#	define _FILE_OFFSET_BITS 64
#	define _LARGEFILE_SOURCE
#	define _LARGEFILE64_SOURCE
#   endif
// use kaio by default
#	if defined(USE_LIBC_AIO) || defined(USE_POSIX_AIO)
#		ifdef HAVE_SYS_TYPES_H
#			include <sys/types.h>
#		endif
#		if defined HAVE_SYS_AIO_H
#			include <sys/aio.h>
#		elif HAVE_AIO_H
#			include <aio.h>
#		endif
#		ifdef HAVE_FEATURES_H
#			include <features.h>
#		endif
typedef struct aiocb aiocb_t;
#	endif
#	ifdef USE_LINUX_AIO
#		define _FILE_OFFSET_BITS 64
#		define _LARGEFILE_SOURCE
#		define _LARGEFILE64_SOURCE
#		include <libaio.h>
typedef struct iocb iocb_t;
#	endif
#	ifdef HAVE_FEATURES_H
#		include <features.h>
#	endif
#endif

#include <stdlib.h>
#include <stdio.h>


#define CURRENT_POS -1
typedef unsigned long long offset_t;

// is this a struct to store RECORDS to write out?
// the way it's used looks like it's not intended
// for generic data...
typedef struct _aIORec {
#if defined USE_POSIX_AIO
// posix aio uses the aiocb_t type to describe aio requests
	aiocb_t *aio_cb;
#elif defined USE_LINUX_AIO
	iocb_t* aio_cb;
#elif defined USE_LIBC
#elif defined USE_WIN32
    // win32-specific data.
    // this is a pointer because windows needs it to
    // be in a fixed spot.  But we have to resize the
    // data structure that contains these, so we need
    // to allocate them separately.
    // unfortunately, this means we need to do linear
    // search to figure out what thing in the queue some
    // completion corresponds to.  Fortunately, this
    // rarely needs to be done.  I think this is The
    // Right Thing, given the tools and our goals.
    OVERLAPPED * w32overlapped;
#endif
    // must do linear search to find specific operations,
    // but no big deal.
    int operation;
    char * buf;
    offset_t size;
    offset_t count;	//what is count for??
    offset_t pos;
    struct _aIORec * next;
    struct _aIORec * last;
    
} aIORec;


// users don't need to concern themselves with this.
typedef struct _aFILE {
#if defined(USE_POSIX_AIO)||defined(USE_LINUX_AIO)
	int file_descriptor;
#elif defined USE_LIBC
    FILE * libchandle;
#elif defined USE_WIN32
    HANDLE w32handle;
#endif
    // read or write (both read and write??)
    int mode;
    // file seek pointer
    unsigned int filep_high;
    unsigned int filep_low;
    // is a read/write operation in progress?
    int busy;
    // operation serial number (to ensure serial operation).
    int op;
    // are we to be closed?
    int toclose;
    // queue of io operations
    aIORec *queuehead, *queuetail;
} aFILE;


enum {
    A_READ,
    A_WRITE
};


// these work just like fopen and fclose
aFILE * aOpen( const char * path, int mode );
// close will block until all operations
// on the file are complete.
int aClose( aFILE * file );

// these allow you to queue reads and writes.
// these return 0 for a failure, or an operation
// code that can be checked for completion with
// a_OperationComplete
int aWrite( void * buffer, offset_t size, offset_t count, aFILE * file, offset_t pos );
int aRead( void * buffer, offset_t size, offset_t count, aFILE * file, offset_t pos );

// returns 1 if the operation was completed, 0 otherwise.
int aOperationComplete( aFILE * file, int operation );

// returns 1 if the file is doing IO, 0 otherwise.
int aFileBusy( aFILE * file );

// blocks and waits for the specified operation to
// complete.
void aWaitComplete( aFILE * file, int operation );
    
// blocks and waits for the file to not be busy
// and for *all* IO operations to complete.
void aWaitNotBusy( aFILE * file );

// polls the aio file to see if anything's completed, and
// starts the next queued up jobs if they are.  does not
// block.
void aUpdateOperations( aFILE * file );


// for files open for writing, ensures that all data is
// safely on disk (flushes buffer cache).
void aFlush( aFILE *file );

// get the size in records of a particular file
// used when skipping the binning phase
unsigned long aStatSize( const char * path );

// get the size in bytes of a particular file
unsigned long long aStatFileSize( const char * path );

#endif /* _asyncio_h_ */
