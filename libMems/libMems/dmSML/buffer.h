#ifndef _buffer_h_
#define _buffer_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/dmSML/asyncio.h"

// forward decl for the benefit of iodevice_t
// (can't be avoided)
typedef struct buffer_s buffer_t;

enum {
    DEV_FREE,
    DEV_BUSY,
};

/*
================
iodevice_t
An IO device represents a physical disk.  This is used to make
sure that we're not doing more than one operation on any disk
at a time.  The reason is if we are, the OS threads that do the
asynchronous IO will contend for the disk and we'll start seeking.
Seeking is bad.
================
*/
typedef struct iodevice_s {
    int         op;     // an operation number used to enforce serial operation on each device.
    int         state;  // either DEV_FREE or DEV_BUSY as above.
    buffer_t    *buf;   // if we're DEV_BUSY, the buffer we're operating on.
} iodevice_t;


/*
================
buffer_state_t
Buffers are used for IO, and the IO is asynchronous.  Buffers
have a bit of state to indicate what their current IO status is.
Ordinarily, buffers are in an OP_NONE state, to indicate no
operation is being performed.  When an operation is initiated,
the buffer transitions to OP_PENDING and changes to a valid operation
code when the operation actually starts.  When the operation completes,
the buffer is transitioned to the OP_FINISHED state, so that the
program may determine when buffers have completed their operations.
Then the app should transition the state back to OP_NONE.
================
*/
enum {
    OP_PENDING = -2,
    OP_FINISHED = -1,
    OP_NONE = 0,
};



#define CompareKeyPtrs( a, b ) \
((int)(a)->key[0]-(int)(b)->key[0] ? (int)(a)->key[0]-(int)(b)->key[0] :    \
(int)(a)->key[1]-(int)(b)->key[1] ? (int)(a)->key[1]-(int)(b)->key[1] :     \
(int)(a)->key[2]-(int)(b)->key[2] ? (int)(a)->key[2]-(int)(b)->key[2] :     \
(int)(a)->key[3]-(int)(b)->key[3] ? (int)(a)->key[3]-(int)(b)->key[3] :     \
(int)(a)->key[4]-(int)(b)->key[4] ? (int)(a)->key[4]-(int)(b)->key[4] :     \
(int)(a)->key[5]-(int)(b)->key[5] ? (int)(a)->key[5]-(int)(b)->key[5] :     \
(int)(a)->key[6]-(int)(b)->key[6] ? (int)(a)->key[6]-(int)(b)->key[6] :     \
(int)(a)->key[7]-(int)(b)->key[7] ? (int)(a)->key[7]-(int)(b)->key[7] :     \
(int)(a)->key[8]-(int)(b)->key[8] ? (int)(a)->key[8]-(int)(b)->key[8] :     \
(int)(a)->key[9]-(int)(b)->key[9] ? (int)(a)->key[9]-(int)(b)->key[9] : 0)

#define COMPARE_KEYS( a, b ) \
((a).key[0]!=(b).key[0] ? (a).key[0]-(b).key[0] :    \
(a).key[1]!=(b).key[1] ? (a).key[1]-(b).key[1] :     \
(a).key[2]!=(b).key[2] ? (a).key[2]-(b).key[2] :     \
(a).key[3]!=(b).key[3] ? (a).key[3]-(b).key[3] :     \
(a).key[4]!=(b).key[4] ? (a).key[4]-(b).key[4] :     \
(a).key[5]!=(b).key[5] ? (a).key[5]-(b).key[5] :     \
(a).key[6]!=(b).key[6] ? (a).key[6]-(b).key[6] :     \
(a).key[7]!=(b).key[7] ? (a).key[7]-(b).key[7] :     \
(a).key[8]!=(b).key[8] ? (a).key[8]-(b).key[8] :     \
(a).key[9]!=(b).key[9] ? (a).key[9]-(b).key[9] : 0)




// this is the record as in the files to be sorted
typedef struct record_s {
    unsigned char key[10];
    unsigned char num[1];
    unsigned char payload[1];
} record_t;





int CompareKeys_qsort_wrapper( const void *r1, const void *r2 );
int CompareKeys( const record_t *r1, const record_t *r2 );





/*
================
buffer_t
This is the unit of information most commonly dealt with.
We read into these, and write these out, and use these for
binning.  A single Working Set of these buffers should be
used for the duration of the program, and they should be
managed with the buffer lists below.
sizeof( buffer_t ) == 32, so the overhead isn't bad.
================
*/
struct buffer_s {
    
    aFILE           *file;      // the file this buffer is attached to for IO ops
    iodevice_t      *device;    // which IO device this is on (for scheduling IO)
    int             operation;  // either OP_NONE, OP_FINISHED, or the op #.
    offset_t        numrecs;    // the number of valid records in this buffer.
    offset_t        totalrecs;  // number of real records in recs
    int             fileop;     // operation number on device to ensure serialized ops.
    record_t        *recs;      // actual record storage
    struct buffer_s *next;      // for chaining lists together.
    struct buffer_s *last;
    offset_t        io_size;	// amount of data for i/o, need not be equal to numrecs
    long long		input_pos;	// the sequence offset that this data was read from, only valid during binning phase
    offset_t		io_pos;		// the file offset for I/O, set to CURRENT_POS to use the current file seek pointer
};


/*
================
buffer_list_t
Buffer lists are used to manage pools, like the free list,
the reading list, the to process list, and a list for each bin.
We use circular lists because they're simpler.
================
*/
typedef struct buffer_list_s { 
    int nitems;
    buffer_t * head;
} buffer_list_t;


/*
================
working_set_t
Working sets are collections of buffers.  They are useful so that
you can use a fixed amount of memory to deal with things.  The
problem then becomes internal working set management.
================
*/
typedef struct working_set_s {
    offset_t    size;       // actual size of working set in bytes
    int         nbufs;
    buffer_t    *bufs;
} working_set_t;

    


// Working Set support.
// returns resulting size of the entire structure.
// goalsize is the desired size of the working set in bytes, minbufsize and maxbufsize
// are the minimum and maximum desired number of records in buffers.  The buffers will
// be allocated with random sizes in this range until the desired goalsize is reached.
// this will return 0 in the case of a malloc error or if the goalsize is too small to
// have any buffers allocated for it.
int MakeWorkingSet( working_set_t * ws, offset_t goalsize, offset_t minrecs, offset_t maxrecs );

// Working Set support.
// Reorganize the working set with a different distribution of buffers.
void ReorganizeWorkingSet( working_set_t * ws, offset_t minrecs, offset_t maxrecs );

    
// this updates all the IO on the working set buffers, querying those that
// are not in OP_FINISHED or OP_NONE and putting those that finish into OP_FINISHED
void UpdateWSIOFinishedState( working_set_t * ws );


// this updates the IO on a particular device.  this routine and the one above
// should probably be called as this one after that one, and in addtion, this
// one called for every device in the system.
void UpdateDeviceIOExecuteState( working_set_t * ws, iodevice_t * dev );


// buffer list manipulations
// returns argument
buffer_list_t * InitList( buffer_list_t * list );
void PushHead( buffer_list_t * list, buffer_t * item );
void PushTail( buffer_list_t * list, buffer_t * item );
buffer_t * PopHead( buffer_list_t * list );
buffer_t * PopTail( buffer_list_t * list );
// returns second argument
buffer_t * RemoveItem( buffer_list_t * list, buffer_t * item );

// read and write to/from disk.
void ReadBuffer( buffer_t * buffer, offset_t num_recs, iodevice_t * dev );
void WriteBuffer( buffer_t * buffer, offset_t num_recs, iodevice_t * dev );


#endif /* _buffer_h_ */

