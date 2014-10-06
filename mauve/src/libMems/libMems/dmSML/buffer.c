#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <time.h>
#include <stddef.h>
#include "libMems/dmSML/buffer.h"
#include <string.h>

// portably fills an int with reasonably random bits.
// one assumption is that MAX_RAND is bigger than 256.
static int BigRandom() {
    static char firsttime = 1;
    int i, result;
    if( firsttime ) {
        firsttime = 0;
        srand( 0 );
        //srand( time( NULL ) );
    }
    
    result = 0;
    for( i = 0; i < sizeof( result ); i++ ) {
        result <<= sizeof( result );
        result ^= rand();
    }
    // the funny test here because if result == INT_MIN on a
    // two's complement machine, -result *also* == INT_MIN.
    return( result < 0 ? (-result < 0 ? 0 : -result) : result );
}


// Working Set support.
// returns resulting size of the entire structure.
int MakeWorkingSet( working_set_t * ws, offset_t goalsize, offset_t minrecs, offset_t maxrecs ) {
// wrap the memory allocation loop with an outer loop
// that will attempt smaller working set sizes if large ones fail to allocate
	while( 1 ){

	    // we incrementally grow the working set to the desired size.
	    // however, we just compute the growth and how big the buffers will be,
	    // then we malloc one single large chunk of memory, and arrange things
	    // such that all of the buffer_ts are contiguous in one chunk, and
	    // all the actual data is contiguous after.
	    offset_t cursize = 0;
	    offset_t overhead = sizeof( ws->bufs[0] );
	    offset_t minsize = overhead + minrecs * sizeof( record_t );
	    offset_t maxsize = overhead + maxrecs * sizeof( record_t );
	    offset_t nbufs = 0;      // number of real buffers pleged to the working set
	    offset_t maxbufs = 256;  // the max number of buffers we track (this grows if necessary)
	    offset_t *buflist = malloc( sizeof( *buflist ) * maxbufs ); // grows when necessary
	    
	    record_t *recordptr;
	    offset_t i;
	    // if we can't possibly do anything useful
	    if( goalsize < minsize || maxrecs < minrecs || !buflist ) {
	    	if( buflist )
		        free( buflist );
	        return( 0 );
	    }

	    // just start allocating buffers until we can't anymore
	    while( goalsize - cursize >= maxsize ) {
	        offset_t randrecs = BigRandom() % (maxrecs - minrecs + 1) + minrecs;
	        if( nbufs == maxbufs ) {
	            // resize the array
	            maxbufs *= 2;
	            buflist = realloc( buflist, sizeof( *buflist ) * maxbufs );
	        }
	        buflist[nbufs++] = randrecs;
	        // update the number of bytes we've currently decided to allocate.
	        cursize += overhead + randrecs * sizeof( record_t );
	    }
		// now we have nbufs buffers, and the number of records they should
		// store is in the buflist list.
		// allocate one big chunk of memory
		printf( "allocating %llu bytes for working set (%llu bufs)\n", cursize, nbufs );

		ws->bufs = malloc( cursize );
		// if it failed to allocate try a smaller size
		if( !ws->bufs ){
			goalsize /= 2;
			continue;
		}

		ws->size = cursize;
		ws->nbufs = nbufs;
		// clear it out
		memset( ws->bufs, 0, cursize );

		// Now fill in the pointers to the records for all the buffers.
		// these all reside after the buffers in the working set.
		// Something convenient from this scheme is that in order to free
		// the working set when we're through, we just free ws->bufs.
		// pointer to first set of records.
		recordptr = (record_t *)( ((ptrdiff_t)ws->bufs) + (ws->nbufs * sizeof( ws->bufs[0] )) );
		for( i = 0; i < nbufs; i++ ) {
		    ws->bufs[i].totalrecs = buflist[i];
		    ws->bufs[i].recs = recordptr;
		    recordptr += ws->bufs[i].totalrecs;
		}

		free( buflist );
	    return( cursize );
	}
    return 0;
}





// Working Set support.
// Reorganize the working set with a different distribution of buffers.
void ReorganizeWorkingSet( working_set_t * ws, offset_t minrecs, offset_t maxrecs ) {
    // we incrementally grow the working set to the desired size.
    // however, we just compute the growth and how big the buffers will be,
    // then we malloc one single large chunk of memory, and arrange things
    // such that all of the buffer_ts are contiguous in one chunk, and
    // all the actual data is contiguous after.
    offset_t goalsize = ws->size;
    offset_t cursize = 0;
    offset_t overhead = sizeof( ws->bufs[0] );
    offset_t minsize = overhead + minrecs * sizeof( record_t );
    offset_t maxsize = overhead + maxrecs * sizeof( record_t );
    offset_t nbufs = 0;      // number of real buffers pledged to the working set
    offset_t maxbufs = 256;  // the max number of buffers we're tracking (this grows if necessary)
    offset_t *buflist = malloc( sizeof( *buflist ) * maxbufs ); // grows when necessary
    offset_t leftovers;
    record_t *recordptr;
    offset_t i;
    
    // if we can't possibly do anything useful
    if( maxrecs < minrecs ) {
        free( buflist );
        return;
    }

    if( goalsize < minsize ) {
        minsize = goalsize;
        minrecs = (minsize-overhead) / sizeof( record_t );
    }
    
    // just start allocating buffers until we can't anymore
    while( goalsize - cursize >= maxsize ) {
        offset_t randrecs = BigRandom() % (maxrecs - minrecs + 1) + minrecs;
        if( nbufs == maxbufs ) {
            // resize the array
            maxbufs *= 2;
            buflist = realloc( buflist, sizeof( *buflist ) * maxbufs );
        }
        buflist[nbufs++] = randrecs;
        // update the number of bytes we've currently decided to allocate.
        cursize += overhead + randrecs * sizeof( record_t );
    }
    
    // clean up the last bit
    if( goalsize - cursize > overhead ) {
        leftovers = (goalsize - cursize - overhead) / sizeof( record_t );
        if( leftovers ) {
            if( nbufs == maxbufs ) {
                // resize the array
                maxbufs *= 2;
                buflist = realloc( buflist, sizeof( *buflist ) * maxbufs );
            }
            buflist[nbufs++] = leftovers;
            cursize += overhead + leftovers * sizeof( record_t );
        }
    }

    // now we have nbufs buffers, and the number of records they should
    // store is in the buflist list.

    ws->nbufs = nbufs;
    // clear it out
    memset( ws->bufs, 0, cursize );
    // Now fill in the pointers to the records for all the buffers.
    // these all reside after the buffers in the working set.
    // Something convenient from this scheme is that in order to free
    // the working set when we're through, we just free ws->bufs.
    // pointer to first set of records.
    recordptr = (record_t *)( ((ptrdiff_t)ws->bufs) + (ws->nbufs * sizeof( ws->bufs[0] )) );
    for( i = 0; i < nbufs; i++ ) {
        ws->bufs[i].totalrecs = buflist[i];
        ws->bufs[i].recs = recordptr;
        recordptr += ws->bufs[i].totalrecs;
    }
    
    free( buflist );
    return;
}









// this updates all the IO on the working set buffers, querying those that
// are not in OP_FINISHED or OP_NONE and putting those that finish into OP_FINISHED
void UpdateWSIOFinishedState( working_set_t * ws ) {
    // gets rid of an indirection in the loop.
    // this method (rather than using an index) 
    // (I also think it's cleaner)
    buffer_t *b;
    // simply walk all of them 
    for( b = ws->bufs; b - ws->bufs < ws->nbufs; b++ ) {
        // real operation #s are whole numbers.
        if( b->operation > OP_NONE ) {
            //printf( "examining operation %d\n", b->operation );
            if( aOperationComplete( b->file, b->operation ) ) {
                //printf( "* Completed operation %d on device %x\n", b->operation, b->device );
                b->operation = OP_FINISHED;
            } else {
                //printf( "operation %d INCOMPLETE IO\n", b->operation );
            }
        }
    }
}



// buffer list manipulations
// returns argument
buffer_list_t * InitList( buffer_list_t * list ) {
    list->head = NULL;
    list->nitems = 0;
    return( list );
}


void PushHead( buffer_list_t * list, buffer_t * item ) {
    // one special case for empty list, because we can't
    // dereference list->head until we assign to it.
    if( list->head == NULL ) {
        list->head = item;
        list->nitems = 1;
        list->head->next = list->head;
        list->head->last = list->head;
        return;
    }
    // other cases are easier, because no more null pointers.
    item->last = list->head->last;
    item->next = list->head;
    list->head->last->next = item;
    list->head->last = item;
    list->head = item;
    // we added an item.
    list->nitems++;
}

void PushTail( buffer_list_t * list, buffer_t * item ) {
    // this is exactly equivalent to doing a PushHead and
    // then backing up the list head one.
    // get the item in there
    PushHead( list, item );
    // back up the head.
    list->head = list->head->last;
}

buffer_t * PopHead( buffer_list_t * list ) {
    buffer_t *ret;
    // just get rid of the head item and return it.
    if( list->head == NULL ) {
        return( NULL );
    }
    list->head->next->last = list->head->last;
    list->head->last->next = list->head->next;
    ret = list->head;
    list->head = list->head->next;
    ret->next = ret->last = NULL;
    list->nitems--;
    if( list->nitems == 0 ) {
        list->head = NULL;
    }
    return( ret );
}

buffer_t * PopTail( buffer_list_t * list ) {
    // just get rid of the tail item and return it.
    if( list->head == NULL ) {
        return( list->head );
    }
    // otherwise, a pop tail is equivalent to moving the
    // head back one and popping head.
    list->head = list->head->last;
    return( PopHead( list ) );
}

// returns second argument
buffer_t * RemoveItem( buffer_list_t * list, buffer_t * item ) {
    // FIXME: handle NULL cases in a reasonable way?
    if( item == list->head ) {
        return( PopHead( list ) );
    }
    item->next->last = item->last;
    item->last->next = item->next;
    item->next = item->last = NULL;
    list->nitems--;
    if( list->nitems == 0 ) {
        list->head = NULL;
    }
    return( item );
}






int CompareKeys_qsort_wrapper( const void *r1, const void *r2 ) {

    return( CompareKeys( (record_t *)r1, (record_t *)r2 ) );

}



int CompareKeys( const record_t *r1, const record_t *r2 ) {

    return( COMPARE_KEYS( *r1, *r2 ) );
    //return( memcmp( r1->key, r2->key, sizeof( r1->key ) ) );

}







// This *must* enforce a serialized order for reading and writing, lest
// we write sorted data out in the wrong order!
void UpdateDeviceIOExecuteState( working_set_t * ws, iodevice_t * dev ) {
    // check to see if the device's IO job completed    
    if( !dev->buf || dev->state == DEV_FREE || dev->buf->operation == OP_FINISHED ) {
        // find another job to take its place and execute it.
        buffer_t *b;
        buffer_t *found_buf = NULL;
        dev->state = DEV_FREE;
        dev->buf = NULL;
        // simply walk all of them, find the operation on this device
        // that has the lowest op number for its file.  This is made "more fair"
        // by picking the first operation that matches the device, then finding
        // all other buffers that operate on the same file
        for( b = ws->bufs; b - ws->bufs < ws->nbufs; b++ ) {
            // is this one that should be executed next?
            
            if( b->operation == OP_PENDING && b->device == dev ) {
                if( !found_buf ) {
                    found_buf = b;
                } else if( (b->file == found_buf->file) && 
                    (b->fileop < found_buf->fileop) ) {
                    found_buf = b;
                }
            }
            
            /*
            if( b->operation == OP_PENDING && b->device == dev ) {
            dev->buf = b;
            b->operation = b->file->mode == A_READ 
            ? aRead( b->recs, sizeof( b->recs[0] ), b->numrecs, b->file )
            : aWrite( b->recs, sizeof( b->recs[0] ), b->numrecs, b->file );
            dev->state = DEV_BUSY;
            //printf( "* Created operation %d on device %x\n", b->operation, b->device );
            // found one, so quit.
            break;
            }
            */
        }
        
        if( found_buf ) {
            dev->buf = found_buf;
            found_buf->operation = found_buf->file->mode == A_READ 
                ? aRead( found_buf->recs, 1, found_buf->io_size, found_buf->file, found_buf->io_pos )
                : aWrite( found_buf->recs, 1, found_buf->io_size, found_buf->file, found_buf->io_pos );
            dev->state = DEV_BUSY;
        }
        
    }
}

// read and write to/from disk.
void ReadBuffer( buffer_t * buffer, offset_t num_recs, iodevice_t * dev ) {
	buffer->io_size = num_recs * sizeof( record_t );
    buffer->numrecs = num_recs;
    buffer->device = dev;
    buffer->fileop = buffer->file->op++;
    buffer->io_pos = CURRENT_POS;
    if( buffer->operation != OP_NONE ) {
        printf( "weird!\n" );
    } else {
        buffer->operation = OP_PENDING;
    }
    //printf( "* Initiated (pending) operation on %x\n", dev );
}

void WriteBuffer( buffer_t * buffer, offset_t num_recs, iodevice_t * dev ) {
    // exactly the same as a read -- the operation is just scheduled.
    // the exact nature (read or write) is determined by the mode
    // of the opened file at the time operation is in fact
    // executed.
    ReadBuffer( buffer, num_recs, dev );
    
}

