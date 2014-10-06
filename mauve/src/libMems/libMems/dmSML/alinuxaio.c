#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/dmSML/alinuxaio.h"
#ifdef USE_LINUX_AIO

#include <libaio.h>

#include "libMems/dmSML/asyncio.h"
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
/*
#define __NR_io_setup		245
#define __NR_io_destroy		246
#define __NR_io_getevents	247
#define __NR_io_submit		248
#define __NR_io_cancel		249
*/

io_context_t ctx_id = NULL;

#ifndef __u64
typedef unsigned long long __u64;
#endif

__u64 current_id = 0;

unsigned event_max = 10000;

// error = sys_io_destroy( ctx_id );


typedef struct completion_id_s {
	__u64 data;
	struct completion_id_s* next;
	struct completion_id_s* last;
} completion_id_t;

typedef struct completion_id_list_s { 
    int nitems;
    completion_id_t * head;
} completion_id_list_t;

// buffer list manipulations
// returns argument
completion_id_list_t * InitListComp( completion_id_list_t * list );
void PushHeadComp( completion_id_list_t * list, completion_id_t * item );
void PushTailComp( completion_id_list_t * list, completion_id_t * item );
completion_id_t * PopHeadComp( completion_id_list_t * list );
completion_id_t * PopTailComp( completion_id_list_t * list );
// returns second argument
completion_id_t * RemoveItemComp( completion_id_list_t * list, completion_id_t * item );


// buffer list manipulations
// returns argument
completion_id_list_t * InitListComp( completion_id_list_t * list ) {
    list->head = NULL;
    list->nitems = 0;
    return( list );
}


void PushHeadComp( completion_id_list_t * list, completion_id_t * item ) {
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

void PushTailComp( completion_id_list_t * list, completion_id_t * item ) {
    // this is exactly equivalent to doing a PushHead and
    // then backing up the list head one.
    // get the item in there
    PushHeadComp( list, item );
    // back up the head.
    list->head = list->head->last;
}

completion_id_t * PopHeadComp( completion_id_list_t * list ) {
    completion_id_t *ret;
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

completion_id_t * PopTailComp( completion_id_list_t * list ) {
    // just get rid of the tail item and return it.
    if( list->head == NULL ) {
        return( list->head );
    }
    // otherwise, a pop tail is equivalent to moving the
    // head back one and popping head.
    list->head = list->head->last;
    return( PopHeadComp( list ) );
}

// returns second argument
completion_id_t * RemoveItemComp( completion_id_list_t * list, completion_id_t * item ) {
    // FIXME: handle NULL cases in a reasonable way?
    if( item == list->head ) {
        return( PopHeadComp( list ) );
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


completion_id_list_t *completion_list = NULL;

int OpenLinux( aFILE * file, const char *path, int mode ){
	long error;
	if( ctx_id == 0 ){
		error = io_queue_init( event_max, &ctx_id );
		if( error != 0 )
			perror( "io_setup" );
	}
	if( completion_list == NULL ){
		completion_list = (completion_id_list_t*)malloc( sizeof( completion_id_list_t ) );
		completion_list = InitListComp( completion_list );
	}
		
	if(mode == A_READ){
		file->file_descriptor = open(path, O_LARGEFILE | O_RDONLY, S_IREAD | S_IWRITE | S_IRGRP | S_IWGRP );
	}else{
		file->file_descriptor = open(path, O_RDWR | O_CREAT | O_TRUNC | O_LARGEFILE,  S_IREAD | S_IWRITE | S_IRGRP | S_IWGRP);
	}
	if(file->file_descriptor < 0){
		
		perror(path);
	}
	return file->file_descriptor >= 0;
}

int CloseLinux( aFILE * file ){	
	return close( file->file_descriptor ) == 0;
}

void CleanupLinux(){
	// free the completion list
	free( completion_list );
	completion_list = NULL;
	ctx_id = NULL;
}

int FillAIOStruct( aFILE * file, aIORec * rec ){
// fill the request data structure
	rec->aio_cb = (iocb_t*) malloc( sizeof(iocb_t));
	if(rec->aio_cb == 0)
		return 0;

	memset(rec->aio_cb, 0, sizeof(iocb_t));
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

//	rec->aio_cb->aio_data = current_id++;
	rec->aio_cb->aio_fildes = file->file_descriptor;
	rec->aio_cb->u.c.offset = file->filep_high;
	rec->aio_cb->u.c.offset <<= 32;
	rec->aio_cb->u.c.offset |= file->filep_low;
	rec->aio_cb->u.c.buf = rec->buf;
	rec->aio_cb->u.c.nbytes = rec->size * rec->count;
	
	return 1;
}

int WriteLinux( aFILE * file, aIORec * rec ){
        int req_error;
	struct iocb *request_array[] = { rec->aio_cb };
	if( FillAIOStruct( file, rec ) ){
		// request the io
		rec->aio_cb->aio_lio_opcode = IO_CMD_PWRITE;
		req_error = io_submit( ctx_id, 1, &rec->aio_cb );
		if(req_error != 1){
			printf("write_submit: io_submit res=%d [%s]\n", req_error, strerror(-req_error));
            printf( "aiocb->aio_fildes = %d\n", rec->aio_cb->aio_fildes );
            printf( "aiocb->u.c.offset = %llu\n", rec->aio_cb->u.c.offset );
            printf( "aiocb->u.c.buf = %lx\n", rec->aio_cb->u.c.buf );
            printf( "aiocb->u.c.nbytes = %llu\n", rec->aio_cb->u.c.nbytes );
            printf( "aiocb->aio_reqprio = %d\n", rec->aio_cb->aio_reqprio );
		}
		return req_error == 1;
	}
	return 0;
}

int ReadLinux( aFILE * file, aIORec * rec ){
	int req_error;
	struct iocb *request_array[] = { rec->aio_cb };
// fill the request data structure
	if( FillAIOStruct( file, rec ) ){
	// request the io
		rec->aio_cb->aio_lio_opcode = IO_CMD_PREAD;
		req_error = io_submit( ctx_id, 1, &rec->aio_cb );
        if(req_error != 1){
			printf("read_submit: io_submit res=%d [%s]\n", req_error, strerror(-req_error));
//                printf( "aiocb->aio_filedes = %d\n", rec->aio_cb->aio_filedes );
//                printf( "aiocb->aio_offset = %llu\n", rec->aio_cb->aio_offset );
//                printf( "aiocb->aio_buf = %lx\n", rec->aio_cb->aio_buf );
//                printf( "aiocb->aio_nbytes = %llu\n", rec->aio_cb->aio_nbytes );
                printf( "aiocb->aio_reqprio = %d\n", rec->aio_cb->aio_reqprio );
        }
		return req_error == 1;
	}
	return 0;
}


// PRECONDITION:  file->queuetail is not null
// simply queries wether the first request submitted to the file has
// completed yet.
int QueryLastCompleteLinux( aFILE * file ){
	int rval;
	int compI;
	completion_id_t *comp;
	struct io_event ioe;
	struct timespec zero_wait;

	zero_wait.tv_sec = 0;
	zero_wait.tv_nsec = 10000000;
	
	rval = io_getevents( ctx_id, 0, 1, &ioe, &zero_wait );
	if( rval == 1 ){
		completion_id_t *completion = (completion_id_t*)malloc( sizeof(completion_id_t) );
		completion->data = ioe.data;
		PushTailComp( completion_list, completion );
	}
	comp = completion_list->head;
	for( compI = 0; compI < completion_list->nitems; compI++ ){
		if( comp->data == ioe.data )
			break;
	}
	if( compI != completion_list->nitems ){
		RemoveItemComp( completion_list, comp );
		return 1; // success
	}
	return 0;	// hasn't completed yet
}

#endif
