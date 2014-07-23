#ifndef _sorting_h_
#define _sorting_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/dmSML/buffer.h"


// define this if you want to use the qsort only version
// of dmsort.
#define USE_QSORT_ONLY




// START configurable values

// the number of bits in each radix.
#define RADIX_BITS 12

// the number of bins to qsort during each call to RadixSort during the qsort phase
#define SORT_BINS_SIZE 1000
// the number of records to copy into sorted order during each call to CopySortedData
#define COPY_CHUNK_SIZE 50000
#define HISTOGRAM_CHUNK_SIZE 50000
#define PTR_COPY_CHUNK_SIZE 50000

// END configurable values

// sorting states -- this is for the second phase, after binning
#define WAIT_WRITE      (-100)
#define WAIT_READ       (-200)
#define SORTING         (-300)
#define BUSY_READ       (-400)
#define BUSY_WRITE      (-500)
#define SORTING_SCRATCH (-600)
#define WRITE_RESTRUCTURE (-700)

enum{
	CalculateHistogram = 0,	// At this stage we compute a histogram on the current radix
	CopyPointers = 1,		// This stage copies the pointers into (more) sorted order
	QsortPointers = 2,	// This stage qsorts the pointers
	CopyData = 3		// This stage copies the data into totally sorted order
};

typedef struct sort_buf_s {
    int state;          // WAIT_READ, WAIT_WRITE, SORTING, BUSY
    int bin;            // what bin this buffer holds right now.
    iodevice_t *dev;
    buffer_t *buf;		// the buffer where records live
    buffer_t *radix_tmp;		// temp space for the radix sort copy
    record_t **rec_ptrs;		// array of pointers to records
	
	unsigned base_number;
	unsigned divisor;
	unsigned histogram_size;
    unsigned *histogram;		// the histogram of bins
    unsigned *cur_ptr_offsets;	// the locations to copy data in each histogram bucket
	unsigned cur_position;	// the current record or bin position in the current stage.
	int sort_state;			// current state of the sort algorithm
} sort_buf_t;

// Need NumBins so that we can compute the amount already sorted
extern int NumBins;

//typedef unsigned long long uint64;

/* Fills and returns a new sort_buf_t with the appropriate
 * data.
 */
void InitRadixSort( sort_buf_t* sortbuf, buffer_t* scratch_buffer );

/* Checks the current state of the radix sort and performs a fixed
 * amount of sorting computation before returning.
 * call until state is set to WriteData
 */
void RadixSort( sort_buf_t* sortbuffer );

#endif /* _sorting_h_ */
