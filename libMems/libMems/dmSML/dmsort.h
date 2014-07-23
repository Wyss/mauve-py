#ifndef __DMSORT_H__
#define __DMSORT_H__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "libMems/dmSML/util.h"
#include "libMems/dmSML/timing.h"
#include "libMems/dmSML/asyncio.h"
#include "libMems/dmSML/buffer.h"
#include "libMems/dmSML/sorting.h"
#include "libMems/dmSML/sml.h"

// define this if you're using the ASCII sortgen data.
// don't define if you're using random data (dmsortgen)
//#define ASCII_KEYBYTES

// define this if using dmSML with sequences that have large
// stretches of NNNNN...  such as an unfinished eukaryote
//#define NNNNN_KEYBYTES

// define this if you want to measure the overlapping
// of your sorting with I/O in the sorting phase --
// this makes the sort routine do nothing.
//#define NO_SORT_PERF_TEST

// define the following if you don't want to write
// data during the sort phase in order to get timings
//#define NO_WRITE_PERF_TEST

// define this to skip the binning phase in order to
// perform measurements on the sort phase.  The bin
// files to use during sorting must already exist (duh!)
//#define NO_BINNING_PERF_TEST

// define this to test the performance of binning and
// restructuring without bin writing
//#define NO_BIN_WRITE_PERF_TEST

// define this to test the performance without restructuring
// each SML bin
//#define NO_RESTRUCTURE_PERF_TEST

#ifndef NELEMS
#define NELEMS(x) \
    ( sizeof((x)) / sizeof((x)[0]) )
#endif

#define MIN(x,y)    ((x)<(y)?(x):(y))
#define MINRECS     (1311)
#define MAXRECS     (1311)


// this is somewhat less appealing than a config file,
// but speed is critical and parsing a config file at
// startup is just inconvenient.  Besides, specifying
// what we care about is easy enough this way.
typedef struct device_s {
    const char      *devname;
    const char      *path;
    iodevice_t      dev;
} device_t;


// ugly hack
#define BIN_SPECIAL     (-10000)



// what we use to represent a bin.
typedef struct bin_s {
    aFILE               *file;      // File we write/read on.
    int                 dev;        // This is an index into the Devices table.
    offset_t            nrecs;      // Number of records written to bin.
    buffer_list_t       bufs;       // Our list of buffers that holds our data.
    char*				fname;		/**< The file name of this bin */
} bin_t;

typedef struct seqbuf_s {
	aFILE				*file;		// Output file
	int					dev;		// device table index for output file
	offset_t			bufpos;		// position in current buffer
	uint64				seq_pos;	// position in sequence that is next to translate
	buffer_list_t		bufs;		// list of buffers for data
} seqbuf_t;

enum dm_errors {
	SUCCESS,
	TOO_FEW_BINS,
	TOO_MANY_BINS,
	INPUT_NOT_OPENED,
	INVALID_WS_SIZE,
	SEQUENCE_TOO_SHORT,
	OUTPUT_NOT_OPENED,
	INVALID_NUMRECS,
	NO_FREE_BUFFERS,
	BIN_NOT_OPENED,
};


void print_usage( const char* pname );


static buffer_t * AllocateFree( void );

static int ComputeBinNumber( const unsigned char key[10] );

// just like ComputeBinNumber except we reserve one bin for zero keys.
static int ComputeNNNNNBinNumber( const unsigned char key[10] );

static int ComputeAsciiBinNumber( const unsigned char key[10] );

static void DoBinning( void );

void FinishBinning();

offset_t CalculateDataReadSize( buffer_t* b );

static void DoReading( void );

static void HandleBinWriteCompletions( void );

static void HandleSeqbufWriteCompletions( void );

#define ALPHA_BITS 2

static void Translate32(uint32* dest, const char* src, const unsigned len);

void RestructureReadSMLBins( void );

static void HandleReadingCompletions( void );

int InitdmSML( long working_mb, long buffer_size, const char* input_filename, const char* output_filename, const char* const* scratch_paths, uint64 seed );

void DisplayStatusHeader( void );

void DisplayStatus( void );

void UpdateIOState( void );

void EnsureAllOperationsComplete( void );

void BinningPhase( void );

void SortReading( void );

#ifdef USE_QSORT_ONLY

int comp_keys( record_t a, record_t b );

void QBrute( record_t a[], int lo, int hi );

void QSort( record_t a[], int lo0, int hi0 );

void RecSort( record_t a[], int nelems );

int SortBuffer( buffer_t * buf );

void SortSorting( void );

#elif defined NO_SORT_PERF_TEST

void SortSorting( void );

#else 

sort_buf_t* CurrentSortBuf;
buffer_t* SortScratchBuffer;

void SortSorting( void );

#endif

void RestructureSMLBinsForWrite( void );

int CalculateSortWriteSize( int sortI );

void SortWriting( void );

void SortHandleCompletions( void );

void SortUpdateIOState();

void SortingEnsureAllOperationsComplete();

void SortingPhase( void );

int dmsort( void );

int dmSML( const char* input_file, const char* output_file, const char* const* scratch_paths, uint64 seed );


#endif // __DMSORT_H__
