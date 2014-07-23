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
#include "libMems/dmSML/dmsort.h"

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
// #define NO_BINNING_PERF_TEST

// define this to test the performance of binning and
// restructuring without bin writing
//#define NO_BIN_WRITE_PERF_TEST

// define this to test the performance without restructuring
// each SML bin
//#define NO_RESTRUCTURE_PERF_TEST

/*

#define NELEMS(x) \
    ( sizeof((x)) / sizeof((x)[0]) )


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
*/

device_t *Devices;
int NumDevices;



// ugly hack
//#define BIN_SPECIAL     (-10000)

int NSortBufs;
sort_buf_t *SortBufs;


// how the working set is allocated originally.
offset_t BufferSizeMin;
offset_t BufferSizeMax;



/*
// what we use to represent a bin.
typedef struct bin_s {
    aFILE               *file;      // File we write/read on.
    int                 dev;        // This is an index into the Devices table.
    offset_t            nrecs;      // Number of records written to bin.
    buffer_list_t       bufs;       // Our list of buffers that holds our data.
} bin_t;

*/

// number specified by a cmdline param at runtime.
bin_t   *Bins;
int     NumBins;
int     NumBinDevs;  // number of binning devices

/*
typedef struct seqbuf_s {
	aFILE				*file;		// Output file
	int					dev;		// device table index for output file
	offset_t			bufpos;		// position in current buffer
	uint64				seq_pos;	// position in sequence that is next to translate
	buffer_list_t		bufs;		// list of buffers for data
} seqbuf_t;
*/

seqbuf_t Seqbuf;
    
aFILE   *Data;       // the data to sort
int     DataDev;     // the device the data file is on.


const char *OutFileName = "unset";     // the output file name.
aFILE   *Output;     // the output file (sorted data goes here)
int     OutputDev;   // the device the output goes on.



int BinToRead, BinToWrite, BinToSort;



working_set_t   WS;      // the Working Set we use to do our sorting.

offset_t NumRecs;        // the total number of blocks to process
offset_t RecsProcessed;  // number of blocks processed (put in bins to write out)
offset_t RecsRead;       // number of records fully read in.
offset_t RecsUnread;     // number of blocks on disk (not yet had 'read' called)
offset_t RecsCommitted;  // number of records committed to be written.
offset_t RecsWritten;    // number of records actually written on disk.


// timers
double RunningTime;
dmtimer_t *RunningTimer;
double BinningTime;
dmtimer_t *BinningTimer;
double SortingTime;
dmtimer_t *SortingTimer;

double QSortTime;
dmtimer_t *QSortTimer;

double ReadIdleTime;
dmtimer_t *ReadIdleTimer;
double SortIdleTime;
dmtimer_t *SortIdleTimer;
double WriteIdleTime;
dmtimer_t *WriteIdleTimer;

// buffer lists
buffer_list_t   Free;           // the free list
buffer_list_t   ToProcess;      // list read and to be processed
buffer_list_t   Reading;        // the list that's waiting on stuff to read.
buffer_list_t   Restructure;	// buffers that need post-read and pre-binning processing


static buffer_t * AllocateFree( void ) {
    buffer_t * ret;
    if( Free.nitems ) {
        ret = PopHead( &Free );
    } else {
        printf( "error: called AllocateFree but free list is empty\n" );
        return( NULL );
    }
    ret->device = NULL;
    ret->file = NULL;
    ret->last = ret->next = NULL;
    ret->numrecs = 0;
    ret->operation = OP_NONE;
    return( ret );
}


static unsigned int divisor = 0;

static int ComputeBinNumber( const unsigned char key[10] ) {
    int i;
    unsigned int keyval = 0;
    // how many bits can we use for the binning number?
    // first time through, compute divisor
    // assume even distribution
    // strange constant is 256^3, because we're dealing
    // with effectively a base 256 number here, and we can
    // only handle 3 places without overflowing.
    if( divisor == 0 ) {
        divisor = (unsigned)16777216 / (unsigned)NumBins;
        // need ceiling of this
        divisor += (unsigned)16777216 % (unsigned)NumBins ? 1 : 0;
        printf( "Divisor is: %u\n", divisor );
    }
    // now we compute the number represented by the first 3
    // characters of the key, and divide it by divisor, the
    // integral part gives the bin number.
    for( i = 0; i < 3; i++ ) {
        keyval <<= 8;
        keyval += key[i];
    }
//    printf( "Key is %.2x %.2x %.2x \n", key[0],key[1], key[2] );
    
//    printf( "Keyval is: %u\n", keyval );
//    printf( "Bin is: %u\n", keyval / divisor );
    return( keyval / divisor );
}

// just like ComputeBinNumber except we reserve one bin for zero keys.
static int ComputeNNNNNBinNumber( const unsigned char key[10] ) {
    int i;
    unsigned int keyval = 0;
    if( divisor == 0 ) {
        divisor = (unsigned)16777216 / ((unsigned)NumBins - 1);
        // need ceiling of this
        divisor += (unsigned)16777216 % ((unsigned)NumBins - 1) ? 1 : 0;
        printf( "Divisor is: %u\n", divisor );
    }
    // now we compute the number represented by the first 3
    // characters of the key, and divide it by divisor, the
    // integral part gives the bin number.
    for( i = 0; i < 3; i++ ) {
        keyval <<= 8;
        keyval += key[i];
    }
//    printf( "Key is %.2x %.2x %.2x \n", key[0],key[1], key[2] );
    
//    printf( "Keyval is: %u\n", keyval );
//    printf( "Bin is: %u\n", keyval / divisor );
	if( keyval == 0 )
		return 0;
    return ( keyval / divisor ) + 1;
}



static int ComputeAsciiBinNumber( const unsigned char key[10] ) {
    int i;
    unsigned int keyval = 0;
    // how many bits can we use for the binning number?
    // first time through, compute divisor
    if( divisor == 0 ) {
        // strange constant is 95^4 -- the max possible value
        // of the first five key characters + 1.
        divisor = 81450625 / NumBins;
        // need ceiling of this
        divisor += 81450625 % NumBins ? 1 : 0;
    }
    // now we compute the number represented by the first 4
    // characters of the key, and divide it by divisor, the
    // integral part gives the bin number.
    
    for( i = 0; i < 4; i++ ) {
        keyval *= 95;
        keyval += key[i] - ' ';
    }
    return( keyval / divisor );
}



static offset_t         consumed_recs = 0;
static buffer_t         *toprocess = NULL;

static void DoBinning( void ) {
    //printf( "---------------  do binning -------------\n" );
    while( 1 ) {
        int bin = -1;
        // if we don't already have a buffer to process, see if we
        // can get one.
        if( toprocess == NULL ) {
            //printf( "toprocess == null -- no currently processing buffer\n" );
            if( ToProcess.nitems ) {
                //printf( "getting one off ToProcess list\n" );
                toprocess = PopHead( &(ToProcess) );
                consumed_recs = 0;
            } else {
                // we can't get anything to process
                //printf( "nothing to process\n" );
                return;
            }
        }
        //printf( "processing records in current toprocess buffer\n" );
        // try to process all the records in the toprocess buffer.
        //printf( "for( ; consumed_recs (%d) < toprocess->numrecs (%d); ... ) {\n", consumed_recs, toprocess->numrecs );
        for( ; consumed_recs < toprocess->numrecs; consumed_recs++, RecsProcessed++ ) {
            
            buffer_t *headbuf;
            record_t *rec = &(toprocess->recs[consumed_recs]);
            
            // find what bin this next record belongs in.
#ifdef ASCII_KEYBYTES
            bin = ComputeAsciiBinNumber( rec->key );
#else
#ifdef NNNNN_KEYBYTES
			bin = ComputeNNNNNBinNumber( rec->key );
#else
            bin = ComputeBinNumber( rec->key );
#endif
#endif
            if( (bin >= NumBins) || (bin < 0) ) {
                printf( "error: invalid bin from ComputeBinNumber: %d\n", bin );
            }

            //printf( "record bound for bin %d\n", bin );

            // now, let's see what the situation is with that bin and its
            // buffers.  In particular, do we have a spot to put this record?
            headbuf = Bins[bin].bufs.head;
            // if we have a buffer, and the buffer is full or executing or
            // if there's no buffer at all, let's try to get one
            if( !headbuf || 
                headbuf->numrecs == headbuf->totalrecs || 
                headbuf->operation != OP_NONE ) {
                //printf( "headbuf busy or full -- op: %d, numrecs: %d, totalrecs: %d\n",
                    //headbuf->operation, headbuf->numrecs, headbuf->totalrecs );
                // first see if this is our 'special' buffer and if we can use it
                if( headbuf->operation == BIN_SPECIAL ) {
                    //printf( "headbuf is only one left and finished so reclaiming for use\n" );
                    headbuf->numrecs = 0;
                    headbuf->operation = OP_NONE;
                } else {
                    //printf( "trying to get buffer from free list\n" );
                    if( Free.nitems ) {
                        //printf( "got one from freelist\n" );
                        PushHead( &(Bins[bin].bufs), AllocateFree() );
                        headbuf = Bins[bin].bufs.head;
                    } else {
//                        printf( "no free buffers to use for bin -- binning BLOCKS!\n" );
                        return;
                    }
                }
            }
            // now headbuf must exist, and it must be non-full so we can
            // add our item.
            headbuf->recs[headbuf->numrecs++] = *rec;
            Bins[bin].nrecs++;
            //printf( "added rec to bin\n" );
            // if we made it full, write the thing
            if( headbuf->numrecs >= headbuf->totalrecs ) {
                //printf( "writing bin buffer because full\n" );
                headbuf->file = Bins[bin].file;
                headbuf->device = &(Devices[Bins[bin].dev].dev);
                RecsCommitted += headbuf->numrecs;
#ifdef NO_BIN_WRITE_PERF_TEST
				// just put it in the finished state
				headbuf->operation = OP_FINISHED;
#else
                WriteBuffer( headbuf, headbuf->numrecs, headbuf->device );
#endif
                headbuf = NULL;
            }
            
        }
        
        // if we hit the end of this buffer,
        // put it back on the free list, and start the loop over
        if( consumed_recs >= toprocess->numrecs ) {
            //printf( "finished with this block\n" );
            PushTail( &Free, toprocess );
            toprocess = NULL;
        }

        //printf( "going back for more\n" );
        
    }
}





void FinishBinning() {
    int i;
    buffer_t *b;
    offset_t recs = 0;
    // be sure to finish off the write operations.
    for( i = 0; i < NumBins; i++ ) {
        //printf( "bin: %d, nrecs: %d, operation: %d\n", i, Bins[i].nrecs, Bins[i].operation );
        while( Bins[i].bufs.nitems ) {
            // walk through the buffers, and if they haven't been executed,
            // execute them.
            b = PopHead( &(Bins[i].bufs) );
            if( b->operation == OP_NONE && b->numrecs ) {
                recs += b->numrecs;
                b->file = Bins[i].file;
                b->device = &(Devices[Bins[i].dev].dev);
#ifdef NO_BIN_WRITE_PERF_TEST
				// just put it in the finished state
				b->operation = OP_FINISHED;
#else
                WriteBuffer( b, b->numrecs, b->device );
#endif
            }
        }
    }
    RecsCommitted += recs;
}



offset_t CalculateDataReadSize( buffer_t* b ){
// commented version is for traditional dmsort
//	return MIN(b->totalrecs, RecsUnread) * sizeof( record_t );
	return MIN(b->totalrecs + mask_length - 1, RecsUnread + mask_length - 1 );
}

static void DoReading( void ) {
    buffer_t * b;
    //printf( "do reading\n" );
    if( RecsUnread && Free.nitems ) {
        // allocate a buffer
        b = AllocateFree();
        
        // start reading into it.
        b->file = Data;
        ReadBuffer( b, MIN(b->totalrecs, RecsUnread), &(Devices[DataDev].dev) );

        b->input_pos = NumRecs - RecsUnread;
       	// need to step back mask_length - 1 characters to get the complete sequence!!
//		if( b->input_pos >= mask_length - 1 )
//			b->input_pos -= mask_length - 1;
        b->io_pos = b->input_pos;
//		printf( "Reading offset %llu\n", b->io_pos );
        b->io_size = CalculateDataReadSize( b );
        // decrement recsunread appropriately
        RecsUnread -= MIN(MIN(b->totalrecs,RecsUnread),RecsUnread);
        
        // put the thing on the Reading list.
        //printf( "new buffer on reading list\n" );
        PushTail( &Reading, b );
    }
}





static void HandleBinWriteCompletions( void ) {
    int i;
    buffer_t *b, *tmpnext;
    //printf( "handle bin write completions\n" );
    for( i = 0; i < NumBins; i++ ) {
        b = Bins[i].bufs.head;
        do {
            if( !b ) {
                break;
            }
            tmpnext = b->next;
            if( b->operation == OP_FINISHED ) {
                RecsWritten += b->numrecs;
                if( Bins[i].bufs.nitems > 1 ) {
                    b->operation = OP_NONE;
                    PushHead( &Free, RemoveItem( &(Bins[i].bufs), b ) );
                } else {
                    b->operation = BIN_SPECIAL;
                }
            }
            b = tmpnext;
        } while( b != Bins[i].bufs.head && Bins[i].bufs.nitems > 1 );
    }
}

static void HandleSeqbufWriteCompletions( void ) {
    buffer_t *b, *tmpnext;
    //printf( "handle bin write completions\n" );
    b = Seqbuf.bufs.head;
    do {
        if( !b ) {
            break;
        }
        tmpnext = b->next;
        if( b->operation == OP_FINISHED ) {
            if( Seqbuf.bufs.nitems > 1 ) {
                b->operation = OP_NONE;
                PushHead( &Free, RemoveItem( &(Seqbuf.bufs), b ) );
            } 
        }
        b = tmpnext;
    } while( b != Seqbuf.bufs.head && Seqbuf.bufs.nitems > 1 );
}

#define ALPHA_BITS 2

static void Translate32(uint32* dest, const char* src, const unsigned len){
	uint8 start_bit = 0;
	unsigned cur_word = 0;
	uint32 word_mer = 0;
	uint32 i = 0;
	if( len == 0 )
		return;
	for(i=0; i < len; i++){
//		uint32 tmp = DNA_TABLE[src[i]];
		if(start_bit + ALPHA_BITS <= 32){
			word_mer <<= ALPHA_BITS;
			word_mer |= DNA_TABLE[src[i]];
			dest[cur_word] = word_mer;
			start_bit += ALPHA_BITS;
			if(start_bit >= 32 && i < len - 1){
				word_mer = 0;
				start_bit %= 32;
				cur_word++;
			}
		}else{
			printf("Error, this should never happen with DNA sequence\n" );
/*			uint8 over_bits = (start_bit + ALPHA_BITS) % 32;
			uint32 tmp2 = tmp;
			tmp2 <<= 32 - over_bits;
			tmp >>= over_bits;
			dest[cur_word] |= tmp;
			cur_word++;
			dest[cur_word] = 0;
			dest[cur_word] |= tmp2;
			start_bit = over_bits;
*/		}
	}
	if( start_bit != 0 ){
		dest[cur_word] <<= 32 - start_bit;
	}
}


void RestructureReadSMLBins( void ) {
	char little_endian = 1;
	mask_t bit;
	mask_t mer, rc_mer;
	record_t forward, reverse;
	record_t begin[6];	// the first six records could potentially overwrite the sequence
	int i;
	offset_t seqI, extras, weight;
	char* sequence;
	sml_t *sml;

    buffer_t *b, *tmpnext;
	
	// variables for translation to 2-bit per base
    buffer_t *headbuf;
	int seq_bit;
	int seq_word;
	int word_remainder;
	offset_t translate_length;
	int config_value = 4554307;
//	int seq_offset;
    // go through and see if any have completed.
    b = Restructure.head;
    do {
        if( !b ) {
            break;
        }
		// is this the buffer we need to translate next?
		if( b->input_pos != Seqbuf.seq_pos ){
			b = b->next;
			continue;
		}
		
        tmpnext = b->next;
		sequence = (char *)b->recs;
		sml = (sml_t*)b->recs;

		// translate the sequence that was just read and write it out
        headbuf = Seqbuf.bufs.head;
        // if we have a buffer, and the buffer is full or executing or
        // if there's no buffer at all, let's try to get one
        if( !headbuf || 
            headbuf->operation != OP_NONE ) {
            //printf( "headbuf busy or full -- op: %d, numrecs: %d, totalrecs: %d\n",
                //headbuf->operation, headbuf->numrecs, headbuf->totalrecs );
            // first see if this is our 'special' buffer and if we can use it
            if( headbuf->operation == OP_FINISHED ) {
                //printf( "headbuf is only one left and finished so reclaiming for use\n" );
                headbuf->numrecs = 0;
                headbuf->operation = OP_NONE;
	            Seqbuf.bufpos = 0;
            } else {
                //printf( "trying to get buffer from free list\n" );
                if( Free.nitems ) {
//                    printf( "got one from freelist\n" );
                    PushHead( &(Seqbuf.bufs), AllocateFree() );
                    headbuf = Seqbuf.bufs.head;
		            Seqbuf.bufpos = 0;
                } else {
//                    printf( "no free buffers to use for Seqbuf -- restructuring BLOCKS!\n" );
                    return;
                }
            }
        }

		seq_bit = Seqbuf.bufpos * 2;
		seq_word = seq_bit / 32;
		word_remainder = seq_bit % 32;
		if( word_remainder != 0 ){
			seq_word++;
		}
		
//			int end_bit = 2 * (Seqbuf->bufpos + b->io_size - mask_length + 1);
//			int end_remainder = end_bit % 32;
		translate_length = b->io_size - mask_length + 1 - (word_remainder / 2);
		if( b->io_size + b->input_pos >= NumRecs ){
			// this is the last I/O, translate the whole thing
			translate_length += mask_length - 1;
		}
//			translate_length -= end_remainder / 2;
		
		// The number of bytes in headbuf->recs must ALWAYS be divisible by 4 when using
		// Translate32, otherwise corruption will result
#ifndef NO_RESTRUCTURE_PERF_TEST
		Translate32( (uint32*)(headbuf->recs) + seq_word, ((char*)b->recs) + (word_remainder / 2), translate_length );
#endif
		
		// need to fill in beginning
		if( word_remainder != 0 ){
			int begin_mer = 0;
			for( seqI = 0; seqI < word_remainder / 2; seqI++ ){
				begin_mer <<= 2;
				begin_mer |= DNA_TABLE[ sequence[ seqI ] ];
			}
//			((uint32*)headbuf->recs)[ seq_word - 1 ] <<= 32 - word_remainder;
			((uint32*)headbuf->recs)[ seq_word - 1 ] |= begin_mer;
		}
		
		Seqbuf.bufpos += translate_length + (word_remainder / 2);
		Seqbuf.seq_pos += translate_length + (word_remainder / 2);

        // if we made it full, write the thing
        // each buf will consume headbuf->totalrecs / 4 bytes.
        // there are headbuf->totalrecs * sizeof( record_t ) bytes available in the Seqbuf.
        // thus we can fit 4 * sizeof( record_t ) bufs in each Seqbuf
        if( Seqbuf.bufpos == headbuf->totalrecs * sizeof( record_t ) * 4 ||
			b->io_size + b->input_pos >= NumRecs ) {
            //printf( "writing bin buffer because full\n" );
            headbuf->file = Seqbuf.file;
            headbuf->device = &(Devices[Seqbuf.dev].dev);
            WriteBuffer( headbuf, headbuf->totalrecs, headbuf->device );
            headbuf->io_size = Seqbuf.bufpos / 4;
            if( b->io_size + b->input_pos >= NumRecs ){
            	offset_t offI = 0;
            	offI = headbuf->io_size % 4;
            	if( offI != 0 )
	            	headbuf->io_size += 4 - offI;
            	for( offI = 0; offI < 8; offI++ )
            		((char*)headbuf->recs)[ headbuf->io_size + offI ] = 0;
            	headbuf->io_size += 8;
            }
            headbuf = NULL;
        }else if( Seqbuf.bufpos > headbuf->totalrecs * sizeof( record_t ) * 4 ){
        	printf( "Error.  Over filled Seqbuf\n" );
        }


		// translate the sequence according to the current sequence mask
#ifndef NO_RESTRUCTURE_PERF_TEST
		for( seqI = b->io_size - mask_length + 1; seqI > 0; seqI-- ){
			bit = 1;
			bit <<= mask_length - 1;
			mer = 0;
			weight = 0;
			for( i = 0; i < mask_length; i++ ){
				if( bit & seed_mask ){
					mer <<= 2;
					mer |= DNA_TABLE[ sequence[ seqI + i - 1 ] ];
				}
				bit >>= 1;
			}
			// copy the mer from the 64-bit integer based on the endian-ness of the system
			// copy mer to forward key
			mer <<= 64 - (2 * mask_weight);
//			if( seqI + b->input_pos == config_value )
//				__asm( nop );
			if( little_endian ){
				for( i = 0; i < MASK_T_BYTES; i++ )
					forward.key[i] = ((char*)(&mer))[ sizeof( mer ) - i - 1 ];

			}else{
				for( i = 0; i < MASK_T_BYTES; i++ )
					forward.key[i] = ((char*)(&mer))[ i ];
			}

			// reverse complement the mer
			mer = ~mer;
			for( i = 0; i < 64; i += 2 ){
				rc_mer <<= 2;
				rc_mer |= mer & 3;
				mer >>= 2;
			}
			rc_mer <<= 64 - (2 * mask_weight);
			// copy mer to reverse key
			if( little_endian ){
				for( i = 0; i < MASK_T_BYTES; i++ )
					reverse.key[i] = ((char*)(&rc_mer))[ sizeof( mer ) - i - 1 ];
			}else{
				for( i = 0; i < MASK_T_BYTES; i++ )
					reverse.key[i] = (((char*)(&rc_mer))[i]);
			}
			// put the lesser key in forward
			if( COMPARE_KEYS( forward, reverse ) > 0)
				forward = reverse;
			
			// watch out for the last 6 records
			if( seqI <= 6 ){
				begin[ seqI - 1] = forward;
			}else{
				b->recs[ seqI - 1 ] = forward;
				// set the position
				sml[ seqI - 1 ].pos = b->input_pos + seqI - 1;
			}
		}

		extras = b->io_size - mask_length + 1 < 6 ? b->io_size - mask_length + 1 : 6;
		
		// fill in the first six records
		for(; seqI < extras; seqI++ ){
			b->recs[ seqI ] = begin[ seqI ];
			// set the position
			sml[ seqI ].pos = b->input_pos + seqI;
		}
#else
	if(1){	// define a new scope so the variables can be local
	// simulate random data in each bin
    int i;
    unsigned int keyval = 0;
	unsigned int tmpval = 0;
    if( divisor == 0 ) {
        divisor = (unsigned)16777216 / (unsigned)NumBins;
        // need ceiling of this
        divisor += (unsigned)16777216 % (unsigned)NumBins ? 1 : 0;
        printf( "Divisor is: %u\n", divisor );
    }
	for( seqI = 0; seqI < b->numrecs; seqI++ ){
		tmpval = keyval;
	    for( i = 3; i > 0; i-- ) {
			b->recs[ seqI ].key[ i - 1 ] = (tmpval & 0xFF);
			b->recs[ seqI ].key[ i - 1 ] = 0;
			tmpval >>= 8;
		}
		keyval += divisor;
	}
	}
#endif
		
		// b has been restructured, add it to the ToProcess list
        PushTail( &ToProcess, RemoveItem( &Restructure, b ) );

		
        b = tmpnext;
    } while( b != Restructure.head && Restructure.nitems );
}

static void HandleReadingCompletions( void ) {
    buffer_t *b, *tmpnext;
    // just go through and see if any have completed.
    b = Reading.head;
    do {
        if( !b ) {
            break;
        }
        tmpnext = b->next;
        if( b->operation == OP_FINISHED ) {
            // migrate this to the toprocess list
            b->operation = OP_NONE;
            PushTail( &Restructure, RemoveItem( &Reading, b ) );
            // bookkeeping
            RecsRead += b->numrecs;
        }
        b = tmpnext;
    } while( b != Reading.head && Reading.nitems );
}



void print_usage( const char* pname ){
	printf( "Usage: %s <-m Working set size in MB> <-b buffer size> <-i input file> <-o output file> [-n number of records] <bin directory> <num bins> ... [bin directory] [num bins]\n", pname );
}

int InitdmSML( long working_mb, long buffer_size, const char* input_filename, const char* output_filename, const char* const* scratch_paths, uint64 seed ) {
    int i, j;
    offset_t desired_ws_size, actual_ws_size;
    SMLHeader_t header;
    struct {
        const char * bin_dev;
        int devnum;
        int nbins;
    } bins[8];

	char *bin_name;
	int scratchI = 0;

    // initialize the timing stuff
    InitTime();

    // start the running timer now.
    RunningTime = 0;
    RunningTimer = StartTimer();

	if( working_mb != 0 ){

	desired_ws_size = working_mb;
	desired_ws_size *= 1024 * 1024;	// convert to bytes

	}else{
	// set desired working set size to half of physical memory...
#ifdef WIN32
	{
/*	MEMORYSTATUSEX ms;
	memset( &ms, 0, sizeof( MEMORYSTATUSEX ) );
	GlobalMemoryStatusEx( &ms );
	desired_ws_size = ms.ullTotalPhys / 2;
*/
	MEMORYSTATUS ms;
	memset( &ms, 0, sizeof( MEMORYSTATUS ) );
	GlobalMemoryStatus( &ms );
	desired_ws_size = ms.dwTotalPhys / 2;
	}
#else

    {
	// get it from /proc/meminfo
	FILE *fp = fopen("/proc/meminfo", "r");
	if ( fp )
	{
		long memTotal;

		char buf[1024];
		if ( fgets(buf, sizeof(buf), fp) )
		{
			sscanf(buf, "MemTotal: %ld kB", &memTotal);
			fprintf( stderr, buf );
		}
		fclose(fp);
		// allocate about 6/10 of physical memory
		// leave the rest for buffer cache
		desired_ws_size = memTotal * 512;
	}
	}

#endif
	// never allocate more than 2GB
	if( desired_ws_size / 1024  > 2048 * 1024 ){
		desired_ws_size = 1024 * 1024;
		desired_ws_size *= 2048;
	}
//	desired_ws_size /= sizeof( record_t ); // get working set size in records
	}
	
	if( buffer_size == 0 ){
		buffer_size = 1;
		while( desired_ws_size / (buffer_size*sizeof(record_t)) > 2048 ){
			buffer_size *= 2;
		}
	}

	BufferSizeMin = BufferSizeMax = buffer_size;
	OutFileName = output_filename;
	
	// find out how many scratch paths were given before the null terminator
	for( ; ; scratchI++ ){
		if( !scratch_paths || scratch_paths[ scratchI ] == NULL )
			break;
	}
	
    

	NumBinDevs = scratchI;
	NumDevices = 2 + NumBinDevs;
	Devices = (device_t*)malloc( NumDevices * sizeof(device_t) );
	DataDev = 0;
	OutputDev = 1;
	Devices[DataDev].devname = "Input device";
	Devices[DataDev].path = input_filename;
	Devices[DataDev].dev.buf = NULL;
	Devices[OutputDev].devname = "Output device";
	Devices[OutputDev].path = OutFileName;
	Devices[OutputDev].dev.buf = NULL;
    
    
    if( NumBinDevs == 0 ) {
    	return TOO_FEW_BINS;
    } else if( NumBinDevs > 8 ) {
    	return TOO_MANY_BINS;
    }
	
	NumRecs = aStatFileSize( input_filename );

	// calculate number of bins using nrecs and ws_size
	NumBins = desired_ws_size / (200 * NumBinDevs);
	NumBins = NumRecs / NumBins;
	NumBins = NumBins < 5 * NumBinDevs ? 5 * NumBinDevs : NumBins;	// don't allow fewer than 5 bins per dev
	// round for equal number of bins per dev
	if( NumBins % NumBinDevs != 0 )
		NumBins = ( (NumBins / NumBinDevs) + 1 ) * NumBinDevs;
	printf( "Creating %d bin files\n", NumBins );
	for( i = 2; i < NumDevices; i++ ){
		bin_name = (char*)malloc( 10 );
		strcpy( bin_name, "bin dev__" );
		bin_name[8] = 0x40 + i - 2;
		Devices[i].devname = bin_name;
		Devices[i].path = scratch_paths[ i - 2 ];
		Devices[i].dev.buf = NULL;
		bins[i - 2].bin_dev = bin_name;
		bins[i - 2].nbins = NumBins / NumBinDevs;	// allocate even an portion of bins per device
		bins[i - 2].devnum = i;
	}
	
    // get buffer size.
    if( BufferSizeMin == 0 ) {
        BufferSizeMin = MINRECS;
        BufferSizeMax = MAXRECS;
    }


    // open the input file
    Data = aOpen( input_filename, A_READ );
	if( Data == NULL ) {
	        printf( "couldn't open data file\n" );
		return INPUT_NOT_OPENED;
	}
   
    // get working set size
    if( desired_ws_size == 0 ) {
        printf( "invalid working set size (%llu) -- must be at least 0\n", desired_ws_size );
    	return INVALID_WS_SIZE;
    }
	
	// init translation table
	DNA_TABLE = CreateBasicDNATable();

    // open the output file
    Output = aOpen( OutFileName, A_WRITE );
    if( !Output ) {
        printf( "couldn't open output file!\n" );
    	return OUTPUT_NOT_OPENED;
    }
	
	header = InitSML( Output, NumRecs, seed );
	seed_mask = header.seed;
	mask_length = header.seed_length;
	mask_weight = header.seed_weight;
	
	if( NumRecs <= mask_length - 1 ){
	        printf( "Sequence must be at least %d characters in length\n", mask_length );
		return SEQUENCE_TOO_SHORT;
	}

	NumRecs -= mask_length - 1;
	printf( "NumRecs is: %llu \n", NumRecs );
    // get the number of records we should process
    RecsProcessed = 0;
    RecsUnread = NumRecs;
    if( NumRecs <= 0 ) {
    	return INVALID_NUMRECS;
        printf( "invalid NumRecs: %llu\n", NumRecs );
    }
    
    
    
    // go ahead and create the working set.
    actual_ws_size = MakeWorkingSet( &WS, desired_ws_size, BufferSizeMin, BufferSizeMax );
    printf( "desired working set: %llu, actual working set: %llu\n", 
        desired_ws_size, actual_ws_size );

    // initialize the Free list -- just put all the buffers on it.
    for( i = 0; i < WS.nbufs; i++ ) {
        PushHead( &Free, &(WS.bufs[i]) );
    }

    printf( "working set size        : %llu\n", actual_ws_size );
    printf( "total buffers           : %d\n", WS.nbufs );
    // FIXME: can any touching of the memory here help us?
    // toprocess and reading list empty to start
    ToProcess.nitems = Reading.nitems = 0;
    ToProcess.head = Reading.head = NULL;
	Restructure.nitems = 0;
	Restructure.head = NULL;
		
	// allocate Seqbuf
	Seqbuf.file = Output;
	Seqbuf.dev = OutputDev;
	Seqbuf.bufpos = 0;
	Seqbuf.seq_pos = 0;
    if( Free.nitems ) {
        PushHead( &(Seqbuf.bufs), AllocateFree() );
    } else {
        printf( "error: could not give a buffer to Seqbuf\n" );
        return NO_FREE_BUFFERS;
    }

    // allocate the bins.
    Bins = malloc( sizeof( *Bins ) * NumBins );
    memset( Bins, 0, sizeof( *Bins ) * NumBins );

    // allocate the bins in a round-robin fashion, so when we read
    // things back for sorting, we're not swamping one device at a time --
    // instead, things are spread out.
    printf( "opening %d bins\n", NumBins );
    j = -1;
    for( i = 0; i < NumBins; i++ ) {
        // find a bin on the next device.
        while( 1 ) {
            j = (j+1) % NumBinDevs;
            if( bins[j].nbins ) {
                // make this bin on that device, and
                // round-robin switch to the next device.
                const char *fname = Fmt("%sout%05d.binned",Devices[bins[j].devnum].path,i);
                Bins[i].dev = bins[j].devnum;
                Bins[i].fname = malloc( strlen( fname ) + 1 );
                strcpy( Bins[i].fname, fname );

#ifndef NO_BINNING_PERF_TEST
                Bins[i].file = aOpen( fname, A_WRITE );
                //printf( "opened '%s' on device '%s'\n", fname, Devices[bins[j].devnum].devname );
                if( Bins[i].file == NULL ) {
                    printf( "couldn't open output bin file '%s'\n", fname );
					return BIN_NOT_OPENED;
                }
#else
                Bins[i].nrecs = aStatSize( fname );
		if( Bins[i].nrecs == 0 ){
			// just make sure the file exists
	                Bins[i].file = aOpen( fname, A_WRITE );
			aClose( Bins[i].file );
			Bins[i].file = NULL;
		}
#endif // NO_BINNING_PERF_TEST
                bins[j].nbins--;
                break;
            }
        }
    }

    // now we allocate one buffer for each bin
    // and each bin will hold onto at least one buffer
    // so that we can guarantee no locking cases
    for( i = 0; i < NumBins; i++ ) {
        if( Free.nitems ) {
            PushHead( &(Bins[i].bufs), AllocateFree() );
        } else {
            printf( "error: could not give one buffer to each bin\n" );
	        return NO_FREE_BUFFERS;
        }
    }
	
	// all went well
	return 0;
}



void DisplayStatusHeader( void ) {
    printf( "time recs_read recs_processed recs_committed recs_written binning_rate free reading toprocess bins restructure\n" );
}


void DisplayStatus( void ) {

    printf( "%f %llu %llu %llu %llu %f %d %d %d %d %d\n",
        RunningTime, RecsRead, RecsProcessed, RecsCommitted, RecsWritten, 
        RecsProcessed/RunningTime, Free.nitems, Reading.nitems, ToProcess.nitems, 
        WS.nbufs - Free.nitems - Reading.nitems - ToProcess.nitems - Restructure.nitems, Restructure.nitems );

    /*
    int i;
    printf( "-----------------------------------------------------------\n" );
    printf( "Records Processed : %d/%d\n", RecsProcessed, NumRecs );
    printf( "Records Committed : %d\n", RecsCommitted );
    printf( "Records Written   : %d\n", RecsWritten );
    printf( "Records Read      : %d\n", RecsRead );
    printf( "Running Time      : %f seconds\n", RunningTime );
    printf( "Binning Rate      : %f records/sec  (%f bytes/sec)\n",
        RecsProcessed / RunningTime, RecsProcessed * sizeof(record_t) / RunningTime );
    printf( "Freelist entries  : %d\n", Free.nitems );
    printf( "Reading entries   : %d\n", Reading.nitems );
    printf( "ToProcess entries : %d\n", ToProcess.nitems );
    printf( "Bin entries:\n" );
    for( i = 0; i < NumBins; i++ ) {
        printf( "  %4d : %4d\n", i, Bins[i].bufs.nitems );
    }    
    printf( "Device status:\n" );
    for( i = 0; i < NumDevices; i++ ) {
        printf( "  %d : '%16s' : '%16s' : %s\n", i, Devices[i].devname,
            Devices[i].path, Devices[i].dev.state == DEV_FREE ? "FREE" : "BUSY" );
    }
    */
}


void UpdateIOState( void ) {
    int i;
    //printf( "update io state\n" );
    
    // first update aio ops on the data file
    aUpdateOperations( Data );
    // next update aio ops on the bin files
    for( i = 0; i < NumBins; i++ ) {
        aUpdateOperations( Bins[i].file );
    }
    // update aio ops on the output file
    aUpdateOperations( Output );
    // next, let the working set adjust operation states and such
    UpdateWSIOFinishedState( &WS );
    // finally, let the devices start new operations if possible.
    for( i = 0; i < NumDevices; i++ ) {
        UpdateDeviceIOExecuteState( &WS, &(Devices[i].dev) );
    }
    
}


void EnsureAllOperationsComplete( void ) {
    int i;
    int not_complete = 1;
    dmtimer_t *wait;
    wait = StartTimer();
    while( not_complete ) {
        UpdateIOState();
        // see if we're done
        not_complete = 0;
        for( i = 0; i < WS.nbufs; i++ ) {
            if( WS.bufs[i].device &&
                WS.bufs[i].file &&
                (WS.bufs[i].operation == OP_PENDING || WS.bufs[i].operation > OP_NONE) ) {
                not_complete = 1;
                break;
            }
        }
    }
    printf( "Ensure All Operations Complete: %d msec\n", ReadTimer( wait ) );
    StopTimer( wait );
}





static double lasttime = 0;

void BinningPhase( void ) {

    int i;
    // for progress output
    int iter;
    int timeaccum;

    // the main loop.
    printf( "----------------- Starting -----------------\n" );
    printf( "working set buffers : %d\n", WS.nbufs );
    printf( "number of bins      : %d\n", NumBins );
    timeaccum = 0;
    iter = 0;
    DisplayStatusHeader();
    while( RecsProcessed < NumRecs ) {

        // print status every few seconds or so.
        // not until timing gets fixed
        //if( RunningTime - lasttime >= 5.0f ) {
        if( (RunningTime - lasttime) >= 2.0f ) {
            DisplayStatus();
            lasttime = RunningTime;
        }
        
        // keep the async io running
        // first update the operations on all our files.
        UpdateIOState();
        
        // Handle read and write completions
        // (transition reads to ToProcess, writes to Free)
        HandleReadingCompletions();
		HandleSeqbufWriteCompletions();
		RestructureReadSMLBins();
        HandleBinWriteCompletions();

        // do reading and binning
        DoReading();
        DoBinning();
        
        // finish up the loop.
        iter++;

        RunningTime = (double)ReadTimer( RunningTimer ) / 1000.0;

    }
    
    printf( "total iters: %d\n", iter );
    // now, we *must* take care to make sure all writes have completed
    // We can't simply call aClose on a file.  It's true that that will
    // wait until all the currently scheduled operations on that file
    // complete, but with the device method, we only allow one operation
    // on any device at a time.  Thus, we must ask the device managers
    // to complete their own IO.
    // FIXME: this could potentially be moved into the buffer stuff for
    // a DeviceClose type of call, but then if there is lots of stuff
    // pending, unless DeviceClose could know about more than one device
    // at a time, we would get effectively synchronous IO here, so we
    // have the ugly hack for now.
    FinishBinning();
    EnsureAllOperationsComplete();

    // close the input file.
    aClose( Data );
    Data = NULL;
    // Finally, close all the bin files
    for( i = 0; i < NumBins; i++ ) {
        aClose( Bins[i].file );
        Bins[i].file = NULL;
    }
    printf( "Finally, RecsCommitted: %llu\n", RecsCommitted );

    DisplayStatus();

}





void SortReading( void ) {

    int i;

    // if anything is in WAIT_READ, and we have crap to read yet,
    // start reading it in.

    for( i = 0; i < NSortBufs; i++ ) {
        // quick out if we're done reading.
        if( BinToRead >= NumBins ) {
            return;
        }
        if( SortBufs[i].state == WAIT_READ ) {
            // schedule a read here.
            const char *fname = Fmt("%sout%05d.binned",Devices[Bins[BinToRead].dev].path,BinToRead);
            aFILE *in = aOpen( fname, A_READ );
            if( !in ) {
                printf( "couldn't open '%s' to read!\n", fname );
            }
            if( Bins[BinToRead].nrecs > SortBufs[i].buf->totalrecs ) {
                printf( "buffer not big enough to hold bin!\n" );
            }
            SortBufs[i].bin = BinToRead;
            SortBufs[i].dev = &(Devices[Bins[BinToRead].dev].dev);
            SortBufs[i].state = BUSY_READ;
            SortBufs[i].buf->file = in;
            ReadBuffer( SortBufs[i].buf, Bins[BinToRead].nrecs, SortBufs[i].dev );
            printf( "scheduled read of bin %d\n", BinToRead );
            BinToRead++;
            return;
        }
    }

}



#ifdef USE_QSORT_ONLY

int comp_keys( record_t a, record_t b ){
	int compval;
	sml_t *mer_a, *mer_b;
	mer_a = (sml_t*)&a;
	mer_b = (sml_t*)&b;
/*	if( ( mer_a->pos == 4554307 &&
		mer_b->pos == 4407600 ) ||
		( mer_a->pos == 4407600 &&
		mer_b->pos == 4554307 ) )
		__asm( nop );
*/	compval = COMPARE_KEYS( a, b );
	return compval;
}

void QBrute( record_t a[], int lo, int hi ) {
    if ((hi-lo) == 1) {
        if( comp_keys( a[hi], a[lo] ) < 0 ) {
            record_t T = a[lo];
            a[lo] = a[hi];
            a[hi] = T;
        }
    }
    if ((hi-lo) == 2) {
        int pmin = comp_keys( a[lo], a[lo+1] ) < 0 ? lo : lo+1;
        pmin = comp_keys( a[pmin], a[lo+2] ) < 0 ? pmin : lo+2;
        if (pmin != lo) {
            record_t T = a[lo];
            a[lo] = a[pmin];
            a[pmin] = T;
        }
        QBrute(a, lo+1, hi);
    }
    if ((hi-lo) == 3) {
        int pmin, pmax;
        pmin = comp_keys( a[lo], a[lo+1] ) < 0 ? lo : lo+1;
        pmin = comp_keys( a[pmin], a[lo+2] ) < 0 ? pmin : lo+2;
        pmin = comp_keys( a[pmin], a[lo+3] ) < 0 ? pmin : lo+3;
        if (pmin != lo) {
            record_t T = a[lo];
            a[lo] = a[pmin];
            a[pmin] = T;
        }
        pmax = comp_keys( a[hi], a[hi-1] ) > 0 ? hi : hi-1;
        pmax = comp_keys( a[pmax], a[hi-2] ) > 0 ? pmax : hi-2;
        if (pmax != hi) {
            record_t T = a[hi];
            a[hi] = a[pmax];
            a[pmax] = T;
        }
        QBrute(a, lo+1, hi-1);
    }
}



void QSort( record_t a[], int lo0, int hi0 ) {
    
    int lo = lo0;
    int hi = hi0;
    
    record_t pivot;

    if ((hi-lo) <= 3) {
        QBrute(a, lo, hi);
        return;
    }
    
    // Pick a pivot and move it out of the way
    pivot = a[(lo + hi) / 2];
    a[(lo + hi) / 2] = a[hi];
    a[hi] = pivot;
    
    while( lo < hi ) {

    // Search forward from a[lo] until an element is found that
    // is greater than the pivot or lo >= hi 
        //while( a[lo] <= pivot && lo < hi ) {
        while( (comp_keys( a[lo], pivot ) <= 0) && lo < hi ) {
            lo++;
        }
        
        //
        //  Search backward from a[hi] until element is found that
        //  is less than the pivot, or hi <= lo 
        //
        //while (pivot <= a[hi] && lo < hi ) {
        while( (comp_keys( pivot, a[hi] ) <= 0) && lo < hi ) {
            hi--;
        }
        
        //
        //  Swap elements a[lo] and a[hi]
        //
        if( lo < hi ) {
            record_t T = a[lo];
            a[lo] = a[hi];
            a[hi] = T;
        }
    }
    
    //
    // Put the median in the "center" of the list
    //
    a[hi0] = a[hi];
    a[hi] = pivot;
    
    //
    // Recursive calls, elements a[lo0] to a[lo-1] are less than or
    // equal to pivot, elements a[hi+1] to a[hi0] are greater than
    // pivot.
    //
    QSort( a, lo0, lo-1 );
    QSort( a, hi+1, hi0 );
}





void RecSort( record_t a[], int nelems ) {

    QSort( a, 0, nelems-1 );

}


int SortBuffer( buffer_t * buf ) {

    RecSort( buf->recs, buf->numrecs );
    return( 1 );

}


void SortSorting( void ) {

    int i, finished;
    int lowest = -1;
    QSortTimer = StartTimer();

    for( i = 0; i < NSortBufs; i++ ) {
        if( SortBufs[i].state == SORTING ) {
            if( lowest == -1 || SortBufs[i].bin < SortBufs[lowest].bin ) {
                lowest = i;
            }
        }
    }

    if( lowest != -1 ) {
        printf( "sorting bin %d\n", SortBufs[lowest].bin );
        finished = SortBuffer( SortBufs[lowest].buf );
        if( finished ) {
            SortBufs[lowest].state = WRITE_RESTRUCTURE;
//            SortBufs[lowest].state = WAIT_WRITE;
        }
    }

    QSortTime += ReadTimer( QSortTimer ) / 1000.0;
    StopTimer( QSortTimer );

}

#elif defined NO_SORT_PERF_TEST



void SortSorting( void ) {
    
    int i;

    QSortTimer = StartTimer();

    for( i = 0; i < NSortBufs; i++ ) {
        if( SortBufs[i].state == SORTING ) {
            SortBufs[i].state = WAIT_WRITE;
        }
    }

    QSortTime += ReadTimer( QSortTimer ) / 1000.0;
    StopTimer( QSortTimer );

}



#else 

sort_buf_t* CurrentSortBuf;
buffer_t* SortScratchBuffer;

void SortSorting( void ) {

    int i;

    QSortTimer = StartTimer();

    // SortData -- sort everything in SORTING -- if it finishes, transition
    // to WAIT_WRITE.
	if( CurrentSortBuf == NULL ){
	    for( i = 0; i < NSortBufs; i++ ) {
	        // if this one is ready to sort, and it's the bin we're looking for...
	        if( SortBufs[i].state == SORTING && SortBufs[i].bin == BinToSort ) {
	        	CurrentSortBuf = &SortBufs[i];
	        	InitRadixSort( CurrentSortBuf, SortScratchBuffer );
	            printf( "scheduling sort of bin %d\n", BinToSort );
	            break;
	        }
	    }
	}
	
	// if there is something to sort right now then try to sort it.
	if( CurrentSortBuf != NULL ){
		if( CurrentSortBuf->state != WRITE_RESTRUCTURE ){

			// automatically transitions to WAIT_WRITE when done.
			RadixSort( CurrentSortBuf );

			// prepare this bin for writing and setup to sort the next
			if( CurrentSortBuf->state == WRITE_RESTRUCTURE ){
				CurrentSortBuf = NULL;
	            BinToSort++;
	        }
		}
    }

    QSortTime += ReadTimer( QSortTimer ) / 1000.0;
    StopTimer( QSortTimer );

}

#endif

void RestructureSMLBinsForWrite( void ) {
    int i;
    offset_t j;
	position_t* positions;
	sml_t *sml;

    for( i = 0; i < NSortBufs; i++ ) {
        // if this one is ready to be restructured...
        if( SortBufs[i].state == WRITE_RESTRUCTURE ) {
            printf( "restructuring bin %d\n", SortBufs[i].bin );
            positions = (position_t*)SortBufs[i].buf->recs;
            sml = (sml_t*)SortBufs[i].buf->recs;
            for( j = 0; j < Bins[SortBufs[i].bin].nrecs; j++ ){
            	positions[ j ] = sml[ j ].pos;
            }
            
            // set its state for writing
            SortBufs[i].state = WAIT_WRITE;
        }
    }
}

// use this version if no pre-write modifications are required
/*
void RestructureSMLBinsForWrite( void ) {
    int i;

    for( i = 0; i < NSortBufs; i++ ) {
        // if this one is ready to be restructured...
        if( SortBufs[i].state == WRITE_RESTRUCTURE ) {
            // set its state for writing
            SortBufs[i].state = WAIT_WRITE;
        }
    }
}
*/

int CalculateSortWriteSize( int sortI ){
     return Bins[SortBufs[sortI].bin].nrecs * sizeof( position_t );
}

void SortWriting( void ) {

    int i;

    for( i = 0; i < NSortBufs; i++ ) {
        // if this one is ready to write, and it's the bin we're looking for...
        if( SortBufs[i].state == WAIT_WRITE && SortBufs[i].bin == BinToWrite ) {
#ifdef NO_WRITE_PERF_TEST
			// skip writing by setting the state to wait_read
            SortBufs[i].state = WAIT_READ;
#else
            printf( "scheduling write of bin %d\n", BinToWrite );
            // write it out.
            SortBufs[i].dev = &(Devices[OutputDev].dev);
            SortBufs[i].state = BUSY_WRITE;
            SortBufs[i].buf->file = Output;
            WriteBuffer( SortBufs[i].buf, Bins[SortBufs[i].bin].nrecs, &(Devices[OutputDev].dev) );
			SortBufs[i].buf->io_size = CalculateSortWriteSize( i );
#endif // NO_WRITE_PERF_TEST
            BinToWrite++;
        }
    }

}





void SortHandleCompletions( void ) {

    int i;

    // transition states of those that finished.
    for( i = 0; i < NSortBufs; i++ ) {
        if( SortBufs[i].state == BUSY_READ || SortBufs[i].state == BUSY_WRITE ) {
            if( SortBufs[i].buf->operation == OP_FINISHED ) {
                //printf( "operation finished on buf %d\n", i );
                SortBufs[i].buf->operation = OP_NONE;
                SortBufs[i].state = SortBufs[i].state == BUSY_READ ? SORTING : WAIT_READ;
#ifdef NNNNN_KEYBYTES
				// bin 0 doesn't need to be sorted
				if( SortBufs[i].bin == 0 && SortBufs[i].state == SORTING )
					SortBufs[i].state = WAIT_WRITE;
#endif
            }
        }
    }

}





void SortUpdateIOState() {

    int i;
    //printf( "update io state\n" );
    
    // first update aio ops on the data file
    aUpdateOperations( Output );
    // next update aio ops on the sortbuf files
    for( i = 0; i < NSortBufs; i++ ) {
        if( SortBufs[i].buf->file ) {
            aUpdateOperations( SortBufs[i].buf->file );
        }
    }
    // next, let the working set adjust operation states and such
    UpdateWSIOFinishedState( &WS );
    // finally, let the devices start new operations if possible.
    for( i = 0; i < NumDevices; i++ ) {
        UpdateDeviceIOExecuteState( &WS, &(Devices[i].dev) );
    }

}






void SortingEnsureAllOperationsComplete() {
    int i;
    int not_complete = 1;
    dmtimer_t *wait;
    wait = StartTimer();
    while( not_complete ) {
        SortUpdateIOState();
        // see if we're done
        not_complete = 0;
        for( i = 0; i < WS.nbufs; i++ ) {
            if( WS.bufs[i].device &&
                WS.bufs[i].file &&
                (WS.bufs[i].operation == OP_PENDING || WS.bufs[i].operation > OP_NONE) ) {
                not_complete = 1;
                break;
            }
        }
    }

    // flush the output file to disk.
    aFlush( Output );

    printf( "Sort Ensure All Operations Complete: %d msec\n", ReadTimer( wait ) );
    StopTimer( wait );
}







void SortingPhase( void ) {

    // now reorganize the working set, and start up the sort procedure.

    // we need to have the ability to read from N bin files at a time, where
    // N is the number of bin devices.

    // We read entire bin files at a time into each slot.  We wait for the
    // first one to finish, and then we sort it.  We can start sorting the
    // others too, as they finish.  When the first sort is done, we write it
    // out to the sorted output file, similarly we write everything out in
    // order.  When one is confirmed finished writing, we can start reading
    // the next bin file from that device in.

    int i;
    offset_t recs_per_buffer;
    offset_t biggest_bin = 0;
    offset_t biggest_nrecs = 0;
    
    NSortBufs = NumBinDevs;

    for( i = 0; i < NumBins; i++ ) {
        if( Bins[i].nrecs > biggest_nrecs ) {
            biggest_nrecs = Bins[i].nrecs;
            biggest_bin = i;
        }
    }
    
    //recs_per_buffer = (WS.size / sizeof( record_t )) / NSortBufs;

    recs_per_buffer = biggest_nrecs;

    if( (WS.size / sizeof( record_t )) < (unsigned)recs_per_buffer ) {
        printf( "working set holds %llu recs, but we need %llu\n", 
            (WS.size / sizeof( record_t )), recs_per_buffer );
    }

    NSortBufs = (WS.size / sizeof( record_t )) / recs_per_buffer;

    printf( "NSortBufs = %d\n", NSortBufs );

    // this goes from 0 to NumBins-1 as we read stuff.
    BinToRead = 0;
    BinToWrite = 0;
	BinToSort = 0;
	
    printf( "reorganizing working set: %llu recs per buffer, %d sort bufs\n", recs_per_buffer, NSortBufs );
    ReorganizeWorkingSet( &WS, recs_per_buffer, recs_per_buffer );

#if !defined USE_QSORT_ONLY && !defined NO_SORT_PERF_TEST
	// steal the last buffer for scratch space
    NSortBufs--;
	SortScratchBuffer = &(WS.bufs[NSortBufs]);
	SortScratchBuffer->operation = SORTING_SCRATCH;
#endif
    
    // nbufs should be same as NumBinDevs
    printf( "reorganized working set has %d buffers of %llu bytes\n", WS.nbufs, recs_per_buffer * sizeof(record_t) );
    SortBufs = malloc( sizeof( *SortBufs ) * NSortBufs );
    memset( SortBufs, 0, sizeof( *SortBufs ) * NSortBufs );

    // put everything in WAIT_READ;
    
    for( i = 0; i < NSortBufs; i++ ) {

        SortBufs[i].state = WAIT_READ;
        SortBufs[i].buf = &(WS.bufs[i]);
        SortBufs[i].dev = NULL;

    }
#ifdef NNNNN_KEYBYTES
	// process the first bin then restructure the working set again
	
    while( BinToWrite < 1 ) {
        SortReading();
        SortSorting();
        RestructureSMLBinsForWrite();
        SortWriting();
        SortUpdateIOState();
        SortHandleCompletions();
    }
    SortingEnsureAllOperationsComplete();

    for( i = 1; i < NumBins; i++ ) {
        if( Bins[i].nrecs > biggest_nrecs ) {
            biggest_nrecs = Bins[i].nrecs;
            biggest_bin = i;
        }
    }
    recs_per_buffer = biggest_nrecs;
    if( (WS.size / sizeof( record_t )) < (unsigned)recs_per_buffer ) {
        printf( "working set holds %llu recs, but we need %llu\n", 
            (WS.size / sizeof( record_t )), recs_per_buffer );
    }
    NSortBufs = (WS.size / sizeof( record_t )) / recs_per_buffer;
    printf( "NSortBufs = %d\n", NSortBufs );
    // this goes from 0 to NumBins-1 as we read stuff.
    BinToRead = 1;
    BinToWrite = 1;
	BinToSort = 1;
	
    printf( "reorganizing working set: %llu recs per buffer, %d sort bufs\n", recs_per_buffer, NSortBufs );
    ReorganizeWorkingSet( &WS, recs_per_buffer, recs_per_buffer );

#if !defined USE_QSORT_ONLY && !defined NO_SORT_PERF_TEST
	// steal the last buffer for scratch space
    NSortBufs--;
	SortScratchBuffer = &(WS.bufs[NSortBufs]);
	SortScratchBuffer->operation = SORTING_SCRATCH;
#endif
    
    // nbufs should be same as NumBinDevs
    printf( "reorganized working set has %d buffers of %llu bytes\n", WS.nbufs, recs_per_buffer * sizeof(record_t) );
    SortBufs = malloc( sizeof( *SortBufs ) * NSortBufs );
    memset( SortBufs, 0, sizeof( *SortBufs ) * NSortBufs );

    // put everything in WAIT_READ;
    for( i = 0; i < NSortBufs; i++ ) {
        SortBufs[i].state = WAIT_READ;
        SortBufs[i].buf = &(WS.bufs[i]);
        SortBufs[i].dev = NULL;
    }
#endif    
    

    while( BinToWrite < NumBins ) {
        
        // ReadFiles -- schedule reading operations if we can (are any buffers
        // in WAIT_READ?)
        //printf( "sortreading\n" );
        SortReading();
        
        // SortData -- sort everything in SORTING -- if it finishes, transition
        // to WAIT_WRITE.
        //printf( "sortsorting\n" );
        SortSorting();
        
        // Perform any necessary post-sort processing on the data to prepare it for
        // writing out to the sorted file
        RestructureSMLBinsForWrite();
        // WriteFiles -- schedule writing operations for everything in WAIT_WRITE, if
        // it is the next file we need to write (make sure to schedule in order).
        //printf( "sortwriting\n" );
        SortWriting();
        
        // update io state
        //printf( "sortupdateiostate\n" );
        SortUpdateIOState();


        // HandleCompletions -- if something finishes,
        // if it was reading, transition to SORTING
        // if it was writing, transition to WAIT_READ.
        //printf( "sorthandlecompletions\n" );
        SortHandleCompletions();
        
        
    }

    SortingEnsureAllOperationsComplete();

    printf( "QSort took %f seconds\n", QSortTime );

}







int dmsort() {


    // Do the first pass binning stuff
    BinningTimer = StartTimer();
#ifndef NO_BINNING_PERF_TEST
    BinningPhase();
    BinningTime = ReadTimer( BinningTimer ) / 1000.0;
#endif // NO_BINNING_PERF_TEST
    StopTimer( BinningTimer );


    // Do the second pass sort
    SortingTimer = StartTimer();
    SortingPhase();
    SortingTime = ReadTimer( SortingTimer ) / 1000.0;
    StopTimer( SortingTimer );


    RunningTime = ReadTimer( RunningTimer ) / 1000.0;
    StopTimer( RunningTimer );

    printf( "total time      : %f sec\n", RunningTime );
    printf( "binning time    : %f sec (%f%%)\n", BinningTime, BinningTime/RunningTime * sizeof(record_t) );
    printf( "sorting time    : %f sec (%f%%)\n", SortingTime, SortingTime/RunningTime * sizeof(record_t) );
    
    printf( "total rate      : %f MB/sec\n", (((double)NumRecs)/10485.760)/RunningTime );
    printf( "total bin rate  : %f MB/sec\n", (((double)NumRecs)/10485.760)/BinningTime );
    printf( "total sort rate : %f MB/sec\n", (((double)NumRecs)/10485.760)/SortingTime );

    return 0;
}


int dmSML( const char* input_file, const char* output_file, const char* const* scratch_paths, uint64 seed ) {
	long working_mb = 300;
	long buffer_size = 1000;
	int rval = 0;
	int i = 0;
	rval = InitdmSML( 0, 0, input_file, output_file, scratch_paths, seed );
	if( rval != 0 )
		return rval;
	rval = dmsort();
	
	// Hey slob!  cleanup after yourself!
	for( i = 0; i < NumBins; i++ ){
		removeFile( Bins[ i ].fname, FALSE );
		free( Bins[ i ].fname );
	}
	if( Bins )
		free( Bins );
	Bins = NULL;
	NumBins = 0;
//	for( i = 0; i < NumDevices; i++ )
//		free( Devices[i].devname );
	NumDevices = 0;
	if( Devices )
		free( Devices );
	Devices = NULL;
	if( SortBufs )
		free( SortBufs );
	SortBufs = NULL;

	NSortBufs = 0;

	BufferSizeMin = 0;
	BufferSizeMax = 0;
	
    memset( &Seqbuf, 0, sizeof( seqbuf_t ) );

	DataDev = 0;

	OutFileName = "unset";

	// close the sorted file
    aClose( Output );
    Output = NULL;
	OutputDev = 0;

	BinToRead = 0;
	BinToWrite = 0;
	BinToSort = 0;

	free( WS.bufs );
	memset( &WS, 0, sizeof( working_set_t ) );

	NumRecs = 0;
	RecsProcessed = 0;
	RecsRead = 0;
	RecsUnread = 0;
	RecsCommitted = 0;
	RecsWritten = 0;
	

// timers
	RunningTime = 0;
	RunningTimer= NULL;
	BinningTime = 0;
	BinningTimer= NULL;
	SortingTime = 0;
	SortingTimer = NULL;

	QSortTime = 0;
	QSortTimer = NULL;

	ReadIdleTime = 0;
	ReadIdleTimer = NULL;
	SortIdleTime = 0;
	SortIdleTimer = NULL;
	WriteIdleTime = 0;
	WriteIdleTimer = NULL;
	
	
	memset( &Free, 0, sizeof( buffer_list_t ) );
	memset( &ToProcess, 0, sizeof( buffer_list_t ) );
	memset( &Reading, 0, sizeof( buffer_list_t ) );
	memset( &Restructure, 0, sizeof( buffer_list_t ) );

	// static variables
	divisor = 0;
	consumed_recs = 0;
	toprocess = NULL;
	lasttime = 0;
	
	// from asyncio.c
//	OperationNumber = 0;
	
	return rval;
}

