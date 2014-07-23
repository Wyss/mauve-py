#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/dmSML/sml.h"
#include "libMems/SeedMasks.h"


SMLHeader_t InitSML( aFILE* file, uint64 file_size, uint64 seed ){
	SMLHeader_t header;
	int retcode;
	
	header.version = 5;
	header.alphabet_bits = 2;
	header.seed = seed;
	header.seed_length = getSeedLength( seed );
	header.seed_weight = getSeedWeight( seed );
	header.length = file_size;
	header.unique_mers = -1;
	header.word_size = 32;
	header.little_endian = 1;
	header.id = 0;
	header.circular = 0;
	memcpy(header.translation_table, CreateBasicDNATable(), UINT8_MAX);
	header.description[ 0 ] = 0;
	
	retcode = aWrite( (void*)&header, sizeof( header ), 1, file, 0 );
	if( retcode == 0 )
		printf( "Error writing to SML\n" );
	aWaitComplete( file, retcode );
	return header;
}

/*
// use this version of RestructureReadSMLBins when no restructuring is necessary
void RestructureReadSMLBins( void ) {
    buffer_t *b, *tmpnext;
    // go through and see if any have completed.
    b = Restructure.head;
    do {
        if( !b ) {
            break;
        }

        tmpnext = b->next;
		
		// b has been restructured, add it to the ToProcess list
        PushTail( &ToProcess, RemoveItem( &Restructure, b ) );
        // bookkeeping
        RecsRead += b->numrecs;
		
        b = tmpnext;
    } while( b != Restructure.head && Restructure.nitems );
}
*/
