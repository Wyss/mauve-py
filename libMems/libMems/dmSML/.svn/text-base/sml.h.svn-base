#ifndef _sml_h_
#define _sml_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/dmSML/asyncio.h"
#include <string.h>
#include "libGenome/gnDefs.h"

#ifndef UINT8_MAX

#define UINT8_MAX 256
typedef unsigned char uint8;
typedef unsigned uint32;
typedef unsigned long long uint64;

#endif

static uint8* CreateBasicDNATable(){
	uint8* bdt = (uint8*)malloc( sizeof(uint8) * UINT8_MAX );
	memset(bdt, 0, UINT8_MAX);
	bdt['c'] = 1;
	bdt['C'] = 1;
	bdt['b'] = 1;
	bdt['B'] = 1;
	bdt['y'] = 1;
	bdt['Y'] = 1;
	bdt['g'] = 2;
	bdt['G'] = 2;
	bdt['s'] = 2;
	bdt['S'] = 2;
	bdt['k'] = 2;
	bdt['K'] = 2;
	bdt['t'] = 3;
	bdt['T'] = 3;
	return bdt;
}

static uint8* DNA_TABLE;
typedef unsigned position_t;
typedef unsigned long long mask_t;
#define MASK_T_BYTES 8
static mask_t seed_mask = 0x7FFFFFFF;
static int mask_length = 31;
static int mask_weight = 31;

#define DESCRIPTION_SIZE 2048	/**< Number of bytes for the freeform text description of an SML */


typedef signed short sarID_t;

typedef struct SMLHeader_s{
	uint32 version;						/**< Format version - 4 bytes */
	uint32 alphabet_bits;				/**< Bits per character in the alphabet - 4 bytes */
	uint64 seed;						/**< The pattern used in each seed */
	uint32 seed_length;					/**< The length of the seed mask */
	uint32 seed_weight;					/**< The weight of the seed mask */
	uint64 length;						/**< length of the sequence before circularity - 8 bytes */
	uint32 unique_mers;					/**< Number of unique mers in the sequence 4 bytes */
	uint32 word_size;					/**< Word size on the machine the sequence was translated */
	boolean little_endian;				/**< Is the byte order little endian?  0==no, !0==yes */
	signed short id;					/**< Obsolete ID value - 1 byte, eaten by alignment? */
	boolean circular;					/**< Circularity of sequence - 1 byte */
	uint8 translation_table[UINT8_MAX];	/**< Translation table for ascii characters to binary values -- 256 bytes */
	char description[DESCRIPTION_SIZE]; /**< Freeform text description of sequence data -- 2048 bytes */
} SMLHeader_t;


typedef struct sml_s {
		char key[8];
		position_t pos;
} sml_t;

SMLHeader_t InitSML( aFILE* file, uint64 file_size, uint64 seed );


#endif /* _sml_h_ */
