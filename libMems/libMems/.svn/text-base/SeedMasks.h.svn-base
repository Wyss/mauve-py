/*******************************************************************************
 * $Id: SortedMerList.h,v 1.13 2004/02/27 23:08:55 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _SeedMasks_h_
#define _SeedMasks_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef __cplusplus
#include <cmath>
#else
#include <math.h>
#endif
#include "libGenome/gnDefs.h"

/* Seed patterns taken from: AE Darling, T Treangen, L Zhang, C Kuiken, X Messeguer, NT Perna
 * "Procrastination leads to efficient match filtration for local multiple alignment" 
 * Lecture Notes in Bioinformatics 4175:126-137 Springer-Verlag 2006
 */

/**
 * returns the array of default seed mask patterns
 * Each seed is represented by a pair of 32 bit integers
 */
#ifdef __cplusplus
static
#endif
uint32** seedMasks();

/**
 * the first three seed masks in each of the following are
 * 'good' seeds according to Louxin Zhang
 */
#ifdef __cplusplus
inline static
#endif
uint32** seedMasks(){
	static uint32 seed_masks_3[] = 
	{
		0,0xb, //0b1011
		0,0,		0,0,		0,0,		0,0,		0,0,	};
	static uint32 seed_masks_4[] = 
	{
		0,0x3b, //0b101011,
		0,0,		0,0,		0,0,		0,0,		0,0,	};
	static uint32 seed_masks_5[] = 
	{
		0,0x6b,		//0b1101011,
		0,0x139,	//0b100111001,
		0,0x193,	//0b110010011,
		0,0x6b,		//0b1101011,
		0,0,		0,0,	};
	static uint32 seed_masks_6[] = 
	{
		0,0x58D, //0b10110001101,
		0,0x653, //0b11001010011,
		0,0x1AB, //0b110101011,
		0,0xdb,	//0b11011011,
		0,0,		0,0,	};
	static uint32 seed_masks_7[] = 
	{
		0,0x1953,	//0b1100101010011
		0,0x588d,	//0b101100010001101
		0,0x688b,	//0b110100010001011
		0,0x17d,	//0b101111101,
		0,0x164d,	//0b1011001001101,
		0,0,		0,0,	};
	static uint32 seed_masks_8[] = 
	{
		0,0x3927, //0b11100100100111,
		0,0x1CA7, //0b1110010100111,
		0,0x6553, //0b110010101010011,
		0,0xb6d,	//0b101101101101,
		0,0,		0,0,	};
	static uint32 seed_masks_9[] = 
	{
		0,0x7497,	//0b111010010010111,
		0,0x1c927,	//0b11100100100100111,
		0,0x72a7,	//0b111001010100111,
		0,0x6fb,	//0b11011111011,
		0,0x16ed,	//0b1011011101101,
		0,0,	};
	static uint32 seed_masks_10[] = 
	{
		0,0x1d297,	//		0,0b11101001010010111,
		0,0x3A497,  //		0,0b111010010010010111,
		0,0xE997,  //		0,0b1110100110010111,
		0,0x6D5B,  //		0,0b110110101011011,
		0,0,		0,0,	};
	static uint32 seed_masks_11[] = 
	{
		0,0x7954f,	//0b11110010101001111,
		0,0x75257,	//0b1110101001001010111,
		0,0x1c9527,	//0b111001001010100100111,
		0,0x5bed,	//0b101101111101101,  // third b.p. coding pattern
		0,0x5b26d,	//0b1011011001001101101,
		0,0,	};
	static uint32 seed_masks_12[] = 
	{
		0,0x7954f,	//		0,0b1111001010101001111,
		0,0x3D32F,  //		0,0b111101001100101111,
		0,0x768B7,  //		0,0b1110110100010110111,
		0,0x5B56D,  //		0,0b1011011010101101101,
		0,0,		0,0,	};
	static uint32 seed_masks_13[] = 
	{
		0,0x792a4f,	//0b11110010010101001001111,
		0,0x1d64d7,	//0b111010110010011010111,
		0,0x1d3597,	//0b111010011010110010111,
		0,0x1b7db,	//0b11011011111011011,  // third b.p. coding pattern
		0,0x75ad7,	//0b1110101101011010111,
		0,0,	};
	static uint32 seed_masks_14[] = 
	{
		0,0x1e6acf,  //		0,0b111100110101011001111,
		0,0xF59AF,   //		0,0b11110101100110101111,
		0,0x3D4CAF,  //		0,0b1111010100110010101111,
		0,0x35AD6B,  //		0,0b1101011010110101101011,
		0,0,		0,0,	};
	static uint32 seed_masks_15[] = 
	{
		0,0x7ac9af,	//0b11110101100100110101111,
		0,0x7b2a6f,	//0b11110110010101001101111,
		0,0x79aacf,	//0b11110011010101011001111,
		0,0x16df6d,	//0b101101101111101101101,	// third b.p. coding pattern
		0,0x6b5d6b,	//0b11010110101110101101011,
		0,0,	};
	static uint32 seed_masks_16[] = 
	{
		0,0xf599af,  //		0,0b111101011001100110101111,
		0,0xEE5A77,  //		0,0b111011100101101001110111,
		0,0x7CD59F,  //		0,0b11111001101010110011111,
		0,0xEB5AD7,  //		0,0b111010110101101011010111,
		0,0,		0,0,	};
	static uint32 seed_masks_17[] =
	{
		0,0x6dbedb,	//0b11011011011111011011011,	// third b.p. coding pattern
		0,0,		0,0,		0,0,		0,0,		0,0,	};
	static uint32 seed_masks_18[] =
	{
		0,0x3E6B59F,//		0,0b11111001101011010110011111,
		0,0x3EB335F,//		0,0b11111010110011001101011111,
		0,0x7B3566F,//		0,0b111101100110101011001101111,
		0,0,		0,0,		0,0,	};

	static uint32 seed_masks_19[] =
	{
		0,0x7b974ef,	//0b111101110010111010011101111
		0,0x7d6735f,	//0b111110101100111001101011111
		0,0x1edd74f,	//0b1111011011101011101101111
		0,0,		0,0,		0,0,	};
	static uint32 seed_masks_20[] =
	{
		0,0x1F59B35F,	//0b11111010110011011001101011111,
		0,0x3EDCEDF,	//0b11111011011100111011011111,
		0,0xFAE675F,	//0b1111101011100110011101011111,
		0,0,		0,0,		0,0,	};
	static uint32 seed_masks_21[] =
	{
		0,0x7ddaddf,	//0b111110111011010110111011111,
		0,0xaeb3f,		//0b11111100110101110101100111111,
		0,0x7eb76bf,	//0b111111010110111011010111111,
		0,0,		0,0,		0,0,	};
	// default to solid seeds for weight 22+
	static uint32 seed_masks_22[] =
	{
		0,0x003fffff,
		0,0,		0,0,		0,0,		0,0,		0,0,	};
	static uint32 seed_masks_23[] =
	{
		0,0x007fffff,
		0,0,		0,0,		0,0,		0,0,		0,0,	};
	static uint32 seed_masks_24[] =
	{
		0,0x00ffffff,
		0,0,		0,0,		0,0,		0,0,		0,0,	};
	static uint32 seed_masks_25[] =
	{
		0,0x01ffffff,
		0,0,		0,0,		0,0,		0,0,		0,0,	};
	static uint32 seed_masks_26[] =
	{
		0,0x03ffffff,
		0,0,		0,0,		0,0,		0,0,		0,0,	};
	static uint32 seed_masks_27[] =
	{
		0,0x07ffffff,
		0,0,		0,0,		0,0,		0,0,		0,0,	};
	static uint32 seed_masks_28[] =
	{
		0,0x0fffffff,
		0,0,		0,0,		0,0,		0,0,		0,0,	};
	static uint32 seed_masks_29[] =
	{
		0,0x1fffffff,
		0,0,		0,0,		0,0,		0,0,		0,0,	};
	static uint32 seed_masks_30[] =
	{
		0,0x3fffffff,
		0,0,		0,0,		0,0,		0,0,		0,0,	};
	static uint32 seed_masks_31[] =
	{
		0,0x7fffffff,
		0,0,		0,0,		0,0,		0,0,		0,0,	};

	static uint32 no_seeds[] = 
	{
		0,0,	
		0,0,	
		0,0,	
		0,0,	
		0,0,	
		0,0,	
	};
	
	static uint32* seed_masks[] = 
	{
	no_seeds,
	no_seeds,
	no_seeds,
	seed_masks_3,
	seed_masks_4,
	seed_masks_5,
	seed_masks_6,
	seed_masks_7,
	seed_masks_8,
	seed_masks_9,
	seed_masks_10,
	seed_masks_11,
	seed_masks_12,
	seed_masks_13,
	seed_masks_14,
	seed_masks_15,
	seed_masks_16,
	seed_masks_17,
	seed_masks_18,
	seed_masks_19,
	seed_masks_20,
	seed_masks_21,
	seed_masks_22,
	seed_masks_23,
	seed_masks_24,
	seed_masks_25,
	seed_masks_26,
	seed_masks_27,
	seed_masks_28,
	seed_masks_29,
	seed_masks_30,
	seed_masks_31,
	};
	
	return seed_masks;
}

static const int CODING_SEED = 3;
static const int SOLID_SEED = INT_MAX;

/**
 * Returns a solid seed of a given weight.
 */
#ifdef __cplusplus
static
#endif
int64 getSolidSeed( int weight );

#ifdef __cplusplus
inline static
#endif
int64 getSolidSeed( int weight ){
	int64 seed = 1;
	seed <<= weight;
	seed--;
	return seed;
};



/**
 * returns a seed of a given weight.  Setting seed_rank > 0 will select a seed
 * of a lower sensitivity rank according to Choi et. al. 2004
 */
#ifdef __cplusplus
static int64 getSeed( int weight, int seed_rank = 0 );
#else
int64 getSeed( int weight, int seed_rank );
#endif

#ifdef __cplusplus
inline static
#endif
int64 getSeed( int weight, int seed_rank ){
	uint32** masks;
	int high;
	int low;
	int i = 1;
	int64 seed = 0;
	if( seed_rank == SOLID_SEED )
		return getSolidSeed( weight );

	masks = seedMasks();
	if(weight > 31)
		return getSolidSeed(32);
	if( seed_rank > 5 )
		return getSolidSeed(weight);
	if( masks[weight][seed_rank*2+1] == 0 )
		return getSolidSeed(weight);
	high = masks[ weight ][ seed_rank*2 ];
	low = masks[ weight ][ seed_rank*2 + 1 ];
	
	seed |= high;
	seed <<= 32;
	seed |= low;
	return seed;
};


/**
 * calculates the length of a seed pattern
 */
#ifdef __cplusplus
static
#endif
int getSeedLength( int64 seed );

#ifdef __cplusplus
inline static
#endif
int getSeedLength( int64 seed ){
	int right_bit = -1;
	int left_bit = -1;
	uint bitI = 0;
	for( ; bitI < 64; ++bitI ){
		if( (seed & 1) == 1 ){
			left_bit = bitI;
			if( right_bit == -1 )
				right_bit = bitI;
		}
		seed >>= 1;
	}
	if( left_bit != -1 )
		return left_bit - right_bit + 1;
	return 0;
}

/**
 * calculates the weight of a seed pattern
 */
#ifdef __cplusplus
static
#endif
int getSeedWeight( int64 seed );

#ifdef __cplusplus
inline static
#endif
int getSeedWeight( int64 seed ){
	int weight = 0;
	uint bitI = 0;
	for( ; bitI < 64; ++bitI ){
		if( (seed & 1) == 1 ){
			++weight;
		}
		seed >>= 1;
	}
	return weight;
}

const uint MIN_DNA_SEED_WEIGHT = 5;
const uint MAX_DNA_SEED_WEIGHT = 31;

/**
 * Calculate the default seed weight based on sequence length
 */
#ifdef __cplusplus
static
#endif
uint getDefaultSeedWeight( gnSeqI avg_sequence_length );

#ifdef __cplusplus
inline static
#endif
uint getDefaultSeedWeight( gnSeqI avg_sequence_length ){
	uint mer_size = (uint)ceil((log( (double)avg_sequence_length ) / log( 2.0 ))/1.5);
	// don't allow even weights-- they can be palindromic
	if( !(mer_size & 0x1 ) )
		++mer_size;
	mer_size = mer_size < MIN_DNA_SEED_WEIGHT ? 0 : mer_size;
	if( avg_sequence_length == 0 )
		mer_size = 0;

	// 31 is the maximum DNA seed weight
	mer_size = mer_size > MAX_DNA_SEED_WEIGHT ? MAX_DNA_SEED_WEIGHT : mer_size;
	return mer_size;
}


#endif // _SeedMasks_h_
