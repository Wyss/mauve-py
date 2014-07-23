/*******************************************************************************
 * $Id: HybridAbstractMatch.h,v 1.8 2004/02/27 23:08:55 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef __HybridAbstractMatch_h__
#define __HybridAbstractMatch_h__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnClone.h"
#include "libGenome/gnDefs.h"
#include "libMems/AbstractMatch.h"
#include <vector>
#include <limits>
#include <cstring>
namespace mems {

/**
 * The HybridAbstractMatch implements the AbstractMatch interface in a way
 * that allows matches with a large SeqCount and low Multiplicity to be stored efficiently
 */
template< unsigned FIXED_SEQ_COUNT=2, class int64Alloc=std::allocator<int64>, class uintAlloc=std::allocator<uint> >
class HybridAbstractMatch : public AbstractMatch {
public:
	HybridAbstractMatch() : m_seq_count(0) 
	{
		memset(fixed_seq_ids, 0xFF, sizeof(fixed_seq_ids));
		memset(fixed_starts, 0, sizeof(fixed_starts));
	}
	/**
	 * Creates a new HybridAbstractMatch.
	 * @param seq_count The total number of sequences in the alignment
	 */
	HybridAbstractMatch(const uint seq_count )
		: m_seq_count(seq_count)
	{
		memset(fixed_seq_ids, 0xFF, sizeof(fixed_seq_ids));
		memset(fixed_starts, 0, sizeof(fixed_starts));
	}


	// use compiler-generated copy constructor, assignment operator, and destructor

	// see AbstractMatch base class documentation for these functions

	int64 Start(uint seqI) const;
	void SetStart(uint seqI, int64 startI);
	uint Multiplicity() const
	{
		uint mult = 0;
		for( size_t fI = 0; fI < FIXED_SEQ_COUNT; ++fI )
			mult += fixed_seq_ids[fI] != NO_SEQ ? 1 : 0;
		return mult + (uint)seq_ids.size();
	}
	uint SeqCount() const{return m_seq_count;}
	uint FirstStart() const;
	virtual void Invert();

	gnSeqI LeftEnd(uint seqI) const;
	orientation Orientation(uint seqI) const;
	void SetLeftEnd(uint seqI, gnSeqI position);
	void SetOrientation(uint seqI, orientation o);
	
	// these functions manipulate the start coordinates quickly
	virtual void MoveStart(int64 move_amount);
	virtual void MoveEnd(int64 move_amount);

	virtual boolean operator==( const HybridAbstractMatch& ham ) const;

	virtual uint UsedSeq( uint seqI ) const { 
		if(seqI < FIXED_SEQ_COUNT) return fixed_seq_ids[seqI];
		return seq_ids[seqI];
	}

protected:
	uint m_seq_count;

	static const uint NO_SEQ = UINT_MAX;

	// storage for a fixed number of seqs
	uint fixed_seq_ids[FIXED_SEQ_COUNT];
	int64 fixed_starts[FIXED_SEQ_COUNT];

	// storage for any number of seqs
	std::vector<uint, uintAlloc > seq_ids;
	std::vector<int64, int64Alloc > starts;

	uint SeqToIndex( uint seqI ) const;

	// for use by derived classes in order to swap contents
	void swap( HybridAbstractMatch* other );
};


template< unsigned FIXED_SEQ_COUNT, class gnSeqIAlloc, class uintAlloc >
void HybridAbstractMatch< FIXED_SEQ_COUNT, gnSeqIAlloc, uintAlloc >::swap( HybridAbstractMatch* other )
{
	std::swap( m_seq_count, other->m_seq_count );

	uint tmp_ids[FIXED_SEQ_COUNT];
	for( int i = 0; i < FIXED_SEQ_COUNT; i++ ) tmp_ids[i] = other->fixed_seq_ids[i];
	for( int i = 0; i < FIXED_SEQ_COUNT; i++ ) other->fixed_seq_ids[i] = fixed_seq_ids[i];
	for( int i = 0; i < FIXED_SEQ_COUNT; i++ ) fixed_seq_ids[i] = tmp_ids[i];

	int64 tmp_starts[FIXED_SEQ_COUNT];
	for( int i = 0; i < FIXED_SEQ_COUNT; i++ ) tmp_starts[i] = other->fixed_starts[i];
	for( int i = 0; i < FIXED_SEQ_COUNT; i++ ) other->fixed_starts[i] = fixed_starts[i];
	for( int i = 0; i < FIXED_SEQ_COUNT; i++ ) fixed_starts[i] = tmp_starts[i];

	std::swap( seq_ids, other->seq_ids );
	std::swap( starts, other->starts );
}

template< unsigned FIXED_SEQ_COUNT, class gnSeqIAlloc, class uintAlloc >
uint HybridAbstractMatch< FIXED_SEQ_COUNT, gnSeqIAlloc, uintAlloc >::FirstStart() const
{
	uint minI = NO_SEQ;
	std::size_t i = 0;
	for( ; i < FIXED_SEQ_COUNT; ++i )
		minI = fixed_seq_ids[i] < minI ? fixed_seq_ids[i] : minI;
	for( i = 0; i < seq_ids.size(); ++i )
		minI = seq_ids[i] < minI ? seq_ids[i] : minI;
	return minI;
}

template< unsigned FIXED_SEQ_COUNT, class gnSeqIAlloc, class uintAlloc >
uint HybridAbstractMatch< FIXED_SEQ_COUNT, gnSeqIAlloc, uintAlloc >::SeqToIndex( uint seqI ) const
{
	uint posI = 0;
	for( ; posI < FIXED_SEQ_COUNT; ++posI )
		if( fixed_seq_ids[posI] == seqI )
			break;
	if(posI < FIXED_SEQ_COUNT)
		return posI;
	for( posI = 0; posI < seq_ids.size(); ++posI )
		if( seq_ids[posI] == seqI )
			break;
	if( posI == seq_ids.size() )
		return NO_SEQ;
	return posI + FIXED_SEQ_COUNT;
}


template< unsigned FIXED_SEQ_COUNT, class gnSeqIAlloc, class uintAlloc >
int64 HybridAbstractMatch< FIXED_SEQ_COUNT, gnSeqIAlloc, uintAlloc >::Start(uint seqI) const
{
	uint posI = SeqToIndex( seqI );
	if( posI == NO_SEQ )
		return NO_MATCH;
	if( posI < FIXED_SEQ_COUNT )
		return fixed_starts[posI];
	return starts[posI-FIXED_SEQ_COUNT];
}


template< unsigned FIXED_SEQ_COUNT, class gnSeqIAlloc, class uintAlloc >
void HybridAbstractMatch< FIXED_SEQ_COUNT, gnSeqIAlloc, uintAlloc >::SetStart(uint seqI, int64 startI)
{
	uint posI = SeqToIndex( seqI );
	if( startI == NO_MATCH && posI == NO_SEQ )
		return;
	if( posI == NO_SEQ )
	{
		for( size_t i = 0; i < FIXED_SEQ_COUNT; ++i )
			if( fixed_seq_ids[i] == NO_SEQ )
			{
				posI = i;
				break;
			}
	}
	if( posI < FIXED_SEQ_COUNT )
	{
		if( startI == NO_MATCH )
			fixed_seq_ids[posI] = NO_SEQ;
		else
			fixed_seq_ids[posI] = seqI;
		fixed_starts[posI] = startI;
	}
	else
	{
		posI -= FIXED_SEQ_COUNT;
		if( startI == NO_MATCH )
		{
			seq_ids.erase( seq_ids.begin() + posI );
			starts.erase( starts.begin() + posI );
			return;
		}
		if( posI >= seq_ids.size() )
		{
			seq_ids.push_back(seqI);
			starts.push_back(startI);
		}else{
			starts[posI] = startI; 
		}
	}
}


template< unsigned FIXED_SEQ_COUNT, class gnSeqIAlloc, class uintAlloc >
void HybridAbstractMatch< FIXED_SEQ_COUNT, gnSeqIAlloc, uintAlloc >::Invert()
{
	for( size_t i = 0; i < FIXED_SEQ_COUNT; ++i )
		fixed_starts[i] = -fixed_starts[i];
	for( size_t i = 0; i < starts.size(); ++i )
		starts[i] = -starts[i];
}



template< unsigned FIXED_SEQ_COUNT, class gnSeqIAlloc, class uintAlloc >
gnSeqI HybridAbstractMatch< FIXED_SEQ_COUNT, gnSeqIAlloc, uintAlloc >::LeftEnd(uint seqI) const
{ 
	uint posI = SeqToIndex( seqI );
	if( posI == NO_SEQ )
		return NO_MATCH;
	if( posI < FIXED_SEQ_COUNT )
		return genome::absolut(fixed_starts[posI]);
	return genome::absolut(starts[posI-FIXED_SEQ_COUNT]);
}


template< unsigned FIXED_SEQ_COUNT, class gnSeqIAlloc, class uintAlloc >
AbstractMatch::orientation HybridAbstractMatch< FIXED_SEQ_COUNT, gnSeqIAlloc, uintAlloc >::Orientation(uint seqI) const
{ 
	uint posI = SeqToIndex( seqI );
	if( posI == NO_SEQ )
		return undefined;
	if( posI < FIXED_SEQ_COUNT )
		return fixed_starts[posI] < 0 ? reverse : forward;
	return starts[posI-FIXED_SEQ_COUNT] < 0 ? reverse : forward;
}


template< unsigned FIXED_SEQ_COUNT, class gnSeqIAlloc, class uintAlloc >
void HybridAbstractMatch< FIXED_SEQ_COUNT, gnSeqIAlloc, uintAlloc >::SetLeftEnd(uint seqI, gnSeqI position)
{ 
	uint posI = SeqToIndex( seqI );
	orientation o = posI == NO_SEQ || position == NO_MATCH ? undefined : Orientation( seqI );
	SetStart(seqI,position);
	if( o != undefined )
		SetOrientation(seqI, o);
}

template< unsigned FIXED_SEQ_COUNT, class gnSeqIAlloc, class uintAlloc >
void HybridAbstractMatch< FIXED_SEQ_COUNT, gnSeqIAlloc, uintAlloc >::SetOrientation(uint seqI, orientation o)
{ 
	if( o == undefined )
	{
		SetStart(seqI, NO_MATCH);
		return;
	}
	uint posI = SeqToIndex( seqI );
	if( posI == NO_SEQ )
		throw "ArrayIndexOutOfBounds!\n";
	int oi = o == reverse ? -1 : 1;
	if( posI < FIXED_SEQ_COUNT )
	{
		fixed_starts[posI] = genome::absolut(fixed_starts[posI]) * oi;
		return;
	}
	starts[posI-FIXED_SEQ_COUNT] = genome::absolut(starts[posI-FIXED_SEQ_COUNT]) * oi;
}

template< unsigned FIXED_SEQ_COUNT, class gnSeqIAlloc, class uintAlloc >
void HybridAbstractMatch< FIXED_SEQ_COUNT, gnSeqIAlloc, uintAlloc >::MoveStart(int64 move_amount)
{
	for( size_t i=0; i < FIXED_SEQ_COUNT; ++i )
		if( fixed_starts[i] > 0 )
			fixed_starts[i] += move_amount;
	for( size_t i=0; i < starts.size(); ++i )
		if( starts[i] > 0 )
			starts[i] += move_amount;
}

template< unsigned FIXED_SEQ_COUNT, class gnSeqIAlloc, class uintAlloc >
void HybridAbstractMatch< FIXED_SEQ_COUNT, gnSeqIAlloc, uintAlloc >::MoveEnd(int64 move_amount)
{
	for( size_t i=0; i < FIXED_SEQ_COUNT; ++i )
		if( fixed_starts[i] < 0 )
			fixed_starts[i] -= move_amount;
	for( size_t i=0; i < starts.size(); ++i )
		if( starts[i] < 0 )
			starts[i] -= move_amount;
}

template< unsigned FIXED_SEQ_COUNT, class gnSeqIAlloc, class uintAlloc >
boolean HybridAbstractMatch< FIXED_SEQ_COUNT, gnSeqIAlloc, uintAlloc >::operator==( const HybridAbstractMatch< FIXED_SEQ_COUNT, gnSeqIAlloc, uintAlloc >& sam ) const
{
	for( size_t i = 0; i < FIXED_SEQ_COUNT; ++i )
	{
		if( fixed_seq_ids[i] == NO_SEQ )
			continue;
		if( Start(fixed_seq_ids[i]) !=  sam.Start(fixed_seq_ids[i]) )
			return false;
	}
	for( size_t i = 0; i < seq_ids.size(); ++i )
	{
		if( seq_ids[i] == NO_SEQ )
			continue;
		if( Start(seq_ids[i]) !=  sam.Start(seq_ids[i]) )
			return false;
	}
	return Multiplicity() == sam.Multiplicity();
}


}

#endif // __HybridAbstractMatch_h__
