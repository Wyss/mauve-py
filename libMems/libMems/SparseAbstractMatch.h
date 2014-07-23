/*******************************************************************************
 * $Id: SparseAbstractMatch.h,v 1.8 2004/02/27 23:08:55 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef __SparseAbstractMatch_h__
#define __SparseAbstractMatch_h__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnClone.h"
#include "libGenome/gnDefs.h"
#include "libMems/AbstractMatch.h"
#include <vector>
#include <limits>

namespace mems {

//template< class gnSeqIAlloc=boost::pool_allocator<gnSeqI>, class uintAlloc=boost::pool_allocator<uint> >
/**
 * The SparseAbstractMatch implements the AbstractMatch interface in a way
 * that allows matches with a large SeqCount and low Multiplicity to be stored efficiently
 */
template< class gnSeqIAlloc=std::allocator<gnSeqI>, class uintAlloc=std::allocator<uint> >
class SparseAbstractMatch : public AbstractMatch {
public:
	SparseAbstractMatch() : m_seq_count(0) {}
	/**
	 * Creates a new SparseAbstractMatch.
	 * @param seq_count The total number of sequences in the alignment
	 */
	SparseAbstractMatch(const uint seq_count );

	// use compiler-generated copy constructor, assignment operator, and destructor

	// see AbstractMatch base class documentation for these functions

	int64 Start(uint seqI) const;
	void SetStart(uint seqI, int64 startI);
	uint Multiplicity() const{return (uint)seq_ids.size();}
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

	virtual boolean operator==( const SparseAbstractMatch& sam ) const;

	virtual uint UsedSeq( uint seqI ) const;
protected:

	std::vector<uint, uintAlloc > seq_ids;
	uint m_seq_count;
	std::vector<gnSeqI, gnSeqIAlloc > leftend;
	bitset_t orient;	// bitset_t has its own allocator
	uint SeqToIndex( uint seqI ) const;

	// for use by derived classes in order to swap contents
	void swap( SparseAbstractMatch* other );	
};


template< class gnSeqIAlloc, class uintAlloc >
SparseAbstractMatch< gnSeqIAlloc, uintAlloc >::SparseAbstractMatch(const uint seq_count ) :
m_seq_count(seq_count)
{}

template< class gnSeqIAlloc, class uintAlloc >
void SparseAbstractMatch< gnSeqIAlloc, uintAlloc >::swap( SparseAbstractMatch* other )
{
	std::swap(seq_ids, other->seq_ids);
	std::swap(m_seq_count, other->m_seq_count);
	std::swap(leftend, other->leftend);
	std::swap(orient, other->orient);
}

template< class gnSeqIAlloc, class uintAlloc >
uint SparseAbstractMatch< gnSeqIAlloc, uintAlloc >::FirstStart() const
{
	uint minI = (std::numeric_limits<uint>::max)();
	for( std::size_t i = 0; i < seq_ids.size(); ++i )
		minI = seq_ids[i] < minI ? seq_ids[i] : minI;
	return minI;
}

template< class gnSeqIAlloc, class uintAlloc >
uint SparseAbstractMatch< gnSeqIAlloc, uintAlloc >::SeqToIndex( uint seqI ) const
{
	uint posI = 0;
	for( ; posI < seq_ids.size(); ++posI )
		if( seq_ids[posI] == seqI )
			break;
	return posI;
}


template< class gnSeqIAlloc, class uintAlloc >
int64 SparseAbstractMatch< gnSeqIAlloc, uintAlloc >::Start(uint seqI) const
{
	uint posI = SeqToIndex( seqI );
	if( posI >= seq_ids.size() )
		return NO_MATCH;
	int64 s = leftend[posI];
	return orient.test(posI)? -s : s;
}


template< class gnSeqIAlloc, class uintAlloc >
void SparseAbstractMatch< gnSeqIAlloc, uintAlloc >::SetStart(uint seqI, int64 startI)
{
	uint posI = SeqToIndex( seqI );
	if( startI == NO_MATCH && posI >= seq_ids.size() )
		return;
	if( startI == NO_MATCH )
	{
		seq_ids.erase( seq_ids.begin() + posI );
		leftend.erase( leftend.begin() + posI );
		for( size_t i = posI; i + 1 < orient.size(); ++i )
			orient.set( i, orient.test( i + 1 ) );
		orient.resize( orient.size()-1 );
		return;
	}
	if( posI >= seq_ids.size() )
	{
		seq_ids.push_back(seqI);
		leftend.push_back(genome::absolut(startI));
		orient.resize( orient.size() + 1, (startI < 0) );
	}else{
		leftend[posI] = genome::absolut(startI); 
		orient.set(posI, startI < 0);
	}
}


template< class gnSeqIAlloc, class uintAlloc >
void SparseAbstractMatch< gnSeqIAlloc, uintAlloc >::Invert()
{
	orient.flip();
}



template< class gnSeqIAlloc, class uintAlloc >
gnSeqI SparseAbstractMatch< gnSeqIAlloc, uintAlloc >::LeftEnd(uint seqI) const
{ 
	uint posI = SeqToIndex( seqI );
	return posI < leftend.size() ? leftend[posI] : 0;
}


template< class gnSeqIAlloc, class uintAlloc >
AbstractMatch::orientation SparseAbstractMatch< gnSeqIAlloc, uintAlloc >::Orientation(uint seqI) const
{ 
	uint posI = SeqToIndex( seqI );
	if( posI < leftend.size() && leftend[posI] != NO_MATCH )
		return orient.test(posI) ? reverse : forward; 
	return undefined;
}


template< class gnSeqIAlloc, class uintAlloc >
void SparseAbstractMatch< gnSeqIAlloc, uintAlloc >::SetLeftEnd(uint seqI, gnSeqI position)
{ 
	uint posI = SeqToIndex( seqI );
	if( position == NO_MATCH && posI >= seq_ids.size() )
		return;
	if( posI >= leftend.size() )
	{
		seq_ids.push_back(seqI);
		leftend.push_back(position);
		orient.resize( orient.size() + 1 );	// defaults to false
	}else if( position == NO_MATCH )
	{
		seq_ids.erase( seq_ids.begin() + posI );
		leftend.erase( leftend.begin() + posI );
		for( size_t i = posI; i + 1 < orient.size(); ++i )
			orient.set( i, orient.test( i + 1 ) );
		orient.resize( orient.size()-1 );
		return;
	}

	leftend[posI]=position; 
}

template< class gnSeqIAlloc, class uintAlloc >
void SparseAbstractMatch< gnSeqIAlloc, uintAlloc >::SetOrientation(uint seqI, orientation o)
{ 
	uint posI = SeqToIndex( seqI );
	// just assume that posI is in-bounds... if not throw an exception!
	if( posI >= orient.size() )
		throw "ArrayIndexOutOfBounds!\n";
	orient.set(posI, o == reverse);
}

template< class gnSeqIAlloc, class uintAlloc >
void SparseAbstractMatch< gnSeqIAlloc, uintAlloc >::MoveStart(int64 move_amount)
{
	for( uint i=0; i < leftend.size(); ++i )
		if( orient.test(i) == false && leftend[i] != NO_MATCH )
			leftend[i] += move_amount;
}

template< class gnSeqIAlloc, class uintAlloc >
void SparseAbstractMatch< gnSeqIAlloc, uintAlloc >::MoveEnd(int64 move_amount)
{
	for( uint i=0; i < leftend.size(); ++i )
		if( orient.test(i) && leftend[i] != NO_MATCH )
			leftend[i] += move_amount;
}

template< class gnSeqIAlloc, class uintAlloc >
boolean SparseAbstractMatch< gnSeqIAlloc, uintAlloc >::operator==( const SparseAbstractMatch< gnSeqIAlloc, uintAlloc >& sam ) const
{
	for( uint i=0; i < leftend.size(); ++i ){
		if( leftend[i] != sam.leftend[i] ||
			(leftend[i] != 0 && orient.test(i) != sam.orient.test(i)))
			return false;
	}
	return true;
}

template< class gnSeqIAlloc, class uintAlloc >
uint SparseAbstractMatch< gnSeqIAlloc, uintAlloc >::UsedSeq( uint seqI ) const
{
	uint count = 0;
	for( uint i = 0; i < leftend.size(); i++ )
	{
		if(leftend[i] != 0)
			count++;
		if( count > seqI )
			return i;
	}
	return (std::numeric_limits<uint>::max)();
}

}

#endif // __SparseAbstractMatch_h__
