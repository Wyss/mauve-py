/*******************************************************************************
 * $Id: DenseAbstractMatch.h,v 1.8 2004/02/27 23:08:55 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef __DenseAbstractMatch_h__
#define __DenseAbstractMatch_h__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnClone.h"
#include "libMems/AbstractMatch.h"
#include <limits>

namespace mems {

/**
 * The DenseAbstractMatch implements the AbstractMatch interface in a way
 * that is most efficient when Multiplicity and SeqCount are identical or 
 * nearly so.  It stores all data inline in a fixed size arrays, affording it
 * storage in a contiguous chunk of memory.
 */
template< unsigned int MAX_SEQS >
class DenseAbstractMatch : public AbstractMatch 
{
public:
	DenseAbstractMatch();
	/**
	 * Creates a new AbstractMatch.
	 * @param seq_count The total number of sequences in the alignment
	 */
	DenseAbstractMatch(const uint seq_count );
	// use the compiler generated copy constructor, assignment operator, and destructor

	virtual AbstractMatch* Clone() const = 0;
	
	// see AbstractMatch base class documentation for these functions

	int64 Start(uint seqI) const{
		int64 s = leftend[seqI];
		return orient[seqI]? -s : s;
	}
	void SetStart(uint seqI, int64 startI)
	{
		SetLeftEnd( seqI, genome::absolut(startI) );
		orient[seqI] = startI < 0;
	}
	uint Multiplicity() const{return m_multiplicity;}
	uint SeqCount() const{return m_seq_count;}
	virtual uint FirstStart() const;
	virtual void Invert();

	virtual gnSeqI LeftEnd(uint seqI) const{ return leftend[seqI]; }
	virtual orientation Orientation(uint seqI) const;
	virtual void SetLeftEnd(uint seqI, gnSeqI position)
	{ 
		if( position == NO_MATCH && leftend[seqI] != NO_MATCH )
			--m_multiplicity;
		else if( position != NO_MATCH && leftend[seqI] == NO_MATCH )
			++m_multiplicity;
		leftend[seqI]=position; 
	}
	virtual void SetOrientation(uint seqI, orientation o){ orient[seqI]= (o == reverse); }
	
	virtual boolean operator==( const DenseAbstractMatch& dam ) const;
	
	void MoveStart(int64 move_amount);

	void MoveEnd(int64 move_amount);

	virtual uint UsedSeq( uint seqI ) const {
		return seqI;
	}

protected:

	uint m_seq_count;
	gnSeqI leftend[ MAX_SEQS ];
	bool orient[ MAX_SEQS ];
	uint m_multiplicity;
};

template< unsigned int MAX_SEQS >
DenseAbstractMatch<MAX_SEQS>::DenseAbstractMatch() :
m_seq_count(0),
m_multiplicity(0)
{
	memset( leftend, 0, MAX_SEQS * sizeof(gnSeqI) );
	memset( orient, 0, sizeof( orient ) );
}

template< unsigned int MAX_SEQS >
DenseAbstractMatch<MAX_SEQS>::DenseAbstractMatch(const uint seq_count ) :
m_seq_count(seq_count),
m_multiplicity(0)
{
	memset( leftend, 0, MAX_SEQS * sizeof(gnSeqI) );
	memset( orient, 0, sizeof( orient ) );
}

template< unsigned int MAX_SEQS >
boolean DenseAbstractMatch<MAX_SEQS>::operator==( const DenseAbstractMatch<MAX_SEQS>& dam ) const
{
	for( uint seqI = 0; seqI < m_seq_count; ++seqI )
	{
		if( leftend[seqI] != dam.leftend[seqI] ||
			(leftend[seqI] != 0 && orient[seqI] != orient[seqI]))
			return false;
	}
	return true;
}

template< unsigned int MAX_SEQS >
AbstractMatch::orientation DenseAbstractMatch<MAX_SEQS>::Orientation(uint seqI) const
{ 
	if( leftend[seqI] != NO_MATCH && seqI < m_seq_count )
		return orient[seqI] ? reverse : forward; 
	return undefined;
}

template< unsigned int MAX_SEQS >
void DenseAbstractMatch<MAX_SEQS>::Invert()
{
	for( uint seqI = 0; seqI < MAX_SEQS; ++seqI )
		orient[seqI] = !orient[seqI];
}

template< unsigned int MAX_SEQS >
uint DenseAbstractMatch<MAX_SEQS>::FirstStart() const
{
	for( uint m_firstStart = 0; m_firstStart < SeqCount(); ++m_firstStart )
		if( leftend[m_firstStart] != NO_MATCH )
			return m_firstStart;
	return (std::numeric_limits<uint>::max)();
}

template< unsigned int MAX_SEQS >
void DenseAbstractMatch<MAX_SEQS>::MoveStart(int64 move_amount)
{
	for( uint i=0; i < m_seq_count; ++i )
		if( leftend[i] != NO_MATCH && orient[i] == false )
			leftend[i] += move_amount;
}

template< unsigned int MAX_SEQS >
void DenseAbstractMatch<MAX_SEQS>::MoveEnd(int64 move_amount)
{
	for( uint i=0; i < m_seq_count; ++i )
		if( leftend[i] != NO_MATCH && orient[i] )
			leftend[i] += move_amount;
}


typedef DenseAbstractMatch<2> DenseAbstractMatch2;
typedef DenseAbstractMatch<4> DenseAbstractMatch4;
typedef DenseAbstractMatch<8> DenseAbstractMatch8;
typedef DenseAbstractMatch<16> DenseAbstractMatch16;
typedef DenseAbstractMatch<32> DenseAbstractMatch32;
typedef DenseAbstractMatch<64> DenseAbstractMatch64;
typedef DenseAbstractMatch<128> DenseAbstractMatch128;

}

#endif // _DenseAbstractMatch_h_
