/*******************************************************************************
 * $Id: AbstractMatch.h,v 1.8 2004/02/27 23:08:55 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef __AbstractMatch_h__
#define __AbstractMatch_h__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnClone.h"
#include <vector>
#include <algorithm>
#include <boost/type_traits/remove_pointer.hpp>
#include <boost/type_traits/add_pointer.hpp>
#include <boost/dynamic_bitset.hpp>
#include <libMems/SlotAllocator.h>
#include <libMems/configuration.h>

namespace mems {

static const gnSeqI NO_MATCH = 0;


#ifdef WIN32 
/** define this to force all matches to use boost allocators instead of new/delete */
//#define _USE_BOOST_MATCH_ALLOCATOR
//typedef boost::dynamic_bitset<unsigned, boost::pool_allocator<unsigned> > bitset_t;

// slot allocator turns out to have the fastest new/free implementation for single object allocations
#define _USE_SLOT_ALLOCATOR
#else
#define _USE_SLOT_ALLOCATOR
#endif
typedef boost::dynamic_bitset<> bitset_t;

#ifdef _USE_SLOT_ALLOCATOR
#include "libMems/SlotAllocator.h"
#elif defined(_USE_BOOST_MATCH_ALLOCATOR)
#include <boost/pool/pool_alloc.hpp>
#endif

template< typename T >
T* m_allocateAndCopy( const T& t )
{
#ifdef _USE_SLOT_ALLOCATOR
	SlotAllocator<T>& sat = SlotAllocator<T>::GetSlotAllocator();
	T* newt = sat.Allocate();
	newt = new(newt) T(t);	// construct a new T at the address given by newt
//	*newt = t;
	return newt;
#elif defined(_USE_BOOST_MATCH_ALLOCATOR)
	boost::fast_pool_allocator< T > fpa;
	T* newt = boost::fast_pool_allocator< T >::allocate();
	fpa.construct(newt, t);
	return newt;
#else
	return new T(t);
#endif
}

template< typename T >
void m_free( T* t )
{
#ifdef _USE_SLOT_ALLOCATOR
	SlotAllocator<T>& sat = SlotAllocator<T>::GetSlotAllocator();
	sat.Free(t);
#elif defined(_USE_BOOST_MATCH_ALLOCATOR)
	boost::fast_pool_allocator< T > fpa;
	fpa.destroy(t);
	boost::fast_pool_allocator< T >::deallocate(t);
#else
	delete t;
#endif
}

/**
 * AbstractMatch is a pure virtual base class that defines an interface for 
 * both gapped and ungapped alignments among several sequences or several regions
 * of the same sequence 
 */
class AbstractMatch : public genome::gnClone {
public:
	
	enum orientation {
		forward,	/**< the alignment is on the forward strand */
		reverse,	/**< alignment on the reverse strand */
		undefined	/**< there is no alignment on either strand */
	};

	/** creates a copy of this using a boost::pool::fast_pool_allocator */
	virtual AbstractMatch* Copy() const = 0;

	/** frees storage used by this object in a boost::pool::fast_pool_allocator */
	virtual void Free() = 0;
	
	/** Returns the length of this match */
	virtual gnSeqI Length( uint seqI ) const = 0;

	/** Sets the length of this match to @param len */
	virtual void SetLength( gnSeqI len, uint seqI ) = 0;

	/** Deprecated:  use LeftEnd() and Orientation() instead.
	 * Returns the start coordinate of this match in sequence @param startI */
	virtual int64 Start(uint startI) const = 0;

	/** Deprecated: use SetLeftEnd() and SetOrientation instead
	 * Sets the start in sequence @param seqI of this match to @param start */
	virtual void SetStart(uint seqI, int64 start) = 0;

	/** Deprecated: use LeftEnd() instead
	 * Returns the start coordinate of this match in sequence @param seqI */
	int64 operator[](uint seqI) const{return Start(seqI);}	// this is a synonym for Start()

	/** Deprecated: use RightEnd() instead
	 * Returns the last coordinate of this match in sequence @param seqI */
	virtual int64 End(uint seqI) const;

	/** Returns the left end coordinate of this match at the seqI'th matching position/sequence */
	virtual gnSeqI LeftEnd(uint seqI) const = 0;

	/** Returns the right-end coordinate of this match at the seqI'th matching position/sequence 
	    (equal to LeftEnd(seqI) + Length(seqI) - 1) */
	virtual gnSeqI RightEnd(uint seqI) const{ return LeftEnd(seqI) + Length( seqI ) - 1; };

	/** Returns the orientation of this match at the startI'th matching position/sequence, 
	 *  either AbstractMatch::forward or AbstractMatch::reverse 
	 */
	virtual orientation Orientation(uint seqI) const = 0;

	/** sets the left end coordinate of this match in the seqI'th matching position/sequence */
	virtual void SetLeftEnd(uint seqI, gnSeqI start) = 0;

	/** sets the relative orientation of this match in the seqI'th matching position/sequence */
	virtual void SetOrientation(uint seqI, orientation o) = 0;

	/** Shift the left-end coordinates in forward oriented positions by a given amount */
	virtual void MoveStart(int64 move_amount) = 0;
	/** Shift the left-end coordinates  in reverse oriented positions by a given amount */
	virtual void MoveEnd(int64 move_amount) = 0;

	/** Returns the multiplicity of the match.  e.g. the number of sequences this match occurs in */
	virtual uint Multiplicity() const = 0;

	/** Returns the number of sequences in the alignment which contains this match */
	virtual uint SeqCount() const = 0;

	/** Returns the index of the first sequence this match occurs in */
	virtual uint FirstStart() const = 0;
	
	/** Returns the total length of this alignment in columns */
	virtual gnSeqI AlignmentLength() const = 0;

	/** Inverts the coordinates of this match */
	virtual void Invert() = 0;
	
	//warning:  none of the following do bounds checking.
	/** 
	 * Deprecated:  Use CropLeft and CropRight instead
	 * Removes the first <code>crop_amount</code> base pairs from the beginning of the match.
	 */
	virtual void CropStart(gnSeqI crop_amount) = 0;
	/** 
	 * Deprecated:  Use CropLeft and CropRight instead
	 * Removes the last <code>crop_amount</code> base pairs from the end of the match.
	 */
	virtual void CropEnd(gnSeqI crop_amount) = 0;

	/**
	 * Crop this match from the left
	 * Removes the first <code>crop_amount</code> positions from the left side of the match.
	 */
	virtual void CropLeft(gnSeqI crop_amount, uint seqI) = 0;
	/**
	 * Crop this match from the right
	 * Removes the last <code>crop_amount</code> positions from the right side of the match.
	 */
	virtual void CropRight(gnSeqI crop_amount, uint seqI) = 0;
	
//	virtual AbstractMatch* Split( gnSeqI before_column ) = 0;

	/**
	 * Gets a copy of the alignment as an array of dynamic_bitsets
	 */
	virtual void GetAlignment( std::vector< bitset_t >& align_matrix ) const = 0;

	/** Given an alignment column index, this function returns the corresponding sequence coordinates
	 *  and whether each sequence is aligned in that column 
	 *  If a given sequence is not represented in the requested column, the position returned 
	 *  in pos should be that of the first nucleotide to the left of the requested column.  If no
	 *  nucleotides exist to the left of the requested column, then a NO_MATCH is returned in pos
	 *  for that sequence.
	 */
	virtual void GetColumn( gnSeqI col, std::vector<gnSeqI>& pos, std::vector<bool>& column ) const = 0;

//	gnSeqI SeqPosToColumn( uint seq, int64 pos) const = 0;
	/** returns true if the given row,column of the alignment has a gap character */
	virtual bool IsGap( uint seq, gnSeqI col ) const = 0;
	/** Returns the id of the i-th defined sequence in this match */ 
	virtual uint UsedSeq( uint seqI ) const = 0;
};

inline
int64 AbstractMatch::End(uint endI) const
{
	if( Start(endI) > 0 )
		return Start(endI) + Length(endI) - 1;
	return Start(endI);
}


template< typename MatchType >
class AbstractMatchStartComparator {
public:
	AbstractMatchStartComparator( unsigned seq = 0 ){
		m_seq = seq;
	}
	AbstractMatchStartComparator( const AbstractMatchStartComparator& msc ){
		m_seq = msc.m_seq;
	}
	AbstractMatchStartComparator<MatchType>& operator=( const AbstractMatchStartComparator<MatchType>& msc )
	{
		m_seq = msc.m_seq;
	}
	// TODO??  make this do a wraparound comparison if all is equal?
	boolean operator()(const MatchType& a, const MatchType& b) const{
		int start_diff = std::max( a.FirstStart(), m_seq ) - std::max( a.FirstStart(), m_seq );
		if(start_diff == 0){
			uint m_count = a.SeqCount();
			m_count = m_count <= b.SeqCount() ? m_count : b.SeqCount();
			for(uint seqI = m_seq; seqI < m_count; seqI++){
				gnSeqI a_start = a.Orientation(seqI) == AbstractMatch::forward ? a.LeftEnd( seqI ) : a.RightEnd( seqI );
				gnSeqI b_start = b.Orientation(seqI) == AbstractMatch::forward ? b.LeftEnd( seqI ) : b.RightEnd( seqI );
				int64 diff = a_start - b_start;
				if(a_start == NO_MATCH || b_start == NO_MATCH)
					continue;
				else if(a_start == b_start)
					continue;
				else
					return a_start < b_start;
			}
		}
		return start_diff < 0;
	}
private:
	unsigned m_seq;
};

template< typename MatchType >
class AbstractMatchSingleStartComparator {
public:
	AbstractMatchSingleStartComparator( unsigned seq = 0 ){
		m_seq = seq;
	}
	AbstractMatchSingleStartComparator( const AbstractMatchSingleStartComparator& msc ){
		m_seq = msc.m_seq;
	}
	AbstractMatchSingleStartComparator<MatchType>& operator=( const AbstractMatchSingleStartComparator<MatchType>& msc )
	{
		m_seq = msc.m_seq;
	}
	/**
	 * Compare on only one sequence.  Undefined matches are less than defined matches
	 */
	boolean operator()(const MatchType& a, const MatchType& b) const{
		int64 a_start = a.LeftEnd( m_seq ), b_start = b.LeftEnd( m_seq );
		if( a_start == NO_MATCH || b_start == NO_MATCH ){
			if( b_start != NO_MATCH )
				return true;
			return false;
		}

		return a_start < b_start;
	}
private:
	unsigned m_seq;
};



template< typename MatchType >
class MatchStartComparator {
public:
	MatchStartComparator( unsigned seq = 0 ){
		m_seq = seq;
	}
	MatchStartComparator( const MatchStartComparator& msc ){
		m_seq = msc.m_seq;
	}
	MatchStartComparator<MatchType>& operator=( const MatchStartComparator<MatchType>& msc )
	{
		m_seq = msc.m_seq;
	}
	// TODO??  make this do a wraparound comparison if all is equal?
	boolean operator()(const MatchType* a, const MatchType* b) const{
		int start_diff = std::max( a->FirstStart(), m_seq ) - std::max( a->FirstStart(), m_seq );
		if(start_diff == 0){
			uint m_count = a->SeqCount();
			m_count = m_count <= b->SeqCount() ? m_count : b->SeqCount();
			for(uint seqI = m_seq; seqI < m_count; seqI++){
				gnSeqI a_start = a->Orientation(seqI) == AbstractMatch::forward ? a->LeftEnd( seqI ) : a->RightEnd( seqI );
				gnSeqI b_start = b->Orientation(seqI) == AbstractMatch::forward ? b->LeftEnd( seqI ) : b->RightEnd( seqI );
				int64 diff = a_start - b_start;
				if(a_start == NO_MATCH || b_start == NO_MATCH)
					continue;
				else if(a_start == b_start)
					continue;
				else
					return a_start < b_start;
			}
		}
		return start_diff < 0;
	}
private:
	unsigned m_seq;
};

template< typename MatchType >
class SingleStartComparator {
public:
	SingleStartComparator( unsigned seq = 0 ){
		m_seq = seq;
	}
	SingleStartComparator( const SingleStartComparator& msc ){
		m_seq = msc.m_seq;
	}
	SingleStartComparator<MatchType>& operator=( const SingleStartComparator<MatchType>& msc )
	{
		m_seq = msc.m_seq;
	}
	/**
	 * Compare on only one sequence.  Undefined matches are less than defined matches
	 */
	boolean operator()(const MatchType* a, const MatchType* b) const{
		int64 a_start = a->LeftEnd( m_seq ), b_start = b->LeftEnd( m_seq );
		if( a_start == NO_MATCH || b_start == NO_MATCH ){
			if( b_start != NO_MATCH )
				return true;
			return false;
		}

		return a_start < b_start;
	}
private:
	unsigned m_seq;
};


template< typename MatchType >
class SSC {
public:
	SSC( unsigned seq = 0 ){
		m_seq = seq;
	}
	SSC( const SSC<MatchType>& msc ){
		m_seq = msc.m_seq;
	}
	SSC<MatchType>& operator=( const SSC<MatchType>& msc )
	{
		m_seq = msc.m_seq;
	}
	boolean operator()( const typename boost::add_pointer<MatchType>::type& a, 
		const typename boost::add_pointer<MatchType>::type& b) const
	{
		return operator()(*a,*b);
	}
	/**
	 * Compare on only one sequence.  Undefined matches are less than defined matches
	 */
	boolean operator()(const typename boost::remove_pointer<MatchType>::type& a, 
		const typename boost::remove_pointer<MatchType>::type& b) const{
		int64 a_start = a.LeftEnd( m_seq ), b_start = b.LeftEnd( m_seq );
		if( a_start == NO_MATCH || b_start == NO_MATCH ){
			if( b_start != NO_MATCH )
				return true;
			return false;
		}

		return a_start < b_start;
	}
private:
	unsigned m_seq;
};

}

#endif // __AbstractMatch_h__
