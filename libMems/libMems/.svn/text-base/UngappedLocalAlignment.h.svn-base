/*******************************************************************************
 * $Id: UngappedLocalAlignment.h,v 1.10 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef __UngappedLocalAlignment_h__
#define __UngappedLocalAlignment_h__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnClone.h"
#include "libGenome/gnException.h"
#include "libMems/AbstractMatch.h"

namespace mems {

/**
 * The UngappedLocalAlignment class stores the location of an <b>equal size</b> (inexact or exactly) 
 * matching region between several sequences.  This class can use one of several storage schemes
 * such as DenseAbstractMatch or SparseAbstractMatch
 */
template< class AbstractMatchImpl >
class UngappedLocalAlignment : public AbstractMatchImpl 
{

public:
	UngappedLocalAlignment();
	/**
	 * Creates a new UngappedLocalAlignment.
	 * @param seq_count The total number of sequences in the alignment
	 */
	UngappedLocalAlignment( const uint seq_count );

	// use trivial copy constructor, destructor, and operator =

	UngappedLocalAlignment* Clone() const;
	UngappedLocalAlignment* Copy() const;
	virtual void Free();

	/** comparison operator, compares two UngappedLocalAlignmentes to see if they are the same */
	boolean operator==(const UngappedLocalAlignment& mhe) const;

	gnSeqI Length( uint seqI = (std::numeric_limits<uint>::max)() ) const
	{
		if( seqI == (std::numeric_limits<uint>::max)() ) 
			return m_length;
		if( this->LeftEnd(seqI) == NO_MATCH )
			return 0;
		return m_length;
	}
	void SetLength(gnSeqI len, uint seqI = 0){m_length = len;}
	gnSeqI AlignmentLength() const{return m_length;}
	
	//warning:  none of the following do bounds checking.
	virtual void Move( int64 distance );
	virtual void CropStart(gnSeqI crop_amount);
	virtual void CropEnd(gnSeqI crop_amount);
	virtual void ExtendStart(gnSeqI extend_amount);
	virtual void ExtendEnd(gnSeqI extend_amount);

	virtual void CropLeft(gnSeqI crop_amount, uint seqI);
	virtual void CropRight(gnSeqI crop_amount, uint seqI);

	void GetAlignment( std::vector< bitset_t >& align_matrix ) const;

	void GetColumn( gnSeqI col, std::vector<gnSeqI>& pos, std::vector<bool>& column ) const;

	/**
	 * Writes the location of this UngappedLocalAlignment to the specified output stream (e.g. cout).
	 */
	template<typename AMImpl> friend std::ostream& operator<<(std::ostream& os, const UngappedLocalAlignment<AMImpl>& ula); //write to source.

	bool IsGap( uint seqI, gnSeqI col ) const {
		return (this->LeftEnd(seqI) != NO_MATCH && col < m_length);
	}

protected:

	gnSeqI m_length;
};


template< class AbstractMatchImpl >
UngappedLocalAlignment< AbstractMatchImpl >::UngappedLocalAlignment() : AbstractMatchImpl()
{
}


template< class AbstractMatchImpl >
UngappedLocalAlignment< AbstractMatchImpl >::UngappedLocalAlignment(uint seq_count)
 : AbstractMatchImpl( seq_count )
{
}


template< class AbstractMatchImpl >
UngappedLocalAlignment< AbstractMatchImpl >* 
UngappedLocalAlignment< AbstractMatchImpl >::Clone() const
{
	return new UngappedLocalAlignment(*this);
}

template< class AbstractMatchImpl >
UngappedLocalAlignment<AbstractMatchImpl>* UngappedLocalAlignment<AbstractMatchImpl>::Copy() const
{
	return m_allocateAndCopy( *this );
}
template< class AbstractMatchImpl >
void UngappedLocalAlignment<AbstractMatchImpl>::Free()
{
	m_free(this);
}

template< class AbstractMatchImpl >
boolean UngappedLocalAlignment<AbstractMatchImpl>::operator==(const UngappedLocalAlignment& ula) const
{
	if(m_length != ula.m_length)
		return false;
	return AbstractMatchImpl::operator==(ula);
}

template< class AbstractMatchImpl >
void UngappedLocalAlignment<AbstractMatchImpl>::Move( int64 distance )
{
	for( uint32 i=0; i < AbstractMatchImpl::SeqCount(); i++ ){
		int64 start = AbstractMatchImpl::Start(i);
		if( start != NO_MATCH )
			AbstractMatchImpl::SetStart(i, start + distance );
	}
}

template< class AbstractMatchImpl >
void UngappedLocalAlignment<AbstractMatchImpl>::CropStart(gnSeqI crop_amount)
{
	if( crop_amount > m_length )
		Throw_gnEx( genome::SeqIndexOutOfBounds() );
	m_length -= crop_amount;
	AbstractMatchImpl::MoveStart(crop_amount);
}

template< class AbstractMatchImpl >
void UngappedLocalAlignment<AbstractMatchImpl>::CropEnd(gnSeqI crop_amount){
	if( crop_amount > m_length )
		Throw_gnEx( genome::SeqIndexOutOfBounds() );
	m_length -= crop_amount;
	AbstractMatchImpl::MoveEnd(crop_amount);
}

template< class AbstractMatchImpl >
void UngappedLocalAlignment<AbstractMatchImpl>::ExtendStart(gnSeqI extend_amount){
	m_length += extend_amount;
	int64 amt = extend_amount;
	AbstractMatchImpl::MoveStart(-amt);
}

template< class AbstractMatchImpl >
void UngappedLocalAlignment<AbstractMatchImpl>::ExtendEnd(gnSeqI extend_amount){
	m_length += extend_amount;
	int64 amt = extend_amount;
	AbstractMatchImpl::MoveEnd(-amt);
}

template< class AbstractMatchImpl >
void UngappedLocalAlignment<AbstractMatchImpl>::CropLeft(gnSeqI crop_amount, uint seqI)
{
	if(AbstractMatchImpl::Orientation(seqI) == AbstractMatch::forward)
		CropStart(crop_amount);
	else
		CropEnd(crop_amount);
}

template< class AbstractMatchImpl >
void UngappedLocalAlignment<AbstractMatchImpl>::CropRight(gnSeqI crop_amount, uint seqI)
{
	if(AbstractMatchImpl::Orientation(seqI) == AbstractMatch::forward)
		CropEnd(crop_amount);
	else
		CropStart(crop_amount);
}

template< class AbstractMatchImpl >
void UngappedLocalAlignment< AbstractMatchImpl >::GetAlignment( std::vector< bitset_t >& align_matrix ) const
{
	align_matrix = std::vector< bitset_t >(this->SeqCount(), bitset_t( this->AlignmentLength(), false ) );
	for( uint seqI = 0; seqI < this->SeqCount(); seqI++ )
	{
		if( this->LeftEnd(seqI) != NO_MATCH )
			align_matrix[seqI].flip();
	}
}

template< typename AbstractMatchImpl >
std::ostream& operator<<(std::ostream& os, const UngappedLocalAlignment< AbstractMatchImpl >& ula);

template< typename AbstractMatchImpl >
std::ostream& operator<<(std::ostream& os, const UngappedLocalAlignment< AbstractMatchImpl >& ula){ //write to stream.
	os << ula.m_length;
	for(uint i=0; i < ula.SeqCount(); i++)
		os << '\t' << ula.Start(i);
	return os;
}


template< class AbstractMatchImpl >
void UngappedLocalAlignment< AbstractMatchImpl >::GetColumn( gnSeqI col, std::vector<gnSeqI>& pos, std::vector<bool>& column ) const
{
	pos = std::vector<gnSeqI>(this->SeqCount(), NO_MATCH);
	column = std::vector<bool>(this->SeqCount(), true);
	for( uint seqI = 0; seqI < this->SeqCount(); seqI++ )
	{
		if( this->Orientation(seqI) == AbstractMatch::forward )
			pos[seqI] = this->LeftEnd(seqI) + col;
		else if( this->Orientation(seqI) == AbstractMatch::reverse )
			pos[seqI] = this->RightEnd(seqI) - col;
		else
			column[seqI] = false;
	}
}

}

#endif // _UngappedLocalAlignment_h_
