/*******************************************************************************
 * $Id: MatchProjectionAdapter.h,v 1.8 2004/02/27 23:08:55 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef __MatchProjectionAdapter_h__
#define __MatchProjectionAdapter_h__

#include "libMems/AbstractMatch.h"
#include <vector>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

namespace mems {

/**
 * MatchProjectionAdapter is a wrapper around an AbstractMatch that effectively projects a multi-match to a
 * subset match.  The adapter class forwards most function calls to the original match
 * class, to which it stores a pointer.  Use of non-const functions results in undefined state.
 */
class MatchProjectionAdapter : public mems::AbstractMatch
{
public:
	MatchProjectionAdapter() : m(NULL){};
	MatchProjectionAdapter( mems::AbstractMatch* match, const std::vector< size_t >& projection ) :
	  seq(projection)
	{
		m = match->Copy();
	}

	MatchProjectionAdapter( const MatchProjectionAdapter& mpa ) : 
	seq( mpa.seq )
	{
		if( mpa.m != NULL )
			m = mpa.m->Copy();
		else
			m = NULL;
	}

	~MatchProjectionAdapter()
	{
		if( m != NULL )
			m->Free();
	}

	MatchProjectionAdapter* Clone() const { return new MatchProjectionAdapter( *this ); }

	inline
	MatchProjectionAdapter* Copy() const
	{
		return m_allocateAndCopy( *this );
	}

	void Free()
	{
		m_free(this);
	}

	MatchProjectionAdapter& operator=( const MatchProjectionAdapter& mpa )
	{
		if( m != NULL )
			m->Free();
		m = mpa.m->Copy();
		seq = mpa.seq;
		return *this;
	}

	//
	// forward all function calls to match
	//
	gnSeqI Length( uint seqI ) const { return m->Length(seq[seqI]); }
	void SetLength( gnSeqI len, uint seqI ) { m->SetLength(len, seq[seqI]); }
	int64 Start(uint startI) const { return m->Start(seq[startI]); }
	void SetStart(uint seqI, int64 start) { m->SetStart(seq[seqI],start); }
	gnSeqI LeftEnd(uint seqI) const { return m->LeftEnd(seq[seqI]); }
	orientation Orientation(uint seqI) const { return m->Orientation(seq[seqI]); }
	void SetLeftEnd(uint seqI, gnSeqI start) { m->SetLeftEnd(seq[seqI],start); }
	void SetOrientation(uint seqI, orientation o) { m->SetOrientation(seq[seqI],o); }
	void MoveStart(int64 move_amount) { m->MoveStart(move_amount); }
	void MoveEnd(int64 move_amount) { m->MoveEnd(move_amount); }
	uint Multiplicity() const 
	{ 
		size_t mult = 0;
		for( size_t projI = 0; projI < seq.size(); projI++ )
			if( m->LeftEnd(projI) != mems::NO_MATCH )
				++mult;
		return mult; 
	}
	uint SeqCount() const { return seq.size(); }
	uint FirstStart() const { return 0; }	
	gnSeqI AlignmentLength() const { return m->AlignmentLength(); }
	void Invert() { m->Invert(); }
	void CropStart(gnSeqI crop_amount) { m->CropStart(crop_amount); }
	void CropEnd(gnSeqI crop_amount) { m->CropEnd(crop_amount); }
	void CropLeft(gnSeqI crop_amount, uint seqI) { m->CropLeft(crop_amount, seq[seqI]); }
	void CropRight(gnSeqI crop_amount, uint seqI) { m->CropRight(crop_amount, seq[seqI]); }
	void GetAlignment( std::vector< mems::bitset_t >& align_matrix ) const 
	{
		std::vector< mems::bitset_t > aln_mat;
		m->GetAlignment(aln_mat);
		align_matrix.clear();
		for( size_t seqI = 0; seqI < seq.size(); ++seqI )
			align_matrix.push_back(aln_mat[seq[seqI]]);
	}
	void GetColumn( gnSeqI col, std::vector<gnSeqI>& pos, std::vector<bool>& column ) const 
	{
		std::vector<gnSeqI> m_pos; 
		std::vector<bool> m_column;
		m->GetColumn(col,m_pos,m_column);
		pos.clear();
		for( size_t seqI = 0; seqI < seq.size(); ++seqI )
		{
			pos.push_back(m_pos[seq[seqI]]);
			column.push_back(m_column[seq[seqI]]);
		}
	}
	bool IsGap( uint seqI, gnSeqI col ) const { return m->IsGap( seq[seqI],col ); }
	uint UsedSeq( uint seqI ) const 
	{
		uint c = 0;
		for( uint i = 0; i < seq.size(); i++ )
		{
			if(m->Start(seq[i]) != 0)
				c++;
			if(c>seqI)
				return i;
		}
		return (std::numeric_limits<uint>::max)();
	};

	mems::AbstractMatch* m;
	std::vector< size_t > seq;
};

}

#endif // __MatchProjectionAdapter_h__
