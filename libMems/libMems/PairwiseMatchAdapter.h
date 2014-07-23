#ifndef __PairwiseMatchAdapter_h__
#define __PairwiseMatchAdapter_h__

#include "libMems/AbstractMatch.h"
#include "libMems/ProgressiveAligner.h"
#include <vector>

namespace mems {

/**
 * PairwiseMatchAdapter is a wrapper around an AbstractMatch that effectively projects a multi-match to a
 * pairwise match.  The adapter class forwards most function calls to the original match
 * class, to which it stores a pointer.  Use of non-const functions results in undefined state.
 */
class PairwiseMatchAdapter : public mems::AbstractMatch
{
public:
	PairwiseMatchAdapter() : m(NULL) {}
	PairwiseMatchAdapter( AbstractMatch* match, uint seq1, uint seq2 ) :
	  m(match)
	{
		seq[0] = seq1;
		seq[1] = seq2;
		inverted = false;
	}

	PairwiseMatchAdapter* Clone() const { return new PairwiseMatchAdapter( *this ); }

	PairwiseMatchAdapter* Copy() const
	{
		return m_allocateAndCopy( *this );
	}

	void Free()
	{
		m_free(this);
	}

	//
	// forward all function calls to match
	//
	gnSeqI Length( uint seqI ) const { return m->Length(seq[seqI]); }
	void SetLength( gnSeqI len, uint seqI ) { m->SetLength(len, seq[seqI]); }
	int64 Start(uint startI) const { 
		if(inverted)
			return -m->Start(seq[startI]);
		return m->Start(seq[startI]); 
	}
	void SetStart(uint seqI, int64 start) { m->SetStart(seq[seqI],start); }
	gnSeqI LeftEnd(uint seqI) const { return m->LeftEnd(seq[seqI]); }
	orientation Orientation(uint seqI) const { 
		orientation o = m->Orientation(seq[seqI]);
		if(inverted && o != AbstractMatch::undefined )
			o = o == AbstractMatch::forward ? AbstractMatch::reverse : AbstractMatch::forward; 
		return o; 
	}
	void SetLeftEnd(uint seqI, gnSeqI start) { m->SetLeftEnd(seq[seqI],start); }
	void SetOrientation(uint seqI, orientation o) { m->SetOrientation(seq[seqI],o); }
	void MoveStart(int64 move_amount) { m->MoveStart(move_amount); }
	void MoveEnd(int64 move_amount) { m->MoveEnd(move_amount); }
	uint Multiplicity() const { return 2; }
	uint SeqCount() const { return 2; }
	uint FirstStart() const { return 0; }	
	gnSeqI AlignmentLength() const { return m->AlignmentLength(); }
	void Invert() { inverted = !inverted; }
	void CropStart(gnSeqI crop_amount) { m->CropStart(crop_amount); }
	void CropEnd(gnSeqI crop_amount) { m->CropEnd(crop_amount); }
	void CropLeft(gnSeqI crop_amount, uint seqI) { m->CropLeft(crop_amount, seq[seqI]); }
	void CropRight(gnSeqI crop_amount, uint seqI) { m->CropRight(crop_amount, seq[seqI]); }
	void GetAlignment( std::vector< mems::bitset_t >& align_matrix ) const 
	{
		if( inverted )
			m->Invert();
		std::vector< mems::bitset_t > aln_mat;
		m->GetAlignment(aln_mat);
		align_matrix.clear();
		align_matrix.push_back(aln_mat[seq[0]]);
		align_matrix.push_back(aln_mat[seq[1]]);
		if( inverted )
			m->Invert();
	}
	void GetColumn( gnSeqI col, std::vector<gnSeqI>& pos, std::vector<bool>& column ) const 
	{
		if( inverted )
			m->Invert();
		std::vector<gnSeqI> m_pos; 
		std::vector<bool> m_column;
		m->GetColumn(col,m_pos,m_column);
		pos.clear();
		pos.push_back(m_pos[seq[0]]);
		pos.push_back(m_pos[seq[1]]);
		column.push_back(m_column[seq[0]]);
		column.push_back(m_column[seq[1]]);
		if( inverted )
			m->Invert();
	}

	bool IsGap( uint seqI, gnSeqI col ) const { return m->IsGap( seq[seqI],col ); }
	uint UsedSeq( uint seqI ) const 
	{
		if(m->Start(seq[0]) != NO_MATCH)
			return 0;
		if(m->Start(seq[1]) != NO_MATCH)
			return 1;
		return (std::numeric_limits<uint>::max)();
	};

	AbstractMatch* m;
	TrackingMatch* tm;
	uint seq[2];
	bool inverted;
};

}

#endif	// __PairwiseMatchAdapter_h__  

