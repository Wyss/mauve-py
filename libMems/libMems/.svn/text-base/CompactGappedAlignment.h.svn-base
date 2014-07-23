/*******************************************************************************
 * $Id: CompactGappedAlignment.h,v 1.12 2004/04/19 23:10:50 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef __CompactGappedAlignment_h__
#define __CompactGappedAlignment_h__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnDebug.h"
#include "libGenome/gnFilter.h"
#include "libGenome/gnSequence.h"
#include "libMems/SparseAbstractMatch.h"
#include "libMems/HybridAbstractMatch.h"
#include "libMems/AbstractGappedAlignment.h"
#include "libMems/UngappedLocalAlignment.h"

#include <algorithm>

#ifdef WIN32
#include "windows.h"
#endif

namespace mems {

/**
 * The CompactGappedAlignment stores a gapped alignment as a bit-vector
 * Rather than using one byte per aligned position, this class uses one bit, making
 * particularly space efficient
 */
template< class BaseType = AbstractGappedAlignment< HybridAbstractMatch<> > >
class CompactGappedAlignment : public BaseType
{
public:
	CompactGappedAlignment() : BaseType(){};
	CompactGappedAlignment( uint seq_count, gnSeqI align_length );
	CompactGappedAlignment( std::vector< bitset_t >& aln_mat, gnSeqI alignment_length );
	
	template< class MatchType >
	CompactGappedAlignment( MatchType& m ) : 
		BaseType( m.SeqCount(), m.AlignmentLength() ),
		bcount( std::vector< std::vector< size_t > >( m.SeqCount() ) )
	{
		m.GetAlignment(align_matrix);

		for( uint seqI = 0; seqI < this->SeqCount(); seqI++ )
		{
			this->SetStart(seqI, m.Start(seqI));
			if( m.Start(seqI) != NO_MATCH )
				this->SetLength(m.Length(seqI), seqI);
			else
				this->SetLength(0, seqI);
		}

		this->create_bitcount();

		if( !this->validate() )
			std::cerr << "kahnstruct error\n";
	}

	CompactGappedAlignment* Clone() const { return new CompactGappedAlignment( *this ); }
	CompactGappedAlignment* Copy() const;
	virtual void Free();
	
	void SetAlignment( const std::vector< std::string >& seq_align );

	void SetAlignment( std::vector< bitset_t >& seq_align );

	// Inherited methods from AbstractMatch:
	virtual void Invert();
	virtual void CropStart(gnSeqI crop_amount);
	virtual void CropEnd(gnSeqI crop_amount);

	virtual void CropLeft(gnSeqI crop_amount, uint seqI);
	virtual void CropRight(gnSeqI crop_amount, uint seqI);

	void GetAlignment( std::vector< bitset_t >& align_matrix ) const;

	/** allows a peek at the data inside this alignment.  don't change it or the CompactGappedAlignment will become corrupt */
	const std::vector< bitset_t >& GetAlignment() const{ return align_matrix; }

//	friend void GetAlignment( const CompactGappedAlignment& ga, const std::vector< genome::gnSequence* >& seq_table, std::vector<std::string>& alignment );
	
	void GetColumn( gnSeqI col, std::vector<gnSeqI>& pos, std::vector<bool>& column ) const;

	/** returns true if the given row,column of the alignment has a gap character */
	virtual bool IsGap( uint seq, gnSeqI col ) const;
	/** translate a cga to a new coordinate system */
	void translate( CompactGappedAlignment& cga, uint cga_seq, uint my_seq, bool add_bits = true );

	bool validate() const;
	bool validate_bitcount() const;

	void copyRange( CompactGappedAlignment& dest, gnSeqI left_column, gnSeqI length );
	gnSeqI SeqPosToColumn( uint seq, int64 pos);

	/** Eliminates any columns that contain only gap characters */
	void CondenseGapColumns();

	void swap( CompactGappedAlignment& other ){ swap(&other); }

protected:
	// for use by derived classes in order to swap contents
	void swap( CompactGappedAlignment* other ){
		std::swap( align_matrix, other->align_matrix );
		std::swap( bcount, other->bcount );
		BaseType::swap( other );
	}

	std::vector< bitset_t > align_matrix;		/**< aligned positions have true values, gaps are false */
	std::vector< std::vector< size_t > > bcount;

	void create_bitcount();
	gnSeqI SeqPosToColumn( gnSeqI pos, const bitset_t& bvec, const std::vector< size_t >& index ) const;

};

static bool debug_cga = false;

template< class BaseType >
CompactGappedAlignment<BaseType>* CompactGappedAlignment<BaseType>::Copy() const
{
	return m_allocateAndCopy( *this );
}

template< class BaseType >
void CompactGappedAlignment<BaseType>::Free()
{
	m_free(this);
}

template< class BaseType >
bool CompactGappedAlignment<BaseType>::validate() const
{
	if( !debug_cga )
		return true;
	bool good = true;
	for( uint seqI = 0; seqI < this->SeqCount(); seqI++ )
	{
		if( this->AlignmentLength() != align_matrix[seqI].size() )
		{
			good = false;
			std::cerr << "vanishing pig trick\n";
			genome::breakHere();
		}
		gnSeqI count = align_matrix[seqI].count();
		if( count > 0 && this->LeftEnd(seqI) == 0 )
		{
			good = false;
			std::cerr << "boner_McHoserknob\n";
			genome::breakHere();
		}
		if( (count == 0 || this->Length(seqI) == 0) && this->LeftEnd(seqI) != 0 )
		{
			good = false;
			std::cerr << "Length(" << seqI << "): " << this->Length(seqI) << std::endl;
			std::cerr << "LeftEnd(seqI): " << this->LeftEnd(seqI) << std::endl;
			std::cerr << "spumante explosion\n";
			genome::breakHere();
		}
		if( count != this->Length(seqI) )
		{
			std::cerr << "seqI: " << seqI << " count: " << count << "  Length(seqI): " << this->Length(seqI) << std::endl;
			std::cerr << "LeftEnd(seqI): " << this->LeftEnd(seqI) << std::endl;
			std::cerr << "lendo mismatcho\n";
			genome::breakHere();
			return false;
		}
//		std::vector< std::vector< size_t > > tmp_bcount = bcount;
//		create_bitcount();
//		if( !tmp_bcount == bcount )
//		{
//			good = false;
//			std::cerr << "bcount mismatch!!!\n";
//		}
//		bcount = tmp_bcount;

	}
	if( good )	// check for all gap cols
	{
/*	allow gap cols...
		for( size_t colI = 0; colI < this->AlignmentLength(); ++colI )
		{
			bool aa = false;
			for( uint seqI = 0; seqI < this->SeqCount(); seqI++ )
				aa = aa || align_matrix[seqI].test(colI);
			if( aa == false )
			{
				std::cerr << "gap col at " << colI << std::endl;
				genome::breakHere();
			}
		}
		*/
	}
	return 	validate_bitcount() && good;
}


template< class BaseType >
CompactGappedAlignment<BaseType>::CompactGappedAlignment( std::vector< bitset_t >& aln_mat, gnSeqI alignment_length ) :
BaseType( aln_mat.size(), alignment_length ),
align_matrix( aln_mat ),
bcount( std::vector< std::vector< size_t > >( aln_mat.size() ) )
{
	this->create_bitcount();
	this->validate_bitcount();
}

template< class BaseType >
CompactGappedAlignment<BaseType>::CompactGappedAlignment( uint seq_count, gnSeqI align_length ) : 
BaseType( seq_count, align_length )
{}


template< class BaseType >
void CompactGappedAlignment<BaseType>::SetAlignment( const std::vector< std::string >& seq_align ){
	if( seq_align.size() == 0 )
	{
		this->SetAlignmentLength(0);
		return;
	}
	this->SetAlignmentLength(seq_align[0].size());
	align_matrix = std::vector< bitset_t >( seq_align.size(), bitset_t( seq_align[0].size(), false ) );
	bcount = std::vector< std::vector<size_t> >( seq_align.size() );
	for( size_t seqI = 0; seqI < seq_align.size(); seqI++ )
	{
		bool nonzero = false;
		for( size_t charI = 0; charI < seq_align[seqI].size(); charI++ )
			if( seq_align[seqI][charI] != '-' )
			{
				align_matrix[seqI].set(charI);
				nonzero = true;
			}
	}
	this->create_bitcount();
}

template< class BaseType >
void CompactGappedAlignment<BaseType>::SetAlignment( std::vector< bitset_t >& seq_align )
{
	std::swap( align_matrix, seq_align );
	seq_align.clear();
	if( align_matrix.size() > 0 )
		this->SetAlignmentLength( align_matrix[0].size() );
	else
		this->SetAlignmentLength(0);
	bcount = std::vector< std::vector<size_t> >(align_matrix.size());
	this->create_bitcount();
	this->validate_bitcount();
}

template< class BaseType >
void CompactGappedAlignment<BaseType>::GetAlignment( std::vector< bitset_t >& align_matrix ) const
{
	align_matrix = this->align_matrix;
}

template< class BaseType >
bool CompactGappedAlignment<BaseType>::IsGap( uint seq, gnSeqI col ) const
{
	return !align_matrix[seq][col];
}

static const unsigned INDEX_INTERVAL = 512;

template< class BaseType >
void CompactGappedAlignment<BaseType>::create_bitcount()
{
	for( uint seqI = 0; seqI < this->SeqCount(); seqI++ )
	{
//		if( this->LeftEnd(seqI) == NO_MATCH )
//			continue;
		bitset_t& bvec = align_matrix[seqI];
		bcount[seqI].clear();
		bcount[seqI].push_back(0);
		for( size_t indie = 0; indie + INDEX_INTERVAL <= bvec.size(); indie += INDEX_INTERVAL )
		{
			size_t end = indie + INDEX_INTERVAL;
			size_t ct = 0;
			for( size_t i = indie; i < end; ++i )
				ct += bvec.test(i);
			bcount[seqI].push_back( ct + bcount[seqI].back() );
		}
	}
}

template< class BaseType >
bool CompactGappedAlignment<BaseType>::validate_bitcount() const
{
	if( !debug_cga )
		return true;
	bool valid = true;	// innocent until proven guilty
	for( uint seqI = 0; seqI < this->SeqCount(); seqI++ )
	{
		gnSeqI count = align_matrix[seqI].count();
		size_t bc_len = align_matrix[seqI].size() / INDEX_INTERVAL;
		if( count < INDEX_INTERVAL && bcount[seqI].size() == 0 )
			continue;	// a-ok here
		if( bc_len + 1 != bcount[seqI].size() && (bcount[seqI].back() % INDEX_INTERVAL != 0) )
		{
			std::cerr << "bitcount problem, bc_len + 1: " << bc_len + 1 << " and bcount[seqI].size(): " << bcount[seqI].size() << std::endl;
			std::cerr << "count: " << count << " and bcount[seqI].back(): " << bcount[seqI].back() << std::endl;
			valid = false;
		}
		if( count - bcount[seqI].back() > INDEX_INTERVAL )
		{
			std::cerr << "bitcount problem, count: " << count << " and bcount[seqI].back(): " << bcount[seqI].back() << std::endl;
			valid = false;
		}
	}
	return valid;
}

template< class BaseType > 
gnSeqI CompactGappedAlignment<BaseType>::SeqPosToColumn( uint seq, int64 pos )
{
	if( this->Orientation(seq) == AbstractMatch::forward )
		pos = genome::absolut(pos) - this->LeftEnd(seq) + 1;
	else
		pos = this->RightEnd(seq)-genome::absolut(pos) + 1;	// is this right?
	return SeqPosToColumn(pos, align_matrix[seq], bcount[seq]);
}

template< class BaseType > 
gnSeqI CompactGappedAlignment<BaseType>::SeqPosToColumn( gnSeqI pos, const bitset_t& bvec, const std::vector< size_t >& index ) const
{
	std::vector<size_t>::const_iterator iter = std::lower_bound(index.begin(), index.end(), pos);
	--iter;
	size_t cur_pos = *iter;
	size_t col = iter - index.begin();
	col *= INDEX_INTERVAL;
	if( col == 0 )
		col = bvec.find_first();
	else
		col = bvec.find_next(col-1);
	for( ++cur_pos; cur_pos < pos; ++cur_pos )
		col = bvec.find_next(col);
	return col;
}

template< class BaseType >
void CompactGappedAlignment<BaseType>::translate( CompactGappedAlignment& cga, uint cga_seq, uint my_seq, bool add_bits ) // const
{
	AbstractMatch::orientation my_orient = this->Orientation(my_seq);

	if( cga.Length(cga_seq)  > this->Length(my_seq) )
	{
		std::cerr << "Oh scheisskopf.  What are you trying to do to me??\n";
		std::cerr << "cga.Length(" << cga_seq << "): " << cga.Length(cga_seq) << std::endl;
		std::cerr << "Length(" << my_seq << "): " << this->Length(my_seq) << std::endl;
		genome::breakHere();
	}

	gnSeqI prev_lend = cga.LeftEnd(cga_seq);
	gnSeqI prev_len = cga.Length(cga_seq);
	gnSeqI my_lend = this->LeftEnd(my_seq);
	gnSeqI my_len = this->Length(my_seq);
	gnSeqI my_count = 0;
	uint seqI = 0;

	// what assumptions should be made about cga?
	// does it already have the correct left-end relative to this?
	// no, it needs to have a left-end relative to the first aligned char in this
	size_t cur_bit = 0;

	// determine left_bit
	size_t left_bit = this->SeqPosToColumn(cga.LeftEnd(cga_seq), align_matrix[my_seq], bcount[my_seq]);
	// determine right_bit
	size_t right_bit = this->SeqPosToColumn(cga.RightEnd(cga_seq), align_matrix[my_seq], bcount[my_seq]);
	if( right_bit > 4000000000u )
	{
		std::cerr << "cga doesn't fit\n";
		std::cerr << "cga.RightEnd(cga_seq) " << cga.RightEnd(cga_seq) << std::endl;
		std::cerr << "RightEnd(my_seq): " << this->RightEnd(my_seq) << std::endl;
		std::cerr << "cga.LeftEnd(cga_seq) " << cga.LeftEnd(cga_seq) << std::endl;
		std::cerr << "LeftEnd(my_seq): " << this->LeftEnd(my_seq) << std::endl;
		std::cerr << "cga.AlignmentLength(): " << cga.AlignmentLength() << std::endl;
		std::cerr << "AlignmentLength(): " << this->AlignmentLength() << std::endl;
		genome::breakHere();
	}
	right_bit++;
	if( right_bit == 0 )
		right_bit = this->AlignmentLength();

	cga.SetLeftEnd(cga_seq,left_bit+1);

	// add on length of unaligned left and right sides
	size_t cga_left = cga.align_matrix[cga_seq].find_first();

	size_t somesize = (right_bit - left_bit) - cga.Length(cga_seq) + cga.AlignmentLength();

	size_t cga_bit = cga_left;
	size_t my_bit = left_bit;
	size_t xlat_bit = cga_left;
	size_t added_bits = 0;
	// copy in everything up to cga_left
	std::vector< bitset_t > xrated( cga.SeqCount(), bitset_t( somesize, false ) );
	for( size_t seqI = 0; seqI < xrated.size(); ++seqI )
		for( size_t asdf = cga.align_matrix[seqI].find_first(); asdf < cga_left; asdf = cga.align_matrix[seqI].find_next(asdf) )
			xrated[seqI].set(asdf);
	
	while(xlat_bit < somesize)
	{
		// assume that align_matrix[my_seq][my_bit] is set
		if( !align_matrix[my_seq].test(my_bit) )
		{
			std::cerr << "ohhhhhhzheiss!\n";
			genome::breakHere();
		}
		// copy the column in cga
		for( size_t seqI = 0; seqI < xrated.size(); ++seqI )
			xrated[seqI].set( xlat_bit, cga.align_matrix[seqI].test(cga_bit) );

		++cga_bit;
		++xlat_bit;

		if( xlat_bit >= somesize )
			break;

		// TODO: should this condition be replaced by cropping xlat_bit + diff - 1 down to < somesize?
		if( cga.align_matrix[cga_seq].test(cga_bit) )
		{
			size_t next_bit = align_matrix[my_seq].find_next(my_bit);
			if( next_bit > 4000000000u )
				genome::breakHere();
			size_t diff = next_bit - my_bit;
			if( diff > 1 && add_bits )
			{
				if( xlat_bit + diff - 1 >= somesize )
				{
					std::cerr << "ERRRORRR porker!!\n";
					genome::breakHere();
				}
				for( size_t i = xlat_bit; i < xlat_bit + diff - 1; ++i )
					xrated[cga_seq].set(i);
				added_bits += diff-1;
			}
			my_bit = next_bit;
			xlat_bit += diff - 1;
		}
	}

	cga.align_matrix = xrated;
	cga.create_bitcount();
	cga.SetLength(cga.Length(cga_seq)+added_bits,cga_seq);
	cga.SetAlignmentLength(somesize);
	if( !cga.validate() )
	{
		std::cerr << "prev_lend: " << prev_lend << std::endl;
		std::cerr << "prev_len: " << prev_len << std::endl;
		std::cerr << "translate error\n";
		genome::breakHere();
	}
}


template< class BaseType >
void CompactGappedAlignment<BaseType>::Invert(){
	for(uint seqI = 0; seqI < this->SeqCount(); seqI++)
	{
		if( this->LeftEnd(seqI) == NO_MATCH )
			continue;
		bitset_t& fwd = align_matrix[seqI];
		bitset_t rev(this->AlignmentLength());
		size_t r = this->AlignmentLength();
		for( size_t i = 0; i < fwd.size(); ++i )
			rev.set( --r, fwd.test(i) );
		fwd.swap(rev);
	}
	this->create_bitcount();
	BaseType::Invert();
	if( !this->validate() )
	{
		std::cerr << "invert error\n";
	}
}

template< class BaseType >
void CompactGappedAlignment<BaseType>::CropStart(gnSeqI crop_amount){
	if( crop_amount > this->AlignmentLength() )
		Throw_gnEx( genome::SeqIndexOutOfBounds() );
	if( crop_amount == 0 )
		return;

	gnSeqI pre_alignlen = this->AlignmentLength();
	gnSeqI pre_lend0 = this->LeftEnd(0);

	std::vector<gnSeqI> pos;
	std::vector<bool> column;
	GetColumn( crop_amount-1, pos, column );

	for( uint i=0; i < this->SeqCount(); i++ ){
		if( this->LeftEnd(i) == NO_MATCH )
		{
			align_matrix[i].resize(this->AlignmentLength()-crop_amount);
			align_matrix[i] = align_matrix[i];	// force reallocation on "optimized" windows builds
			continue;
		}

		align_matrix[i] >>= crop_amount;	// why not shift left?  is this a bug in boost::dynamic_bitset?
		align_matrix[i].resize(this->AlignmentLength()-crop_amount);
		align_matrix[i] = align_matrix[i];	// force reallocation on "optimized" windows builds
		size_t char_count = this->Orientation(i) == AbstractMatch::forward ? pos[i] - this->LeftEnd(i) + 1 : this->RightEnd(i) - pos[i] + 1;

		if( pos[i] > 0 && char_count > 0 )
		{
			this->SetLength(this->Length(i)-char_count, i);
			if( this->Length(i) == 0 )
				this->SetStart(i, NO_MATCH);
			if( this->Orientation(i) == AbstractMatch::forward )
				this->SetStart(i, this->Start(i) + char_count);
		}else if( pos[i] == 0 && this->Orientation(i) == AbstractMatch::reverse )
		{
			// this sequence was completely obliterated by the crop
			this->SetLength(0, i);
			this->SetStart(i, NO_MATCH);
		}
	}

	this->SetAlignmentLength( this->AlignmentLength() - crop_amount );
	this->create_bitcount();
	if( !this->validate() )
	{
		std::cerr << "pre_lend0: " << pre_lend0 << std::endl;
		std::cerr << "pre_alignlen: " << pre_alignlen << std::endl;
		std::cerr << "CropStart error\n";
	}
}

template< class BaseType >
void CompactGappedAlignment<BaseType>::CropEnd(gnSeqI crop_amount){
	if( crop_amount > this->AlignmentLength() )
		Throw_gnEx( genome::SeqIndexOutOfBounds() );
	if( crop_amount == 0 )
		return;

	std::vector<gnSeqI> pos;
	std::vector<bool> column;
	this->GetColumn( this->AlignmentLength()-crop_amount, pos, column );

	for( uint i=0; i < this->SeqCount(); i++ ){
		align_matrix[i].resize( this->AlignmentLength() - crop_amount );
		align_matrix[i] = align_matrix[i];	// force reallocation on "optimized" windows builds
		if( this->LeftEnd(i) == NO_MATCH )
			continue;
		AbstractMatch::orientation orient = this->Orientation(i);
		if( pos[i] > 0 )
		{
			gnSeqI char_count = pos[i] - (orient == AbstractMatch::forward ? (column[i] ? 1 : 0 ) : (column[i] ? 0 : 1 ) );
			char_count = orient == AbstractMatch::forward ? char_count - this->LeftEnd(i) + 1 : this->RightEnd(i) - char_count;
			if( char_count == 0 && align_matrix[i].count() > 0)
			{
				std::cerr << "orienatation: " << (orient == AbstractMatch::forward ? "forward\n" : (orient == AbstractMatch::reverse ? "reverse\n" : "undef\n"));
				std::cerr << "lend: " << this->LeftEnd(i) << std::endl;
				std::cerr << "length: " << this->Length(i) << std::endl;
				std::cerr << "count: " << align_matrix[i].count() << std::endl;
			}
			gnSeqI deleted = this->Length(i) - char_count;
			this->SetLength(char_count, i);
			if( this->Length(i) == 0 )
				this->SetStart(i, 0);
			if( this->Start(i) < 0 )
				this->SetStart(i, this->Start(i)-deleted);
		}else if( orient == AbstractMatch::forward ){
			this->SetLength(0, i);
			this->SetStart(i, 0);
		}
	}
	SetAlignmentLength( this->AlignmentLength() - crop_amount );
	this->create_bitcount();
	if( !this->validate() )
		std::cerr << "CropEnd error\n";
}

template< class BaseType >
void CompactGappedAlignment<BaseType>::CropLeft(gnSeqI crop_amount, uint seqI)
{
	if( crop_amount == 0 )
		return;

	gnSeqI pre_len = this->Length(seqI);
	// count "crop_amount" characters into seqI and crop there
	if( this->Orientation(seqI) == AbstractMatch::forward )
	{
		size_t left_col = this->SeqPosToColumn(crop_amount, align_matrix[seqI], bcount[seqI]) + 1;
		this->CropStart(left_col);
	}else{
		size_t left_col = this->SeqPosToColumn(this->Length(seqI) - crop_amount + 1, align_matrix[seqI], bcount[seqI]);
		if( left_col > 4000000000u )
		{
			std::cerr << this->LeftEnd(seqI) << std::endl;
			std::cerr << this->LeftEnd(0) << std::endl;
			std::cerr << "bogus cropper cga\n";
		}
		this->CropEnd(this->AlignmentLength()-left_col);
	}
	if( this->Length(seqI) != pre_len - crop_amount )
	{
		std::cerr << this->LeftEnd(seqI) << std::endl;
		std::cerr << this->LeftEnd(0) << std::endl;
		std::cerr << "bad cropperLeftie\n";
	}
	if( !this->validate() )
		std::cerr << "CropLeft error\n";
}

template< class BaseType >
void CompactGappedAlignment<BaseType>::CropRight(gnSeqI crop_amount, uint seqI)
{
	if( crop_amount == 0 )
		return;

	gnSeqI pre_len = this->Length(seqI);
	gnSeqI pre_lend = this->LeftEnd(seqI);
	gnSeqI pre_lend0 = this->LeftEnd(0);
	if( this->Orientation(seqI) == AbstractMatch::forward )
	{
		// count "crop_amount" characters into seqI and crop there
		size_t right_col = this->SeqPosToColumn(this->Length(seqI) - crop_amount + 1, align_matrix[seqI], bcount[seqI]);
		this->CropEnd( this->AlignmentLength()-right_col );
	}else
	{
		size_t right_col = this->SeqPosToColumn(crop_amount, align_matrix[seqI], bcount[seqI]) + 1;
		if( right_col > 4000000000u )
		{
			std::cerr << this->LeftEnd(seqI) << std::endl;
			std::cerr << this->LeftEnd(0) << std::endl;
			std::cerr << "bogus cropper cga\n";
		}
		this->CropStart( right_col );
	}
	if( this->Length(seqI) != pre_len - crop_amount )
	{
		std::cerr << this->LeftEnd(seqI) << std::endl;
		std::cerr << this->LeftEnd(0) << std::endl;
		std::cerr << "bad cropperight\n";
	}
	if( !this->validate() )
		std::cerr << "CropRight error\n";
}

template< class BaseType >
void CompactGappedAlignment<BaseType>::GetColumn( gnSeqI col, std::vector<gnSeqI>& pos, std::vector<bool>& column ) const
{
	pos = std::vector<gnSeqI>(this->SeqCount(), NO_MATCH);
	column = std::vector<bool>(this->SeqCount(), false);
	for( uint seqI = 0; seqI < this->SeqCount(); seqI++ )
	{
		if( align_matrix[seqI][col] )
			column[seqI] = true;

		gnSeqI count = 0;
		if( this->LeftEnd(seqI) != NO_MATCH )
		{
			size_t col_index = col / INDEX_INTERVAL;
			for( size_t i = col_index * INDEX_INTERVAL; i <= col; i++ )
				count += align_matrix[seqI].test(i);
			count += bcount[seqI][col_index];
		}

		if( count > 0 && this->Orientation(seqI) == AbstractMatch::forward )
			pos[seqI] = this->LeftEnd(seqI) + count - 1;
		else if( this->Orientation(seqI) == AbstractMatch::reverse && !(count == this->Length(seqI) && !column[seqI]) )
			pos[seqI] = this->RightEnd(seqI) - count + 1;
	}
}

template< class BaseType >
void CompactGappedAlignment<BaseType>::copyRange( CompactGappedAlignment& dest, gnSeqI left_column, gnSeqI length )
{
	if( left_column + length > this->AlignmentLength() )
		Throw_gnEx( genome::SeqIndexOutOfBounds() );
//	if( length == 0 )
//		return;

	// first copy the coordinates
	dest = CompactGappedAlignment(this->SeqCount(), length);
	for( uint i=0; i < this->SeqCount(); i++ ){
		dest.SetStart(i, this->Start(i));
		if( this->Orientation(i) != AbstractMatch::undefined )
			dest.SetLength(this->Length(i), i);
	}
	// then trim the coordinates appropriately

	gnSeqI pre_alignlen = this->AlignmentLength();
	gnSeqI pre_lend0 = this->LeftEnd(0);

	std::vector< bitset_t > dest_mat(this->SeqCount(), bitset_t(length));
	std::vector<gnSeqI> pos;
	std::vector<bool> column;
	std::vector<gnSeqI> left_cc(this->SeqCount(), 0);
	if( left_column > 0 )
	{
		this->GetColumn( left_column-1, pos, column );
		for( uint i=0; i < this->SeqCount(); i++ ){
			if( this->LeftEnd(i) == NO_MATCH )
				continue;

			size_t char_count = this->Orientation(i) == AbstractMatch::forward ? pos[i] - this->LeftEnd(i) + 1 : this->RightEnd(i) - pos[i] + 1;
			if( pos[i] > 0 && char_count > 0 )
			{
				left_cc[i] = char_count;
				if( dest.Orientation(i) == AbstractMatch::forward )
					dest.SetStart(i, dest.Start(i) + char_count);
			}else if( pos[i] == 0 && dest.Orientation(i) == AbstractMatch::reverse )
			{
				// this sequence was completely obliterated by the crop
				dest.SetStart(i, NO_MATCH);
			}
		}
	}

// now trim up the right side...
	gnSeqI right_trim = this->AlignmentLength() - left_column - length;

	if( right_trim > 0 )
	{
		this->GetColumn( this->AlignmentLength()-right_trim, pos, column );

		for( uint i=0; i < this->SeqCount(); i++ ){
			if( this->LeftEnd(i) == NO_MATCH )
				continue;
			AbstractMatch::orientation orient = this->Orientation(i);
			if( pos[i] > 0 )
			{
				gnSeqI char_count = pos[i] - (orient == AbstractMatch::forward ? (column[i] ? 1 : 0 ) : (column[i] ? 0 : 1 ) );
				char_count = orient == AbstractMatch::forward ? char_count - this->LeftEnd(i) + 1 : this->RightEnd(i) - char_count;
				char_count -= left_cc[i];
				gnSeqI deleted = this->Length(i) - char_count;
				if( dest.Start(i) < 0 )
					dest.SetStart(i, dest.Start(i)-deleted+left_cc[i]);	// fixme: is this off-by-one?
			}else if( orient == AbstractMatch::forward ){
				dest.SetStart(i, NO_MATCH);
			}
		}
	}

	for( size_t i = 0; i < dest_mat.size(); ++i )
	{
		size_t count = 0;
		for( size_t j = 0; j < length; ++j )
		{
			if(align_matrix[i].test(j+left_column))
			{
				dest_mat[i].set(j, true);
				++count;
			}
		}
		dest.SetLength(count, i);
		if( count == 0 )
			dest.SetStart(i, NO_MATCH);
	}
	dest.SetAlignment(dest_mat);

	dest.create_bitcount();
	if( !dest.validate() )
	{
		std::cerr << "pre_lend0: " << pre_lend0 << std::endl;
		std::cerr << "pre_alignlen: " << pre_alignlen << std::endl;
		std::cerr << "CropStart error\n";
	}

}

template< class BaseType >
void CompactGappedAlignment<BaseType>::CondenseGapColumns()
{
	const size_t len = this->AlignmentLength();
	size_t d = 0;	// destination index
	for( size_t i = 0; i < len; ++i )
	{
		size_t seqI = 0;
		// check whether this is a gap col
		for( ; seqI < align_matrix.size(); ++seqI )
			if( this->LeftEnd(seqI) != 0 && align_matrix[seqI].test(i) )
				break;

		// copy if not a gap col (and i != d )
		if( seqI < align_matrix.size() )
		{
			if( i != d )
			{
				for( seqI = 0; seqI < align_matrix.size(); ++seqI )
					align_matrix[seqI].set( d, align_matrix[seqI].test(i)  );
			}
			d++;
		}
		else
			std::cout << "";
	}
	this->SetAlignmentLength(d);
	for( size_t seqI = 0; seqI < align_matrix.size(); ++seqI )
	{
		align_matrix[seqI].resize(d);
		align_matrix[seqI] = align_matrix[seqI];	// force reallocation on "optimized" windows builds
	}
	this->create_bitcount();
}


}

namespace std {
template<> inline
void swap( mems::CompactGappedAlignment<>& a, mems::CompactGappedAlignment<>& b )
{
	a.swap(b);
}
}


#endif // __CompactGappedAlignment_h__

