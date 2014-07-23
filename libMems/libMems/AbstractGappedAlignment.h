/*******************************************************************************
 * $Id: AbstractGappedAlignment.h,v 1.12 2004/04/19 23:10:50 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef __AbstractGappedAlignment_h__
#define __AbstractGappedAlignment_h__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/AbstractMatch.h"
#include "libGenome/gnFilter.h"

namespace mems {

template<class AbstractMatchImpl>
class AbstractGappedAlignment : public AbstractMatchImpl
{
public:
	AbstractGappedAlignment();
	AbstractGappedAlignment( uint seq_count, gnSeqI align_length );
	
	/**
	 * Sets the alignment 
	 * @param seq_align	should be in row/column format, e.g. one string per sequence (row)
	 */
	virtual void SetAlignment( const std::vector< std::string >& seq_align ) = 0;

	// Inherited methods from AbstractMatch:
	gnSeqI Length( uint seqI = UINT_MAX ) const; 
	virtual void SetLength( gnSeqI len, uint seqI ) { length[ seqI ] = len; }

	gnSeqI AlignmentLength() const {return align_length;}
	void SetAlignmentLength(gnSeqI len){ align_length = len; }

protected:
	// for use by derived classes in order to swap contents
	void swap( AbstractGappedAlignment* other );	
private:
	std::vector< gnSeqI > length;
	gnSeqI align_length;
};


template<class AbstractMatchImpl>
AbstractGappedAlignment<AbstractMatchImpl>::AbstractGappedAlignment() : AbstractMatchImpl()
{
	align_length = 0;
}

template<class AbstractMatchImpl>
AbstractGappedAlignment<AbstractMatchImpl>::AbstractGappedAlignment( uint seq_count, gnSeqI align_length ) : AbstractMatchImpl( seq_count )
{
	length = std::vector< gnSeqI >( seq_count, 0 );
	this->align_length = align_length;
}

template<class AbstractMatchImpl>
void AbstractGappedAlignment<AbstractMatchImpl>::swap( AbstractGappedAlignment* other )
{
	std::swap( length, other->length );
	std::swap( align_length, other->align_length );
	AbstractMatchImpl::swap( other );
}

template<class AbstractMatchImpl>
gnSeqI AbstractGappedAlignment<AbstractMatchImpl>::Length( uint seqI ) const 
{
	if( seqI == UINT_MAX )
		return align_length;
	return length[ seqI ]; 
}

//template<class AbstractGappedAlignmentImpl>
void GetAlignment( const AbstractMatch& ga, const std::vector< genome::gnSequence* >& seq_table, std::vector<std::string>& alignment );

//template<class AbstractGappedAlignmentImpl>
inline
void GetAlignment( const AbstractMatch& ga, const std::vector< genome::gnSequence* >& seq_table, std::vector<std::string>& alignment )
{
	std::vector< bitset_t > aln_mat;
	ga.GetAlignment(aln_mat);
	alignment = std::vector<std::string>( aln_mat.size() );
	const genome::gnFilter* comp_filter = genome::gnFilter::DNAComplementFilter();
	for( std::size_t seqI = 0; seqI < alignment.size(); seqI++ )
	{
		alignment[seqI] = std::string( aln_mat[0].size(), '-' );
		if( ga.LeftEnd(seqI) == NO_MATCH )
			continue;
		std::string cur_seq;
		seq_table[seqI]->ToString( cur_seq, ga.Length(seqI), ga.LeftEnd(seqI) );
		if( ga.Orientation(seqI) == AbstractMatch::reverse )
			comp_filter->ReverseFilter(cur_seq);
		std::size_t cI = 0; 
		for( std::size_t gI = 0; gI < alignment[seqI].size(); gI++ )
			if( aln_mat[seqI][gI] )
				alignment[seqI][gI] = cur_seq[cI++];
	}
}

}

#endif // __AbstractGappedAlignment_h__

