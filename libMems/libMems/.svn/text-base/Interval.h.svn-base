/*******************************************************************************
 * $Id: GenericInterval.h,v 1.4 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef __Interval_h__
#define __Interval_h__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnClone.h"
#include "libGenome/gnDebug.h"
#include "libMems/SparseAbstractMatch.h"
#include "libMems/gnAlignedSequences.h"
#include "libMems/AbstractGappedAlignment.h"
#include "libMems/Match.h"
#include "libMems/GappedAlignment.h"
#include <iostream>
#include <vector>
#include "libMems/twister.h"

//#include "boost/pool/object_pool.hpp"

namespace mems {

// adapter function to allow inserts on reverse iterators
template< typename ListType, typename RanIt, typename Ty >
void insert( ListType& the_list, std::reverse_iterator<RanIt>& riter, Ty& val )
{
	the_list.insert( riter.base(), val );
	++riter;	// need to shift riter
}
template< typename ListType, typename Ty >
void insert( ListType& the_list, const typename ListType::iterator& iter, Ty& val )
{
	the_list.insert( iter, val );
}


template< class GappedBaseImpl = AbstractGappedAlignment< SparseAbstractMatch<> > >
class GenericInterval : public GappedBaseImpl
{
public:
	GenericInterval(){};

//	GenericInterval( uint seq_count, gnSeqI aln_length) : GappedBaseImpl( seq_count, aln_length ){};

	/** construct from a MatchList or a vector of pointers to AbstractMatches */
	template<typename BidIt>
	GenericInterval( BidIt it_begin, const BidIt& it_end ) : GappedBaseImpl( (*it_begin)->SeqCount(), 0 )
	{
		std::vector<gnSeqI> pos((*it_begin)->SeqCount(), NO_MATCH);
		for( ; it_begin != it_end; ++it_begin )
			this->matches.push_back( (*it_begin)->Copy() );
		CalculateOffset();
		addUnalignedRegions();
		CalculateAlignmentLength();
		ValidateMatches();
	}

	GenericInterval( const GenericInterval& iv );
	~GenericInterval();
	GenericInterval& operator=( const GenericInterval& iv );
	
	GenericInterval* Clone() const;
	GenericInterval* Copy() const;
	virtual void Free();
	
	/** Set the matches in this interval *without* making a copy.  The GenericInterval takes ownership of matches */
	template< class MatchVector >
	void SetMatches( MatchVector& matches )
	{
		// Set the SeqCount and other bits
		Match m( matches[0]->SeqCount() );
		std::vector<AbstractMatch*> tmp(1, &m);
		*this = GenericInterval( tmp.begin(), tmp.end() );

		// then delete the allocated dummy match
		for( std::size_t mI = 0; mI < this->matches.size(); mI++ )
			this->matches[mI]->Free();
		
		// now set the matches and update the interval data
		this->matches.resize(matches.size());
		std::copy(matches.begin(), matches.end(), this->matches.begin());
//		this->matches.insert( this->matches.end(), matches.begin(), matches.end() );
		CalculateOffset();
 	    addUnalignedRegions();
		CalculateAlignmentLength();
		ValidateMatches();

		// finally, clear the user supplied matches to indicate that we own the memory
		matches.clear();
	}

	/** Set the matches in this interval *without* cloberring the interval.*/
	template< class MatchVector >
	void SetMatchesTemp( MatchVector& matches )
	{
		// Set the SeqCount and other bits
		Match m( matches[0]->SeqCount() );
		std::vector<AbstractMatch*> tmp(1, &m);
		*this = GenericInterval( tmp.begin(), tmp.end() );

		// then delete the allocated dummy match
		for( std::size_t mI = 0; mI < this->matches.size(); mI++ )
			this->matches[mI]->Free();
		
		// now set the matches and update the interval data
		this->matches.resize(matches.size());
		std::copy(matches.begin(), matches.end(), this->matches.begin());
		CalculateOffset();
		CalculateAlignmentLength();
		ValidateMatches();

		// finally, clear the user supplied matches to indicate that we own the memory
		matches.clear();
	}
	/**
	 * Writes this GenericInterval to the specified output stream (e.g. cout).
	 */
	template<typename BaseImpl> friend std::ostream& operator<<(std::ostream& os, const GenericInterval<BaseImpl>& iv); //write to source.

	/**
	 * Reads a GenericInterval from the specified input stream (e.g. cin).
	 */
	template<typename BaseImpl> friend std::istream& operator>>(std::istream& is, const GenericInterval<BaseImpl>& iv); //read from source

	// Inherited methods from AbstractMatch:
	void Invert();
	void CropStart(gnSeqI crop_amount);
	void CropEnd(gnSeqI crop_amount);
	void MoveStart(int64 move_amount);
	void MoveEnd(int64 move_amount);

	virtual void CalculateOffset();

	void add( AbstractMatch* am ){ matches.push_back( am->Copy() ); }
	
	/** 
	 * Get a gnAlignedSequences object
	 * TODO: get rid of this
	 */
	virtual void GetAlignedSequences( gnAlignedSequences& gnas, const std::vector< genome::gnSequence* >& seq_table ) const;

	void GetAlignment( std::vector< bitset_t >& align_matrix ) const;

	void CropLeft( gnSeqI amount, uint seqI );
	void CropRight( gnSeqI amount, uint seqI );

	void SetAlignment( const std::vector< std::string >& seq_align );

	// TODO: get rid of code that uses this hack...
	const std::vector<AbstractMatch*>& GetMatches() const{ return matches; }
	void StealMatches( std::vector<AbstractMatch*>& matches );

	/** marbles the gaps so that no sequence has more than "size" contiguous gaps */
	void Marble( gnSeqI size );

	void GetColumn( gnSeqI col, std::vector<gnSeqI>& pos, std::vector<bool>& column ) const;
	
	bool IsGap( uint seq, gnSeqI col ) const;

	/** self test code */
	void ValidateMatches() const;

	void swap( GenericInterval& other ){ swap(&other); }

protected:
	// for use by derived classes in order to swap contents
	void swap( GenericInterval* other ){
		std::swap( matches, other->matches );
		GappedBaseImpl::swap( other );
	}
	std::vector< AbstractMatch* > matches;
private:
	void addUnalignedRegions();
	void FindMatchPos( uint seqI, gnSeqI pos, size_t& matchI, gnSeqI& match_pos );
	void GetColumnAndMatch( gnSeqI col, std::vector<gnSeqI>& pos, std::vector<bool>& column, size_t& matchI, gnSeqI& match_col ) const;
	void CalculateAlignmentLength();
};

typedef GenericInterval<> Interval;


template<class GappedBaseImpl>
GenericInterval<GappedBaseImpl>* GenericInterval<GappedBaseImpl>::Copy() const
{
	return m_allocateAndCopy( *this );
}
template<class GappedBaseImpl>
void GenericInterval<GappedBaseImpl>::Free()
{
	m_free(this);
}

template<class GappedBaseImpl>
GenericInterval<GappedBaseImpl>::~GenericInterval()
{
	for( std::size_t mI = 0; mI < matches.size(); mI++ )
		matches[mI]->Free();
}

template<class GappedBaseImpl>
void GenericInterval<GappedBaseImpl>::StealMatches( std::vector<AbstractMatch*>& matches ){
	matches = this->matches;
	this->matches.clear();
	for( uint seqI = 0; seqI < this->SeqCount(); seqI++ )
	{
		this->SetLeftEnd( seqI, NO_MATCH );
		this->SetLength( 0, seqI );
	}
	this->SetAlignmentLength(0);
}

template<class GappedBaseImpl>
GenericInterval<GappedBaseImpl>::GenericInterval( const GenericInterval<GappedBaseImpl>& iv )
{
	*this = iv;
}

template<class GappedBaseImpl>
GenericInterval<GappedBaseImpl>& GenericInterval<GappedBaseImpl>::operator=( const GenericInterval& iv )
{
	GappedBaseImpl::operator=( iv );
	for( std::size_t mI = 0; mI < matches.size(); mI++ )
		matches[mI]->Free();
	matches.clear();
	for( std::size_t mI = 0; mI < iv.matches.size(); mI++ )
		matches.push_back( iv.matches[mI]->Copy() );
	return *this;
}

template<class GappedBaseImpl>
GenericInterval<GappedBaseImpl>* GenericInterval<GappedBaseImpl>::Clone() const 
{
	return new GenericInterval( *this );
}


static bool debug_interval = false;

template<class GappedBaseImpl>
void GenericInterval<GappedBaseImpl>::ValidateMatches() const
{
	if( !debug_interval )
		return;
	if( matches.size() == 0 )
	{
//		genome::breakHere();
//		std::cerr << "iv has no matches\n";
		return;
	}
	for( uint seqI = 0; seqI < matches[0]->SeqCount(); ++seqI )
	{
		gnSeqI prev_rend = this->LeftEnd(seqI);
		if( this->Orientation(seqI) == AbstractMatch::forward )
		{
			for( size_t mI = 0; mI < matches.size(); ++mI )
			{
				if( matches[mI]->LeftEnd(seqI) != NO_MATCH )
				{
					if( prev_rend != matches[mI]->LeftEnd(seqI) )
					{
						std::cerr << "iv broken\n";
						std::cerr << "seqI: " << seqI << "\t prev_rend: " << prev_rend << std::endl;
						std::cerr << "mI: " << mI << "\tlend: " << matches[mI]->LeftEnd(seqI) << std::endl;
						genome::breakHere();
					}
					prev_rend = matches[mI]->RightEnd(seqI) + 1;
				}
			}
		}else if( this->Orientation(seqI) == AbstractMatch::reverse )
		{
			for( size_t mI = matches.size(); mI > 0; mI-- )
			{
				if( matches[mI-1]->LeftEnd(seqI) != NO_MATCH )
				{
					if( prev_rend != matches[mI-1]->LeftEnd(seqI) )
					{
						std::cerr << "iv broken 2\n";
						genome::breakHere();
					}
					prev_rend = matches[mI-1]->RightEnd(seqI) + 1;
				}
			}
		}

		if( this->Orientation(seqI) != AbstractMatch::undefined && this->Length(seqI) == 0 )
		{
			genome::breakHere();
			std::cerr << "ERROR: confused interval\n";
		}
	}
}

template<class GappedBaseImpl>
void GenericInterval<GappedBaseImpl>::GetColumnAndMatch( gnSeqI col, std::vector<gnSeqI>& pos, std::vector<bool>& column, size_t& matchI, gnSeqI& match_col ) const
{
	// bail when the appropriate match is found
	gnSeqI col_pos = 0;
	size_t mI = 0;
	pos.clear();
	for( uint seqI = 0; seqI < this->SeqCount(); ++seqI )
	{
		if( this->LeftEnd(seqI) == NO_MATCH )
			pos.push_back(NO_MATCH);
		else if( this->Orientation(seqI) == AbstractMatch::forward )
			pos.push_back(this->LeftEnd(seqI));
		else
			pos.push_back(this->RightEnd(seqI)+1);
	}

	column = std::vector<bool>(this->SeqCount(), false);

	for( ; mI < matches.size(); ++mI )
	{
		uint seqI = 0;

		gnSeqI diff = matches[mI]->AlignmentLength();
		diff = col_pos + diff <= col ? diff : col - col_pos;

		for( seqI = 0; seqI < this->SeqCount(); ++seqI )
			if( this->Orientation(seqI) == AbstractMatch::forward )
				pos[seqI] += diff;
			else if( this->Orientation(seqI) == AbstractMatch::reverse )
				pos[seqI] -= diff;

		col_pos += diff;

		if( col_pos >= col && diff < matches[mI]->AlignmentLength() )
		{
			std::vector<gnSeqI> m_pos;
			matches[mI]->GetColumn( diff, m_pos, column );
			for( uint seqI = 0; seqI < this->SeqCount(); ++seqI )
				if( m_pos[seqI] != NO_MATCH )
					pos[seqI] = m_pos[seqI];
			matchI = mI;
			match_col = diff;
			break;
		}
	}
}

template<typename ListType, typename Iter>
void AddGapMatches( ListType& the_list, const Iter& first, const Iter& last, 
				   uint seqI, int64 left_end, int64 right_end, 
				   AbstractMatch::orientation seq_orient, uint seq_count )
{
	Iter iter = first;
	int64 pos = left_end-1;
    //MatchList& tmp_list;
    std::vector< std::pair<Match*,Iter> > insert_pos;
	for( ; iter != last; ++iter )
	{
		if( (*iter)->LeftEnd(seqI) != NO_MATCH )
		{
			gnSeqI len = (*iter)->LeftEnd(seqI)-pos-1;

            //tjt: there are perfectly valid chains that blow up when this is enabled
            //i.e:      
            //                         <----c1----><----d1---->
            //          <--a1---><---b1--->
            // pos would get set to b1->RightEnd() since diff between a1 & b1 == 0
            // but then c1->LeftEnd < pos, so genome::breakHere() gets called
            // this is because SetMatches() gets called before finalize(), but should it??

            if( len > 4000000000u )
			{
				std::cerr << "triplebogus interval data\n";
				std::cerr << "(*iter)->LeftEnd(" << seqI << "): " << (*iter)->LeftEnd(seqI) << std::endl;
				std::cerr << "pos: " << pos << std::endl;
				genome::breakHere();
			}

			if( len > 0 )
			{
				Match tmp(seq_count);
				Match* new_m = tmp.Copy();
				new_m->SetLeftEnd(seqI, pos + 1);
				new_m->SetOrientation(seqI, seq_orient);
				new_m->SetLength(len);
				pos = (*iter)->RightEnd(seqI);
				//insert(the_list, iter, new_m);	// this may move iter
                //tmp_list.push_back(new_m);
                insert_pos.push_back(make_pair(new_m,iter));
			}
            else
				pos = (*iter)->RightEnd(seqI);
		}
	}
    for ( uint i = 0; i < insert_pos.size(); i++)
    {
        insert(the_list, insert_pos.at(i).second, insert_pos.at(i).first);
    }
	if( right_end != pos )
	{
		Match tmp(seq_count);
		Match* new_m = tmp.Copy();
		new_m->SetLeftEnd(seqI, pos+1);
		new_m->SetLength(right_end-pos-1);
		insert(the_list, iter, new_m);
	}
}

// The best steaks are well marbled
template<class GappedBaseImpl>
void GenericInterval<GappedBaseImpl>::Marble( gnSeqI size )
{
	if( this->SeqCount() > 2 )
		throw "I can't handle that many at once\n";
	if( this->Multiplicity() < 2 )
		return;	// can't marble unless there are at least two seqs

	// first break up all the pieces
	std::list<AbstractMatch*> mlist;
	mlist.insert( mlist.end(), matches.begin(), matches.end() );
	std::list<AbstractMatch*>::iterator m_iter = mlist.begin();
	for(; m_iter != mlist.end(); ++m_iter )
	{
		if( (*m_iter)->Multiplicity() != 1 || (*m_iter)->AlignmentLength() <= size )
			continue;
		// which seq are we working with?
		uint seqI = 0;
		for( ; seqI < (*m_iter)->SeqCount(); seqI++ )
			if( (*m_iter)->LeftEnd(seqI) != NO_MATCH )
				break;
		AbstractMatch* left_iv = (*m_iter)->Copy();
		left_iv->CropEnd( left_iv->AlignmentLength() - size );
		(*m_iter)->CropStart( size );
		m_iter = mlist.insert( m_iter, left_iv );
	}
	matches.clear();
	matches.insert( matches.end(), mlist.begin(), mlist.end() );
	this->ValidateMatches();

	// now interleave the gaps
	std::vector< std::vector<AbstractMatch*>::iterator > seq_iter( this->SeqCount(), matches.begin() );
	std::vector< AbstractMatch* > interleaved(matches.size());
	std::vector<AbstractMatch*>::iterator anchor = matches.begin();
	for( uint seqI = 0; seqI < this->SeqCount(); seqI++ )
	{
		if( this->LeftEnd(seqI) == NO_MATCH )
			continue;
		for( ; seq_iter[seqI] != matches.end() && (*seq_iter[seqI])->LeftEnd(seqI) == NO_MATCH; ++seq_iter[seqI] );
	}
	for( ; anchor != matches.end() && (*anchor)->Multiplicity() < this->SeqCount(); ++anchor );
	size_t cur = 0;
	while(true)
	{
		// increment anchor if an iter has caught up to it...
		uint seqI = 0;
		do{
			for( seqI = 0; seqI < this->SeqCount(); seqI++ )
			{
				if( seq_iter[seqI] == anchor && anchor != matches.end() )
				{
					for( uint seqJ = 0; seqJ < this->SeqCount(); seqJ++ )
					{
						// add anything in seq_iter[seqJ]
						while( seq_iter[seqJ] != anchor )
						{
							interleaved[cur++] = *(seq_iter[seqJ]);
							for( ++seq_iter[seqJ]; seq_iter[seqJ] != matches.end() && (*seq_iter[seqJ])->LeftEnd(seqJ) == NO_MATCH; ++seq_iter[seqJ] );
						}
						// don't end on an anchor
						for( ++seq_iter[seqJ]; seq_iter[seqJ] != matches.end() && (*seq_iter[seqJ])->LeftEnd(seqJ) == NO_MATCH; ++seq_iter[seqJ] );
					}
					// increment anchor
					interleaved[cur++] = *anchor;
					for( ++anchor; anchor != matches.end() && (*anchor)->Multiplicity() < this->SeqCount(); ++anchor );

					break;
				}
			}
		}while( seqI < this->SeqCount() );

		size_t diff1 = anchor - seq_iter[0];
		size_t diff2 = anchor - seq_iter[1];
		if( diff1 == 0 && diff2 == 0 )
			break;
		// sample from a binomial with p(success) = diff1 / diff1+diff2
//		double samp = ((double)rand())/((double)RAND_MAX);
		double samp = RandTwisterDouble();
		// add one of the intervals and move on to the next...
		if( diff2 == 0 || (samp < .5 && diff1 > 0) )
		{
			interleaved[cur++] = *(seq_iter[0]);
			for( ++seq_iter[0]; seq_iter[0] != matches.end() && (*seq_iter[0])->LeftEnd(0) == NO_MATCH; ++seq_iter[0] );
		}else{
			interleaved[cur++] = *(seq_iter[1]);
			for( ++seq_iter[1]; seq_iter[1] != matches.end() && (*seq_iter[1])->LeftEnd(1) == NO_MATCH; ++seq_iter[1] );
		}
	}
	matches = interleaved;
	this->ValidateMatches();
}


template<class GappedBaseImpl>
void GenericInterval<GappedBaseImpl>::CropStart(gnSeqI crop_amount)
{
	if( crop_amount > this->AlignmentLength() )
		Throw_gnEx( genome::SeqIndexOutOfBounds() );
	if( crop_amount == 0 )
		return;
	
	std::vector<bool> col;
	std::vector<gnSeqI> pos;
	size_t matchI = 0;
	gnSeqI match_col;
	this->GetColumnAndMatch( crop_amount, pos, col, matchI, match_col );

	// delete everything before matchI
	for( size_t mI = 0; mI < matchI; ++mI )
		matches[mI]->Free();
	matches.erase(matches.begin(), matches.begin()+matchI);

	// crop from within matchI
	matches[0]->CropStart(match_col);

	this->CalculateOffset();
	this->ValidateMatches();
}

template<class GappedBaseImpl>
void GenericInterval<GappedBaseImpl>::CropEnd(gnSeqI crop_amount)
{
	if( crop_amount > this->AlignmentLength() )
		Throw_gnEx( genome::SeqIndexOutOfBounds() );
	if( crop_amount == 0 )
		return;
	std::vector<bool> col;
	std::vector<gnSeqI> pos;
	size_t matchI = 0;
	gnSeqI match_col;
	this->GetColumnAndMatch( this->AlignmentLength()-crop_amount, pos, col, matchI, match_col );

	// delete everything after matchI
	size_t plusmatch = match_col == 0 ? 0 : 1;
	for( size_t mI = matchI+plusmatch; mI < matches.size(); ++mI )
		matches[mI]->Free();
	matches.erase(matches.begin()+matchI+plusmatch, matches.end());

	// crop from within matchI
	if( matches.size() > 0 && plusmatch == 1 )
		matches.back()->CropEnd(matches.back()->AlignmentLength() - match_col);

	this->CalculateOffset();
	this->ValidateMatches();
}

template<class GappedBaseImpl>
void GenericInterval<GappedBaseImpl>::GetAlignment( std::vector< bitset_t >& align_matrix ) const
{
	gnSeqI cur_col = 0;
	align_matrix = std::vector< bitset_t >( this->SeqCount(), bitset_t(this->AlignmentLength(),false) );
	for( uint matchI = 0; matchI < matches.size(); ++matchI ){
		std::vector< bitset_t > aln_mat;
		matches[matchI]->GetAlignment( aln_mat );
		for( uint seqI = 0; seqI < this->SeqCount(); ++seqI )
		{
			if( matches[matchI]->LeftEnd(seqI) == NO_MATCH || matches[matchI]->Length(seqI) == 0 )
				continue;

			size_t ct = 0;
			gnSeqI len = matches[matchI]->Length(seqI);
			for( bitset_t::size_type pos = aln_mat[seqI].find_first(); ct < len; pos = aln_mat[seqI].find_next(pos) )
			{
				align_matrix[seqI].set( cur_col + pos );
				ct++;
			}
		}
		cur_col += matches[matchI]->AlignmentLength();
	}
}

template<class GappedBaseImpl>
void GenericInterval<GappedBaseImpl>::CropLeft( gnSeqI amount, uint seqI )
{
	if( amount > this->Length(seqI) )
		Throw_gnEx( genome::SeqIndexOutOfBounds() );
	if( this->LeftEnd(seqI) == NO_MATCH || amount == 0 )
		return;

	// for debugging
	gnSeqI pre_len = this->Length(seqI);
	gnSeqI pre_lend = this->LeftEnd(seqI);

	gnSeqI match_pos;
	size_t mI;
	this->FindMatchPos(seqI, amount, mI, match_pos);
	if( matches[mI]->Orientation(seqI) == this->Orientation(seqI) )
		matches[mI]->CropLeft(match_pos, seqI);
	else
		matches[mI]->CropRight(match_pos, seqI);

	if( matches[mI]->Length(seqI) == 0 )
		std::cerr << "Big fat zero 1\n";

	// get rid of everything to the left of mI
	if( this->Orientation(seqI) == AbstractMatch::forward )
	{
		for( size_t m = 0; m < mI; m++ )
			matches[m]->Free();
		matches.erase(matches.begin(), matches.begin()+mI);
	}else{
		for( size_t m = mI+1; m < matches.size(); m++ )
			matches[m]->Free();
		matches.erase(matches.begin()+mI+1, matches.end());
	}

	this->CalculateOffset();
	this->ValidateMatches();

	if( this->Length(seqI) != pre_len - amount )
	{
		std::cerr << "Error intercroplef\n";
		std::cerr << "pre len: " << pre_len << std::endl;
		std::cerr << "pre lend: " << pre_lend << std::endl;
		std::cerr << "amount: " << amount << std::endl;
		std::cerr << "LeftEnd(seqI) " << this->LeftEnd(seqI) << std::endl;
		std::cerr << "Length(seqI) " << this->Length(seqI) << std::endl;
		std::cerr << "AlignmentLength() " << this->AlignmentLength() << std::endl;
		genome::breakHere();
	}
}

template<class GappedBaseImpl>
void GenericInterval<GappedBaseImpl>::CropRight( gnSeqI amount, uint seqI )
{
	if( amount > this->Length(seqI) )
		Throw_gnEx( genome::SeqIndexOutOfBounds() );

	if( this->LeftEnd(seqI) == NO_MATCH || amount == 0 )
		return;

	// for debugging
	gnSeqI pre_len = this->Length(seqI);
	gnSeqI pre_lend = this->LeftEnd(seqI);

	gnSeqI left_amount = this->Length(seqI) - amount;
	gnSeqI match_pos;
	size_t mI;
	this->FindMatchPos(seqI, left_amount, mI, match_pos);
	if( matches[mI]->Orientation(seqI) == this->Orientation(seqI) )
		matches[mI]->CropRight(matches[mI]->Length(seqI)-match_pos, seqI);
	else
		matches[mI]->CropLeft(matches[mI]->Length(seqI)-match_pos, seqI);

	if( matches[mI]->Length(seqI) == 0 )
		mI += this->Orientation(seqI) == AbstractMatch::forward ? -1 : 1;	// delete this match too

	// get rid of everything to the left of mI
	if( this->Orientation(seqI) == AbstractMatch::forward )
	{
		for( size_t m = mI+1; m < matches.size(); m++ )
			matches[m]->Free();
		matches.erase(matches.begin()+(mI+1), matches.end());
	}else{
		for( size_t m = 0; m < mI; m++ )
			matches[m]->Free();
		matches.erase(matches.begin(), matches.begin()+mI);
	}

	this->CalculateOffset();
	this->ValidateMatches();

	if( this->Length(seqI) != pre_len - amount )
	{
		std::cerr << "Error intercropright\n";
		std::cerr << "pre len: " << pre_len << std::endl;
		std::cerr << "pre lend: " << pre_lend << std::endl;
		std::cerr << "amount: " << amount << std::endl;
		std::cerr << "LeftEnd(seqI) " << this->LeftEnd(seqI) << std::endl;
		std::cerr << "Length(seqI) " << this->Length(seqI) << std::endl;
		std::cerr << "AlignmentLength() " << this->AlignmentLength() << std::endl;
		genome::breakHere();
	}
}

template<class GappedBaseImpl>
void GenericInterval<GappedBaseImpl>::MoveStart(int64 move_amount)
{
	GappedBaseImpl::MoveStart(move_amount);
	for( size_t mI = 0; mI < matches.size(); mI++ )
		matches[mI]->MoveStart(move_amount);
}

template<class GappedBaseImpl>
void GenericInterval<GappedBaseImpl>::MoveEnd(int64 move_amount)
{
	GappedBaseImpl::MoveEnd(move_amount);
	for( size_t mI = 0; mI < matches.size(); mI++ )
		matches[mI]->MoveEnd(move_amount);
}


template< class MatchVector >
void FindBoundaries( const MatchVector& matches, std::vector<gnSeqI>& left_ends, std::vector<gnSeqI>& lengths, std::vector<bool>& orientations )
{
	uint seqI;
	boolean zero_exists = false;
	uint seq_count = matches.front()->SeqCount();
	left_ends = std::vector<gnSeqI>( seq_count, NO_MATCH );
	lengths = std::vector<gnSeqI>( seq_count, 0 );
	orientations = std::vector<bool>( seq_count, false );

	// find leftend in each forward sequence
	uint matchI = 0;
	for(; matchI != matches.size(); ++matchI )
	{
		zero_exists = false;
		for( seqI = 0; seqI < seq_count; ++seqI )
		{
			if( left_ends[seqI] == NO_MATCH && matches[matchI]->Orientation(seqI) == AbstractMatch::forward )
			{
				left_ends[seqI] = matches[ matchI ]->LeftEnd(seqI);
				orientations[seqI] = true;
			}
			else if( left_ends[seqI] == NO_MATCH )
				zero_exists = true;
		}
		if( !zero_exists )
			break;
	}

	// find end in each forward sequence
	for( matchI = matches.size(); matchI > 0; matchI-- )
	{
		zero_exists = false;
		for( seqI = 0; seqI < seq_count; ++seqI )
		{
			if( lengths[seqI] == 0 &&
				matches[ matchI - 1 ]->Orientation(seqI) == AbstractMatch::forward )
			{
					lengths[seqI] = matches[matchI - 1]->LeftEnd(seqI) + matches[matchI - 1]->Length(seqI) - left_ends[seqI];
			}
			if( left_ends[seqI] != NO_MATCH && lengths[seqI] == 0 )
				zero_exists = true;
		}
		if( !zero_exists )
			break;
	}

	// find start in each reverse sequence
	for( matchI = matches.size(); matchI > 0; matchI-- )
	{
		zero_exists = false;
		for( seqI = 0; seqI < seq_count; ++seqI )
		{
			if( left_ends[seqI] == NO_MATCH && matches[ matchI - 1 ]->Orientation(seqI) == AbstractMatch::reverse )
				left_ends[seqI] = matches[matchI - 1]->LeftEnd(seqI);
			if( left_ends[seqI] == NO_MATCH )
				zero_exists = true;
		}
		if( !zero_exists )
			break;
	}

	// find end in each reverse sequence
	for( matchI = 0; matchI != matches.size(); ++matchI )
	{
		zero_exists = false;
		for( seqI = 0; seqI < seq_count; ++seqI )
		{
			if( lengths[seqI] == 0 &&
				matches[matchI]->Orientation(seqI) == AbstractMatch::reverse )
			{
					lengths[seqI] = matches[matchI]->Length(seqI)+(matches[matchI]->LeftEnd(seqI) - left_ends[seqI]);
			}
			if( lengths[seqI] == 0 )
				zero_exists = true;
		}
		if( !zero_exists )
			break;
	}
}


template<class GappedBaseImpl>
void GenericInterval<GappedBaseImpl>::addUnalignedRegions()
{
	std::list<AbstractMatch*> new_matches(matches.begin(), matches.end());

	for( uint seqI = 0; seqI < this->SeqCount(); ++seqI )
	{
		if( this->LeftEnd(seqI) == NO_MATCH )
			continue;
		if(this->Orientation(seqI) == AbstractMatch::forward)
			AddGapMatches( new_matches, new_matches.begin(), new_matches.end(), seqI, this->LeftEnd(seqI), this->RightEnd(seqI), this->Orientation(seqI), this->SeqCount() );
		else
			AddGapMatches( new_matches, new_matches.rbegin(), new_matches.rend(), seqI, this->LeftEnd(seqI), this->RightEnd(seqI), this->Orientation(seqI), this->SeqCount() );
	}
	matches.clear();
	matches.insert(matches.end(), new_matches.begin(), new_matches.end());
}

template<class GappedBaseImpl>
void GenericInterval<GappedBaseImpl>::Invert(){
	GappedBaseImpl::Invert();
	for( uint matchI = 0; matchI < matches.size(); ++matchI )
		matches[ matchI ]->Invert();

	std::reverse( matches.begin(), matches.end() );
}

template<class GappedBaseImpl>
void GenericInterval<GappedBaseImpl>::GetColumn( gnSeqI col, std::vector<gnSeqI>& pos, std::vector<bool>& column ) const
{
	size_t matchI;
	gnSeqI match_col;
	this->GetColumnAndMatch( col, pos, column, matchI, match_col );
}

template<class GappedBaseImpl>
void GenericInterval<GappedBaseImpl>::FindMatchPos( uint seqI, gnSeqI pos, size_t& matchI, gnSeqI& match_pos )
{
	match_pos = pos;
	int diff_amt = 0;
	int incr = 1;
	matchI = 0;
	size_t end_mI = matches.size();
	if( this->Orientation(seqI) == AbstractMatch::reverse )
	{
		diff_amt = -1;
		incr = -1;
		matchI = matches.size();
		end_mI = 0;
	}

	for( ; matchI != end_mI; matchI+=incr )
	{
		if( matches[matchI+diff_amt]->LeftEnd(seqI) == NO_MATCH )
			continue;
		if( matches[matchI+diff_amt]->Length(seqI) <= match_pos )
			match_pos -= matches[matchI+diff_amt]->Length(seqI);
		else
			break;
	}

	if( this->Orientation(seqI) == AbstractMatch::reverse )
		matchI--;
}

template<class GappedBaseImpl>
void GenericInterval<GappedBaseImpl>::CalculateOffset(){
	std::vector<gnSeqI> left_end( this->SeqCount(), NO_MATCH );
	std::vector<gnSeqI> length( this->SeqCount(), 0 );
	std::vector<bool> orientation;
	if( this->matches.size() > 0 )
		FindBoundaries( this->matches, left_end, length, orientation );
	for( uint seqI = 0; seqI < this->SeqCount(); seqI++ )
	{
		if( left_end[seqI] != 0 )
		{
			this->SetLeftEnd(seqI, left_end[seqI]);
			this->SetLength(length[seqI], seqI);
			if( orientation[seqI] )
				this->SetOrientation(seqI, AbstractMatch::forward);
			else
				this->SetOrientation(seqI, AbstractMatch::reverse);
		}else if( this->LeftEnd(seqI) != NO_MATCH )
		{
			this->SetLength(0, seqI);
			this->SetLeftEnd(seqI, NO_MATCH);
		}

	}

	this->CalculateAlignmentLength();
}

template<class GappedBaseImpl>
void GenericInterval<GappedBaseImpl>::SetAlignment( const std::vector< std::string >& seq_align )
{
	GappedAlignment* ga = new GappedAlignment(seq_align.size(), seq_align[0].size());
	matches.clear();
	matches.push_back(ga);
	ga->SetAlignment(seq_align);
	for( uint seqI = 0; seqI < this->SeqCount(); ++seqI )
	{
		ga->SetStart(seqI, this->Start(seqI));
		ga->SetLength(this->Length(seqI), seqI);
	}
}


/**
 * Writes this GenericInterval to the specified output stream (e.g. cout).
 */
template<class GappedBaseImpl>
std::ostream& operator<<(std::ostream& os, const GenericInterval<GappedBaseImpl>& cr){
	try{
	for( uint matchI = 0; matchI < cr.matches.size(); ++matchI ){
		const AbstractMatch* m = cr.matches[ matchI ];
		const GappedAlignment* clust = dynamic_cast< const GappedAlignment* >( m );
		if( clust != NULL )
			os << *clust;
		const Match* match = dynamic_cast< const Match* >( m );
		if( match != NULL )
			os << *match;
		os << std::endl;
	}
	}catch(...){
		std::cerr << "Exceptional handler\n";
	}
	return os;
}

template<class GappedBaseImpl>
void GenericInterval<GappedBaseImpl>::CalculateAlignmentLength()
{
	gnSeqI aln_len = 0;
	// count each match's alignment length
	for( size_t mI = 0; mI < matches.size(); ++mI )
		aln_len += matches[mI]->AlignmentLength();
	this->SetAlignmentLength(aln_len);
}

template<class GappedBaseImpl>
void GenericInterval<GappedBaseImpl>::GetAlignedSequences( gnAlignedSequences& gnas, const std::vector< genome::gnSequence* >& seq_table ) const 
{
	gnas.names.clear();
	for( uint seqI = 0; seqI < seq_table.size(); ++seqI ){
		std::string name;
		if( seq_table[ seqI ]->contigListSize() > 0 )
			name = seq_table[ seqI ]->contigName( 0 );
		gnas.names.push_back( name );
		gnas.positions.push_back(this->Start(seqI));
	}
	mems::GetAlignment( *this, seq_table, gnas.sequences );
}

template<class GappedBaseImpl>
bool GenericInterval<GappedBaseImpl>::IsGap( uint seq, gnSeqI col ) const
{
	std::vector<gnSeqI> pos;
	std::vector<bool> column;
	GetColumn(col, pos, column);
	return column[seq];
}

}

namespace std {
template<> inline
void swap( mems::Interval& a, mems::Interval& b )
{
	a.swap(b);
}
}

#endif	// __Interval_h__
