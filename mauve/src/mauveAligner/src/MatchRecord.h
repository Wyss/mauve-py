#ifndef __MatchRecord_h__
#define __MatchRecord_h__

#include "libMems/MuscleInterface.h"
#include "libMems/AbstractMatch.h"
#include "libMems/SparseAbstractMatch.h"
#include "libMems/AbstractGappedAlignment.h"
#include "libMems/Interval.h"
#include "libMems/CompactGappedAlignment.h"
#include "libMems/MatchProjectionAdapter.h"
#include <iostream>
#include <set>
#include <vector>
//#include <boost/variant.hpp>

// forward declaration
class MatchLink;
class MatchRecord;
class CompactMatchRecord;
class GappedMatchRecord;
class UngappedMatchRecord;
class CompactUngappedMatchRecord;

/** stores a link between a subset and a superset match */
class MatchLink
{
public:
	MatchLink() : superset(NULL), subset(NULL) {};
	MatchLink( MatchRecord* super, MatchRecord* sub, boost::dynamic_bitset<>& comp_list, std::vector< size_t > comp_map ) :
		superset( super ), subset( sub ), super_component_list( comp_list ), sub_to_super_map( comp_map ) {};
	void clear()
	{
		superset = NULL;
		subset = NULL;
		super_component_list.clear();
		sub_to_super_map.clear();
	}
	MatchRecord* superset;	/**< The superset match connected by this link */
	MatchRecord* subset;  /**< The subset match connected by this link */
	boost::dynamic_bitset<> super_component_list;	/**< this gets sized to be equal to superset->Multiplicity() and tracks which components of the superset are linked */
	std::vector< size_t > sub_to_super_map;	/**< mapping of subset components to superset components */
};

class MatchRecord : public mems::SparseAbstractMatch<>
{
public:
	MatchRecord() : mems::SparseAbstractMatch<>() { clear(); }
	MatchRecord( uint seq_count ): mems::SparseAbstractMatch<>( seq_count ){ clear(); }
	GappedMatchRecord* subsuming_match;
	std::vector< size_t > subsumption_component_map;
	std::vector< MatchLink > left_subset_links;			/**< Links to nearby subset matches on the left side */
	std::vector< MatchLink > right_subset_links;		/**< Links to nearby subset matches on the right side */
	MatchLink left_superset;							/**< The left-side superset, if one exists */
	MatchLink right_superset;							/**< The right-side superset, if one exists */
	std::vector< MatchLink > extra_left_subsets;		/**< left-side subsets that were further away than the first linked subset on the left side */
	std::vector< MatchLink > extra_right_subsets;		/**< right-side subsets that were further away than the first linked subset on the right side */
	std::vector< MatchRecord* > chained_matches;
	std::vector< std::vector< size_t > > chained_component_maps;	/**< maps components in this match to those in chained matches */
	bool tandem;			/**< set to true if components of the match are chainable to each other (tandem repeats)*/
	bool extended;			/**< set to false prior to extending this match */
	bool is_novel_subset;
    bool dont_extend;
    
	void clear()
	{
		subsuming_match = NULL;
		left_superset.clear();
		right_superset.clear();
		tandem = false;
		extended = false;
		dont_extend = false;
        is_novel_subset = false;
	}
};

class CompactMatchRecord : public mems::DenseAbstractMatch<1>
{
public:
	CompactMatchRecord() : mems::DenseAbstractMatch<1>() { clear(); }
	CompactMatchRecord( uint seq_count ): mems::DenseAbstractMatch<1>( seq_count ){ clear(); }
	GappedMatchRecord* subsuming_match;

	void clear()
	{
		subsuming_match = NULL;
	}
};

/**
 * An ungapped alignment that also stores a match record
 */
class CompactUngappedMatchRecord : public mems::UngappedLocalAlignment< CompactMatchRecord >
{
public:

	CompactUngappedMatchRecord(){};

	/** always set seq_count, don't worry about align_length */
	CompactUngappedMatchRecord( uint seq_count, gnSeqI align_length ) : mems::UngappedLocalAlignment< CompactMatchRecord >( seq_count )
	{
		subsuming_match = NULL;
	}

	CompactUngappedMatchRecord* Clone() const { return new CompactUngappedMatchRecord( *this ); }
	CompactUngappedMatchRecord* Copy() const;
	virtual void Free();
};

inline
CompactUngappedMatchRecord* CompactUngappedMatchRecord::Copy() const
{
	return m_allocateAndCopy( *this );
}
inline
void CompactUngappedMatchRecord::Free()
{
	m_free(this);
}


/**
 * An ungapped alignment that also stores a match record
 */
class UngappedMatchRecord : public mems::UngappedLocalAlignment< MatchRecord >
{
public:

	UngappedMatchRecord(){};

	/** always set seq_count, don't worry about align_length */
	UngappedMatchRecord( uint seq_count, gnSeqI align_length ) : mems::UngappedLocalAlignment< MatchRecord >( seq_count )
	{
		subsuming_match = NULL;
	}

	UngappedMatchRecord* Clone() const { return new UngappedMatchRecord( *this ); }
	UngappedMatchRecord* Copy() const;
	virtual void Free();

	friend std::ostream& operator<<(std::ostream& os, const UngappedMatchRecord& mr); //write to source.
};

inline
UngappedMatchRecord* UngappedMatchRecord::Copy() const
{
	return m_allocateAndCopy( *this );
}
inline
void UngappedMatchRecord::Free()
{
	m_free(this);
}


/**
 * The gapped match record class.  Abuses the Interval class to store a chain of other matches
 */
class GappedMatchRecord : public mems::GenericInterval< mems::AbstractGappedAlignment< MatchRecord > >
{
public:

	/** always set seq_count, don't worry about align_length */
	GappedMatchRecord() : 
	  mems::GenericInterval< mems::AbstractGappedAlignment< MatchRecord > >()
	{
        //tjt: initialize spscore to 0
        spscore = 0;  
    }

	GappedMatchRecord( UngappedMatchRecord& umr )
	{
        //tjt: initialize spscore to 0
        spscore = 0;
		std::vector<UngappedMatchRecord*> asdf(1, &umr);
		mems::GenericInterval< mems::AbstractGappedAlignment< MatchRecord > > iv( asdf.begin(), asdf.end() );
		mems::GenericInterval< mems::AbstractGappedAlignment< MatchRecord > >::operator=( iv );
		MatchRecord::operator=( umr );
	}

	/** 
	 * Call to indicate that all matches have been placed in the chained_matches list and can be 
	 * converted to a gapped alignment
	 */
	void finalize(std::vector<genome::gnSequence *> seq_table);
    //tjt: should this go somewhere else?
    mems::score_t spscore;
// methods inherited from AbstractGappedAlignment
public:
	GappedMatchRecord* Clone() const { return new GappedMatchRecord( *this ); }
	GappedMatchRecord* Copy() const;
	virtual void Free();

	friend std::ostream& operator<<(std::ostream& os, const GappedMatchRecord& mr); //write to source.
};

inline
GappedMatchRecord* GappedMatchRecord::Copy() const
{
	return m_allocateAndCopy( *this );
}
inline
void GappedMatchRecord::Free()
{
	m_free(this);
}


/** orders on increasing multiplicity */
typedef std::pair< MatchRecord*, std::vector< size_t >* > MatchSortEntry;
class MatchSortEntryCompare
{
public:
	bool operator()( const MatchSortEntry& a, const MatchSortEntry& b )
	{
		return a.first->Multiplicity() < b.first->Multiplicity();
	}
};

template< typename T >
class IsNullPtr
{
public:
	bool operator()( const T* a ){ return a == NULL; }
};

void GappedMatchRecord::finalize( std::vector<genome::gnSequence *> seq_table)
{
	std::vector< mems::AbstractMatch* > iv_matches;
	MatchSortEntryCompare msec;
	std::vector< MatchSortEntry > mse_list( chained_matches.size() );

	for( size_t cI = 0; cI < chained_matches.size(); ++cI )
	{
		mse_list[cI].first = chained_matches[cI];
		mse_list[cI].second = &chained_component_maps[cI];
	}
	std::sort( mse_list.begin(), mse_list.end(), msec );
	// add lowest multiplicity matches first, progressively add higher mult. matches
	std::vector< mems::AbstractMatch* > chain;
	for( size_t cI = 0; cI < mse_list.size(); ++cI )
	{
		mems::MatchProjectionAdapter mpaa( mse_list[cI].first, *(mse_list[cI].second) );
		// clobber any region that overlaps with this mpaa
		for( size_t seqI = 0; seqI < mpaa.SeqCount(); seqI++ )
		{
			size_t csize = chain.size();
			for( size_t mI = 0; mI < csize; mI++ )
			{
				mems::AbstractMatch* m = chain[mI];
				if( m == NULL )
					continue;
                if (m->LeftEnd(seqI) == 0 && m->Length( seqI ) == 0)
                    continue; //should we throw error here?
				if( m->RightEnd(seqI) < mpaa.LeftEnd(seqI) )
					continue;	// no overlap here!
				if( m->LeftEnd(seqI) > mpaa.RightEnd(seqI) )
					continue;	// no overlap, woohoo!
				if( m->LeftEnd(seqI) < mpaa.LeftEnd(seqI) &&
					m->RightEnd(seqI) >= mpaa.LeftEnd(seqI) )
				{
	                // take the part of m to the left of mpaa and put it at the end of our chain
					mems::AbstractMatch* m_left = m->Copy();
					m_left->CropRight( m_left->RightEnd(seqI) - mpaa.LeftEnd(seqI) + 1, seqI );
					m->CropLeft( m_left->Length(seqI), seqI );
                    chain.push_back(m_left);
				}
				// now m is guaranteed to have left-end >= mpaa
				if( m->RightEnd(seqI) <= mpaa.RightEnd(seqI) )
				{
					// m is completely contained inside mpaa, so get rid of it
					m->Free();
					chain[mI] = NULL;
					continue;
				}

    			m->CropLeft( mpaa.RightEnd(seqI) - m->LeftEnd(seqI) + 1, seqI );		
			}
		}
		// get rid of any null entries in the chain
		std::vector< mems::AbstractMatch* >::iterator end_iter = std::remove( chain.begin(), chain.end(), (AbstractMatch*)NULL );
		chain.erase( end_iter, chain.end() );
		chain.push_back( mpaa.Copy() );
		if( chain.back()->Orientation(0) == AbstractMatch::reverse )
			chain.back()->Invert();
	}
	
	if( chain.size() == 0 )
	{
		*this = GappedMatchRecord();
		return;
	}
	mems::MatchStartComparator< mems::AbstractMatch > asc(0);
	std::sort( chain.begin(), chain.end(), asc );
	// aed: At this point the matches in chain are in sorted order, so the region betweeen each of them is what should get fed to muscle
	//      will need to feed AbstractMatch instead of Match to MuscleInterface::Align though
	std::vector< mems::AbstractMatch* >::iterator chain_begin = chain.begin();
	uint chainsize = chain.size()-1;

	if (1)
	{
	try{
	for( uint i = 0; i < chainsize; i++ )
	{
		mems::GappedAlignment* cr = NULL;
		boolean align_success = false;
		// attempt a muscle alignment
		cr = new mems::GappedAlignment();
		mems::AbstractMatch* m1 = chain.at(i);
		mems::AbstractMatch* m2 = chain.at(i+1);
	
		align_success = mems::MuscleInterface::getMuscleInterface().Align( *cr,  m1 , m2,  seq_table );
		if( align_success )
		{
            //cerr << "muscle alignment success!!" << endl;
			iv_matches.push_back( cr );
			// aed: just insert the resulting GappedAlignment objects into chain
			chain.insert(chain.begin()+(i+1), cr);
			chainsize++;
			// tjt: skip over newly inserted item
			i++;		
		}
        else
            continue;
		
	}
	
	}catch( genome::gnException& gne ){
		std::cerr << gne << std::endl;
	}catch(std::exception& e){
		std::cerr << e.what() << std::endl;
		std::cerr << chain.size() << std::endl;
	}catch(...){
		std::cerr << "matrix exception?\n";
	}
	}

	MatchRecord* mr = this->Copy();
	SetMatches( chain );
	//tjt: now chain should be empty
	// don't keep a potentially huge tree of GappedMatchRecords.  instead, flatten to a single cga
	mems::CompactGappedAlignment<> tmpcga(*this);
	chain.push_back(tmpcga.Copy());
	SetMatches( chain );
	//tjt: assign this to slot allocated & copied MatchRecord
	MatchRecord::operator=(*mr);
	mr->Free();
}

std::ostream& operator<<(std::ostream& os, const UngappedMatchRecord& ula);
std::ostream& operator<<(std::ostream& os, const UngappedMatchRecord& ula){ //write to stream.
	os << ula.AlignmentLength();
	for(uint i=0; i < ula.SeqCount(); i++)
		os << '\t' << ula.Start(i);
	return os;
}

std::ostream& operator<<(std::ostream& os, const GappedMatchRecord& ula);
std::ostream& operator<<(std::ostream& os, const GappedMatchRecord& ula){ //write to stream.
	os << ula.AlignmentLength();
	for(uint i=0; i < ula.SeqCount(); i++)
		os << '\t' << ula.Start(i);
	os << "\nlens:";
	for(uint i=0; i < ula.SeqCount(); i++)
		os << '\t' << ula.Length(i);

	return os;
}

#endif // __MatchRecord_h__
