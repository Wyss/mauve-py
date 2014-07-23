/*******************************************************************************
 * $Id: Match.h,v 1.10 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef __MatchHashEntry_h__
#define __MatchHashEntry_h__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnClone.h"
#include <iostream>
#include <set>
#include "libMems/Match.h"

namespace mems {

/**
 * The Match class stores the location of an <b>equal size</b> (inexact or exactly) 
 * matching region
 * between several sequences.  There are numerous functions in this
 * class which can be used to compare and manipulate this match.
 */

class MatchHashEntry : public Match
{
public:
	enum MemType
	{
		seed,
		extended
	};

public:
	MatchHashEntry();
	/**
	 * Creates a new Match.
	 * @param seq_count The total number of sequences in the alignment
	 * @param mersize The size of the mers used in the sorted mer lists.
	 * @param m_type The type of mem to create, can either be a seed or already extended.
	 * @see MemType
	 */
	MatchHashEntry( const uint seq_count, const gnSeqI mersize, const MemType m_type = seed );
	MatchHashEntry* Clone() const;
	MatchHashEntry* Copy() const;
	virtual void Free();
	MatchHashEntry( const MatchHashEntry& mhe ){ *this = mhe; }
	MatchHashEntry& operator=(const MatchHashEntry& mhe);

	/** comparison operator, compares two matches to see if they are the same */
	boolean operator==(const MatchHashEntry& mhe) const;


	/** @return true if this match has already been extended */
	boolean Extended() const{return m_extended;}
	/** Sets this match to be extended if the value passed in "extended" is true */
	void SetExtended(boolean extended){m_extended = extended;}
	/** @return the mer size of the sorted mer lists used to find this match */
	uint MerSize() const{return m_mersize;}

	/**
	 * Calculates the generalized offset and other bookkeeping information
	 * for this mem.  This should <b>always</b> be called after changing the start
	 * positions of the mem.
	 */
	virtual void CalculateOffset();
	
	/** Returns the generalized offset of this match */
	int64 Offset() const{return m_offset;};

	/** Sets the generalized offset of this match to "offset" */
	void SetOffset(int64 offset){m_offset = offset;};		

	static boolean offset_lessthan(const MatchHashEntry& a, const MatchHashEntry& b);
	static boolean start_lessthan_ptr(const MatchHashEntry* a, const MatchHashEntry* b);
	static bool start_lessthan(const MatchHashEntry& a, const MatchHashEntry& b);
	static boolean strict_start_lessthan_ptr(const MatchHashEntry* a, const MatchHashEntry* b);
	/** compare the end of a to the start of b 
	 */
	static int64 end_to_start_compare(const MatchHashEntry& a, const MatchHashEntry& b);
	static int64 start_compare(const MatchHashEntry& a, const MatchHashEntry& b);

	/**
	 *	Will return true if this match contains mhe
	 *  Containment implies that a match has a length >= the contained
	 *  match, it has coordinates in every genome the contained match has,
	 *  the difference in start positions in each genome is the same.
	 * @param mhe The match to check for containment.
	 * @return True if this match contains mhe.
	 */
	boolean Contains(const MatchHashEntry& mhe) const;

private:

	boolean m_extended;
	gnSeqI m_mersize;
	int64 m_offset;
};

inline
MatchHashEntry* MatchHashEntry::Copy() const
{
	return m_allocateAndCopy(*this);
}
inline
void MatchHashEntry::Free()
{
	m_free(this);
}

inline
bool MatchHashEntry::start_lessthan(const MatchHashEntry& a, const MatchHashEntry& b){
	return start_lessthan_ptr(&a, &b);
}

class MheCompare {
public:
	bool operator()(const MatchHashEntry* a, const MatchHashEntry* b) const{
		if( a->FirstStart() > b->FirstStart() ){
			return true;
		}else if( a->FirstStart() == b->FirstStart() ){
			// check that the matches hit the same genomes
			for( size_t i = a->FirstStart(); i < a->SeqCount(); i++ )
			{
				if( a->LeftEnd(i) == NO_MATCH && b->LeftEnd(i) != NO_MATCH )
					return true;
				else if( a->LeftEnd(i) != NO_MATCH && b->LeftEnd(i) == NO_MATCH )
					return false;
			}
			//offsets are the same, check for containment...
			if(a->Contains(*b) || b->Contains(*a)){
				return false;
			}else
				return MatchHashEntry::strict_start_lessthan_ptr(a, b);
		}
		return false;
	}
};

}

#endif // __MatchHashEntry_h__
