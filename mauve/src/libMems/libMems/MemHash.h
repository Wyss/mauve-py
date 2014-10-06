/*******************************************************************************
 * $Id: MemHash.h,v 1.23 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _MemHash_h_
#define _MemHash_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <set>
#include <map>
#include <iostream>

#include "libMems/MatchFinder.h"
#include "libMems/Match.h"
#include "libGenome/gnException.h"
#include "libMems/MatchList.h"
#include "libMems/MatchHashEntry.h"
#include "libMems/SlotAllocator.h"
#include "boost/pool/object_pool.hpp"

namespace mems {

static const uint32 DEFAULT_MEM_TABLE_SIZE = 40000;
static const uint32 DEFAULT_REPEAT_TOLERANCE = 0;
static const uint32 DEFAULT_ENUMERATION_TOLERANCE = 1;

/**
 * MemHash implements an algorithm for finding exact matches of a certain minimal
 * length in several sequences.
 */
class MemHash : public MatchFinder{

   

public:
	MemHash();
	~MemHash();
	MemHash(const MemHash& mh);
	MemHash& operator=( const MemHash& mh );
	virtual MemHash* Clone() const;
	virtual void Clear();
	virtual void ClearSequences();
	
	/**
	 * Finds all maximal exact matches in the sequences contained by "match_list"
	 * The resulting list of matches is stored within "match_list"
	 */
	virtual void FindMatches( MatchList& match_list );
	virtual void FindMatchesFromPosition( MatchList& match_list, const std::vector<gnSeqI>& start_points );

	/**
	 * Generates exact matches for the sequences loaded into this MemHash 
	 */
	virtual boolean CreateMatches();

	/**
	 * Returns the size of the hash table being used. 
	 * @return the size of the hash table being used. 
	 */
	virtual uint32 TableSize() const {return table_size;};
	/**
	 * Sets the size of the hash table to new_table_size.
	 * @param new_table_size The new hash table size
	 */
	virtual void SetTableSize(uint32 new_table_size);
	/**
	 * Creates a new MatchList instance which contains all the matches found by calling Create().
	 */
	virtual MatchList GetMatchList() const;
	/**
	 * Places pointers to the mems that have been found into the vector mem_list
	 * @param mem_list an empty vector.
	 */
	//virtual void GetMatchList( std::vector<Match*>& mem_list ) const;
	
    /**
	* Use this to convert MatchHashEntry mem list to a generic match list type
    * converts the mem_list into the type specified by MatchListType
	*/
	template< class MatchListType >
	void GetMatchList( MatchListType& mem_list ) const;
	
	/**
	 * Returns the number of mems found 
	 * @return The number of mems found 
	 */
	virtual uint32 MemCount(){return m_mem_count;}
	/**
	 * Returns the number of mers thrown out because they were contained in an existing mem 
	 * @return The number of mers thrown out because they were contained in an existing mem 
	 */
	virtual uint32 MemCollisionCount(){return m_collision_count;}
	virtual void MemTableCount(std::vector<uint32>& table_count){table_count = mem_table_count;}
	/**
	 * Prints the number of matches in each hash table bucket to the ostream os.
	 * @param os The stream to print to.
	 */
	virtual void PrintDistribution(std::ostream& os) const;
	
	/**
	 * Reads in a list of mems from an input stream
	 * @throws A InvalidFileFormat exception if the file format is unknown or the file is corrupt
	 */
	virtual void LoadFile(std::istream& mem_file);
	/**
	 * Writes the matches stored in this MemHash out to the ostream @param mem_file.
	 */
	virtual void WriteFile(std::ostream& mem_file) const;

	/**
	 * Sets the permitted repetitivity of match seeds.  
	 * Set @param repeat_tolerance to 0 to generate MUMs, any higher setting will generate MEMs
	 * Many possible combinations of repetitive seed matches may be ignored, depending on the 
	 * setting of the repeat enumeration tolerance.
	 * @see SetEnumerationTolerance
	 * @param repeat_tolerance the permitted repetitivity of match seeds
	 */
	virtual void SetRepeatTolerance(uint32 repeat_tolerance){m_repeat_tolerance = repeat_tolerance;}
	/**
	 * @return the permitted repetitivity of match seeds.  
	 * @see SetRepeatTolerance
	 */
	virtual uint32 GetRepeatTolerance() const{return m_repeat_tolerance;}
	/**
	 * Sets the match seed repeat enumeration tolerance.
	 * When matching mers are found across sequences which also occur several times in any particular
	 * sequence there are several possible match seeds which could be generated.
	 * The enumeration tolerance controls how many of these possibilities are actually used as match
	 * seeds and extended into full matches.  The selection of actual seeds from the realm of possibilities
	 * is essentially arbitrary, though not explicitly randomized.
	 */
	virtual void SetEnumerationTolerance(uint32 enumeration_tolerance){m_enumeration_tolerance = enumeration_tolerance;}
	/**
	 * @return  the match seed repeat enumeration tolerance.
	 * @see SetEnumerationTolerance
	 */
	virtual uint32 GetEnumerationTolerance() const{return m_enumeration_tolerance;}
	
	/**
	 * Setting this to a non-null value causes matches to be logged as they are found
	 */
	void SetMatchLog( std::ostream* match_log ){ this->match_log = match_log; }

	

	//end void GetMatchList( std::vector<MatchListType*>& mem_list );

protected:
	virtual boolean EnumerateMatches( IdmerList& match_list );
	virtual boolean HashMatch(IdmerList& match_list);
	virtual void SetDirection(MatchHashEntry& mhe);
	virtual MatchHashEntry* AddHashEntry(MatchHashEntry& mhe);
	virtual uint32 quadratic_li(uint32 listI){return (listI*(listI+1))/2;}
		
	uint32 table_size;
	std::vector< std::vector<MatchHashEntry*> > mem_table;
	uint32 m_repeat_tolerance;
	uint32 m_enumeration_tolerance;
	uint64 m_mem_count;
	uint64 m_collision_count;
	std::vector<uint32> mem_table_count;

	std::ostream* match_log;
	SlotAllocator<MatchHashEntry>& allocator;
	std::vector<MatchHashEntry*> allocated;	// used to track what needs to get explicitly destroyed later...
//	boost::object_pool<MatchHashEntry> allocator;
	MheCompare mhecomp;
};


/**
 * Use this to convert MatchHashEntry mem list to a generic match list type
 * converts the mem_list into the type specified by MatchListType
 */
template< class MatchListType >
void MemHash::GetMatchList( MatchListType& mem_list ) const {
	
	mem_list.clear();
	typedef typename MatchListType::value_type MatchType;
   
	//Boost to the rescue! use remove_pointer to get at MatchListType's type
	typedef typename boost::remove_pointer<MatchType>::type SinPtrMatchType;
	SinPtrMatchType mm;

	for(uint32 i=0; i < table_size; ++i)
	{
		std::vector<MatchHashEntry*>::const_iterator iter = mem_table[i].begin();
		for(; iter != mem_table[i].end(); iter++ )
		{
			MatchType m = mm.Copy();
			*m = **iter;
			mem_list.push_back( m );
		}
	}

}


}

#endif //_MemHash_h_
