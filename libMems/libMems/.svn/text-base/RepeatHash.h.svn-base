/*******************************************************************************
 * $Id: RepeatHash.h,v 1.8 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _RepeatHash_h_
#define _RepeatHash_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/MemHash.h"

namespace mems {

/**
 * Finds repeats within a single sequence.
 * This class extends the functionality of MemHash to search for repetitive
 * matches within a single sequence.
 */
class RepeatHash : public MemHash{
public:
	virtual RepeatHash* Clone() const;
	virtual boolean CreateMatches();
protected:

	virtual boolean EnumerateMatches( IdmerList& match_list );
	virtual boolean HashMatch(IdmerList& match_list);
	virtual SortedMerList* GetSar(uint32 sarI) const;
};


inline
SortedMerList* RepeatHash::GetSar(uint32 sarI) const{
	return sar_table[0];
}

inline
bool idmer_greaterthan(idmer& a_v, idmer& m_v){
	return (a_v.mer < m_v.mer);// ? true : false;
};

inline
bool idmer_position_lessthan(idmer& a_v, idmer& m_v){
	return (a_v.position < m_v.position);// ? true : false;
};

}

#endif //_RepeatHash_h_
