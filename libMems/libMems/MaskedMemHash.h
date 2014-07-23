/*******************************************************************************
 * $Id: MaskedMemHash.h,v 1.3 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _MaskedMemHash_h_
#define _MaskedMemHash_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/MemHash.h"

namespace mems {

/**
 * Finds matches that meet a particular sequence mask, e.g. 0b11111 for 5-way matches
 * Doesn't filter anything unless a mask is set using SetMask().  The
 * filter can be cleared by calling SetMask(0)
 */
class MaskedMemHash : public MemHash{
public:
	MaskedMemHash();
	~MaskedMemHash(){};
	MaskedMemHash(const MaskedMemHash& mh);
	MaskedMemHash& operator=( const MaskedMemHash& mh );
	virtual MaskedMemHash* Clone() const;
	virtual void SetMask( uint64 seq_mask ){ this->seq_mask = seq_mask; }
protected:
	/**
	 * Can't find subsets when there is only one permitted sequence mask!
	 */
	virtual void FindSubsets(const Match& mhe, std::vector<Match>& subset_matches){};
	virtual boolean HashMatch(std::list<idmer>& match_list);
	uint64 seq_mask;
};

}

#endif //_MaskedMemHash_h_
