/*******************************************************************************
 * $Id: Match.h,v 1.10 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _RepeatMatch_h_
#define _RepeatMatch_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnClone.h"
#include <iostream>
#include <vector>
#include <set>
#include "libMems/MatchHashEntry.h"

namespace mems {

/**
 * The Match class stores the location of an <b>equal size</b> (inexact or exactly) 
 * matching region
 * between several sequences.  There are numerous functions in this
 * class which can be used to compare and manipulate this match.
 */
class RepeatMatch : public MatchHashEntry {

public:
	RepeatMatch();
	RepeatMatch( const uint32 seq_count, const gnSeqI mersize, const MemType m_type = seed );
	RepeatMatch(const RepeatMatch& mhe);
	~RepeatMatch();
	void FromSeq( uint32 match_id, uint32 seq_id );
	uint32 SeqId( uint32 match_id );
protected:
	std::vector<uint32> m_seq_id;

private:


};
std::ostream& operator<<(std::ostream& os, const RepeatMatch& mhe); //write to source.

}	// namespace mems

#endif // _RepeatMatch_h_

