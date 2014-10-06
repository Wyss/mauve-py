/*******************************************************************************
 * $Id: Match.cpp,v 1.9 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing, copying, and modification
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/RepeatMatch.h"
#include "libGenome/gnException.h"
#include "libGenome/gnDebug.h"

namespace mems {

RepeatMatch::RepeatMatch() : MatchHashEntry()
{
}

RepeatMatch::~RepeatMatch(){

}

void RepeatMatch::FromSeq( uint32 match_id, uint32 seq_id )
{
// unsure what to do with this:  (it doesn't compile)
//	this->m_seq_id.insert( match_id, seq_id);

}

uint32 RepeatMatch::SeqId( uint32 match_id )
{
	return this->m_seq_id.at(match_id);

}

std::ostream& operator<<(std::ostream& os, const RepeatMatch& mhe){ //write to stream.
	os << mhe.Length();
	for(uint32 i=0; i < mhe.SeqCount(); i++)
	{
		
		//if ( mhe.Start(i) < 
		os << '\t' << mhe.Start(i);
	}
	return os;
}

}	// namespace mems
