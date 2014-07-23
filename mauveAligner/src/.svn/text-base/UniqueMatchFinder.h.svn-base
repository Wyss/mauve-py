/*******************************************************************************
 * $Id: UniqueMatchFinder.h,v 1.8 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
 ******************************************************************************/

#ifndef _UniqueMatchFinder_h_
#define _UniqueMatchFinder_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/MemHash.h"

/**
 * Finds all pairwise matches with unique seeds among a group of sequences
 */
class UniqueMatchFinder : public mems::MemHash
{
public:
	UniqueMatchFinder();
	~UniqueMatchFinder();

	UniqueMatchFinder(const UniqueMatchFinder& mh);
	virtual UniqueMatchFinder* Clone() const;
protected:

	virtual boolean EnumerateMatches( mems::IdmerList& match_list );
};

#endif //_UniqueMatchFinder_h_
