/*******************************************************************************
 * $Id: UniqueMatchFinder.cpp,v 1.13 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "UniqueMatchFinder.h"
#include <list>

using namespace std;
using namespace genome;
using namespace mems;

UniqueMatchFinder::UniqueMatchFinder(){
}

UniqueMatchFinder::~UniqueMatchFinder(){
}

UniqueMatchFinder::UniqueMatchFinder(const UniqueMatchFinder& mh) : MemHash(mh){

}

UniqueMatchFinder* UniqueMatchFinder::Clone() const{
	return new UniqueMatchFinder(*this);
}


// enumerate out every pairwise match
boolean UniqueMatchFinder::EnumerateMatches( IdmerList& match_list ){

	match_list.sort(&idmer_id_lessthan);
	IdmerList::iterator iter = match_list.begin();
	IdmerList::iterator iter2 = match_list.begin();
	uint cur_id_count = 1;
	IdmerList unique_list;
	// identify all of the unique seeds and add them to unique_list
	while(iter2 != match_list.end()){
		++iter2;
		if(iter2 == match_list.end() || iter->id != iter2->id){
			if( cur_id_count == 1 )
				unique_list.push_back( *iter );
			else
				cur_id_count = 1;
		}else
			cur_id_count++;
		++iter;
	}
	// hash all unique seeds
	boolean success = true;
	if( unique_list.size() >= 2 )
		success = HashMatch(unique_list);
	return success;
}
