/*******************************************************************************
 * $Id: PairwiseMatchFinder.cpp,v 1.13 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing, copying, and modification
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/PairwiseMatchFinder.h"
#include <list>

using namespace std;
using namespace genome;

namespace mems {

PairwiseMatchFinder::PairwiseMatchFinder(){
}

PairwiseMatchFinder::~PairwiseMatchFinder(){
}

PairwiseMatchFinder::PairwiseMatchFinder(const PairwiseMatchFinder& mh) : MemHash(mh){

}

PairwiseMatchFinder* PairwiseMatchFinder::Clone() const{
	return new PairwiseMatchFinder(*this);
}


// enumerate out every pairwise match
boolean PairwiseMatchFinder::EnumerateMatches( IdmerList& match_list ){

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
	// hash each pair of unique seeds
	boolean success = true;
	for( iter = unique_list.begin(); iter != unique_list.end(); ++iter )
	{
		for( iter2 = iter; iter2 != unique_list.end(); ++iter2 )
		{
			if( iter == iter2 )
				continue;
			IdmerList hash_list;
			hash_list.push_back( *iter );
			hash_list.push_back( *iter2 );
			success = success && HashMatch(hash_list);
		}
	}
	return success;
}

}  // namespace mems
