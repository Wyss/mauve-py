/*******************************************************************************
 * $Id: RepeatHash.cpp,v 1.13 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing, copying, and modification
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/RepeatHash.h"
#include <list>

using namespace std;
using namespace genome;
namespace mems {


RepeatHash* RepeatHash::Clone() const{
	return new RepeatHash(*this);
}

boolean RepeatHash::CreateMatches(){
	if(seq_count == 1){
		MatchFinder::FindMatchSeeds();
		return true;
	}

	return false;
}

boolean RepeatHash::EnumerateMatches( IdmerList& match_list ){
	return HashMatch(match_list);
}

//why have separate hash tables?
// RepeatHashEntries use GENETICIST coordinates.  They start at 1, not 0.
boolean RepeatHash::HashMatch(IdmerList& match_list){
	//check that there is at least one forward component
	match_list.sort(&idmer_position_lessthan);
	// initialize the hash entry
	MatchHashEntry mhe = MatchHashEntry( match_list.size(), GetSar(0)->SeedLength());
	mhe.SetLength( GetSar(0)->SeedLength() );
	
	//Fill in the new Match and set direction parity if needed.
	IdmerList::iterator iter = match_list.begin();

	uint32 repeatI = 0;
	for(; iter != match_list.end(); iter++)
		mhe.SetStart(repeatI++, iter->position + 1);

	SetDirection( mhe );
	mhe.CalculateOffset();
	if(mhe.Multiplicity() < 2){
		cout << "red flag " << mhe << "\n";
	}else{
		AddHashEntry(mhe);
	}
	return true;
}

} // namespace mems
