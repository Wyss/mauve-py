/*******************************************************************************
 * $Id: MaskedMemHash.cpp,v 1.3 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing, copying, and modification
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/MaskedMemHash.h"
#include <list>

using namespace std;
using namespace genome;
namespace mems {

MaskedMemHash::MaskedMemHash(){
	seq_mask = 0;
}


MaskedMemHash::MaskedMemHash(const MaskedMemHash& mh) : MemHash(mh){
	*this = mh;
}

MaskedMemHash& MaskedMemHash::operator=( const MaskedMemHash& mh ){
	seq_mask = mh.seq_mask;
	return *this;
}

MaskedMemHash* MaskedMemHash::Clone() const{
	return new MaskedMemHash(*this);
}

boolean MaskedMemHash::HashMatch(list<idmer>& match_list){
	//check that there is at least one forward component
	match_list.sort(&idmer_id_lessthan);
	// initialize the hash entry
	MatchHashEntry mhe = MatchHashEntry(seq_count, GetSar(0)->SeedLength());
	mhe.SetLength(GetSar(0)->SeedLength());
	
	//Fill in the new Match and set direction parity if needed.
	list<idmer>::iterator iter = match_list.begin();
	for(; iter != match_list.end(); iter++)
		mhe.SetStart(iter->id, iter->position + 1);
	SetDirection(mhe);
	mhe.CalculateOffset();
	uint64 match_number = 0;
	// compute "MatchNumber"
	for( uint seqI = 0; seqI < mhe.SeqCount(); seqI++ )
	{
		match_number <<= 1;
		if( mhe.Start(seqI) != NO_MATCH )
			match_number |= 1;
	}
	if( seq_mask == 0 || match_number == seq_mask )
		AddHashEntry(mhe);

	return true;
}

} // namespace mems
