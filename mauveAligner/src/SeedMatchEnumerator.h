#ifndef __SeedMatchEnumerator_h__
#define __SeedMatchEnumerator_h__

#include "libMems/MatchFinder.h"
#include "libMems/RepeatHash.h"
#include "libMems/MemHash.h"
#include "libMems/MatchList.h"
#include "libMems/SortedMerList.h"
#include "libMems/Match.h"

/**
 * Turns every seed match into a full match without extension.
 */
class SeedMatchEnumerator : public mems::MatchFinder 
{
public:
	virtual SeedMatchEnumerator* Clone() const;

	void FindMatches( mems::MatchList& match_list, size_t min_multi = 2, size_t max_multi = 1000, bool direct_repeats_only = false )
	{
        this->max_multiplicity = max_multi;
        this->min_multiplicity = min_multi;
        this->only_direct = direct_repeats_only;
		for( size_t seqI = 0; seqI < match_list.seq_table.size(); ++seqI ){
			if( !AddSequence( match_list.sml_table[ seqI ], match_list.seq_table[ seqI ] ) ){
				genome::ErrorMsg( "Error adding " + match_list.seq_filename[seqI] + "\n");
				return;
			}
		}
		CreateMatches();
		match_list.clear();
		match_list.insert( match_list.end(), mlist.begin(), mlist.end() );
	}

	virtual boolean CreateMatches();
protected:

	virtual boolean EnumerateMatches( mems::IdmerList& match_list );
	virtual boolean HashMatch(mems::IdmerList& match_list);
	virtual mems::SortedMerList* GetSar(uint32 sarI) const;
	mems::MatchList mlist;
	void SetDirection(mems::Match& mhe);
private:
    //used to store rmin, rmax values
    size_t max_multiplicity;
    size_t min_multiplicity;
	bool only_direct;
};

SeedMatchEnumerator* SeedMatchEnumerator::Clone() const{
	return new SeedMatchEnumerator(*this);
}

inline
mems::SortedMerList* SeedMatchEnumerator::GetSar(uint32 sarI) const{
	return sar_table[0];
}

boolean SeedMatchEnumerator::CreateMatches(){
	if(seq_count == 1){
		MatchFinder::FindMatchSeeds();
		return true;
	}
	return false;
}

boolean SeedMatchEnumerator::EnumerateMatches( mems::IdmerList& match_list ){
	return HashMatch(match_list);
}

boolean SeedMatchEnumerator::HashMatch(mems::IdmerList& match_list){
	//check that there is at least one forward component
	match_list.sort(&mems::idmer_position_lessthan);
	// initialize the hash entry
	mems::Match mhe = mems::Match( match_list.size() );
	mhe.SetLength( GetSar(0)->SeedLength() );
	
	//Fill in the new Match and set direction parity if needed.
	mems::IdmerList::iterator iter = match_list.begin();
    
	uint32 repeatI = 0;
	for(; iter != match_list.end(); iter++)
		mhe.SetStart(repeatI++, iter->position + 1);

	SetDirection( mhe );
	bool found_reverse = false;
	vector< size_t > component_map;
	if(this->only_direct)
	{
		for( uint seqI = 0; seqI < mhe.Multiplicity(); seqI++)
		{
			if (mhe.Orientation(seqI) == 0)
				component_map.push_back(seqI);
			else
				found_reverse = true;
		}
	}
	mems::MatchProjectionAdapter mpaa(mhe.Copy(),  component_map);
	if(mhe.Multiplicity() < 2){
		std::cerr << "red flag " << mhe << "\n";
    }
    //use rmin & rmax to discard irrelevant seed matches
    else if(mhe.Multiplicity() > this->max_multiplicity || mhe.Multiplicity() < this->min_multiplicity )
    {
        ;
    }
	else if(this->only_direct && found_reverse)
	{
		if ( mpaa.Multiplicity() > 1)
		{
			mems::Match new_mhe = mems::Match( mpaa.Multiplicity() );
			new_mhe.SetLength( GetSar(0)->SeedLength() );
			for(uint mult = 0; mult < mpaa.Multiplicity(); mult++)
				new_mhe.SetStart(mult, mpaa.Start(mult));
			mlist.push_back(new_mhe.Copy());
		}
	}
    else{
		mlist.push_back(mhe.Copy());
		
	}
	return true;
}

// evil, evil code duplication.

void SeedMatchEnumerator::SetDirection(mems::Match& mhe){
	//get the reference direction
	boolean ref_forward = false;
	uint32 seqI=0;
	for(; seqI < mhe.SeqCount(); ++seqI)
		if(mhe[seqI] != mems::NO_MATCH){
			ref_forward = !(GetSar(seqI)->GetMer(mhe[seqI] - 1) & 0x1);
			break;
		}
	//set directional parity for the rest
	for(++seqI; seqI < mhe.SeqCount(); ++seqI)
		if(mhe[seqI] != mems::NO_MATCH)
			if(ref_forward == (GetSar(seqI)->GetMer(mhe[seqI] - 1) & 0x1))
				mhe.SetStart(seqI, -mhe[seqI]);
}


#endif	// __SeedMatchEnumerator_h__
