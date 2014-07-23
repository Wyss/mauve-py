#include "libGenome/gnSequence.h"
#include "libMems/Interval.h"
#include "libMems/CompactGappedAlignment.h"
#include "libMems/Islands.h"
#include "libMems/Aligner.h"
#include "libMems/MuscleInterface.h"
#include "libGenome/gnFASSource.h"
#include "libMems/Backbone.h"
#include "libMems/ProgressiveAligner.h"
#include "libMems/HomologyHMM/parameters.h"

#include <iomanip>
#include <iostream>
#include <algorithm>
#include <cctype>

#include "MatchRecord.h"
#include "SeedMatchEnumerator.h"
//#include "procrastUtilities.h"

#include <boost/tuple/tuple.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace std;
using namespace genome;
using namespace mems;

bool print_warnings = false;

enum rvalue { OK=0, FAILED=1, DONE=2, NOVEL=3, FIXME=110}; 
int scoredropoff_matrix[10] = {0,0,4756,9144,13471,17981,25302,30945,38361,40754};
int ccount = 0;
/** A Match Position Entry stores a pointer to a match and its component for a given sequence coordinate */
typedef std::pair< MatchRecord*, size_t >	MatchPositionEntry;
/** the Match Position Lookup Table should be sized match the length of the sequence */
typedef vector< MatchPositionEntry > MatchPositionLookupTable;

/** This class stores a single entry in the neighborhood list */
class NeighborhoodListEntry
{
public:
	MatchRecord* match;
	bool relative_orientation;	/** true for identical (1) and false for opposite (-1) */
	size_t Mi_component;	/** the x value in the paper (Matching component of M_i)*/
	size_t distance;		/** the d value from the paper */
	size_t Mj_component;	/** the y value in the paper (Matching component of M_j) */
};

/** Used to sort the neighborhood list using std::sort */
class NeighborhoodListComparator
{
public:
	bool operator()( const NeighborhoodListEntry& a, const NeighborhoodListEntry& b )
	{
		if( a.match != b.match )
			return a.match < b.match;
		if( a.relative_orientation != b.relative_orientation )
			return a.relative_orientation == false;
		if( a.Mi_component != b.Mi_component )
			return a.Mi_component < b.Mi_component;
		return a.distance < b.distance; 
	}
};


bool scorecmp( GappedMatchRecord* a, GappedMatchRecord* b ) 
{
   // sort first by multipicity, then by spscore
   if( a->Multiplicity() > b->Multiplicity())
       return true;
   else if ( a->Multiplicity() < b->Multiplicity())
       return false;
   else
       return a->spscore > b->spscore;
 }

bool score_by_sp( GappedMatchRecord* a, GappedMatchRecord* b ) 
{
   // sort first by multipicity, then by spscore
   if( a->spscore > b->spscore)
       return true;
   else if ( a->spscore < b->spscore)
       return false;
   else
       return a->Multiplicity() > b->Multiplicity();
 }

bool score_by_length( GappedMatchRecord* a, GappedMatchRecord* b ) 
{
   // sort first by multipicity, then by spscore
   if( a->AlignmentLength() > b->AlignmentLength())
       return true;
   else if ( a->AlignmentLength() < b->AlignmentLength())
       return false;
   else
       return a->spscore > b->spscore;
 }
/** The NeighborhoodGroup contains the MatchRecord pointer, the component map to the match being extended (M_i), and a vector of distances to M_i*/
typedef boost::tuple< MatchRecord*, std::vector< size_t >, std::vector< size_t > > NeighborhoodGroup;

class NeighborhoodGroupComponentCompare
{
public:
	bool operator()( const NeighborhoodGroup& a, const NeighborhoodGroup& b ) const
	{
		return compare(a,b) < 0;
	}
	int compare( const NeighborhoodGroup& a, const NeighborhoodGroup& b ) const
	{
	// compare component map vectors
		// todo: make these buffers persistent to avoid reallocation!!
		vector< size_t > ac(a.get<1>());
		vector< size_t > bc(b.get<1>());
		std::sort(ac.begin(), ac.end());
		std::sort(bc.begin(), bc.end());
		size_t i = 0;
		for( ; i < ac.size() && i < bc.size(); ++i )
		{
			if( ac[i] != bc[i] )
				return ac[i] - bc[i];
		}
		if( i < ac.size() && ac[i] != (std::numeric_limits<size_t>::max)())
			return 1;
		else if( i < bc.size() && bc[i] != (std::numeric_limits<size_t>::max)())
			return -1;

		return 0;
	}
};

class NeighborhoodGroupCompare
{
public:
	bool operator()( const NeighborhoodGroup& a, const NeighborhoodGroup& b )
	{
		int cval = srcc.compare(a,b);
		if( cval != 0 )
			return cval < 0;

	// compare distance vectors
		vector< size_t > ad(a.get<2>());
		vector< size_t > bd(b.get<2>());
		std::sort(ad.begin(), ad.end());
		std::sort(bd.begin(), bd.end());
		size_t i = 0;
		for( ; i < ad.size() && i < bd.size(); ++i )
		{
			if( ad[i] != bd[i] )
				return ad[i] < bd[i];
		}
		if( i < ad.size() )
			return false;
		else if( i < bd.size() )
			return true;

		return false;
	}
protected:
	NeighborhoodGroupComponentCompare srcc;
};


//function to test if a chainable match is OK, i.e. none of this stuff:
// |---m1--->    |----c1--->               |---c2---->  |----m2---->
bool testChainableMatch( MatchRecord* M_i, MatchRecord* M_m, const vector< size_t >& component_map )
{
	bool ok = true;
	// set range to cover M_m
    int right_count = 0;
    int left_count = 0;
	for( size_t x = 0; x < M_i->Multiplicity(); ++x )
	{
		size_t z = component_map[x];

        //if there is no match, we've got a problem
		if( M_i->LeftEnd(x) == NO_MATCH || M_m->LeftEnd(z) == NO_MATCH )
			genome::breakHere();

        //should it be allowed to chain with matches with differing orientation
        //if so, how do we align the gap between these two matches?
        if (M_m->Orientation(z) != M_i->Orientation(x) )
        {
            left_count = 1;
            right_count = 1;
            break;
        }
		int64 lend_diff = M_m->LeftEnd(z) -M_i->LeftEnd(x);
		int64 rend_diff = M_m->RightEnd(z) - M_i->RightEnd(x);
        /// <---m1---->   <----b1--->    <----b2---->  <----m2---->
		if(  rend_diff < 0 && lend_diff < 0)
		{
            //component to chain is to the left of current match component
            if (M_m->Orientation(z) == AbstractMatch::forward)
                left_count++;
            else
                right_count++;
            ok = false;
        }
        else if ( rend_diff > 0 && lend_diff > 0)
        {
            if (M_m->Orientation(z) == AbstractMatch::forward)
                right_count++;
            else   
                left_count++;
            //component to chain is to the right of current match component
            ok = false;
        }
        else
        {
            left_count = 1;
            right_count = 1;
            break;

        }
	}

    //if there are components to the left && right, things are not ok with this chained match
    if (left_count != 0 && right_count != 0)
        ok = false;
    else
        ok = true;
	return ok;
}

bool extendRange( MatchRecord* M_i, MatchRecord* M_m, const vector< size_t >& component_map )
{
	bool changed = false;
		// set range to cover M_m
	for( size_t x = 0; x < M_i->Multiplicity(); ++x )
	{
		size_t z = component_map[x];
		if( M_i->LeftEnd(x) == NO_MATCH || M_m->LeftEnd(z) == NO_MATCH )
			genome::breakHere();
		int64 lend_diff = M_i->LeftEnd(x) - M_m->LeftEnd(z);
		if( lend_diff > 0 )
		{
            if ( M_i->LeftEnd(x) - lend_diff == 0)
                cerr << "extendRange debugme" << endl;
			M_i->SetLeftEnd(x, M_i->LeftEnd(x) - lend_diff);
			M_i->SetLength(M_i->Length(x)+lend_diff, x);
			changed = true;
		}

		int64 rend_diff = M_m->RightEnd(z) - M_i->RightEnd(x);
		if( rend_diff > 0 )
		{
			M_i->SetLength( M_i->Length(x)+rend_diff, x );
			changed = true;
		}

	}
	return changed;
}

bool reduceRange( MatchRecord* M_i, MatchRecord* M_m, const vector< size_t >& component_map )
{
	bool changed = false;
		// set range to cover M_m
	for( size_t x = 0; x < M_i->Multiplicity(); ++x )
	{
		size_t z = component_map[x];
		if( M_i->LeftEnd(x) == NO_MATCH || M_m->LeftEnd(z) == NO_MATCH )
			genome::breakHere();
		int64 lend_diff = M_m->LeftEnd(z) - M_i->LeftEnd(x);
		if( lend_diff > 0 )
		{
            if ( M_i->LeftEnd(x) - lend_diff == 0)
                cerr << "reduceRange debugme" << endl;
			M_i->SetLeftEnd(x, M_i->LeftEnd(x) - lend_diff);
			M_i->SetLength(M_i->Length(x)+lend_diff, x);
			changed = true;
		}
		int64 rend_diff = M_i->RightEnd(x) - M_m->RightEnd(z) ;
		if( rend_diff > 0 )
		{
			M_i->SetLength( M_i->Length(x)+rend_diff, x );
			changed = true;
		}
	}
	return changed;
}

void remapComponents(const vector< size_t >& srcmap, size_t mid_multiplicity, const vector< size_t >& destmap, vector< size_t >& newmap )
{
	vector< size_t > super_map( mid_multiplicity, (std::numeric_limits<size_t>::max)() );
	for( size_t mapI = 0; mapI < destmap.size(); ++mapI )
		super_map[destmap[mapI]] = mapI;
	for( size_t mapI = 0; mapI < srcmap.size(); ++mapI )
		newmap[mapI] = super_map[srcmap[mapI]];
}

void classifyMatch( AbstractMatch* M_i, AbstractMatch* M_j, vector< size_t >& ji_component_map, bool& subsumed, bool& partial, bool superset = false )
{
	subsumed = true;
	partial = false;
	for( size_t i = 0; i < ji_component_map.size(); ++i )
	{
		size_t x = ji_component_map[i];
		size_t y = i;
		int64 lend_diff = M_i->LeftEnd(x) - M_j->LeftEnd(y);
		int64 rend_diff = M_j->RightEnd(y) - M_i->RightEnd(x);
		if (superset)
		{
			lend_diff =  M_j->LeftEnd(y) - M_i->LeftEnd(x); 
			rend_diff =  M_i->RightEnd(x)-M_j->RightEnd(y);
		}
		
		if( lend_diff > 0 || rend_diff > 0 )
			subsumed = false;
		if( lend_diff <= 0 && rend_diff <= 0 )
			partial = true;
	}
}
//same as classifySubset, except for supersets
void classifySuperset( MatchRecord* M_i, NeighborhoodGroup& sr, bool& subsumed, bool& partial )
{
	classifyMatch( M_i, sr.get<0>(), sr.get<1>(), subsumed, partial, true );
}

void classifySubset( MatchRecord* M_i, NeighborhoodGroup& sr, bool& subsumed, bool& partial )
{
	classifyMatch( M_i, sr.get<0>(), sr.get<1>(), subsumed, partial, false );
}

void checkLink( MatchRecord*& mr )
{
	while( mr->subsuming_match != NULL )
		mr = mr->subsuming_match;
}


void checkLink( MatchLink& mlink )
{
	while( mlink.subset->subsuming_match != NULL )
	{
		vector< size_t > new_map( mlink.sub_to_super_map.size() );
		for( size_t i = 0; i < mlink.sub_to_super_map.size(); ++i )
			new_map[i] = mlink.sub_to_super_map[ mlink.subset->subsumption_component_map[i] ];
		swap( new_map, mlink.sub_to_super_map );
		mlink.subset = mlink.subset->subsuming_match;
	}
}

void checkLinkAndComponent( MatchRecord*& mr, size_t& component )
{
	while( mr->subsuming_match != NULL )
	{
		component = mr->subsumption_component_map[component];
		mr = mr->subsuming_match;
	}
}

/** returns one of the superset links f
rom a match.  direction is 1 for left, -1 for right */
MatchLink& getSuperset( MatchRecord* mr, int direction )
{
	if( direction == 1 )
		return mr->left_superset;
	return mr->right_superset;
}

/** returns the subset links for a given direction.  direction is 1 for left, -1 for right */
vector<MatchLink>& getSubsets( MatchRecord* mr, int direction )
{
	if( direction == 1 )
		return mr->left_subset_links;
	return mr->right_subset_links;
}

/** returns the extra subsets for a given direction.  direction is 1 for left, -1 for right */
vector<MatchLink>& getExtraSubsets( MatchRecord* mr, int direction )
{
	if( direction == 1 )
		return mr->extra_left_subsets;
	return mr->extra_right_subsets;
}
//inverse of unlinkSuperset
//linkSuperset then unlinkSuperset should exactly  offset each other
void linkSuperset( MatchRecord* mr, MatchRecord* supermatch, boost::dynamic_bitset<>& comp_list, vector< size_t >& comp_map, int direction )
{
	// update superset links
	MatchLink slink = MatchLink( supermatch, mr, comp_list, comp_map );
	if( slink.superset != NULL )
	{
		slink.subset = mr;
		int parity = mr->Orientation(0) == slink.superset->Orientation(slink.sub_to_super_map[0]) ? 1 : -1;
		getSubsets(slink.superset,-direction*parity).push_back(slink);
	}
	vector< MatchLink >& subsets = getSubsets(mr,direction);
	for( size_t subI = 0; subI < subsets.size(); ++subI )
	{
		subsets[subI].superset = mr;
		int parity = mr->Orientation(subsets[subI].sub_to_super_map[0]) == subsets[subI].subset->Orientation(0) ? 1 : -1;
		getSuperset(subsets[subI].subset, -direction*parity).superset = mr;
	}
    //punt: link extra subsets too!
    vector< MatchLink >& extrasubsets = getExtraSubsets(mr,direction);
	for( size_t subI = 0; subI < extrasubsets.size(); ++subI )
	{
		subsets[subI].superset = mr;
		int parity = mr->Orientation(extrasubsets[subI].sub_to_super_map[0]) == extrasubsets[subI].subset->Orientation(0) ? 1 : -1;
		getSuperset(extrasubsets[subI].subset, -direction*parity).superset = mr;
	}
    
    
	
}
void unlinkSuperset( MatchRecord* mr, int direction )
{
	MatchLink& superlink = getSuperset( mr, direction );
	MatchRecord* super = superlink.superset;
	if( super != NULL )
	{
		int parity = mr->Orientation(0) == super->Orientation(superlink.sub_to_super_map[0]) ? 1 : -1;
		vector< MatchLink >& subs = getSubsets( super, -direction*parity );
		for( size_t subI = 0; subI < subs.size(); ++subI )
		{
			if( subs[subI].subset == mr )
			{
				subs.erase( subs.begin() + subI, subs.begin() + subI + 1 );
				subI--;
			}
		}
        //tjt: unlink extrasubsets!
        vector< MatchLink >& extrasubs = getExtraSubsets( super, -direction*parity );
		for( size_t subI = 0; subI < extrasubs.size(); ++subI )
		{
			if( extrasubs[subI].subset == mr )
			{
				extrasubs.erase( extrasubs.begin() + subI, extrasubs.begin() + subI + 1 );
				subI--;
			}
		}
        
		superlink.clear();
	}
}

void unlinkSupersets( MatchRecord* mr )
{
	unlinkSuperset( mr, 1 );
	unlinkSuperset( mr, -1 );
}

template< class MatchRecordPtrType >
void validate( vector< MatchRecordPtrType >& records )
{
	// make sure all matches have non-zero components
	for( size_t recI = 0; recI < records.size(); ++recI )
	{
		size_t seqI = 0;
		for( ; seqI < records[recI]->SeqCount(); ++seqI )
			if( records[recI]->LeftEnd(seqI) == NO_MATCH )
				break;
		if( seqI < records[recI]->SeqCount() )
		{
			cerr << "missing component\n";
			genome::breakHere();
		}
	}

	// make sure all links are consistent
	for( size_t recI = 0; recI < records.size(); ++recI )
	{
		MatchRecord* mr = records[recI];
		for( int direction = 1; direction >-2; direction -= 2 )
		{
			for( size_t subI = 0; subI < getSubsets(mr, direction).size(); subI++ )
			{
				// follow any stale links
				MatchRecord* sub = getSubsets(mr, direction)[subI].subset;
				size_t sub_mult = sub->Multiplicity();
				while( sub->subsuming_match != NULL )
					sub = sub->subsuming_match;
				size_t parity_seq = getSubsets(mr, direction)[subI].sub_to_super_map[0];
				int parity = mr->Orientation(parity_seq) == sub->Orientation(0) ? 1 : -1;
				// make sure that each of the subsets in these points back to this superset in its own link
				if( getSuperset(sub, -direction*parity).superset != mr )
				{
					cerr << "ohno\n";
					genome::breakHere();
				}
				if( sub_mult != sub->Multiplicity() )
				{
					cerr << "unequal mult\n";
					genome::breakHere();
				}
				if( getSubsets(mr,direction)[subI].super_component_list.count() != getSubsets(mr,direction)[subI].sub_to_super_map.size())
				{
					cerr << "broke\n";
					genome::breakHere();
				}
			}

			// make sure the supersets have this subset
			if( getSuperset(mr,direction).superset != NULL )
			{
				MatchRecord* sup = getSuperset(mr,direction).superset;
				int parity = mr->Orientation(0) == sup->Orientation(getSuperset(mr,direction).sub_to_super_map[0]) ? 1 : -1;
				size_t subI = 0;
				for( ; subI < getSubsets(sup,-direction*parity).size(); subI++ )
				{
					if( getSubsets(sup,-direction*parity)[subI].subset == mr )
						break;
				}
				if( subI == getSubsets(sup,-direction*parity).size() )
				{
					cerr << "oh crap!\n";
					genome::breakHere();
				}
				if( getSuperset(mr,direction).super_component_list.count() != getSuperset(mr,direction).sub_to_super_map.size())
				{
					cerr << "broke 3\n";
					genome::breakHere();
				}
			}
		}
	}
}

void createNeighborhoodGroupList( vector< NeighborhoodGroup >& group_list, vector< vector< size_t > >& group_members, vector< NeighborhoodListEntry >& neighborhood_list )
{
	group_list.resize( group_members.size() );
	for( size_t gI = 0; gI < group_members.size(); gI++ )
	{
		// is this subset completely contained--is it subsumed?
		MatchRecord* M_j = neighborhood_list[group_members[gI][0]].match;

		vector< size_t > component_map(M_j->Multiplicity(), (std::numeric_limits<size_t>::max)());
		vector< size_t > distances(M_j->Multiplicity(), (std::numeric_limits<size_t>::max)());
		for( vector< size_t >::iterator rec_iter = group_members[gI].begin(); rec_iter != group_members[gI].end(); ++rec_iter )
		{
			component_map[neighborhood_list[*rec_iter].Mj_component] = neighborhood_list[*rec_iter].Mi_component;
			distances[neighborhood_list[*rec_iter].Mj_component] = neighborhood_list[*rec_iter].distance;
		}
		group_list[gI].get<0>() = M_j;
		swap( group_list[gI].get<1>(), component_map );
		swap( group_list[gI].get<2>(), distances );
	}

	static NeighborhoodGroupCompare src;
	std::sort( group_list.begin(), group_list.end(), src );
}

/**
 * Assigns the superset link from M_j to M_i.  This function should be called when
 * M_j has been chained as part of M_i and M_j has an outgoing superset link.
 */
void inheritSuperset( MatchRecord* M_i, MatchRecord* M_j, int direction, int parity )
{
	// remap superset components
	vector< size_t > comp_map( M_i->Multiplicity() );
	for( size_t ci = 0; ci < comp_map.size(); ci++ )
		comp_map[ci] = getSuperset( M_j, direction*parity ).sub_to_super_map[ M_j->subsumption_component_map[ci] ];
	// rebuild the superset component list
	boost::dynamic_bitset<> comp_list(getSuperset( M_j, direction*parity ).superset->Multiplicity(), false);
	for( size_t compI = 0; compI < comp_map.size(); ++compI )
		comp_list.set(comp_map[compI]);
	MatchLink& slink = getSuperset(M_i, direction);
	slink = MatchLink( getSuperset( M_j, direction*parity ).superset, M_i, comp_list, comp_map );
	unlinkSuperset(M_j,direction*parity);
	int slink_parity = M_i->Orientation(0) == slink.superset->Orientation(slink.sub_to_super_map[0]) ? 1 : -1;
	getSubsets(slink.superset,-direction*slink_parity).push_back(slink);

}

/**
 * returns either the left or right list, depending on the current direction of extension
 */
vector< NeighborhoodGroup >& selectList( vector< NeighborhoodGroup >& left_list, vector< NeighborhoodGroup >& right_list, int direction )
{
	return direction == 1 ? left_list : right_list;
}



/**
 * Performs a superset link extension on M_i
 */
void supersetLinkExtension( GappedMatchRecord*& M_i, int direction, int& last_linked, 
						   vector< NeighborhoodGroup >& left_deferred_subsets, 
						   vector< NeighborhoodGroup >& right_deferred_subsets, bool chain )
{
	// update the left end and look for another superset to chain with
	// then extend all the way to that match
	MatchRecord* M_j = getSuperset(M_i, direction).superset;
	MatchLink ij_link = getSuperset(M_i, direction);	// make a copy for safekeeping
	int ij_parity = M_i->Orientation(0) == M_j->Orientation(ij_link.sub_to_super_map[0]) ? 1 : -1;

	//
	// Link extension part 1: 
	// extend M_i to include M_j, add M_j to the chained matches

       
    bool changed = extendRange( M_i, M_j, ij_link.sub_to_super_map );
    M_i->chained_matches.push_back(M_j);
    M_i->chained_component_maps.push_back(ij_link.sub_to_super_map);



	// Link extension part 2:
	// figure out whether any subsets between M_j and M_i got subsumed
	for( size_t subtypeI = 0; subtypeI < 2; subtypeI++ )
	{
		vector< MatchLink >* mjsubs;
		if( subtypeI == 0 )
			mjsubs = &getSubsets(M_j, -direction*ij_parity);
		else
			mjsubs = &getExtraSubsets(M_j, -direction*ij_parity);
		vector< MatchLink >& mj_otherside_subsets = *mjsubs;

		for( size_t leftI = 0; leftI < mj_otherside_subsets.size(); ++leftI )
		{
			if( subtypeI == 0 )
				checkLink( mj_otherside_subsets[leftI] );
			MatchLink& jk_link = mj_otherside_subsets[leftI];
			boost::dynamic_bitset<> intersect = ij_link.super_component_list & jk_link.super_component_list;
			MatchRecord* M_k = jk_link.subset;
			if( M_k == M_i )
				continue;	// been there, chained that.
			size_t inter_size = intersect.count();
			if( inter_size < 2 )
				continue;	// no match
			if( inter_size >= M_i->Multiplicity() || M_k->Multiplicity() != inter_size )
				continue;

			// has this guy already been subsumed?  if so then just skip him
			if( M_k->subsuming_match != NULL )
			{
				if( subtypeI != 1 )
					breakHere(); // this should only happen with extra subsets
				mj_otherside_subsets.erase(mj_otherside_subsets.begin()+leftI, mj_otherside_subsets.begin()+leftI+1 );
				leftI--;
				continue;
			}

			// M_k is a subset relative to M_i
			int jk_parity = M_k->Orientation(0) == M_j->Orientation(jk_link.sub_to_super_map[0]) ? 1 : -1;
			int ik_parity = ij_parity * jk_parity;

			vector< size_t > component_map( M_k->Multiplicity() );
			remapComponents(jk_link.sub_to_super_map, M_j->Multiplicity(), ij_link.sub_to_super_map, component_map );

			NeighborhoodGroup sr = boost::make_tuple( M_k, component_map, vector<size_t>( M_k->Multiplicity(), 0 ) );
			// defer it until we're done extending
			selectList( left_deferred_subsets, right_deferred_subsets, -direction ).push_back( sr );
		}
	}

	//
	// Link extension part 3:
	// classify outgoing links that share components with M_i
	unlinkSuperset(M_i,direction);
	vector< size_t > supersets;
	vector< size_t > chainable;
	vector< size_t > subsets;
	vector< size_t > novel_subsets;
	vector< MatchLink >& mj_subsets = getSubsets(M_j, direction*ij_parity);
	for( size_t leftI = 0; leftI < mj_subsets.size(); ++leftI )
	{
		checkLink( mj_subsets[leftI] );
		boost::dynamic_bitset<> intersect = ij_link.super_component_list & mj_subsets[leftI].super_component_list;
		MatchRecord* M_k = mj_subsets[leftI].subset;
		if( M_k == M_i )
			continue;	// been there, chained that.
		size_t inter_size = intersect.count();
		if( inter_size < 2 )
			continue;	// no match
			// M_k is a superset relative to M_i
		if( inter_size == M_i->Multiplicity() && M_k->Multiplicity() > inter_size )
			supersets.push_back(leftI);
		else if( inter_size == M_i->Multiplicity() && M_k->Multiplicity() == inter_size )
			chainable.push_back(leftI);
		else if( inter_size < M_i->Multiplicity() && M_k->Multiplicity() == inter_size )
			subsets.push_back(leftI);
		else
			novel_subsets.push_back(leftI);
	}


	if( supersets.size() > 0 && 1)
	{
//#4018
		cerr << "something is wrong, we should never have supersets during link extension!\n";
		genome::breakHere();
	}

	if (chain)
	{
		for( size_t cI = 0; cI < chainable.size(); ++cI )
		{
			if( chainable.size() > 1 && 1)
			{
				cerr << "bad news bruthah\n";
				genome::breakHere();
			}
			// chain with this guy
			MatchLink& jk_link = mj_subsets[chainable[cI]];
			MatchRecord* M_k = jk_link.subset;
			if( M_k->extended )
			{
				cerr << "extensor crap\n";
				breakHere();
			}
			if( M_k == M_i )
			{
				cerr << "crap\n";
				breakHere();
			}

			// update boundary coordinates
			vector< size_t > component_map( M_i->Multiplicity() );
			remapComponents(ij_link.sub_to_super_map, M_j->Multiplicity(), jk_link.sub_to_super_map, component_map );
			bool changed = extendRange( M_i, M_k, component_map );
			if( changed )
				last_linked = 2;

			// unlink from superset
			int jk_parity = M_k->Orientation(0) == M_j->Orientation(jk_link.sub_to_super_map[0]) ? 1 : -1;
			unlinkSuperset(M_k,-direction*ij_parity*jk_parity);
			// set subsuming match ptrs
			M_k->subsuming_match = M_i;
			M_k->subsumption_component_map = component_map;
			M_i->chained_matches.push_back( M_k );
			M_i->chained_component_maps.push_back( component_map );

			// compensate for the deletion in subsets
			for( size_t subI = 0; subI < chainable.size(); subI++ )
				if( chainable[subI] > chainable[cI] )
					chainable[subI]--;
			for( size_t subI = 0; subI < subsets.size(); subI++ )
				if( subsets[subI] > chainable[cI] )
					subsets[subI]--;

			// inherit M_k's outward superset and stop chaining here
			if( getSuperset( M_k, direction*ij_parity*jk_parity ).superset != NULL )
			{
				inheritSuperset( M_i, M_k, direction, ij_parity*jk_parity );
				last_linked = 2;
				break;
			}
		}
	}
	// process subsets
	for( size_t sI = 0; sI < subsets.size(); ++sI )
	{
		// change M_k to point at M_i
		MatchLink& jk_link = mj_subsets[subsets[sI]];
		MatchRecord* M_k = jk_link.subset;
		int jk_parity = M_k->Orientation(0) == M_j->Orientation(jk_link.sub_to_super_map[0]) ? 1 : -1;
		int ik_parity = ij_parity * jk_parity;

		vector< size_t > component_map( M_k->Multiplicity() );
		remapComponents(jk_link.sub_to_super_map, M_j->Multiplicity(), ij_link.sub_to_super_map, component_map );
		// rebuild the superset component list
		boost::dynamic_bitset<> comp_list(M_i->Multiplicity(), false);
		for( size_t compI = 0; compI < component_map.size(); ++compI )
			if(component_map[compI] != (std::numeric_limits<size_t>::max)())
				comp_list.set(component_map[compI]);
		unlinkSuperset(M_k,-1*direction*ik_parity);

		// add to the deferred subsets list
		NeighborhoodGroup sr = boost::make_tuple( M_k, component_map, vector<size_t>( M_k->Multiplicity(), 0 ) );
		vector< NeighborhoodGroup >& subset_list = selectList( left_deferred_subsets, right_deferred_subsets, direction );
		subset_list.push_back( sr );

		// compensate for the deletion in subsets
		for( size_t subI = 0; subI < subsets.size(); subI++ )
			if( subsets[subI] > subsets[sI] )
				subsets[subI]--;
	}
}

/**
 * Temporary buffers that get used every time a neighborhood list lookup is performed.
 * Storing the buffers persistently prevents repeated memory allocations
 */
class NllBuffers
{
public:
	std::vector< std::vector< size_t > > superset_groups;
	std::vector< std::vector< size_t > > chainable_groups;
	std::vector< std::vector< size_t > > subset_groups;
	std::vector< std::vector< size_t > > novel_subset_groups;
	vector< NeighborhoodListEntry > neighborhood_list;
	vector< pair< size_t, size_t > > j_comp_sort_list;
	vector<size_t> group_entries;

	NllBuffers()
	{
		superset_groups.reserve(100);
		chainable_groups.reserve(100);
		subset_groups.reserve(100);
		novel_subset_groups.reserve(100);
		neighborhood_list.reserve(10000);
		j_comp_sort_list.reserve(1000);
		group_entries.reserve(1000);
	};
	void clear()
	{
		superset_groups.resize(0);
		chainable_groups.resize(0);
		subset_groups.resize(0);
		novel_subset_groups.resize(0);
		neighborhood_list.resize(0);
		j_comp_sort_list.resize(0);
		group_entries.resize(0);
	};
};

NllBuffers nllbufs;

/**
 * Performs a neighborhood list lookup to find other matches nearby the match of interest
 * @param	M_i		The primary match which is under extension
 * @param	match_pos_lookup_table
 * @param	M_e		(Optionally NULL) A gapped extension which will be added to M_i after its neighborhood has been searched
 */
void neighborhoodListLookup( GappedMatchRecord* M_i, 
						   MatchPositionLookupTable& match_pos_lookup_table,
						   vector< NeighborhoodGroup >& superset_list, 
						   vector< NeighborhoodGroup >& chainable_list,
						   vector< NeighborhoodGroup >& subset_list, 
						   vector< NeighborhoodGroup >& novel_subset_list,
						   int direction,
						   uint seed_size,
						   uint w,
						   bitset_t& left_lookups,
						   bitset_t& right_lookups,
						   GappedMatchRecord* M_e
						   )
{
	// make sure storage is empty
	nllbufs.clear();
	//
	// construct a neighborhood list and process the neighborhood groups
	//
	vector< NeighborhoodListEntry >& neighborhood_list = nllbufs.neighborhood_list;
	for( size_t x = 0; x < M_i->Multiplicity(); ++x )
	{
		int o_x = M_i->Orientation(x) == AbstractMatch::forward ? 1 : -1;
		int parity = o_x * direction;
		int64 match_end = parity == 1 ? M_i->LeftEnd(x) : M_i->RightEnd(x) - seed_size + 1;

		if( match_end > 0 )
			if( (direction == 1 && left_lookups.test(match_end)) ||
				(direction == -1 && right_lookups.test(match_end)) )
			{
				if(print_warnings)
					cerr << "looking twice in the same place\n";
//							genome::breakHere();
			}else{
				if( direction == 1 )
					left_lookups.set(match_end);
				if( direction == -1 )
					right_lookups.set(match_end);
			}

		int d = 1;
		int w_end = parity == 1 ? w : w + seed_size;
		// are we cleaning up a gapped extension?  if so, adjust d and w_end so
		// we don't search anything twice and also cover all of the extension area
		if(M_e != NULL)
		{
			int64 me_match_end = parity == 1 ? M_e->LeftEnd(x) : M_e->RightEnd(x)-(M_e->Length(x)-1);
			d = w+1;	// need to start at the begining of the window to properly 
                    // classify all matches subsumed by extension and all novel 
                    // matches which may have been discovered
   			w_end = w + me_match_end - match_end;	// search anything new included in M_e
		}
		for( ; d <= w_end; ++d )
		{
			if( match_end <= parity * d )
				continue;	// we're too close to the beginning
			size_t mplt_index = match_end - parity * d;
			if( mplt_index >= match_pos_lookup_table.size() )
				continue;	// we're too close to the end!

			MatchRecord* M_j = match_pos_lookup_table[ mplt_index ].first;
			size_t y = match_pos_lookup_table[ mplt_index ].second;
			if( M_j == NULL )
				continue;	// no match at this position

			NeighborhoodListEntry nle;
			nle.match = M_j;
			nle.Mi_component = x;
			nle.Mj_component = y;
			// update the link if this one was subsumed
			checkLinkAndComponent( M_j, y );
			int o_y = ((AbstractMatch*)M_j)->Orientation(y) == AbstractMatch::forward ? 1 : -1;
			nle.relative_orientation = o_x * o_y == 1 ? true : false;
			nle.distance = d;
			neighborhood_list.push_back( nle );
			
			if( M_j == M_i )
			{
				M_i->tandem = true;
				break;	// not so fast there cowboy!  can't chain beyond ourself!
			}
		}
	}

	//
	// now classify each group of the neighborhood list and act appropriately
	// group types are superset, chainable, subset, novel subset
	//
	NeighborhoodListComparator nlc;
	std::sort( neighborhood_list.begin(), neighborhood_list.end(), nlc );

    //std::reverse(neighborhood_list.begin(), neighborhood_list.end());

	std::vector< std::vector< size_t > >& superset_groups = nllbufs.superset_groups;
	std::vector< std::vector< size_t > >& chainable_groups = nllbufs.chainable_groups;
	std::vector< std::vector< size_t > >& subset_groups = nllbufs.subset_groups;
	std::vector< std::vector< size_t > >& novel_subset_groups = nllbufs.novel_subset_groups;

	size_t group_end = 0;
	for( size_t prev = 0; prev < neighborhood_list.size(); prev = group_end )
	{
		group_end = prev + 1;
		while( group_end < neighborhood_list.size() && 
			neighborhood_list[prev].match == neighborhood_list[group_end].match && 
			neighborhood_list[prev].relative_orientation == neighborhood_list[group_end].relative_orientation )
		{
			++group_end;
		}
		// the group is everything in the range of prev to end-1
		if( prev + 1 == group_end )
			continue;	// can't do anything with groups of size 1 -- there's no match

		// do something about ties here...???
		// this code selects the *furthest* away match (e.g. that with the largest d)
		// because that's what got sorted in last in the comparator
		// it eliminates both duplicate M_i and duplicate M_j components...
		// FIXME:  is this true?  is it safe?
		vector< pair< size_t, size_t > >& j_comp_sort_list = nllbufs.j_comp_sort_list;
		j_comp_sort_list.resize(0);
		for( size_t i = prev + 1; i < group_end; ++i )
		{
            //selects the *furthest* away match 
			if( neighborhood_list[i-1].Mi_component == neighborhood_list[i].Mi_component )
				continue;
			j_comp_sort_list.push_back(make_pair(neighborhood_list[i-1].Mj_component, i-1));
		}
		j_comp_sort_list.push_back(make_pair(neighborhood_list[group_end-1].Mj_component, group_end-1));
		std::sort(j_comp_sort_list.begin(), j_comp_sort_list.end());
		vector<size_t>& group_entries = nllbufs.group_entries;
		group_entries.resize(0);
		for( size_t i = 1; i < j_comp_sort_list.size(); ++i )
		{
            //selects the *furthest* away match
			if( j_comp_sort_list[i-1].first == j_comp_sort_list[i].first )
				continue;
			group_entries.push_back(j_comp_sort_list[i-1].second);
		}
		group_entries.push_back(j_comp_sort_list.back().second);

		// update the links in case something is subsumed
		for( size_t gI = 0; gI < group_entries.size(); ++gI )
			checkLinkAndComponent( neighborhood_list[group_entries[gI]].match, neighborhood_list[group_entries[gI]].Mj_component );

		// finally, classify the match as one of superset, subset, 
		// chainable, novel subset
		MatchRecord* M_j = neighborhood_list[prev].match;

		if( group_entries.size() == M_i->Multiplicity() && 
			M_j->Multiplicity() > M_i->Multiplicity() )
		{
			// superset
			superset_groups.push_back( group_entries );
		}else
		if( group_entries.size() == M_i->Multiplicity() && 
			M_j->Multiplicity() == M_i->Multiplicity() )
		{
            
            // chainable
			chainable_groups.push_back( group_entries );
		}else
		if( group_entries.size() < M_i->Multiplicity() && 
			group_entries.size() == M_j->Multiplicity() )
		{
			// subset
			subset_groups.push_back( group_entries );
		}else
		{
            // novel subset
			novel_subset_groups.push_back( group_entries );
		}

	}	// end loop that splits the neighborhood into groups

	createNeighborhoodGroupList( superset_list, superset_groups, neighborhood_list );
	createNeighborhoodGroupList( chainable_list, chainable_groups, neighborhood_list );
	createNeighborhoodGroupList( subset_list, subset_groups, neighborhood_list );
	createNeighborhoodGroupList( novel_subset_list, novel_subset_groups, neighborhood_list );
}

/**
 * Chains matches onto M_i or subsumes them as appropriate
 */
void processChainableMatches( GappedMatchRecord*& M_i, vector< NeighborhoodGroup >& chainable_list,
				  int direction, int& last_linked, bool find_novel_subsets, bool chain )
{
	// link the closest possible chainable first.
	for( size_t gI = 0; gI < chainable_list.size(); gI++ )
	{
		MatchRecord* M_j = chainable_list[gI].get<0>();

		vector< size_t >& component_map = chainable_list[gI].get<1>();

		if( M_j == M_i )
		{
			// this is an inverted overlapping repeat, skip it.
			continue;
		}
        if( M_j->extended )
		{
            if ( !find_novel_subsets && (M_i->is_novel_subset ))
            {
                //novel subsets have been disabled!! this is why it wasn't swallowed up!
                continue;
            }
            else
            {
                // oh no!  M_i should have been swallowed up already!
                //tjt: claro, work has been wasted, but bypassing the breakHere() will allow
                //the assumed-to-be subsumed M_i to be detected and updated accordingly
                //but the question remains, why wasn't M_i previously subsumed?
                //1)   what if before gapped extension M_j was not in M_i's neighborhood?
                //     but after gapped extension, M_i is found in M_j's neighborhood and classified as chainable?
			    //cerr << "extensor crap 2\n";
			    //breakHere();
            }
		}

		bool subsumed;
		bool partial;
		classifySubset( M_i, chainable_list[gI], subsumed, partial );
		
		vector< size_t >& yx_map = chainable_list[gI].get<1>();
		vector< size_t > xy_map(yx_map.size());
		for( size_t i = 0; i < yx_map.size(); ++i )
			xy_map[ yx_map[i] ] = i;
//		for( vector< size_t >::iterator rec_iter = chainable_groups[gI].begin(); rec_iter !=  chainable_groups[gI].end(); ++rec_iter )
//			xy_map[ neighborhood_list[*rec_iter].Mi_component ] = neighborhood_list[*rec_iter].Mj_component;

		// if M_j isn't extending the boundaries of every component of M_i then
		// it may be inconsistent with already chained matches.  just subsume it without
		// chaining in that case.
		if( !subsumed && !partial && chain)
		{
			bool ok = testChainableMatch(M_i, M_j, xy_map);
            if (ok)
            {
                M_i->chained_matches.push_back( M_j );
			    M_i->chained_component_maps.push_back( component_map );
			    bool changed = extendRange(M_i, M_j, xy_map);
			    if( changed )
                {
			        // update the left-end and right-end coords
				    last_linked = 2;
                }
            }
            else
                break;
		}
        M_j->subsuming_match = M_i;
		M_j->subsumption_component_map = component_map;
		int parity = M_i->Orientation(0) == M_j->Orientation(xy_map[0]) ? 1 : -1;
		if( getSuperset( M_j, -direction*parity ).superset != NULL )
			unlinkSuperset(M_j,-direction*parity);	// won't be needing this anymore...

		// if M_j has a superset then inherit it and stop chaining here
		if( getSuperset( M_j, direction*parity ).superset != NULL )
		{
			inheritSuperset( M_i, M_j, direction, parity );
			last_linked = 2;	// we may do a link extension!
			break;
		}
	}
}
//processes supersets
void processSupersetMatches( GappedMatchRecord*& M_i, vector< NeighborhoodGroup >& superset_list,
				  int direction, int& last_linked, bool gapped_extension = false )
{
	
	// link the closest possible superset first.
	for( size_t gI = 0; gI < superset_list.size(); gI++ )
	{
		MatchRecord* M_j = superset_list[gI].get<0>();

		vector< size_t >& component_map = M_i->chained_component_maps.at(0);
		boost::dynamic_bitset<> comp_list(M_j->Multiplicity(), false);
		for( size_t compI = 0; compI < M_i->Multiplicity(); ++compI )
			comp_list.set(component_map[compI]);
		if( M_j == M_i )
		{
			// this is an inverted overlapping repeat, skip it.
			continue;
		}
        //tjt: shouldn't the superset always be extended when we reach this point during gapped extension?
		if( M_j->extended && !gapped_extension )
		{
           	// oh no!  M_i should have been swallowed up already!
			cerr << "extensor crap 2\n";
			breakHere();
		}

		bool subsumed;
		bool partial;
		//update classifysubset to ClassifySuperset
		classifySuperset( M_i, superset_list[gI], subsumed, partial );

		if( subsumed && !partial )
		{
			// update the left-end and right-end coords
			bool changed = reduceRange(M_i, M_j, component_map);
		}
		if( partial )
			//some of the components of the superset matches are subsumed
			//punt: what should I do differently here?

		linkSuperset( M_i, M_j, comp_list, component_map,  direction);
		last_linked = 1;// stores the group type that was chained.  1 == superset, 2 == chainable, 0 == none
					
	}
}


/**
 * Performs a gapped extension on a match.  The region either left or right of the match is processed by
 * progressive alignment.
 * @param	M_i		The match to extend
 * @param	seq_table	gnSequences which correspond to each match component
 * @param	params		The Homology HMM parameters to use
 * @param	w		The max gap for chaining.  Used to compute extension lengths.
 * @param	direction	The direction of extension
 * @param	M_e		(output) A MatchRecord containing just the extension, or NULL if extension failed
 * @return	FAILED, OK, or FIXME
 */
int ExtendMatch(GappedMatchRecord*& M_i, vector< gnSequence* >& seq_table, Params& hmm_params, unsigned w, int direction, vector<GappedMatchRecord*>& novel_matches, int gap_open, int gap_extend, int extension_window)
{
	ccount +=1;
	static bool debug_extension = false;
//  punt on this for now..
	bool novel_hss_regions_support = false;
	bool danger_zone_active = true;
	int multi = M_i->Multiplicity();
	double e = 2.71828182845904523536;
//	I think this works a little better...
	
	int extend_length = 80*pow(e,-0.01*multi);
	//use user specified window if requested
	if (extension_window >= 0 )
		extend_length = extension_window;
	vector<int> left_extend_vector(multi,0);
	vector<int> right_extend_vector(multi,0);
	int left_extend_length = extend_length;	
	int right_extend_length = extend_length;

	if ( M_i->tandem )
	{		
        if ( debug_extension)
		    cerr << "Sorry, no extension for tandem repeats.." << endl << endl;	
		return FIXME;
	}

//  careful, if M_i->LeftEnd(j) < extend_length, ToString() will be disappointed...
	for( gnSeqI j = 0; j < multi; j++)
	{
//      now put check for curpos+extend_length<startpos of next match component..
		if( M_i->Orientation(j) == AbstractMatch::reverse )
		{
//          if leftend <= 0 set right extension to 0
			if( M_i->LeftEnd(j) <= 0 || M_i->LeftEnd(j) > 4000000000u )
                right_extend_vector[j] = 0;
//          if extend_length goes too far, set to maximum possible
			else if ( M_i->LeftEnd(j) <= extend_length )
				right_extend_vector[j] = M_i->LeftEnd(j)-1;
//          if we run into another match, don't extend into it
			else if ( j > 0 && M_i->LeftEnd(j) - extend_length <= M_i->RightEnd(j-1) )
            {
				int parity = M_i->Orientation(j) ==  M_i->Orientation(j-1) ? 1 : 1;
                right_extend_vector[j] = parity*(M_i->LeftEnd(j)-M_i->RightEnd(j-1)-1);
            }
//          else everything ok to set to preset extend_length
			else
				right_extend_vector[j] = extend_length-1;

			if(M_i->RightEnd(j) <= 0 || M_i->RightEnd(j) > 4000000000u)
				left_extend_vector.push_back(0);
			else if ( M_i->RightEnd(j) + extend_length > seq_table[0]->length() )
				left_extend_vector[j] = seq_table[0]->length()-M_i->RightEnd(j);
			else if ( j+1 < multi && M_i->RightEnd(j) + extend_length >= M_i->LeftEnd(j+1) )
            {
                int parity = M_i->Orientation(j) ==  M_i->Orientation(j+1) ? 1 : 1;
				left_extend_vector[j] =parity*( M_i->LeftEnd(j+1)-M_i->RightEnd(j)-1);
            }
            else
				left_extend_vector[j] = extend_length-1;
		}
        else
        {
            if( M_i->LeftEnd(j) <= 0 || M_i->LeftEnd(j) > 4000000000u )
			    left_extend_vector[j] = 0;
		    else if ( M_i->LeftEnd(j) <= extend_length )
			    left_extend_vector[j] = M_i->LeftEnd(j)-1;
		    else if ( j > 0 && M_i->LeftEnd(j) - extend_length <= M_i->RightEnd(j-1) )
            {
                int parity = M_i->Orientation(j) ==  M_i->Orientation(j-1) ? 1 : 1;
			    left_extend_vector[j] = parity*(M_i->LeftEnd(j)-M_i->RightEnd(j-1)-1);
            }
            else
			    left_extend_vector[j] = extend_length;

		    if(M_i->RightEnd(j) <= 0 || M_i->RightEnd(j) > 4000000000u)
			    right_extend_vector[j] = 0;
		    else if ( M_i->RightEnd(j) + extend_length > seq_table[0]->length() )
			    right_extend_vector[j] = seq_table[0]->length()-M_i->RightEnd(j)-1;
		    else if ( j+1 < multi && M_i->RightEnd(j) + extend_length >= M_i->LeftEnd(j+1) )
            {
                int parity = M_i->Orientation(j) ==  M_i->Orientation(j+1) ? 1 : 1;
			    right_extend_vector[j] = parity*(M_i->LeftEnd(j+1)-M_i->RightEnd(j)-1);
            }
		    else
			    right_extend_vector[j] = extend_length;	
        }
	}
    
	left_extend_length = *(std::min_element( left_extend_vector.begin(), left_extend_vector.end() ));
	right_extend_length = *(std::min_element( right_extend_vector.begin(), right_extend_vector.end() ));
    left_extend_length = left_extend_length < 0 ? 0 : left_extend_length;
    right_extend_length = right_extend_length < 0 ? 0 : right_extend_length;
    extend_length = direction < 0 ? right_extend_length : left_extend_length;
	const gnFilter* rc_filter = gnFilter::DNAComplementFilter();
	std::vector<std::string> leftExtension(multi);
	GappedAlignment leftside(multi,left_extend_length);
	std::vector<std::string> rightExtension(multi);
	GappedAlignment rightside(multi,right_extend_length);
	vector< string > leftExtension_aln;
	vector< string > rightExtension_aln;
	if ( left_extend_length > 0 && direction == 1  )
	{
//      extract sequence data
		for( gnSeqI j = 0; j < multi; j++)
		{			
			if( M_i->Orientation(j) == AbstractMatch::reverse )
			{			
				seq_table[0]->ToString( leftExtension[j], left_extend_length, M_i->RightEnd(j)+1 );
				leftside.SetLeftEnd(j,M_i->RightEnd(j)+1);
				rc_filter->ReverseFilter(leftExtension[j]);
			}else{
				seq_table[0]->ToString( leftExtension[j], left_extend_length, M_i->LeftEnd(j) - left_extend_length );
				leftside.SetLeftEnd(j,M_i->LeftEnd(j) - left_extend_length);
			}
			leftside.SetOrientation(j,M_i->Orientation(j));
			leftside.SetLength(left_extend_length,j);
		}
		bool align_success = false;
		//mems::MuscleInterface::getMuscleInterface().SetMuscleArguments("-stable -quiet -seqtype DNA");
		align_success = mems::MuscleInterface::getMuscleInterface().CallMuscleFast( leftExtension_aln, leftExtension, gap_open, gap_extend );
		if ( align_success ){		
			leftside.SetAlignment(leftExtension_aln);
			leftside.SetAlignmentLength(leftExtension_aln.at(0).size());
		}else{
			cerr << "Extension failed: Muscle error" << endl;
			return FAILED;
		}
	}
	else if ( right_extend_length > 0 && direction == -1 )
	{
		for( gnSeqI j = 0; j < multi; j++)
		{
			if( M_i->Orientation(j) == AbstractMatch::reverse )
			{			
				rightside.SetLeftEnd(j,M_i->LeftEnd(j) - right_extend_length-1);
				seq_table[0]->ToString( rightExtension[j], right_extend_length, M_i->LeftEnd(j) - right_extend_length-1);
				rc_filter->ReverseFilter(rightExtension[j]);
			}else{
				rightside.SetLeftEnd(j,M_i->RightEnd(j)+1 );
				seq_table[0]->ToString( rightExtension[j], right_extend_length, M_i->RightEnd(j)+1 );
			}
			rightside.SetOrientation(j,M_i->Orientation(j));
			rightside.SetLength(right_extend_length,j);
		}
		bool align_success = false;		
		align_success = mems::MuscleInterface::getMuscleInterface().CallMuscleFast( rightExtension_aln, rightExtension, gap_open, gap_extend );
		if ( align_success ){
			rightside.SetAlignment(rightExtension_aln);
			rightside.SetAlignmentLength(rightExtension_aln.at(0).size());
		}else{
            cerr << "Extension failed: Muscle error" << endl;
			return FAILED;
		}
	}else{
		//what are you even doing here?!?
		if(debug_extension)
        {
            cerr << "Extension failed: No room to extend" << endl;
		}
		return FAILED;
	}

//  tjt: don't use original match, only regions to the left/right
//       since for now we won't modify M_i, even if the homology detection method suggests otherwise
	vector< AbstractMatch* > mlist;
	if( direction == 1 )
		mlist.push_back(leftside.Copy());
	if( direction == -1 )
		mlist.push_back(rightside.Copy());
 
//  createIntervald
	Interval iv;
	iv.SetMatches(mlist);
	CompactGappedAlignment<> tmp_cga;
	CompactGappedAlignment<>* cga = tmp_cga.Copy();
	new (cga)CompactGappedAlignment<>(iv);
	vector< CompactGappedAlignment<>* > cga_list;
	CompactGappedAlignment<>* result;
//  detectAndApplyBackbone
	backbone_list_t bb_list;

	detectAndApplyBackbone( cga, seq_table,result,bb_list, hmm_params, direction != 1, direction == 1  );
	cga->Free();

	bool boundaries_improved = false;
	if( bb_list.size() == 0 || bb_list.at(0).size() == 0)
	{
//      no backbone segment found
        if(debug_extension)
            cerr << "Extension failed: no backbones found during extension" << endl;
		result->Free();
		return FAILED;
	}
	AbstractMatch* extension_bb;
//  tjt: was > before, wasn't taking right backbone???
//  remember, direction == -1 rightward, == 1 leftward
    bool isnovel = false;
    for( size_t bbI = 0; bb_list.size() > 0 && bbI < bb_list[0].size(); bbI++ )
    {
        extension_bb = bb_list.at(0).at(bbI);
        if (extension_bb == NULL)
            continue;
        
	    cga_list.push_back( tmp_cga.Copy() );
	    result->copyRange( *(cga_list.back()), extension_bb->LeftEnd(0), extension_bb->AlignmentLength()-1 );
        int cgalen = cga_list.back()->AlignmentLength();
        int resultlen = result->AlignmentLength();

        //why set this to > 5? what about seed weight? or to > 0? seems strange to allow novel matches of length 1...
	    if( cga_list.back()->Multiplicity() > 1 && cga_list.back()->Length() > 0 && cga_list.back()->AlignmentLength() > 0 )
	    {
//          successful extension!!
//          boundaries were improved, current match extended original match
//          create a GappedMatchRecord for the gapped extension
		    vector< AbstractMatch* > matches( 1, cga_list.back());
//          GappedMatchRecord* M_e = M_i->Copy();
            UngappedMatchRecord tmp(  cga_list.back()->Multiplicity(), cga_list.back()->AlignmentLength() );
            MatchRecord* umr = tmp.Copy();
		    GappedMatchRecord* M_e = dynamic_cast<GappedMatchRecord*>(umr); 
            
		    if( M_e == NULL )
		    {
//              create a new gapped match record for M_i
			    GappedMatchRecord gmr( *(UngappedMatchRecord*)umr );
			    M_e = gmr.Copy();
			    umr->subsuming_match = M_e;
//              umr->subsuming_match = M_i;
			    M_e->chained_matches.push_back( umr );
			    vector< size_t > component_map( M_e->SeqCount() );
			    for( size_t i = 0; i < component_map.size(); ++i )
				    component_map[i] = i;
			    M_e->chained_component_maps.push_back(component_map);
			    swap(umr->subsumption_component_map, component_map);	// swap avoids reallocation
//              update superset and subset links
			    for( int dI = 1; dI > -2; dI -= 2 )
			    {
				    MatchLink& ij_link = getSuperset(M_e,dI);
				    if( ij_link.superset != NULL )
				    {
					    ij_link.subset = M_e;
					    unlinkSuperset(umr,dI);
					    int parity = M_e->Orientation(0) == ij_link.superset->Orientation(ij_link.sub_to_super_map[0]) ? 1 : -1;
					    getSubsets(ij_link.superset,-dI*parity).push_back(ij_link);
				    }
				    vector< MatchLink >& subsets = getSubsets(M_e,dI);
				    for( size_t subI = 0; subI < subsets.size(); ++subI )
				    {
					    subsets[subI].superset = M_e;
					    int parity = M_i->Orientation(subsets[subI].sub_to_super_map[0]) == subsets[subI].subset->Orientation(0) ? 1 : -1;
					    getSuperset(subsets[subI].subset, -dI*parity).superset = M_e;
				    }
				    getSubsets(umr,dI).clear();	// so that validate works...
                }
                //          tjt: clobber M_e's GappedMatchRecord data and set boundaries
                //tjt: we call the Temp version since we don't actually want to do anything with the regions between the matches
                M_e->SetMatchesTemp(matches);//,cga_list.back()->Multiplicity() );
            }
            novel_matches.push_back(M_e->Copy());  
        }
      }
    result->Free();
	if (novel_matches.size() > 0)
        return OK;
    else
        return FAILED;
}//tjt: match should be extended!

/**
 * A class to prioritize match records for extension based on their multiplicity
 */
class ProcrastinationQueue
{
public:
	template< typename MrPtrType >
	ProcrastinationQueue( vector< MrPtrType >& match_record_list ) :
	mhc()
	{
		q.resize( match_record_list.size() );
		std::copy(match_record_list.begin(), match_record_list.end(), q.begin() );
		std::make_heap( q.begin(), q.end(), mhc );
		q_end = q.size();
		q_size = q.size();
	}

	/** pops an element from the queue, maintaining heap order */
	MatchRecord* pop()
	{
		std::pop_heap( q.begin(), q.begin()+q_end, mhc );
		q_end--;
		return *(q.begin() + q_end);
	}

	/** Adds an element to the queue and restores heap order */
	void push( MatchRecord* M_n )
	{
		if( q_end < q.size() )
			q[q_end] = M_n;
		else
		{
			q.push_back(M_n);
		}
		q_size++;
		q_end++;
		std::push_heap(q.begin(), q.begin()+q_end, mhc);
	}
	/** gets the total number of elements that have been placed in the queue */
	size_t size() const{ return q_size; }
	/** returns the current number of elements in the queue */
	size_t end() const{ return q_end; }


	/** defines a multiplicity heap ordering */
	class MultiplicityHeapCompare
	{
	public:
		bool operator()( const MatchRecord* a, const MatchRecord* b )
		{
			return a->Multiplicity() < b->Multiplicity();
		}
	};

private:
	const MultiplicityHeapCompare mhc;
	std::vector< MatchRecord* > q;
	size_t q_end;
	size_t q_size;
};

/**
 * Creates novel subset matches where appropriate and adds them to the procrastination queue
 */
void processNovelSubsetMatches( GappedMatchRecord*& M_i, vector< NeighborhoodGroup >& novel_subset_list,
				bool find_novel_subsets, ProcrastinationQueue& procrastination_queue,
				vector< gnSequence* >& seq_table, int direction, uint w, int& last_linked,
				size_t& novel_subset_count )
{
	// finally process novel subset
	bool prev_linked = false;	// we only want to link the closest of a group with the same components
	int created_thisround = 0;
	static NeighborhoodGroupComponentCompare srcc;
	for( size_t gI = 0; gI < novel_subset_list.size(); gI++ )
	{
		if( !find_novel_subsets )
			continue;	// only find novel subsets if we're supposed to
		if( last_linked != 0 )
			continue;	// only generate subsets when last_linked == 0

		// be sure to handle case where:
		// --M_i-->   --M_j--   <--M_i--
		// that may cause an identical novel subset to get spawned but with
		// M_i and M_j swapped as left and right supersets

		bool same_components = false;
		if( gI > 0 )
			same_components = srcc.compare(novel_subset_list[gI], novel_subset_list[gI-1]) == 0;
		prev_linked = same_components? prev_linked : false;

		if( prev_linked )
			continue;	// don't link a subset with the same components...

		// TODO: handle the tandem repeat case
		if( M_i->tandem )
		{
			// step 1. count tandem repeat components
			// step 2. create a new GappedMatchRecord with one component per tandem component
			// add a GappedMatchRecord with the outer component boundaries
			// do ordinary extension
			// when finalize() gets called, something special needs to happen
			continue;
		}

		MatchRecord* M_j = novel_subset_list[gI].get<0>();
		// if M_j hasn't been extended then we don't do anything yet.
		// we may find this novel subset again when M_j gets extended
		if( M_j->extended == false)
			continue;

		size_t mult = 0;	// multiplicity of novel subset
		for( size_t i = 0; i < novel_subset_list[gI].get<1>().size(); ++i )
			if( novel_subset_list[gI].get<1>()[i] != (std::numeric_limits<size_t>::max)() )
				mult++;

		if( mult < 2 )
			continue;	// can't do anything if there's no match!

		UngappedMatchRecord tmper1(mult,0);
		GappedMatchRecord tmper2(tmper1);  // this is lame
		GappedMatchRecord* M_n = tmper2.Copy();

		size_t mnewi = 0;
		vector< size_t > new_to_i_map(mult);
		vector< size_t > new_to_j_map(mult);
		boost::dynamic_bitset<> ni_list(M_i->Multiplicity());
		boost::dynamic_bitset<> nj_list(M_j->Multiplicity());
		for( size_t i = 0; i < novel_subset_list[gI].get<1>().size(); ++i )
		{
			if( novel_subset_list[gI].get<1>()[i] != (std::numeric_limits<size_t>::max)() )
			{
				new_to_i_map[mnewi] = novel_subset_list[gI].get<1>()[i];
				new_to_j_map[mnewi] = i;
				ni_list.set(new_to_i_map[mnewi]);
				nj_list.set(i);
				M_n->SetStart(mnewi, M_j->Start(i));	// sets left-end and orientation
				M_n->SetLength(M_j->Length(i),mnewi);
				mnewi++;
			}
		}
		if( M_n->Orientation(0) == AbstractMatch::reverse )
			M_n->Invert();
		// before we go any further, make sure that the relevant portion of M_i is not 
		// either completely or partially subsumed by the relevant portion of M_j!
		MatchProjectionAdapter mpaa( M_i, new_to_i_map );
		vector< size_t > mpaa_to_Mn_map( new_to_i_map.size() );
		for( size_t i = 0; i < mpaa_to_Mn_map.size(); ++i )
			mpaa_to_Mn_map[i] = i;
		bool subsumed;
		bool partial;
		classifyMatch( M_n, &mpaa, mpaa_to_Mn_map, subsumed, partial );
		if( subsumed )
		{
			M_n->Free();
			continue;	// there's nothing novel about this subset...
		}
		if( partial )
		{
			// FIXME: we should really spawn a novel subset on the non-subsumed components
			M_n->Free();
			continue;
		}
		created_thisround+= M_n->Multiplicity();

		M_n->chained_matches.push_back(M_j);
		M_n->chained_component_maps.push_back(new_to_j_map);
		//tjt: need to send finalize seq_table for muscle alignment
		M_n->finalize(seq_table);	// make this one a legitimate match...
		M_n->chained_matches.clear();
		M_n->chained_component_maps.clear();

        M_n->is_novel_subset = true;//yep, this is a novel subset
		// create links from M_n to M_i and M_j
		int ni_parity = M_n->Orientation(0) == M_i->Orientation(new_to_i_map[0]) ? 1 : -1;
		int nj_parity = M_n->Orientation(0) == M_j->Orientation(new_to_j_map[0]) ? 1 : -1;
		MatchLink& ni_link = getSuperset(M_n,-direction*ni_parity);
		ni_link = MatchLink(M_i,M_n,ni_list,new_to_i_map);
		getSubsets(M_i,direction).push_back(ni_link);
        //getExtraSubsets(M_i,direction).push_back(ni_link);
		MatchLink& nj_link = getSuperset(M_n,direction*ni_parity);
		nj_link = MatchLink(M_j,M_n,nj_list,new_to_j_map);
		getSubsets(M_j,-direction*ni_parity*nj_parity).push_back(nj_link);
        //getExtraSubsets(M_j,-direction*ni_parity*nj_parity).push_back(nj_link);
		// push M_n onto the heap
		novel_subset_list.push_back(M_n);
		//procrastination_queue.push(M_n);
		novel_subset_count++;
	}


}


/**
 * Writes a set of MatchRecords in eXtended Multi-FastA format
 * @param	seedml	A matchlist containing the seq_table of interest
 * @param	extended_matches	A set of matches to write out
 * @param	xmfa_file	The filename to use for output
 */
void writeXmfa( MatchList& seedml, std::vector< GappedMatchRecord* >& extended_matches, const std::string& xmfa_file )
{
	GenericIntervalList<GappedMatchRecord> gmr_list;
	for( size_t gmrI = 0; gmrI < extended_matches.size(); ++gmrI )
		gmr_list.push_back(*extended_matches[gmrI]);

	if( xmfa_file.length() > 0  && xmfa_file != "-")
	{
		gmr_list.seq_filename.push_back( seedml.seq_filename[0] );
		gmr_list.seq_table.push_back( seedml.seq_table[0] );
		if( xmfa_file == "-" )
			gmr_list.WriteStandardAlignment(cout);
		else
		{
			ofstream xmfa_out(xmfa_file.c_str());
			gmr_list.WriteStandardAlignment(xmfa_out);
			xmfa_out.close();
		}
	}
}

/**
 * Writes a set of MatchRecords in XML format
 * @param	seedml	A matchlist containing the seq_table of interest
 * @param	extended_matches	A set of matches to write out
 * @param	xml_file	The filename to use for output
 */
void writeXML( MatchList& seedml, std::vector< GappedMatchRecord* >& extended_matches, const std::string& xml_file )
{
	 
	GenericIntervalList<GappedMatchRecord> gmr_list;
	for( size_t gmrI = 0; gmrI < extended_matches.size(); ++gmrI )
		gmr_list.push_back(*extended_matches[gmrI]);

	if( xml_file.length() > 0  && xml_file != "-")
	{
		gmr_list.seq_filename.push_back( seedml.seq_filename[0] );
		gmr_list.seq_table.push_back( seedml.seq_table[0] );
		if( xml_file == "-" )
			gmr_list.WriteXMLAlignment(cout);
		else
		{
			ofstream xml_out(xml_file.c_str());
			gmr_list.WriteXMLAlignment(xml_out);
			xml_out.close();
		}
	}
}

class ToUPPER
{
public:
	char operator()( char a ){ return toupper(a); }
};
int main( int argc, char* argv[] )
{
//	debug_interval = true;
	// Declare the supported options.
    bool debug_extension = false;

	string sequence_file = "";
	int extension_window = 0;
	int w = 0;
	int kmersize =0;
	int gap_open = 0;
	int gap_extend = 0;
	uint seed_weight = 0;
    uint min_repeat_length = 0;
    score_t min_spscore = 0;
    uint rmin = 0;
    uint rmax = 0;
	string outputfile = "";
	string output2file = "";
	string xmfa_file = "";
    string xml_file = "";
	string stat_file = "";
	string seed_file = "";
	bool only_direct = false;
	bool load_sml = false;
	bool small_repeats = false;
	bool large_repeats = false;
    bool allow_tandem = false;
    bool allow_redundant = false;
	bool find_novel_subsets = false;
    bool use_novel_matches = true; //should procrast use novel matches found during gapped extension ?
	bool solid_seed = false;
	bool extend_chains = true;
	bool chain = true;
	bool two_hits = false;
	bool unalign = true;
	float percent_id = 0.0;
    float pGoHomo = 0.004f;
    float pGoUnrelated = 0.004f;
    bool only_extended = false;

	po::variables_map vm;
	try {

        po::options_description desc("Allowed options");
        desc.add_options()
			("allow-redundant", po::value <bool>(&allow_redundant)->default_value(true), "allow redundant alignments?")
			("chain", po::value<bool>(&chain)->default_value(true), "chain seeds?")
			("extend", po::value<bool>(&extend_chains)->default_value(true), "perform gapped extension on chains?")
			("window", po::value<int>(&extension_window)->default_value(-1), "size of window to use during gapped extension")
			("gapopen",po::value <int>(&gap_open)->default_value(0), "gap open penalty")
			("gapextend",po::value <int>(&gap_extend)->default_value(0), "gap extension penalty")
			("h", po::value<float>(&pGoHomo)->default_value(0.008f), "Transition to Homologous")
			("help", "get help message")
			("highest", po::value<string>(&stat_file)->default_value("procrast.highest"), "file containing highest scoring alignment for each multiplicity ")
            ("l", po::value <unsigned>(&min_repeat_length)->default_value(1), "minimum repeat length")
			("large-repeats", po::value <bool>(&large_repeats)->default_value(false), "optimize for large repeats")
			("load-sml", po::value <bool>(&load_sml)->default_value(false), "try to load existing SML file?")
			("onlydirect",po::value<bool>(&only_direct)->default_value(false), "only process seed matches on same strand?")
			("onlyextended",po::value<bool>(&only_extended)->default_value(false), "only output extended matches?")
			("output", po::value<string>(&outputfile)->default_value(""), "procrastAligner output ")
			("percentid", po::value<float>(&percent_id)->default_value(0.0), "min repeat family % id")
			("novel-subsets", po::value<bool>(&find_novel_subsets)->default_value(false), "find novel subset matches?")
            ("novel-matches", po::value<bool>(&use_novel_matches)->default_value(true), "use novel matches found during gapped extension?")
			("rmax",  po::value<unsigned>(&rmax)->default_value(500), "maximum repeat multiplicity (max copy number)")
			("rmin" , po::value<unsigned>(&rmin)->default_value(2), "minimum repeat multiplicity (min copy number)")
			("seeds", po::value<string>(&seed_file), "seed output file")
            ("sequence", po::value<string>(&sequence_file), "FastA sequence file")
			("small-repeats", po::value <bool>(&small_repeats)->default_value(false), "optimize for small repeats")
			("score-out", po::value<string>(&output2file)->default_value(""), "output with corresponding score and alignment info ")
			("solid", po::value<bool>(&solid_seed)->default_value(0), "use solid/exact seeds?")
			("sp", po::value <score_t>(&min_spscore)->default_value(0), "minimum Sum-of-Pairs alignment score")
			("tandem", po::value <bool>(&allow_tandem)->default_value(true), "allow tandem repeats?")
			("two-hits", po::value<bool>(&two_hits)->default_value(false), "require two hits within w to trigger gapped extension?")
			("u", po::value<float>(&pGoUnrelated)->default_value(0.001f), "Transition to Unrelated")			
			("unalign", po::value<bool>(&unalign)->default_value(true), "unalign non-homologous sequence?")
			("w", po::value<int>(&w)->default_value(0), "max gap width ")
			("xmfa", po::value<string>(&xmfa_file)->default_value(""), "XMFA format output")
            ("xml", po::value<string>(&xml_file)->default_value(""), "XML format output")
			("z", po::value <unsigned>(&seed_weight)->default_value(0), "seed weight")

        ;

		if( argc < 2 )
		{
            cout << desc << "\n";
            return 1;
		}
                
		

        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);    

        if (vm.count("help")) {
            cout << desc << "\n";
            return 1;
        }

		if (large_repeats && small_repeats)
		{
			cout << "which is it? small or large? can't optimize for both!\n";
			return 1;
		}
		if (seed_weight < 3) {
            cout << "Invalid seed weight, minimum size is 3!\n";
            return 1;
        }
        if (vm.count("rmin")) {
            cout << "setting minimum multiplicity to " 
                 << rmin << ".\n";
        } else {
            cout << "Using default minimum multiplicity (2).\n";
        }

        if (vm.count("rmax")) {
            cout << "setting maximimum multiplicity to " 
                 << rmax << ".\n";
        } else {
            cout << "Using default maximum multiplicity (500).\n";
        }

	    if (rmin > rmax) 
        {
            cout << "rmin > rmax, setting rmax == rmin\n";
            rmax = rmin;
        } 
        if (rmin < 2)
        {
            cout << "rmin < 2, setting rmin == 2\n";
            rmin = 2;
        }
        if (rmax < 2)
        {
            cout << "rmax < 2, setting rmax == 2\n"; 
            rmax = 2;
        }
		if (percent_id > 1 )
			percent_id = 1.0;
		if (vm.count("z")) {
            cout << "seed weight set to " 
                 << seed_weight << ".\n";
        } else {
            cout << "Using default seed weight.\n";
        }
    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }
	


	// Roadmap: 
	// 1. identify seed matches using a Sorted Mer List
	// 2. create a "UngappedMatchRecord" for each seed match and put in the match record list
	// 3. create a Match Position Lookup Table
	// 4. create a multiplicity priority queue
	// 5. create an (empty) Novel Subset Match Record list
	// 6. extend all matches!
	// 7. create a procrastination queue for novel subset matches
	// 8. extend all novel subset matches!
	// 9. create a final list of matches
	// 10. score matches
	// 11. report matches


	//
	// part 1, load sequence and find seed matches using SML and a repeat class...
	//
	MatchList seedml;
	seedml.seq_filename = vector< string >( 1, sequence_file );
	seedml.sml_filename = vector< string >( 1, seedml.seq_filename[0] + ".sslist");
	//seedml.LoadSequences( &cout );
	LoadSequences( seedml, &cout );
	if( seed_weight == 0 )
		seed_weight = (int)((double)getDefaultSeedWeight( seedml.seq_table[0]->length() ) * .9);
	
	int seed_rank = 0;
	if ( solid_seed )
	{
		seed_rank = INT_MAX;
		std::cout << "Using solid seed" << std::endl;
	}
	seedml.LoadSMLs( seed_weight, &cout, seed_rank, solid_seed, !load_sml );
	int64 seed = getSeed( seed_weight, seed_rank);
	uint seed_size = getSeedLength( seed );

    if (min_spscore < 0 )
        min_spscore = 0;

	if( w == 0 )
		w = seed_weight * 3;	// default value
	else if( w < 0 )
	{
		w = 0;
		chain = false;
	}

	cout << "Using seed weight: " << seed_weight << " and w: " << w << endl;
	SeedMatchEnumerator sme;
	sme.FindMatches( seedml, rmin, rmax, only_direct );
	
    // need single nuc & kmer frequency
	string sequence = seedml.seq_table.at(0)->ToString();
	string uppercase = sequence;
	ToUPPER tupperware;
	std::transform(sequence.begin(),sequence.end(), uppercase.begin(), tupperware);
    kmersize =1;
	map<string,gnSeqI> polyfreq;
	map<string,gnSeqI> monofreq;
	map<string, gnSeqI>::iterator it;
	for (gnSeqI i = 0; i <= uppercase.size()-kmersize; i++)
	{
	   string kmer = uppercase.substr(i,kmersize);
	   string nucleotide = uppercase.substr(i,1);
	   if( nucleotide[0] != 'A' &&
			nucleotide[0] != 'C' &&
			nucleotide[0] != 'G' &&
			nucleotide[0] != 'T' )
			nucleotide[0] = 'A';
	   for( size_t kI = 0; kI < kmer.size(); kI++ )
		   if( kmer[kI] != 'A' &&
				kmer[kI] != 'C' &&
				kmer[kI] != 'G' &&
				kmer[kI] != 'T' )
				kmer[kI] = 'A';

	   polyfreq[kmer] +=1;	
	   monofreq[nucleotide] +=1;	
	   //insert( const string& val );
	   //it = find( const string& mer );
       //it->second+=1;	   
	}
    
	Params hmm_params = getAdaptedHoxdMatrixParameters( double(monofreq["G"]+monofreq["C"])/double(sequence.size()) );

	if (percent_id > 0 )
		adaptToPercentIdentity(hmm_params, percent_id);
	hmm_params.iGoHomologous = pGoHomo;
	hmm_params.iGoUnrelated = pGoUnrelated;

	//
	// part 2, convert to match records
	//
	vector< UngappedMatchRecord* > match_record_list(seedml.size());
	size_t component_count = 0;
    bool all_components_overlap = false;
    
    bool prev_overlaps = false;
    uint mi_multiplicity = 0;
    uint mi2_multiplicity = 0;
    uint num_components = 0;
    int overlap_size = 1;
    int hit_match =0;
   
	cout << "Total number of seed matches found: " << seedml.size() << endl;
	vector< pair< int64, UngappedMatchRecord* > > seed_sort_list;
    for( size_t mI = 0; mI < seedml.size(); ++mI )
	{
		UngappedMatchRecord tmp( seedml[mI]->SeqCount(), seedml[mI]->AlignmentLength() );
		match_record_list[mI] = tmp.Copy();
		
		for( size_t seqI = 0; seqI < seedml[mI]->SeqCount(); seqI++ )
		{
			match_record_list[mI]->SetStart( seqI, seedml[mI]->Start( seqI ) );
			match_record_list[mI]->SetLength( seedml[mI]->Length( seqI ), seqI );
		}
        seed_sort_list.push_back(make_pair(match_record_list[mI]->LeftEnd(0), match_record_list[mI]));
        component_count += seedml[mI]->SeqCount();
        seedml[mI]->Free();
    }
    std::sort( seed_sort_list.begin(), seed_sort_list.end() );
	// write seeds to file if requested
	
	ofstream seed_out;
	if ( seed_file.size() > 0)
		seed_out.open(seed_file.c_str());
    //
	// part 3, create a match position lookup table
	//
	vector< pair< gnSeqI, MatchPositionEntry > > mplt_sort_list( component_count );
	size_t compI = 0;
	for( size_t mI = 0; mI < seed_sort_list.size(); ++mI )
	{
		UngappedMatchRecord* mr = seed_sort_list[mI].second;
		if ( seed_file.size() > 0)
			seed_out << *mr << endl;
		for( size_t seqI = 0; seqI < mr->Multiplicity(); ++seqI )
			mplt_sort_list[compI++] = make_pair( mr->LeftEnd( seqI ), make_pair( mr, seqI ) );
	}
	// pairs get ordered on the first element by default 
	std::sort( mplt_sort_list.begin(), mplt_sort_list.end() );
	gnSeqI seq_length = seedml.seq_table[0]->length();
	MatchPositionLookupTable match_pos_lookup_table( seq_length+1, make_pair( (UngappedMatchRecord*)NULL, 0 ) );
	for( size_t i = 0; i < mplt_sort_list.size(); ++i )
    {
		//if ( seed_file.size() > 0)
		//	seed_out << (*(UngappedMatchRecord*)mplt_sort_list[i].second.first) << endl;
        //cerr << mplt_sort_list[i].first << endl;
		match_pos_lookup_table[ mplt_sort_list[i].first ] = mplt_sort_list[i].second;

    }

	//
	// part 4, create a procrastination queue
	//
	ProcrastinationQueue procrastination_queue( match_record_list );

	//
	// part 5, create an (empty) Novel Subset Match Record list
	//
	vector< GappedMatchRecord* > novel_subset_list;

	size_t superset_count = 0;
	size_t chainable_count = 0;
	size_t subset_count = 0;
	size_t novel_subset_count = 0;

	boost::dynamic_bitset<> left_lookups(seedml.seq_table[0]->length(), false);
	boost::dynamic_bitset<> right_lookups(seedml.seq_table[0]->length(), false);

	//
	// part 6, extend all matches!
	//
	vector< GappedMatchRecord* > extended_matches;	/**< The extended matches will be chains of UngappedMatchRecords */

	//for extension
	PairwiseScoringScheme pss = PairwiseScoringScheme(hoxd_matrix,-100,-20);
	
	int curI = 0;
	uint curr_extensions = 0;
	uint max_extensions = 2;
	while(  procrastination_queue.end() > 0 )
	{
		int prevI = curI;
		curI +=1;
		if( (curI * 100) / procrastination_queue.size() != (prevI * 100) / procrastination_queue.size() )
		{
			cout << (curI * 100) / procrastination_queue.size() << "%..";
			cout.flush();
		}
		
		// pop the next match off the heap
		MatchRecord* umr = procrastination_queue.pop(); 
		// if the match has been subsumed then skip it
		if( umr->subsuming_match != NULL )
			continue;
		if( umr->dont_extend == true )
			continue;

//		if( umr == (MatchRecord*)0x01335878 )
//			cout << "umr:\n" << *(UngappedMatchRecord*)umr << endl;

		GappedMatchRecord* M_i = dynamic_cast<GappedMatchRecord*>(umr);
		if( M_i == NULL )
		{
			// create a new gapped match record for M_i
			GappedMatchRecord gmr( *(UngappedMatchRecord*)umr );
			M_i = gmr.Copy();
			umr->subsuming_match = M_i;
			M_i->chained_matches.push_back( umr );
			vector< size_t > component_map( M_i->SeqCount() );
			for( size_t i = 0; i < component_map.size(); ++i )
				component_map[i] = i;
			M_i->chained_component_maps.push_back(component_map);
			swap(umr->subsumption_component_map, component_map);	// swap avoids reallocation
			// update superset and subset links
			for( int dI = 1; dI > -2; dI -= 2 )
			{
				MatchLink& ij_link = getSuperset(M_i,dI);
				if( ij_link.superset != NULL )
				{
					ij_link.subset = M_i;
					unlinkSuperset(umr,dI);
					int parity = M_i->Orientation(0) == ij_link.superset->Orientation(ij_link.sub_to_super_map[0]) ? 1 : -1;
					getSubsets(ij_link.superset,-dI*parity).push_back(ij_link);
				}
				vector< MatchLink >& subsets = getSubsets(M_i,dI);
				for( size_t subI = 0; subI < subsets.size(); ++subI )
				{
					subsets[subI].superset = M_i;
					int parity = M_i->Orientation(subsets[subI].sub_to_super_map[0]) == subsets[subI].subset->Orientation(0) ? 1 : -1;
					getSuperset(subsets[subI].subset, -dI*parity).superset = M_i;
				}
				getSubsets(umr,dI).clear();	// so that validate works...
			}
		}
        else
            cerr << "castdebugme!!\n" << endl;
     
        M_i->extended = true;
		extended_matches.push_back( M_i );
		
		// extend the match in each direction 
		// if a superset exists use that first
		// otherwise create a neighborhood list
		int direction = 1;	// leftward == 1, rightward == -1, done == -3
		//int direction = -1;	// leftward == 1, rightward == -1, done == 3
		int last_linked = 0;	// stores the group type that was chained.  1 == superset, 2 == chainable, 0 == none
		vector< NeighborhoodGroup > left_deferred_subsets;
		vector< NeighborhoodGroup > right_deferred_subsets;
		vector< NeighborhoodGroup > left_deferred_novel_subsets;
		vector< NeighborhoodGroup > right_deferred_novel_subsets;

		score_t score = 0;
		vector< gnSequence* > seqtable( M_i->SeqCount(), seedml.seq_table[0] );
		vector< string > alignment;
		vector<score_t> scores;
		bool extended = false;
		while( direction > -2 )
		{
			last_linked = 0;
			
			// check for superset
			if( getSuperset(M_i, direction).superset != NULL  )
				supersetLinkExtension( M_i, direction, last_linked, left_deferred_subsets, right_deferred_subsets, chain );

//          else
//          A hack to allow our chaining to work without novel subsets would be to
//          perform an additional neighborhood list lookup after superset link
//          extension even if no chainables are found during link extension.
			else
			{
				//
				// perform a neighborhood list extension, 
				// looks for neighboring matches in the match position lookup table
				// 
				vector< NeighborhoodGroup > superset_list;
				vector< NeighborhoodGroup > chainable_list;
				vector< NeighborhoodGroup > subset_list;
				vector< NeighborhoodGroup > novel_subset_list;
				//tjt: ok
				neighborhoodListLookup( M_i, match_pos_lookup_table,
								superset_list, chainable_list, subset_list, novel_subset_list,
								direction, seed_size, w, left_lookups, right_lookups, NULL);

				// tallies for debugging
				superset_count += superset_list.size();
				chainable_count += chainable_list.size();
				subset_count += subset_list.size();
				
				// now process each type of neighborhood group
				// supersets are already done.  happy chrismakwanzuhkkah
				// then process chainable
				processChainableMatches( M_i, chainable_list, direction, last_linked, find_novel_subsets, chain );

				// defer subset processing
				for( size_t gI = 0; gI < subset_list.size(); gI++ )
				{
					vector< NeighborhoodGroup >& cur_subset_list = selectList( left_deferred_subsets, right_deferred_subsets, direction );
					cur_subset_list.push_back( subset_list[gI] );
				}
				// defer novel subset processing
				vector< NeighborhoodGroup >& cur_novel_subset_list = selectList( left_deferred_novel_subsets, right_deferred_novel_subsets, direction );
				cur_novel_subset_list.clear();	// we only process novel subsets on the very last extension
				for( size_t gI = 0; gI < novel_subset_list.size(); gI++ )
					cur_novel_subset_list.push_back( novel_subset_list[gI] );

			} // end if no superset was found then do neighborhood list lookup
			//if find_novel_subsets not enabled, we can avoid this hack? is this true?
			if (!find_novel_subsets)
			{
				vector< NeighborhoodGroup > superset_list;
				vector< NeighborhoodGroup > chainable_list;
				vector< NeighborhoodGroup > subset_list;
				vector< NeighborhoodGroup > novel_subset_list;
				neighborhoodListLookup( M_i, match_pos_lookup_table,
									superset_list, chainable_list, subset_list, novel_subset_list,
									direction, seed_size, w, left_lookups, right_lookups, NULL);

				// defer subset processing
				for( size_t gI = 0; gI < subset_list.size(); gI++ )
				{
					vector< NeighborhoodGroup >& cur_subset_list = selectList( left_deferred_subsets, right_deferred_subsets, direction );
					cur_subset_list.push_back( subset_list[gI] );
				}
			}
			// if we didn't do a chaining or superset extension, try a gapped extension
			if( last_linked == 0 )
			{
                double e = 2.71828182845904523536;
				int rcode =FAILED;
                bool extend_it = false;
                 //extend_length = 0;
				vector<GappedMatchRecord*> novel_matches;	// M_e will contain the extension
				// only extend if two matches are chained if two-hits == true
                // its fast enough now that printing to screen actually slows things down...
				if( extend_chains && (!two_hits || (two_hits && M_i->chained_matches.size() > 1 )))
					rcode = ExtendMatch(M_i, seqtable, hmm_params, w, direction, novel_matches, gap_open, gap_extend, extension_window);
                
				if (rcode == FAILED || rcode == FIXME || novel_matches.size() == 0)
				{
					//end gapped extension  whenever extension fails.
					direction -=2;
					//direction +=2;
					continue;
				}
                else
                {
                    for (size_t mI = 0; mI < novel_matches.size(); mI++ )
                    {
                        //if (novel_matches.at(mI)->Multiplicity() != M_i->Multiplicity() )
                        //    continue;
                        GappedMatchRecord* M_e = novel_matches.at(mI);
                        M_e->extended = false;
                        
                        if (M_e->Multiplicity() > M_i->Multiplicity())//what does this mean??
                            continue;
                        else if (M_e->Multiplicity() == M_i->Multiplicity())
                        {
                            //immediately chainable!
                            if (direction > 0 && mI == novel_matches.size()-1)
                            {
                                
				                extend_it = true;
                                continue;
                            }
                            else if (direction < 0 && mI == 0)
                            {
				                extend_it = true;
                                continue;
                            }
                        }
                        vector< pair< gnSeqI, MatchPositionEntry > > mplt_sort_list( M_e->Multiplicity() );
                        vector< pair< gnSeqI, MatchPositionEntry > > final_mplt_sort_list;
	                    size_t compI = 0;
    	                
	                    for( size_t seqI = 0; seqI < M_e->Multiplicity(); ++seqI )
		                    mplt_sort_list[compI++] = make_pair( M_e->LeftEnd( seqI ), make_pair( M_e, seqI ) );
    	                
	                    // pairs get ordered on the first element by default 
	                    std::sort( mplt_sort_list.begin(), mplt_sort_list.end() );

                        //don't use novel match if it clobbers the existing left end in the MPLT
                        if (use_novel_matches )
                        {
                     
                            bool clobbers_existing_match = false;
                            for( size_t i = 0; i < mplt_sort_list.size(); ++i)
                            {
                                if (match_pos_lookup_table[ mplt_sort_list[i].first ].first != NULL )
                                {
                                    clobbers_existing_match = true;
                                    break;
                                }
                            }
                            if (! clobbers_existing_match )
                            {
                                for( size_t i = 0; i < mplt_sort_list.size(); ++i)
                                    match_pos_lookup_table[ mplt_sort_list[i].first ] =  mplt_sort_list[i].second;
                            }
                            
                        }
                        //now, during the subsequent call to neighborhoodListLookup(), we should
                        //find the novel homologous region and process it accordingly...
                    }
                }
                //update links appropriately, and we can take another round
				//through the evil megaloop, possibly discovering additional chainable
				//seeds or superset links.
				
				// need to update links by looking for matches in the region that was just extended over
				vector< NeighborhoodGroup > superset_list;
				vector< NeighborhoodGroup > chainable_list;
				vector< NeighborhoodGroup > subset_list;
				vector< NeighborhoodGroup > novel_subset_list;

                //if extend_it is true, it means that we can immediately extend
                //M_i with the corresponding result from ExtendMatch()
                if (extend_it)
                {
					M_i->extended = true;
                    //build a component map for the new record
				    vector< size_t > component_map( M_i->Multiplicity() );
				    for( size_t i = 0; i < component_map.size(); ++i )
					    component_map[i] = i;

                    GappedMatchRecord* M_t = NULL;
                    //leftward extension
                    if (direction > 0 )
                        M_t = novel_matches.back();
                    else
                        M_t = novel_matches.front();
                    
                    neighborhoodListLookup( M_i, match_pos_lookup_table,
					            superset_list, chainable_list, subset_list, novel_subset_list,
					            direction, seed_size, w, left_lookups, right_lookups, M_t);
                    
                    M_t->subsuming_match = M_i;
		            M_t->subsumption_component_map = component_map;
	                M_i->chained_matches.push_back( M_t );
	                M_i->chained_component_maps.push_back( component_map );
	                bool changed = extendRange(M_i, M_t, component_map);

                    
                }
                else
                {
                    
					if (!M_i->extended)
						M_i->extended = false;
                    GappedMatchRecord* M_t = NULL;
                    if (direction > 0 )
                        M_t = novel_matches.front();
                    else
                        M_t = novel_matches.back();

			        //we can't extend M_i, but we can classify all of the novel
                    //homologous regions with respect to M_i
			        neighborhoodListLookup( M_i, match_pos_lookup_table,
							    superset_list, chainable_list, subset_list, novel_subset_list,
							    direction, seed_size, w, left_lookups, right_lookups,M_t);
                    
                }
                extended = true;
				// now process each type of neighborhood group
				// if we have completely extended through a superset
				//   then we want to replace that part of the alignment with the superset
				// if the superset continues beyond the end of at least one component, then 
				// we want to create a superset link for it, and process it during a link extension
				if ( superset_list.size() > 0 && chain )
					processSupersetMatches( M_i, superset_list, direction, last_linked, true );
			
				// then process chainable
				if ( chainable_list.size() > 0 )
					processChainableMatches( M_i, chainable_list, direction, last_linked, find_novel_subsets, chain );

				// defer subset processing
				for( size_t gI = 0; gI < subset_list.size(); gI++ )
				{
					vector< NeighborhoodGroup >& cur_subset_list = selectList( left_deferred_subsets, right_deferred_subsets, direction );
					cur_subset_list.push_back( subset_list[gI] );
				}
				// defer novel subset processing
				vector< NeighborhoodGroup >& cur_novel_subset_list = selectList( left_deferred_novel_subsets, right_deferred_novel_subsets, direction );
				cur_novel_subset_list.clear();	// only process novel subsets from the very last extension
				for( size_t gI = 0; gI < novel_subset_list.size(); gI++ )
					cur_novel_subset_list.push_back( novel_subset_list[gI] );

                //just as before, if we didn't extend M_i, change directions and continue on
                if (!extend_it  )
                {
                    direction -=2;
                    continue;
                }

                //otherwise, enable another round of gapped extension in this direction.
            }
		}	// end loop over leftward and rightward extension

		//
		// finalize the alignment -- this resolves overlapping components into a single gapped alignment
		//
		//tjt: need to send finalize seq_table for muscle alignment
		if( M_i == (GappedMatchRecord*)0x00d37364 )
			cerr << "debugmult\n";
        
		// finally process novel subset
		for( int direction = 1; direction >-2; direction -= 2 )
		{
			vector< NeighborhoodGroup >& cur_novel_subset_list = selectList( left_deferred_novel_subsets, right_deferred_novel_subsets, direction );
			processNovelSubsetMatches(M_i, cur_novel_subset_list, find_novel_subsets, procrastination_queue, 
				seedml.seq_table, direction, w, last_linked, novel_subset_count );
		}

		//tjt: make sure finalize only gets called once!
		M_i->finalize(seedml.seq_table);
	    
         if( M_i->SeqCount() == 0)//what the hell?
            continue;
		//
		// process deferred subsets
		//
		for( int direction = 1; direction >-2; direction -= 2 )
		{
			vector< NeighborhoodGroup >& subset_list = selectList( left_deferred_subsets, right_deferred_subsets, direction );
			NeighborhoodGroupCompare ngc;
			NeighborhoodGroupComponentCompare ngcc;
			std::sort( subset_list.begin(), subset_list.end(), ngc );
			bool prev_linked = false;
			for( size_t sI = 0; sI < subset_list.size(); ++sI )
			{
				bool same_components = false;
				if( sI > 0 )
					same_components = ngcc.compare(subset_list[sI], subset_list[sI-1]) == 0;
				prev_linked = same_components? prev_linked : false;

				// check whether each of these ended up getting subsumed
				bool subsumed;
				bool partial;
				classifySubset( M_i, subset_list[sI], subsumed, partial );
				MatchRecord* M_j = subset_list[sI].get<0>();

				if( M_j->subsuming_match != NULL )
				{
					// sometimes duplicate MatchRecord pointers can exist in the subset list when a subset gets found
					// during a neighborhood list lookup but was already linked to a neighboring superset
					// in that case, we just skip the second entry...
					if(M_j->subsuming_match != M_i )
						cerr << "Error processing M_i " << M_i << ": match " << M_j << " was already subsumed\n";
					continue;
				}

				if( subsumed )
				{
					M_j->subsuming_match = M_i;
					M_j->subsumption_component_map = subset_list[sI].get<1>();
					unlinkSupersets(M_j);
					continue;
				}
				if( partial )
				{
					// create a novel subset record, mark this one as subsumed
					// just destroy it for now...
					M_j->dont_extend = true;
					
					unlinkSupersets(M_j);
					for( size_t mjI = 0; mjI < M_j->Multiplicity(); ++mjI )
                    {
						if( match_pos_lookup_table[M_j->LeftEnd(mjI)].first == M_j )
							match_pos_lookup_table[M_j->LeftEnd(mjI)] = make_pair((MatchRecord*)NULL,0);
					}
					
                    continue;
				}

				if( prev_linked )
				{
					// the previous subset has the same components as this one and was linked.
					// we may consider this one an 'extra' if all components are further away
					NeighborhoodGroup cur_group = subset_list[sI];
					subset_list.erase(subset_list.begin() + sI, subset_list.begin() + sI + 1);
					sI--;
					size_t dI = 0;
                    if (subset_list[sI].get<2>().size() < cur_group.get<2>().size())
                    {
                        //debugme: why would this happen?
                        //cerr << "subset_list[" << sI << "].get<2>().size() < cur_group.get<2>().size()" << endl;
                        //cerr << subset_list[sI].get<2>().size() << " < " <<  cur_group.get<2>().size() << endl;
                        //genome::breakHere();
                        continue;
                    }
					for( ; dI < cur_group.get<2>().size(); ++dI )
                    {
                        // if cur_group.get<2)()[dI] <= subset_list[sI].get<2>()[dI],
                        // component dI is closer than a component from the current subset
						if( cur_group.get<2>()[dI] <= subset_list[sI].get<2>()[dI] )
							break;
                    }
                    // all components were the same, yet further away, so consider this an 'extra' subset
					if( dI == cur_group.get<2>().size() )
					{
						// include this in a list of extra subsets
						boost::dynamic_bitset<> tmp_bs(M_i->Multiplicity());
						getExtraSubsets( M_i, direction ).push_back( MatchLink( (MatchRecord*)M_i, M_j, tmp_bs, cur_group.get<1>() ) );
						continue;
					}
					// else we've got a subset tie.
					if(print_warnings)
						cerr << "Subset tie, erasing M_j\n";

					//tjt: why do we need to erase the subset? later this will mean that we can't chain the two tied subsets..
					M_j->dont_extend = true;
					unlinkSupersets(M_j);
					
					for( size_t mjI = 0; mjI < M_j->Multiplicity(); ++mjI )
                    {
						if( match_pos_lookup_table[M_j->LeftEnd(mjI)].first == M_j )
							match_pos_lookup_table[M_j->LeftEnd(mjI)] = make_pair((MatchRecord*)NULL,0);
                    }
					
                    continue;
				}

				int parity = M_i->Orientation( subset_list[sI].get<1>()[0] ) == M_j->Orientation(0) ? 1 : -1;
				// if we have the following case:
				// --M_i-->   --M_j--   <--M_i-- ... ... --M_i-->   --M_j--   <--M_i--
				// then M_j may already be linked to M_i but on the other side
				if( getSuperset(M_j,direction*parity).superset == M_i )
					continue;
				unlinkSuperset( M_j, -direction*parity );
				// it's outside, just link it in
				// rebuild the superset component list
				boost::dynamic_bitset<> comp_list(M_i->Multiplicity(), false);

				for( size_t compI = 0; compI < subset_list[sI].get<1>().size(); ++compI )
                {
                    //debugme: why do I need to check this first?
                    if (  subset_list[sI].get<1>()[compI] != (std::numeric_limits<size_t>::max)())
                        comp_list.set(subset_list[sI].get<1>()[compI]);
                }
                getSuperset(M_j,-direction*parity) = MatchLink( M_i, M_j, comp_list, subset_list[sI].get<1>() );
				getSubsets(M_i,direction).push_back( getSuperset(M_j,-direction*parity));
                //getExtraSubsets(M_i,direction).push_back( getSuperset(M_j,-direction*parity));
				prev_linked = true;
			}
			subset_list.clear();
		}
	}
	cout << "\n# of calls to MUSCLE: " << ccount << endl;
	cout << "------------------------------"  << endl;
	cout << "superset count: " << superset_count << endl;
	cout << "chainable count: " << chainable_count << endl;
	cout << "subset count: " << subset_count << endl;
	cout << "novel subset count: " << novel_subset_count << endl;
	cout << "------------------------------"  << endl;
	// 
	// part 9, create a final list of local multiple alignments (already done in extended_matches)
	//
    vector< GappedMatchRecord* > &final = extended_matches;

	// part 10, score matches
	
	//create output stream
	ostream* output;
	ostream* output2;
 	ofstream score_out_file;
	ofstream aln_out_file;
	ofstream stats_out_file;
	if(stat_file != "" && stat_file != "-")
		stats_out_file.open( stat_file.c_str() );

	if(outputfile == "" || outputfile == "-")
		output = &cout;
	else
	{
		aln_out_file.open( outputfile.c_str() );
		output = &aln_out_file;
	}
	if(output2file == "" || output2file == "-")
		output2 = &cout;
	else
	{
		score_out_file.open( output2file.c_str() );
		output2 = &score_out_file;
	}
	vector< GappedMatchRecord* > scored;
	vector<score_t> scores_final;
	score_t score_final = 0;
    double e = 2.71828182845904523536;
    vector< GappedMatchRecord* >  filtered_final;
    int finalsize = final.size();
    uint alignment_count = 0;
    
    cout << "->Computing Sum-of-Pairs score of all lmas..." << endl;
	for( size_t fI = 0; fI < finalsize; fI++ )
	{
	    vector<string> alignment;
		vector< gnSequence* > seq_table( final[fI]->SeqCount(), seedml.seq_table[0] );
		mems::GetAlignment(*final[fI], seq_table, alignment);	// expects one seq_table entry per matching component
		//send temporary output format to file if requested
        if (alignment.at(0).size() >= min_repeat_length)
        {
			if(only_extended)
			{
				//we don't want it..
				if ( alignment.at(0).size() <= seed_size )
					continue;
			}
            score_final = 0;
            computeSPScore( alignment, pss, scores_final, score_final);
		    //*output << "#procrastAlignment " << ++alignment_count << endl << *final.at(fI) << endl;
            final[fI]->spscore = score_final;
            scored.push_back(final[fI]);

        }
        else
            continue;

	}
    
	if (!allow_redundant)
    {
		cout << "->Removing redudant lmas..." << endl;
	}
	//
	// remove overlapping regions
	//
    // 1) create a vector of CompactMatchRecord* with one entry for each nucleotide in the input sequence.
    //tjt: CompactMatchRecord is an attempt to reduce the space requirements for the method currently used to 
    //remove overlapping regions
	vector< CompactMatchRecord* > match_record_nt(sequence.size());
    for( size_t mI = 0; mI < match_record_nt.size(); ++mI )
	{
		CompactUngappedMatchRecord tmp( 1, 1 );
		match_record_nt[mI] = tmp.Copy();
		match_record_nt[mI]->SetStart( 0, mI );
		match_record_nt[mI]->SetLength( 1, 0 );
        match_record_nt[mI]->subsuming_match = NULL;
    }

    // 2) sort the result GappedMatchRecords 
	if (large_repeats)
		std::sort( scored.begin(), scored.end(), score_by_length );
	else if (small_repeats)
		std::sort( scored.begin(), scored.end(), scorecmp );
	else
		std::sort( scored.begin(), scored.end(), score_by_sp );
    for( size_t fI = 0; fI < scored.size(); fI++ )
    {
        //this shouldn't be the case, but let's be safe
	    if (scored.at(fI)->AlignmentLength() < 1)
            continue;

        //if user wants to remove all overlapping regions among lmas, let's do it!
        if (!allow_redundant)
        {
			
            //for each match compontent in M_i
            for ( size_t seqI = 0; seqI < scored.at(fI)->Multiplicity(); seqI++)
            {
                //if there is no match, we can't do a thing
                if( scored.at(fI)->LeftEnd(seqI) == NO_MATCH )
                    continue;

                //if left/right ends are good, set subsuming_match pointers
                if (scored.at(fI)->LeftEnd(seqI) < 4000000000u && scored.at(fI)->RightEnd(seqI) < 4000000000u)
                {
                    gnSeqI endI = scored.at(fI)->RightEnd(seqI);
                    gnSeqI startI = scored.at(fI)->LeftEnd(seqI);
                    for( ; startI < scored.at(fI)->RightEnd(seqI); startI++)
                    {
                        //3) Mark each entry in the MatchRecord* vector which corresponds to nucleotides contained within the current GMR.  
                        //A pointer to the current GMR can be >stored in each entry
                        if ( match_record_nt.at(startI)->subsuming_match == NULL)
                            match_record_nt.at(startI)->subsuming_match = scored.at(fI);
                    }
                }
        
                size_t left_crop_amt = 0;
                size_t right_crop_amt = 0;
                gnSeqI startI = scored.at(fI)->LeftEnd(seqI);
                //4) When a non-null entry is encountered in the vector, crop out that portion of the current GMR
                while(match_record_nt.at(startI)->subsuming_match != NULL && match_record_nt.at(startI)->subsuming_match != scored.at(fI) && startI < scored.at(fI)->RightEnd(seqI) && scored.at(fI)->Length(seqI) < 4000000000u) 
                {
                    startI++;
                    left_crop_amt++;
                }
                if (left_crop_amt > 0)
                {
                    if (left_crop_amt >= scored.at(fI)->Length(seqI))
                        scored.at(fI)->CropLeft( scored.at(fI)->Length(seqI)-1, seqI);
                    else
                        scored.at(fI)->CropLeft( left_crop_amt, seqI);
                }
                if (scored.at(fI)->LeftEnd(seqI) < 4000000000u && scored.at(fI)->RightEnd(seqI) < 4000000000u && scored.at(fI)->Length(seqI) < 4000000000u)
                {
                    startI = scored.at(fI)->RightEnd(seqI)-1;
                    //4) When a non-null entry is encountered in the vector, crop out that portion of the current GMR
                    while(match_record_nt.at(startI)->subsuming_match != NULL && match_record_nt.at(startI)->subsuming_match != scored.at(fI) && startI >= scored.at(fI)->LeftEnd(seqI))
                    {
                        startI--;
                        right_crop_amt++;
                    }
                }
                if (right_crop_amt > 0)
                {
                    
                    if (right_crop_amt >= scored.at(fI)->Length(seqI))
                        scored.at(fI)->CropRight( scored.at(fI)->Length(seqI)-1, seqI);
                    else
                        scored.at(fI)->CropRight( right_crop_amt, seqI);
                }
            }
        }
		//if ( left_crop_amt == 0 && right_crop_amt == 0)
		//	filtered_final.push_back(scored.at(fI));
		
        if (scored.at(fI)->AlignmentLength() >= min_repeat_length )
        {
			if(only_extended)
			{
				//we don't want it..
				if ( scored.at(fI)->AlignmentLength() <= seed_size )
					continue;
			}
            // yuck,recalculating sp score to update after removing overlapping regions.. 
            // couldn't I just subtract from the original score??
            vector<string> alignment;
	        vector< gnSequence* > seq_table( scored[fI]->SeqCount(), seedml.seq_table[0] );
	        mems::GetAlignment(*scored[fI], seq_table, alignment);	// expects one seq_table entry per matching component
            // 5) put all LMAs above min_repeat_length and min_spscore into final list of scored LMAs
            score_final = 0;
            computeSPScore( alignment, pss, scores_final, score_final);
            scored.at(fI)->spscore  = score_final;
            // pass it through a tandem repeat filter, too
            if ((scored.at(fI)->spscore > min_spscore && ( scored.at(fI)->tandem <= allow_tandem)))
                filtered_final.push_back(scored.at(fI));
        }
        

    }
    cout << "->Writing xmfa & xml output..." << endl;
    std::sort( filtered_final.begin(), filtered_final.end(), scorecmp );
    // write the output to xmfa
    writeXmfa( seedml, filtered_final, xmfa_file );

    // write the output to xml
    writeXML( seedml, filtered_final, xml_file );
    
	// 
	// part 11, report matches in scored order, by multiplicity then by spscore
	//
	output->setf(ios::fixed);
	output->precision(0);

    
	for( size_t sI = 0; sI < filtered_final.size(); ++sI )
	{
	    *output << "#procrastAlignment " << sI+1 << endl << *filtered_final[sI] << endl;
	    *output << "Alignment length: " << filtered_final[sI]->AlignmentLength() << endl;
        *output << "Score: " << filtered_final[sI]->spscore << endl;
	}
	
	///report highest scoring lma for each multiplicity
    cout << "->Calculating highest scoring lma for each multiplicity..." << endl;
	stats_out_file.setf(ios::fixed);
	stats_out_file.precision(0);
	int prev_multiplicity = 0;
    uint record_count = 0;
    for( size_t tI = 0; tI < filtered_final.size(); ++tI )
    {   if (filtered_final[tI]->Multiplicity() != prev_multiplicity)
        {    
            stats_out_file << "#" << record_count+1 << ": r= " << filtered_final[tI]->Multiplicity() << " l= " << filtered_final[tI]->AlignmentLength() << " s= " << filtered_final[tI]->spscore << endl;
            prev_multiplicity = filtered_final[tI]->Multiplicity();
            record_count++;
        }
        else
            continue;
    }
	// clean up
    cout << "->Cleaning up..." << endl;
	for( size_t eI = 0; eI < match_record_list.size(); ++eI )
		match_record_list[eI]->Free();
	for( size_t eI = 0; eI < novel_subset_list.size(); ++eI )
		if( novel_subset_list[eI]->subsuming_match != NULL  )
			novel_subset_list[eI]->Free();
	for( size_t eI = 0; eI < extended_matches.size(); ++eI )
		if( extended_matches[eI]->subsuming_match == NULL )
			extended_matches[eI]->Free();

	for( size_t seqI = 0; seqI < seedml.seq_table.size(); ++seqI )
		delete seedml.seq_table[seqI];
	for( size_t seqI = 0; seqI < seedml.sml_table.size(); ++seqI )
		delete seedml.sml_table[seqI];
    
    cout << "->Done!" << endl;
	return 0;
}

