/*******************************************************************************
 * $Id: MatchFinder.h,v 1.23 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _MatchFinder_h_
#define _MatchFinder_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/SortedMerList.h"
#include "libMems/Match.h"
#include "libMems/MatchList.h"
#include <list>
#include <iostream>
#include <boost/pool/pool_alloc.hpp>

namespace mems {

struct idmer{
	gnSeqI	position;	//starting position of this mer in the genome
	uint64 	mer; 		//the actual sequence
	sarID_t	id;			//the sequence identifier.
};

// typedef std::list<idmer, boost::fast_pool_allocator<idmer> > IdmerList;
// using boost::fast_pool_allocator<idmer> results in a significant speedup
// over std::allocator.  testing on a Salmonella vs. Y. pestis comparison shows
// a 30% speedup
typedef std::list<idmer> IdmerList;

const unsigned int PROGRESS_GRANULARITY = 100;

/**
 * This pure virtual class implements a general framework for finding
 * exactly matching mers.  It is extended by the MemHash and MemScorer
 * classes.
 * @see MemHash
 * @see MemScorer
 */
class MatchFinder : public genome::gnClone{
public:
	MatchFinder();
	~MatchFinder();
	MatchFinder(const MatchFinder& mf);
	virtual void Clear();
	/**
	 * Adds a sequence to use when searching for exact matches.
	 * @param sar A pointer to the sorted mer list for the new sequence
	 * @param seq A pointer to the genome::gnSequence corresponding to the new sequence.
	 */
	virtual boolean AddSequence( SortedMerList* sar, genome::gnSequence* seq = NULL );
	/**
	 * Given the index of a sequence and an index into the sorted mer list, this function
	 * will search the other sorted mer lists for the same mer.  This function returns the
	 * position of the mer in each sequence in the breakpoints vector.
	 */
	virtual void GetBreakpoint( uint32 sarI, gnSeqI startI, std::vector<gnSeqI>& breakpoints ) const;
	virtual uint32 Multiplicity(void){return seq_count;};
	/** NOT IMPLEMENTED: Sets the number of ambiguities allowed in a mer match*/
	virtual void SetAmbiguityTolerance(uint32 ambiguity_tol){ambiguity_tolerance = ambiguity_tol;}
	/** @return the number of ambiguities allowed in a mer match */
	virtual uint32 AmbiguityTolerance(){return ambiguity_tolerance;}
	/** @return The progress of the current operation.  Ranges from 0 to 100.  -1 indicates no computation is being performed */
	virtual float GetProgress() const {return m_progress;}

	/** Finds all the matches between the sequences */
	virtual void FindMatchSeeds();
	/** Finds all the matches between the sequences, starting at a particular offset */
	virtual void FindMatchSeeds( const std::vector<gnSeqI>& start_offsets );

	/**
	 * Logs progress to the designated ostream.  Set to null to skip progress logging.
	 */
	virtual void LogProgress( std::ostream* os );
	void SetOffsetLog( std::ostream* offset_stream ){ this->offset_stream = offset_stream; }
protected:
	/** 
	 * Searches for mer matches in a designated range of the sequence's sorted mer lists 
	 * @throws InvalidData thrown if the start_points are bad or if the sorted mer lists were sorted on different mer sizes
	 * @return true if completed searching, false if repetitive mers were encountered and FindMatches must be called again.
	 */
	virtual boolean SearchRange(std::vector<gnSeqI>& start_points, std::vector<gnSeqI>& search_len);
	/** Called whenever a mer match is found */
	virtual boolean HashMatch(IdmerList& match_list) = 0;
	virtual boolean EnumerateMatches(IdmerList& match_list);

	template< class MatchType >
	void FindSubsets(const MatchType& mhe, std::vector<MatchType>& subset_matches);

	template< class UngappedMatchType >
	void ExtendMatch(UngappedMatchType& mhe, std::vector<UngappedMatchType>& subset_matches, gnSeqI max_backward = GNSEQI_END, gnSeqI max_forward = GNSEQI_END);

	virtual SortedMerList* GetSar(uint32 sarI) const;
	std::vector<SortedMerList*> sar_table;
	std::vector<genome::gnSequence*> seq_table;
	
	uint32 mer_size;
	uint32 seq_count;
	uint32 ambiguity_tolerance;
	
	// for subset matches
	std::vector< std::vector< uint32 > > alpha_map;
	uint alpha_map_size;
	uint alphabet_bits;
	
	float m_progress;
	std::ostream* log_stream;

	uint64 mers_processed;	/**< The number of mers processed thus far */
	uint64 total_mers;	/**< The total number of mers to search */
	std::ostream* offset_stream;	/**< log for the current offset in each SML */
};

/** 
 * InvalidData exceptions are thrown when the input to an algorithm is invalid
 */
CREATE_EXCEPTION( InvalidData );

inline
SortedMerList* MatchFinder::GetSar(uint32 sarI) const{
	return sar_table[sarI];
}

inline
bool idmer_lessthan(idmer& a_v, idmer& m_v){
	return (a_v.mer < m_v.mer);// ? true : false;
};

//id less than function for STL sort functions
inline
bool idmer_id_lessthan(idmer& a_v, idmer& m_v){
	return (a_v.id < m_v.id);// ? true : false;
};



// takes as input a fully extended mem and returns the subset matches on the lower side
template< class MatchType >
void MatchFinder::FindSubsets(const MatchType& mhe, std::vector<MatchType>& subset_matches){

	SMLHeader head = GetSar( 0 )->GetHeader();
	uint shift_amt = 64 - head.alphabet_bits;
	uint rshift_amt = head.alphabet_bits * ( GetSar(0)->SeedLength() - 1 );

	uint seqI, alphaI;

	// initialize subset match data structures
	alpha_map_size = 1;
	alpha_map_size <<= alphabet_bits;
	if( alpha_map.size() != alpha_map_size ){
		alpha_map.clear();
		alpha_map.reserve( alpha_map_size );
		std::vector< uint32 > tmp_list;
		tmp_list.reserve( seq_count );
		for( uint alphaI = 0; alphaI < alpha_map_size; ++alphaI )
			alpha_map.push_back( tmp_list );
	}else{
		for( uint alphaI = 0; alphaI < alpha_map_size; ++alphaI )
			alpha_map[ alphaI ].clear();
	}
	
	
	for( seqI = 0; seqI < seq_count; ++seqI ){
		//check that all mers at the new position match
		int64 mer_to_get = mhe[ seqI ];
		if( mer_to_get == NO_MATCH )
			continue;
		if(mer_to_get < 0){
			mer_to_get *= -1;
			mer_to_get += mhe.Length() - GetSar(0)->SeedLength();
		}

		uint64 cur_mer = GetSar( seqI )->GetMer( mer_to_get - 1 );

		boolean parity;
		if( mhe[ seqI ] < 0 )
			parity = cur_mer & 0x1;
		else
			parity = !(cur_mer & 0x1);

		if( parity ){
			cur_mer >>= shift_amt;
		}else{
			cur_mer <<= rshift_amt;
			cur_mer = ~cur_mer;
			cur_mer >>= shift_amt;
		}

		alpha_map[ cur_mer ].push_back( seqI );

	}
	
	for( alphaI = 0; alphaI < alpha_map_size; ++alphaI ){
		if( alpha_map[ alphaI ].size() < 2 ){
			alpha_map[ alphaI ].clear();
			continue;
		}
		// this is a subset
		MatchType cur_subset = mhe;
		cur_subset.SetLength( mhe.Length() );
		for( uint sqI = 0; sqI < mhe.SeqCount(); ++sqI )
			cur_subset.SetStart( sqI, NO_MATCH );	// init everything to NO_MATCH
		for( uint subI = 0; subI < alpha_map[ alphaI ].size(); ++subI )
			cur_subset.SetStart( alpha_map[ alphaI ][ subI ], mhe[ alpha_map[ alphaI ][ subI ] ] );
		subset_matches.push_back( cur_subset );
		alpha_map[ alphaI ].clear();
	}
}

// BUGS:
// matches which span the end-start of a circular sequence will be hashed a second time
template< class UngappedMatchType >
void MatchFinder::ExtendMatch(UngappedMatchType& mhe, std::vector<UngappedMatchType>& subset_matches, gnSeqI max_backward, gnSeqI max_forward){
	uint64 cur_mer;
	uint64 mer_mask = GetSar(0)->GetSeedMask();

	//which sequences are used in this match?
	uint32* cur_seqs = new uint32[mhe.SeqCount()];
	uint32 used_seqs = 0;
	for(uint32 seqI = 0; seqI < mhe.SeqCount(); ++seqI){
		if(mhe[seqI] != NO_MATCH){
			cur_seqs[used_seqs] = seqI;
			++used_seqs;
		}
	}
	//First extend backwards then extend forwards.  The following loop does them both.
	int jump_size = GetSar(0)->SeedLength();
	uint extend_limit = 0;	/**< Tracks the distance to the most distant overlapping matching seed */
	uint extend_attempts = 0;	/**< Counts the total number of overlapping seeds checked */
	boolean extend_again = false;	/**< Set to true if any overlapping seeds matched, the search will be restarted from that point */
	for(uint32 directionI = 0; directionI < 4; ++directionI){
		//how far can we go?	
		//first calculate the maximum amount of traversal
		//then do fewer comparisons.
		int64 maxlen = GNSEQI_END;
		if(directionI == 0)
			maxlen = max_backward;
		else if(directionI == 1)
			maxlen = max_forward;
		else
			maxlen = GetSar(0)->SeedLength();
		for(uint32 maxI = 0; maxI < used_seqs; ++maxI)
			if(GetSar(cur_seqs[maxI])->IsCircular()){
				if(GetSar(cur_seqs[maxI])->Length() < maxlen)
					maxlen = GetSar(cur_seqs[maxI])->Length();
			}else if(mhe[cur_seqs[maxI]] < 0){
				int64 rc_len = GetSar(cur_seqs[maxI])->Length() - mhe.Length() + mhe[cur_seqs[maxI]] + 1;
				if( rc_len < maxlen)
					maxlen = rc_len;
			}else if(mhe[cur_seqs[maxI]] - 1 < maxlen)
				maxlen = mhe[cur_seqs[maxI]] - 1;
		uint32 j=0;
		uint32 i = used_seqs;	// set to used_seqs in case maxlen is already less than jump size.

		extend_limit = 0;
		extend_attempts = 0;

		while(maxlen - jump_size >= 0){
			mhe.SetLength(mhe.Length() + jump_size);
			maxlen -= jump_size;
			for(j=0; j < used_seqs; ++j){
				if(mhe[cur_seqs[j]] > 0){
					mhe.SetStart(cur_seqs[j], mhe[cur_seqs[j]] - jump_size);
					if(mhe[cur_seqs[j]] <= 0)
						mhe.SetStart(cur_seqs[j], mhe[cur_seqs[j]] + GetSar(cur_seqs[j])->Length());
				}
			}
			//check that all mers at the new position match
			int64 mer_to_get = mhe[cur_seqs[0]];
			if(mer_to_get < 0){
				mer_to_get *= -1;
				mer_to_get += mhe.Length() - GetSar(0)->SeedLength();
			}
			cur_mer = GetSar(cur_seqs[0])->GetSeedMer(mer_to_get - 1);
			boolean parity;
			if( mhe[cur_seqs[0]] < 0 )
				parity = cur_mer & 0x1;
			else
				parity = !(cur_mer & 0x1);
			cur_mer &= mer_mask;

			for(i=1; i < used_seqs; ++i){
				mer_to_get = mhe[cur_seqs[i]];
				if(mer_to_get < 0){
					//Convert the cur_seqs[i] entry since negative implies reverse complement
					mer_to_get *= -1;
					mer_to_get += mhe.Length() - GetSar(0)->SeedLength();
				}
				uint64 comp_mer = GetSar(cur_seqs[i])->GetSeedMer(mer_to_get - 1);
				boolean comp_parity;				
				if( mhe[cur_seqs[i]] < 0 )
					comp_parity = comp_mer & 0x1;
				else
					comp_parity = !(comp_mer & 0x1);
				comp_mer &= mer_mask;
				
				if(cur_mer != comp_mer || parity != comp_parity ){
					if( directionI < 2 )
						maxlen = 0;
					break;
				}
			}
			extend_attempts += jump_size;
			if( i == used_seqs )
				extend_limit = extend_attempts;
			if( directionI > 1 && extend_attempts == GetSar(0)->SeedLength() )
				break;
		}
		//this stuff cleans up if there was a mismatch
		if(i < used_seqs){
			mhe.SetLength(mhe.Length() - jump_size);
			for(;j > 0; j--){
				if(mhe[cur_seqs[j - 1]] >= 0)
					mhe.SetStart(cur_seqs[j - 1], mhe[cur_seqs[j - 1]] + jump_size);
			}
		}
		// check whether any of the overlapping seeds matched.
		// if so, set the match to that length and set the flag to start the search again
		if( directionI > 1 && extend_attempts > 0 ){
			if( extend_limit > 0 )
				extend_again = true;
			// minus jump_size because the cleanup above already moved the length back a little
			int unmatched_diff = extend_attempts - extend_limit;
			if( i < used_seqs )
				unmatched_diff -= jump_size;
			if( (unmatched_diff > mhe.Length()) && unmatched_diff >= 0 )
				std::cerr << "oh sheethockey mushrooms\n";
			mhe.SetLength(mhe.Length() - unmatched_diff);
			for(j=0; j < used_seqs; ++j){
				if(mhe[cur_seqs[j]] > 0){
					mhe.SetStart(cur_seqs[j], mhe[cur_seqs[j]] + unmatched_diff);
					if(mhe[cur_seqs[j]] > GetSar(cur_seqs[j])->Length() )
						mhe.SetStart(cur_seqs[j], mhe[cur_seqs[j]] - GetSar(cur_seqs[j])->Length() );
				}
			}
		}
		//Invert the sequence directions so that we extend in the other direction
		//next time through the loop.  The second time we do this we are setting
		//sequence directions back to normal.
		mhe.Invert();

		//if we've already been through twice then decrease the jump size
		if(directionI >= 1)
			jump_size = 1;
		if( directionI == 3 && extend_again ){
			directionI = -1;	// will become 0 on next iteration
			jump_size = GetSar(0)->SeedLength();
			extend_again = false;
		}
	}
	// after the match has been fully extended, search for subset matches
	// this code only works when using SOLID seeds-- so it's been disabled
/*	if( used_seqs > 2 ){
		FindSubsets( mhe, subset_matches );
		mhe.Invert();
		FindSubsets( mhe, subset_matches );
		mhe.Invert();
	}
*/
	// set the subsets so their reference sequence is always positive
	for(uint32 subsetI = 0; subsetI < subset_matches.size(); ++subsetI){
		if( subset_matches[subsetI][subset_matches[subsetI].FirstStart()] < 0 )
			subset_matches[subsetI].Invert();
		subset_matches[subsetI].CalculateOffset();
	}

	delete[] cur_seqs;
}



}

#endif	//_MatchFinder_h_
