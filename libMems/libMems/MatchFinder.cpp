/*******************************************************************************
 * $Id: MatchFinder.cpp,v 1.39 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing, copying, and modification
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/MatchFinder.h"

using namespace std;
using namespace genome;
namespace mems {

MatchFinder::MatchFinder(){
	mer_size = DNA_MER_SIZE;
	seq_count = 0;
	ambiguity_tolerance = 0;
	m_progress = -1;
	log_stream = NULL;
	offset_stream = NULL;
}

//make sure this calls the destructor on each element
MatchFinder::~MatchFinder(){
}

MatchFinder::MatchFinder(const MatchFinder& mf){
	mer_size = mf.mer_size;
	seq_count = mf.seq_count;
	ambiguity_tolerance = mf.ambiguity_tolerance;

	m_progress = mf.m_progress;
	sar_table = mf.sar_table;
	seq_table = mf.seq_table;
	log_stream = mf.log_stream;
	offset_stream = mf.offset_stream;
}

void MatchFinder::Clear(){
	mer_size = DNA_MER_SIZE;
	seq_count = 0;
	ambiguity_tolerance = 0;
	m_progress = -1;
	sar_table.clear();
	seq_table.clear();
	log_stream = NULL;
	offset_stream = NULL;
}

void MatchFinder::LogProgress( ostream* os ){
	log_stream = os;
}

boolean MatchFinder::AddSequence( SortedMerList* sar, gnSequence* seq ){
	if(sar == NULL){
		Throw_gnExMsg( NullPointer(), "Null SortedMerList pointer" );
	}
	if(sar == NULL){
		Throw_gnExMsg( NullPointer(), "Null gnSequence pointer" );
	}
	
	//check for consistency between sequence length and sorted mer list lengths
/*	if(seq != NULL && seq->length() != sar->Length()){
		cerr << "MatchFinder::AddSequence: Error mismatched sml and sequence length.\n";
		cerr << "Seq length: " << seq->length() << "\tSML length: " << sar->Length() << endl;
		DebugMsg("MatchFinder::AddSequence: Error mismatched sml and sequence length.");
		return false;
	}
*/	
	//passed checks, add it to the data structures
	sar_table.push_back(sar);
	++seq_count;
	if(seq != NULL){
		seq_table.push_back(seq);
	}
	
	SMLHeader header = sar->GetHeader();
	alphabet_bits = header.alphabet_bits;
	
	return true;

}

void MatchFinder::GetBreakpoint( uint32 sarI, gnSeqI startI, vector<gnSeqI>& breakpoints ) const{
	breakpoints.clear();
	
	//put the mer to break on in break_mer
	bmer break_mer  = (*GetSar(sarI))[startI];
	uint64 mer_mask = GetSar(sarI)->GetSeedMask();
	bmer prev_mer = break_mer;
	//search backwards for the first index of this mer
	while((prev_mer.mer & mer_mask) == (break_mer.mer & mer_mask)){
		if(startI == 0){
			startI--;
			break;
		}
		startI--;
		prev_mer = (*GetSar(sarI))[startI];
	}
	++startI;

	//find the mer's location in the other sorted mer lists
	for(uint32 i=0; i < seq_count; ++i){
		if(i == sarI){
			breakpoints.push_back(startI);
		}else{
			gnSeqI cur_start;
			if(GetSar(i)->FindMer(break_mer.mer, cur_start)){
				//we found a match, see how far backwards we can go.
				int64 cur_matchI = cur_start;
				bmer matchmer = (*GetSar(i))[cur_start];
				while(cur_matchI >= 0 && ((matchmer.mer & mer_mask) == (break_mer.mer && mer_mask))){
					cur_matchI--;
					matchmer = (*GetSar(i))[cur_start];
				}
				cur_start = cur_matchI+1;
			}
			breakpoints.push_back(cur_start);
		}
	}
}

void MatchFinder::FindMatchSeeds(){
	vector<gnSeqI> start_points;

	for(uint32 i=0; i < sar_table.size(); ++i){
		start_points.push_back(0);
	}
	FindMatchSeeds( start_points );
}

void MatchFinder::FindMatchSeeds( const vector<gnSeqI>& start_offsets ){
	vector<gnSeqI> start_points = start_offsets;
	vector<gnSeqI> search_len;
	// keep track of the number of mers processed and the total for progress reporting
	mers_processed = 0;
	total_mers = 0;
	m_progress = -1;
	for(uint32 i=0; i < sar_table.size(); ++i){
		search_len.push_back(GNSEQI_END);
		total_mers += search_len[i] == GNSEQI_END ? sar_table[i]->Length() : search_len[i];
		mers_processed += start_points[ i ];
	}
	while( !SearchRange(start_points, search_len) ){
		mers_processed = 0;
		for( uint32 seqI = 0; seqI < sar_table.size(); ++seqI ){
			if( offset_stream != NULL ){
				if( seqI > 0 )
					*offset_stream << '\t';
				*offset_stream << start_points[ seqI ];
			}
			mers_processed += start_points[ seqI ];
		}
		if( offset_stream != NULL ){
			*offset_stream << endl;
			offset_stream->flush();
		}
	}
}

#define MER_REPEAT_LIMIT 1000 // The maximum number of matching mers before they are completely
								// ignored.

boolean print_sp = false;
//startI must be 0
//At most search_length mers in any one genome will be checked.
boolean MatchFinder::SearchRange(vector<gnSeqI>& start_points, vector<gnSeqI>& search_len){
	//picked a semi-arbitrary number for buffer size.
	uint32 MER_BUFFER_SIZE = 10000;
	vector<uint32> mer_index;   // stores the indexes of the current mers in mer_vector
	vector<uint32> mer_baseindex;   // stores the index in the SortedMerList of each of the first mers in mer_vector
	IdmerList cur_mers;	// stores the current mers.
	IdmerList cur_match;	// stores the current matching mers.
	list<uint32> sar_hitlist;	// list of sars to replace
	uint32 read_size;
	
	//make sure there is at least one sequence
	if(sar_table.size() < 1)
		return true;
	
	//check for consistency in seed patterns.
	uint64 mer_mask = sar_table[0]->GetSeedMask();
	uint64 seed = sar_table[0]->Seed();
	mer_size = sar_table[0]->SeedWeight();
	for(uint32 maskI = 0; maskI < sar_table.size(); ++maskI){
		if(seed != sar_table[maskI]->Seed()){
			Throw_gnExMsg(InvalidData(), "Different seed patterns.");
		}
	}
	
	//check that start_points and end_points are ok.
	if((start_points.size() != sar_table.size()) || (search_len.size() != sar_table.size())){
		Throw_gnExMsg(InvalidData(), "Inconsistent search range specification.");
	}
	
	//allocate buffer space
	// stores arrays of bmers for each sml.

	vector< vector< bmer > > mer_vector;
	for( uint vecI = 0; vecI < sar_table.size(); ++vecI ){
		vector< bmer > vec;
		mer_vector.push_back( vec );
	}

	//initialize the data structures
	idmer newmer;
	for(uint32 n = 0; n < sar_table.size(); ++n){
		read_size = MER_BUFFER_SIZE < search_len[n] ? MER_BUFFER_SIZE : search_len[n]; 
		mer_vector[n].reserve(read_size);
		sar_table[n]->Read(mer_vector[n], read_size, start_points[n]);
		mer_index.push_back(0);
		mer_baseindex.push_back(0);
		if( mer_vector[n].size() > 0 ){
			newmer.position = mer_vector[n][0].position;
			newmer.mer = mer_vector[n][0].mer & mer_mask;
			newmer.id = n;
			cur_mers.push_back(newmer);  //cur_mers gets the first mer from each sorted mer list
		}
	}
	
	if( print_sp ){
	cerr << "First mers are: " << mer_vector[0][0].mer << endl;
	cerr << "First mers are: " << mer_vector[1][0].mer << endl;
	cerr << "First mers are: " << mer_vector[2][0].mer << endl;
	print_sp = false;
	}	
	//nobody reads these fucking things.  why am i writing this.because my fucking 
	//roomate needs a goddamn roadmap......   ohhh ecstasy.... haptic pumpkins

	//loop while there is data to hash.
	cur_mers.sort(&idmer_lessthan);
	while(cur_mers.size() > 0){
		IdmerList::iterator mer_iter = cur_mers.begin();
		sarID_t cur_id = mer_iter->id;
		//first check for matches across genomes.
		if(cur_match.size() > 0){
			if(mer_iter->mer > cur_match.begin()->mer){
				//we are done with this matching.  hash it.
				if(cur_match.size() > 1)
					EnumerateMatches(cur_match);
				cur_match.clear();
			}else if(mer_iter->mer < cur_match.begin()->mer){
				//error checking stuff.
				ErrorMsg("Horrible error occurred!!\n");
			}
		}

		if( cur_match.size() > MER_REPEAT_LIMIT ){
			// scan past the repetitive mers
			// create the lexicographically next mer
			uint64 next_mer = cur_match.begin()->mer;
			next_mer += ~mer_mask + 1;
//			cerr << "Searching to: " << next_mer << endl;
			gnSeqI next_pos = 0;
			uint seqI = 0;
			for( ; seqI < sar_table.size(); ++seqI ){
				if( !sar_table[ seqI ]->FindMer( next_mer, next_pos ))
					++next_pos;
				if( next_pos < sar_table[ seqI ]->SMLLength() )
					break;
			}
			vector< gnSeqI > old_starts = start_points;
			if( seqI < sar_table.size() )
				GetBreakpoint( seqI, next_pos, start_points );
			for( int spI = 0; spI < start_points.size(); ++spI ){
				// don't allow it to move backwards!
				start_points[ spI ] = start_points[ spI ] < mer_index[ spI ] + mer_baseindex[ spI ] + old_starts[ spI ] ? old_starts[ spI ] + mer_index[ spI ] + mer_baseindex[ spI ] : start_points[ spI ];
				if( spI < seqI )
					start_points[ spI ] = sar_table[ spI ]->SMLLength();
			} 
			return false;
		}
		//check for matches within the same genome
		gnSeqI merI = mer_index[cur_id];
		boolean buffer_exhausted = merI < mer_vector[cur_id].size() ? false : true;
		while(!buffer_exhausted && (mer_iter->mer == (mer_vector[cur_id][merI].mer & mer_mask))){
			newmer.position = mer_vector[cur_id][merI].position;
			newmer.mer = mer_vector[cur_id][merI].mer & mer_mask;
			newmer.id = cur_id;
			cur_match.push_back(newmer);
			++merI;
			++mer_index[cur_id];
			//check if we've exhausted our buffer
			if(merI == mer_vector[cur_id].size())
				buffer_exhausted = true;
		}

		if(buffer_exhausted)
		{
			//if we've exhausted our buffer then refill it
			mer_baseindex[cur_id] += mer_vector[cur_id].size();
			
			// update the mers processed
			mers_processed += mer_vector[cur_id].size();
			float64 m_oldprogress = m_progress;
			m_progress = ((float64)mers_processed / (float64)total_mers) * PROGRESS_GRANULARITY;
			if( log_stream != NULL ){
				if((int)m_oldprogress != (int)m_progress){
					(*log_stream) << (int)((m_progress / PROGRESS_GRANULARITY) * 100) << "%..";
					log_stream->flush();
				}
				if(((int)m_oldprogress / 10) != ((int)m_progress / 10))
					(*log_stream) << std::endl;
			}
			uint32 read_size = MER_BUFFER_SIZE;
			if(MER_BUFFER_SIZE + mer_baseindex[cur_id] > search_len[cur_id])
				read_size = search_len[cur_id] - mer_baseindex[cur_id];

			sar_table[cur_id]->Read(mer_vector[cur_id], read_size, start_points[cur_id] + mer_baseindex[cur_id]);
			mer_index[cur_id] = 0;
			if(mer_vector[cur_id].size() == 0){
				//remove mer_iter so that this sar is forgotten
				cur_mers.erase(mer_iter);
			}
		}else{
			//if we haven't exhausted our buffer then we must have
			//run out of matching mers.
			//remove mer_iter and put in a new idmer with the same id
			cur_mers.erase(mer_iter);
			newmer.position = mer_vector[cur_id][merI].position;
			newmer.mer = mer_vector[cur_id][merI].mer & mer_mask;
			newmer.id = cur_id;
			mer_iter = cur_mers.begin();
			while(mer_iter != cur_mers.end() && mer_iter->mer < newmer.mer )
				++mer_iter;
			cur_mers.insert(mer_iter, newmer);
		}
		
	}
	//very last match in the dataset wasn't getting hashed.
    if(cur_match.size() > 1)
       EnumerateMatches(cur_match);

	return true;
}

boolean MatchFinder::EnumerateMatches( IdmerList& match_list ){
	//this must call HashMatch on every possible combination of matches in the list.
	if(match_list.size() == 2){
		//this is the smallest possible match.  simply hash it.
		return HashMatch(match_list);
	}
	
	match_list.sort(&idmer_id_lessthan);
	vector<uint32> id_start;
	vector<IdmerList::iterator> id_pos;
	vector<IdmerList::iterator> id_end;
	IdmerList::iterator iter = match_list.begin();
	IdmerList::iterator iter2 = match_list.begin();
	++iter2;
	id_start.push_back(0);
	id_pos.push_back(iter);
	for(uint32 i=0; iter2 != match_list.end(); ++i){
		if(iter->id != iter2->id){
			id_start.push_back(i);
			id_pos.push_back(iter2);
		}
		++iter;
		++iter2;
	}
	//the following loop iterates through all possible combinations of idmers with
	//different id's and hashes them.
	id_end = id_pos;
	id_end.push_back(match_list.end());
	while(true){
		IdmerList cur_match;
		for(uint32 k = 0; k < id_pos.size(); ++k){
			cur_match.push_back(*id_pos[k]);
		}
		HashMatch(cur_match);
		cur_match.clear();

		//increment the iterators (like an odometer)
		uint32 m = id_pos.size() - 1;
		while(true){
			++id_pos[m];
			if(id_pos[m] == id_end[m+1]){
				if(m == 0)
					return true;
				id_pos[m] = id_end[m];
				m--;
			}else
				break;
		}
	}

	return true;
}
/*
boolean MatchFinder::MatchAmbiguities(MatchHashEntry& mhe, uint32 match_size){
	if(ambiguity_tolerance == 0)
		return false;
			//check that all mers at the new position match
	//which sequences are used in this match?
	uint32* cur_seqs = new uint32[mhe.SeqCount()];
	uint32 used_seqs = 0;
	for(uint32 seqI = 0; seqI < mhe.SeqCount(); ++seqI){
		if(mhe[seqI] != NO_MATCH){
			cur_seqs[used_seqs] = seqI;
			++used_seqs;
		}
	}
	string cur_mer, mer_i;
	gnSequence mer_seq;
	int64 mer_to_get = mhe[cur_seqs[0]];
	if(mer_to_get < 0){
		mer_to_get *= -1;
		mer_to_get += mhe.Length() - mer_size;
	}
	cur_mer = seq_table[cur_seqs[0]]->subseq(mer_to_get, match_size).ToString();
	
	for(uint32 i=1; i < used_seqs; ++i){
		mer_to_get = mhe[cur_seqs[i]];
		if(mer_to_get < 0){
			//Convert the cur_seqs[i] entry since negative implies reverse complement
			mer_to_get *= -1;
			mer_to_get += mhe.Length() - mer_size;
		}
		mer_seq = seq_table[cur_seqs[i]]->subseq(mer_to_get, match_size);
		if(mer_seq.compare(cur_mer) != 0){
			delete[] cur_seqs;
			return false;
		}
		mer_i = mer_seq.ToString();
		uint32 ambiguity_count = 0;
		for(uint32 baseI = 0; baseI < match_size; ++baseI)
			if(cur_mer[baseI] != mer_i[baseI])
				++ambiguity_count;
		if(ambiguity_count > ambiguity_tolerance){
			delete[] cur_seqs;
			return false;
		}
	}
	delete[] cur_seqs;
	return true;
}
*/

} // namespace mems
