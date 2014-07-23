/*******************************************************************************
 * $Id: MemHash.cpp,v 1.32 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing, copying, and modification
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/MemHash.h"
#include "libGenome/gnFilter.h"
#include <list>
#include <map>
#include <sstream>

using namespace std;
using namespace genome;
namespace mems {

	MemHash::MemHash() : MatchFinder(), allocator( SlotAllocator<MatchHashEntry>::GetSlotAllocator() )

{
	table_size = DEFAULT_MEM_TABLE_SIZE;
	seq_count = 0;
	m_mem_count = 0;
	m_collision_count = 0;
	m_repeat_tolerance = DEFAULT_REPEAT_TOLERANCE;
	m_enumeration_tolerance = DEFAULT_ENUMERATION_TOLERANCE;
	//allocate the hash table
	mem_table.resize(table_size);
	mem_table_count.reserve( table_size );
	for(uint32 i=0; i < table_size; ++i)
		mem_table_count.push_back(0);
	match_log = NULL;
}

//make sure this calls the destructor on each element
MemHash::~MemHash(){
//	allocator.Free(allocated);
}

MemHash::MemHash(const MemHash& mh) : MatchFinder(mh), allocator( SlotAllocator<MatchHashEntry>::GetSlotAllocator() )
{
	*this = mh;
}

MemHash& MemHash::operator=( const MemHash& mh ){
	table_size = mh.table_size;
	mer_size = mh.mer_size;
	seq_count = mh.seq_count;
	m_mem_count = mh.m_mem_count;
	m_collision_count = mh.m_collision_count;
	m_repeat_tolerance = mh.m_repeat_tolerance;
	m_enumeration_tolerance = mh.m_enumeration_tolerance;
	mem_table.resize(table_size);
	for(uint32 i=0; i < table_size; ++i){
		mem_table_count.push_back(mh.mem_table_count[i]);
		mem_table[i] = mh.mem_table[i];
	}
	match_log = mh.match_log;
	return *this;
}

MemHash* MemHash::Clone() const{
	return new MemHash(*this);
}

void MemHash::ClearSequences()
{
	MatchFinder::Clear();
}

void MemHash::Clear()
{
	MatchFinder::Clear();
	m_mem_count = 0;
	m_collision_count = 0;
	m_repeat_tolerance = DEFAULT_REPEAT_TOLERANCE;
	m_enumeration_tolerance = DEFAULT_ENUMERATION_TOLERANCE;
	//clear the hash table
	for(uint32 listI = 0; listI < table_size; ++listI){
		mem_table[listI].clear();
		mem_table_count[ listI ] = 0;
	}
	match_log = NULL;

	allocator.Free(allocated);
	// WARNING! WARNING! WARNING! this will destroy ALL objects since the allocator has static lifetime!!
//	allocator.Purge();
}

void MemHash::SetTableSize(uint32 new_table_size){
	//allocate the hash table
	table_size = new_table_size;
	mem_table.clear();
	mem_table.resize(table_size);
	mem_table_count.clear();
	mem_table_count.resize(table_size,0);
}

boolean MemHash::CreateMatches(){
	MatchFinder::FindMatchSeeds();
	return true;
}

void MemHash::FindMatches( MatchList& ml ) {
	vector<gnSeqI> start_points;
	for( uint32 seqI = 0; seqI < ml.seq_table.size(); ++seqI ){
		start_points.push_back( 0 );
	}
	FindMatchesFromPosition( ml, start_points );
}

void MemHash::FindMatchesFromPosition( MatchList& ml, const vector<gnSeqI>& start_points ){
	for( uint32 seqI = 0; seqI < ml.seq_table.size(); ++seqI ){
		if( !AddSequence( ml.sml_table[ seqI ], ml.seq_table[ seqI ] ) ){
			ErrorMsg( "Error adding " + ml.seq_filename[seqI] + "\n");
			return;
		}
	}
	MatchFinder::FindMatchSeeds( start_points );

	GetMatchList( ml );
}

MatchList MemHash::GetMatchList() const{
	MatchList ml;
	GetMatchList( ml );
	ml.seq_table = seq_table;
	ml.sml_table = sar_table;

	return ml;
}

// an attempt to do this without sorting, which appears to be very slow...
boolean MemHash::EnumerateMatches( IdmerList& match_list )
{
	vector< uint > enum_tally(seq_count, 0);
	IdmerList::iterator iter = match_list.begin();
	IdmerList hash_list;
	for(; iter != match_list.end(); ++iter)
	{
		if( enum_tally[iter->id] < m_enumeration_tolerance )
		{
			hash_list.push_back(*iter);
		}
		if(enum_tally[iter->id] > m_repeat_tolerance)
			return true;
		++enum_tally[iter->id];
	}

	if(hash_list.size() > 1){
		if(m_enumeration_tolerance == 1)
			return HashMatch(hash_list);
		else
			return MatchFinder::EnumerateMatches( hash_list );
	}
	return true;
}

//why have separate hash tables? dunno.  no reason.  what was i thinking
// at that coffeehouse in portland when i wrote this crappy code?
// MemHashEntries use GENETICIST coordinates.  They start at 1, not 0.
boolean MemHash::HashMatch(IdmerList& match_list){
	//check that there is at least one forward component
//	match_list.sort(&idmer_id_lessthan);
	// initialize the hash entry
	MatchHashEntry mhe = MatchHashEntry(seq_count, GetSar(0)->SeedLength());
	mhe.SetLength( GetSar(0)->SeedLength() );
	
	//Fill in the new Match and set direction parity if needed.
	IdmerList::iterator iter = match_list.begin();
	for(; iter != match_list.end(); ++iter)
		mhe.SetStart(iter->id, iter->position + 1);
	SetDirection(mhe);
	mhe.CalculateOffset();
	if(mhe.Multiplicity() < 2){
		cout << "red flag " << mhe << "\n";
		cout << "match_list.size(): " << match_list.size() << endl;
	}else 
		AddHashEntry(mhe);

	return true;
}

void MemHash::SetDirection(MatchHashEntry& mhe){
	//get the reference direction
	boolean ref_forward = false;
	uint32 seqI=0;
	for(; seqI < mhe.SeqCount(); ++seqI)
		if(mhe[seqI] != NO_MATCH){
			ref_forward = !(GetSar(seqI)->GetDnaSeedMer(mhe[seqI] - 1) & 0x1);
			break;
		}
	//set directional parity for the rest
	for(++seqI; seqI < mhe.SeqCount(); ++seqI)
		if(mhe[seqI] != NO_MATCH)
			if(ref_forward == (GetSar(seqI)->GetDnaSeedMer(mhe[seqI] - 1) & 0x1))
				mhe.SetStart(seqI, -mhe[seqI]);
}

// Tries to add a new mem to the mem hash table
// If the mem already exists in the table, a pointer to it
// is returned.  Otherwise mhe is added and a pointer to
// it is returned.
MatchHashEntry* MemHash::AddHashEntry(MatchHashEntry& mhe){
	//first compute which hash table bucket this is going into
	int64 offset = mhe.Offset();

	uint32 bucketI = ((offset % table_size) + table_size) % table_size;
	vector<MatchHashEntry*>::iterator insert_he;
	insert_he = std::lower_bound(mem_table[bucketI].begin(), mem_table[bucketI].end(), &mhe, mhecomp);
//	insert_he = mem_table[bucketI].find(&mhe);
	if( insert_he != mem_table[bucketI].end() && (!mhecomp(*insert_he, &mhe) && !mhecomp(&mhe, *insert_he)) ){
		++m_collision_count;
		return *insert_he;
	}
	
	//if we made it this far there were no collisions
	//extend the mem into the surrounding region.
	vector<MatchHashEntry> subset_matches;
	if( !mhe.Extended() )
		ExtendMatch(mhe, subset_matches);

	MatchHashEntry* new_mhe = allocator.Allocate();
	new_mhe = new(new_mhe) MatchHashEntry(mhe); 
//	*new_mhe = mhe;
	allocated.push_back(new_mhe);
	
	// can't insert until after the extend!!
	insert_he = std::lower_bound(mem_table[bucketI].begin(), mem_table[bucketI].end(), new_mhe, mhecomp);
	mem_table[bucketI].insert(insert_he, new_mhe);

	// log it.
	if( match_log != NULL ){
		(*match_log) << *new_mhe << endl;
		match_log->flush();
	}
	
	// link up the subset matches
	for(uint32 subsetI = 0; subsetI < subset_matches.size(); ++subsetI){
		MatchHashEntry* submem = AddHashEntry( subset_matches[ subsetI ] );
	}
	
	++mem_table_count[bucketI];
	++m_mem_count;
	return new_mhe;
}

void MemHash::PrintDistribution(ostream& os) const{
    vector<MatchHashEntry*>::const_iterator mem_iter;
	gnSeqI base_count;
	for(uint32 i=0; i < mem_table_count.size(); ++i){
		mem_iter = mem_table[i].begin();
		base_count = 0;
		for(; mem_iter != mem_table[i].end(); ++mem_iter){
			base_count += (*mem_iter)->Length();
		}
		os << i << '\t' << mem_table_count[i] << '\t' << base_count << '\n';
	}
}

void MemHash::LoadFile(istream& mem_file){
	string tag;
	gnSeqI len;
	int64 start;
	MatchHashEntry mhe;
	getline( mem_file, tag );
	stringstream first_mum( tag );
	seq_count = 0;
	first_mum >> len;
	while( first_mum >> start ){
		seq_count++;
	}
	mhe = MatchHashEntry(seq_count, mer_size, MatchHashEntry::seed);
	first_mum.str( tag );
	first_mum.clear();
	for(uint32 seqI = 0; seqI < seq_count; seqI++){
		first_mum >> start;
		mhe.SetStart(seqI, start);
	}
	mhe.SetLength( len );
	mhe.CalculateOffset();
	AddHashEntry(mhe);
	
	while(mem_file.good()){
		mem_file >> len;
		if(!mem_file.good())
			break;
		mhe.SetLength(len);
		for(uint32 seqI = 0; seqI < seq_count; seqI++){
			mem_file >> start;
			mhe.SetStart(seqI, start);
		}
		//break if the stream ended
		if(!mem_file.good())
			break;
		mhe.CalculateOffset();
		AddHashEntry(mhe);
	}
}

void MemHash::WriteFile(ostream& mem_file) const{
	mem_file << "FormatVersion" << '\t' << 1 << "\n";
	mem_file << "SequenceCount" << '\t' << sar_table.size() << "\n";
	for(unsigned int seqI = 0; seqI < seq_count; seqI++){
		mem_file << "Sequence" << seqI << "File";
		gnGenomeSpec* specker = seq_table[seqI]->GetSpec();
		string sourcename = specker->GetName();
		if( sourcename == "" )
			sourcename = "null";
		mem_file << '\t' << sourcename << "\n";
		mem_file << "Sequence" << seqI << "Length";
		mem_file << '\t' << seq_table[seqI]->length() << "\n";
	}
	mem_file << "MatchCount" << '\t' << m_mem_count << endl;
	//get all the mems out of the hash table and write them out
    vector<MatchHashEntry*>::const_iterator mem_table_iter;
	for(uint32 i=0; i < table_size; i++){
		mem_table_iter = mem_table[i].begin();
		for(; mem_table_iter != mem_table[i].end(); mem_table_iter++)
			mem_file << **mem_table_iter << "\n";
	}
}


} // namespace mems
