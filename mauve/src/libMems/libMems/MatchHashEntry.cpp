/*******************************************************************************
 * $Id: MatchHashEntry.cpp,v 1.9 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing, copying, and modification
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/MatchHashEntry.h"
#include "libGenome/gnException.h"
#include "libGenome/gnDebug.h"

using namespace std;
using namespace genome;
namespace mems {

boolean MatchHashEntry::offset_lessthan(const MatchHashEntry& a, const MatchHashEntry& b){
	return a.m_offset < b.m_offset;
}

boolean MatchHashEntry::start_lessthan_ptr(const MatchHashEntry* a, const MatchHashEntry* b){
	int32 start_diff = a->FirstStart() - b->FirstStart();
	if(start_diff == 0){
		uint32 m_count = a->SeqCount();
		m_count = m_count <= b->SeqCount() ? m_count : b->SeqCount();
		for(uint32 seqI = seq_compare_start; seqI < m_count; seqI++){
			int64 a_start = a->Start(seqI), b_start = b->Start(seqI);
			if(a_start < 0)
				a_start = -a_start + a->Length() - a->m_mersize;
			if(b_start < 0)
				b_start = -b_start + b->Length() - b->m_mersize;
			int64 diff = a_start - b_start;
			if(a_start == NO_MATCH || b_start == NO_MATCH)
				continue;
			else if(diff == 0)
				continue;
			else
				return diff < 0;
		}
	}
	return start_diff < 0;
}

boolean MatchHashEntry::strict_start_lessthan_ptr(const MatchHashEntry* a, const MatchHashEntry* b){
	int start_diff = a->FirstStart() - b->FirstStart();
	if(start_diff == 0){
		uint m_count = a->SeqCount();
		m_count = m_count <= b->SeqCount() ? m_count : b->SeqCount();
		for(uint seqI = 0; seqI < m_count; seqI++){
			int64 a_start = a->Start(seqI), b_start = b->Start(seqI);
			if(a_start < 0)
				a_start = -a_start + a->Length() - a->m_mersize;
			if(b_start < 0)
				b_start = -b_start + b->Length() - b->m_mersize;
			int64 diff = a_start - b_start;
			if(diff == 0)
				continue;
			else
				return diff < 0;
		}
	}
	return start_diff < 0;
}


//ignores mem_no_matches
int64 MatchHashEntry::start_compare(const MatchHashEntry& a, const MatchHashEntry& b){
	uint m_count = a.SeqCount();
	m_count = m_count <= b.SeqCount() ? m_count : b.SeqCount();
	for(uint seqI = 0; seqI < m_count; seqI++){
		int64 a_start = a.Start(seqI), b_start = b.Start(seqI);
		if(a_start < 0)
			a_start = -a_start + a.Length() - a.m_mersize;
		if(b_start < 0)
			b_start = -b_start + b.Length() - b.m_mersize;
		int64 diff = a_start - b_start;
		if(a_start == NO_MATCH || b_start == NO_MATCH)
			continue;
		else if(diff == 0)
			continue;
		else
			return diff;
	}
	return 0;
}

int64 MatchHashEntry::end_to_start_compare(const MatchHashEntry& a, const MatchHashEntry& b){
	MatchHashEntry tmp_a = a;
	tmp_a.CropStart(tmp_a.Length()-1);
	return MatchHashEntry::start_compare(tmp_a, b);
}


MatchHashEntry::MatchHashEntry() : 
Match(),
m_extended( false ),
m_mersize( 0 )
{
}


MatchHashEntry::MatchHashEntry(uint32 seq_count, const gnSeqI mersize, MemType m_type) : 
 Match( seq_count ),
 m_mersize( mersize )
{
	m_extended = m_type == extended;
}


MatchHashEntry* MatchHashEntry::Clone() const{
	return new MatchHashEntry(*this);
}

MatchHashEntry& MatchHashEntry::operator=(const MatchHashEntry& mhe)
{
	Match::operator=( mhe );
	m_extended = mhe.m_extended;
	m_mersize = 0;
	m_offset = mhe.m_offset;

	return *this;
}

boolean MatchHashEntry::operator==(const MatchHashEntry& mhe) const
{
	if(m_seq_count != mhe.m_seq_count)
		return false;
	if(m_mersize != mhe.m_mersize)
		return false;
	if(m_extended != mhe.m_extended)
		return false;
	if( m_offset != mhe.m_offset )
		return false;
	return Match::operator ==(mhe);
}

void MatchHashEntry::CalculateOffset()
{
	if( SeqCount() == 0 )
		return;

	int64 tmp_off = 0;
	m_offset = 0;

	uint seqI = FirstStart();
	int64 ref_start = Start(seqI);

	for(seqI++; seqI < SeqCount(); seqI++){
		if(Start(seqI) != NO_MATCH){
			tmp_off = Start(seqI) - ref_start;
			if( Start(seqI) < 0 )
				tmp_off -= (int64)Length( seqI );
			m_offset += tmp_off;
		}
	}
}

// checks if mhe is _perfectly_ contained in this match.
// all offsets in all sequences must be aligned to each other
boolean MatchHashEntry::Contains(const MatchHashEntry& mhe) const{
	uint i;
	int64 diff_i;
	int64 diff;
	uint seq_count = mhe.SeqCount();
	//check for a consistent number of genomes and
	//identical generalized offsets
	if(SeqCount() != seq_count || m_offset != mhe.m_offset)
		return false;

	i = mhe.FirstStart();
	diff = mhe.Start(i) - Start(i);
	if(Start(i) == NO_MATCH)
		return false;

	//check for containment properties
	if(diff < 0 || Length() < mhe.Length() + diff)
		return false;

	//everything is ok so far, check for alignment
	int64 diff_rc = (int64)mhe.Length() - (int64)Length() + diff;
	for(i++; i < seq_count; i++){
		//check for consistent alignment between all genomes
		//in the case of revcomp, diff_i must equal diff_rc
		diff_i = mhe.Start(i) - Start(i);

		//it's ok if neither matches in a sequence
		if(mhe.Start(i) == NO_MATCH && Start(i) == NO_MATCH)
			continue;
		else if(mhe.Start(i) < 0 && diff_rc == diff_i)
			continue;
		else if(diff != diff_i )
			return false;
	}
	//it was contained.
	return true;
}


} // namespace mems
