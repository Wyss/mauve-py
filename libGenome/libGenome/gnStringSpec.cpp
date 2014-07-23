/////////////////////////////////////////////////////////////////////////////
// File:            gnStringSpec.cpp
// Purpose:         implements gnContigSpec for strings
// Description:     stores sequence in memory
// Changes:        
// Version:         libGenome 0.5.1 
// Author:          Aaron Darling 
// Modified by:     
// Copyright:       (c) Aaron Darling 
// Licenses:        See COPYING file for details
/////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include "libGenome/gnStringSpec.h"
#include <string>
#include <cstring>


using namespace std;
namespace genome {


gnStringSpec::gnStringSpec()
{
	gnContigSpec::Clear();
}

gnStringSpec::gnStringSpec( const string& m_string, gnSeqI start, gnSeqI endI, boolean revComp)
{
	m_seqString = m_string;
	m_start = start;
	gnSeqI actual_len = m_seqString.length();

	//reverse comp has the end b.p. first, then start.  switch them.
	m_start = revComp ? endI : start;
	gnSeqI actual_end = revComp ? start : endI;
	//trim start and end down if they are too big.
	actual_end = actual_end < actual_len ? actual_end : actual_len -1;
	m_start = m_start < actual_len ? m_start : actual_len - 1;
	if(actual_len == 0)
		m_start = 0;
	//if start is after end, we're using a circular sequence.
	m_circular = m_start > actual_end ? true : false;
	//if circular, the length will be different.
	m_length = m_circular ? (actual_len - m_start) + actual_end : actual_end - m_start + 1;

	m_reverseComplement = revComp;
	m_SourceContigIndex = ALL_CONTIGS;
}

gnStringSpec::gnStringSpec( const gnStringSpec& s )
{
	m_seqString = s.m_seqString;
	m_sourceName = s.m_sourceName;
	m_name = s.m_name;
	m_start = s.m_start;
	m_length = s.m_length;
	m_reverseComplement = s.m_reverseComplement;
	m_circular = s.m_circular;
	m_SourceContigIndex = s.m_SourceContigIndex;
}
gnStringSpec::~gnStringSpec()
{
	Clear();
}
void gnStringSpec::Clear()
{
	gnContigSpec::Clear();
	m_seqString = "";
}

gnStringSpec* gnStringSpec::CloneRange( const gnSeqI startI, const gnSeqI len ) const{
	gnStringSpec* destSpec = new gnStringSpec();
	destSpec->m_seqString = m_seqString.substr(m_start + startI, len);
	destSpec->m_sourceName = m_sourceName;
	destSpec->m_name = m_name;
	destSpec->m_start = 0;
	destSpec->m_length = destSpec->m_seqString.length();
	destSpec->m_reverseComplement = m_reverseComplement;
	destSpec->m_circular = m_circular;
	destSpec->m_SourceContigIndex = m_SourceContigIndex;
	return destSpec;
}

}	// end namespace genome

