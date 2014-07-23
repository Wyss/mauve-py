/////////////////////////////////////////////////////////////////////////////
// File:            gnSourceSpec.cpp
// Purpose:         implements gnContigSpec for source specs
// Description:     
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


#include "libGenome/gnSourceSpec.h"


using namespace std;
namespace genome {


gnSourceSpec::gnSourceSpec()
{
	Clear();
}

gnSourceSpec::gnSourceSpec( const gnSourceSpec& s )
{
	m_pSource = s.m_pSource;
	m_sourceName = string(s.m_sourceName);
	m_name = string(s.m_name);
	m_SourceContigIndex = s.m_SourceContigIndex;
	m_start = s.m_start;
	m_length = s.m_length;
	m_reverseComplement = s.m_reverseComplement;
	m_circular = s.m_circular;
}

gnSourceSpec::gnSourceSpec( gnBaseSource* source, const uint32 m_ContigIndex, const gnSeqI start, const gnSeqI endI, const boolean revComp)
{
	m_pSource = source;
	m_SourceContigIndex = m_ContigIndex;
	m_name = "";
	m_reverseComplement = revComp;
	m_circular = false;
	m_start = start;
	
	gnSeqI actual_len = source->GetContigSeqLength(m_ContigIndex);
	gnSeqI actual_end = endI;
	if(actual_len == 0)	
		return;	//this is a bogus gnSourceSpec

	//trim start and end down if they are too big.
	m_start = m_start < actual_len ? m_start : actual_len - 1;
	actual_end = actual_end < actual_len ? actual_end : actual_len - 1;
	//set the circularity and length
	if(revComp){
		m_circular = m_start < actual_end ? true : false;
		m_length = ((m_start - actual_end + actual_len) % actual_len);
	}else{
		m_circular = m_start > actual_end ? true : false;
		m_length = ((actual_end - m_start + actual_len) % actual_len);
	}
	if(actual_len != 0)
		m_length++;

}

gnSourceSpec::~gnSourceSpec()
{
}

void gnSourceSpec::Clear()
{
	gnContigSpec::Clear();
	m_SourceContigIndex = 0;
	m_pSource = NULL;
}

gnSourceSpec* gnSourceSpec::CloneRange( const gnSeqI startI, const gnSeqI len ) const{
	gnSourceSpec* destSpec = new gnSourceSpec();
	destSpec->m_pSource = m_pSource;
	destSpec->m_sourceName = m_sourceName;
	destSpec->m_name = m_name;
	destSpec->m_SourceContigIndex = m_SourceContigIndex;
	destSpec->m_start = m_start + startI;
	destSpec->m_length = len < m_length - startI ? len : m_length - startI;
	destSpec->m_reverseComplement = m_reverseComplement;
	destSpec->m_circular = m_circular;
	return destSpec;
}

}	// end namespace genome

