/////////////////////////////////////////////////////////////////////////////
// File:            gnFileContig.h
// Purpose:         File Position holder.
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

#include "libGenome/gnFileContig.h"
#include <iostream>


using namespace std;
namespace genome {


gnFileContig::gnFileContig()
{
	Clear();
}
gnFileContig::gnFileContig( string nameStr, const uint64 s, const uint64 e )
{
	Clear();
	m_name = nameStr;
	m_fileStartEnd.first = s;
	m_fileStartEnd.second = e;
}
gnFileContig::gnFileContig( const gnFileContig& fc )
{
	m_name = fc.m_name;
	m_seqLength = fc.m_seqLength;
	m_fileStartEnd = fc.m_fileStartEnd;
	for( uint32 i=0; i < CONTIG_SECTION_SIZE; ++i )
		m_startEndArray[i] = fc.m_startEndArray[i];
	m_repeatSeqGap = fc.m_repeatSeqGap;
	m_repeatSeqGapSize = fc.m_repeatSeqGapSize;
}
gnFileContig::~gnFileContig()
{
}
void gnFileContig::Clear()
{
	m_name = "";
	m_seqLength = 0;
	m_fileStartEnd = pair<uint64,uint64>(0,0);
	for( uint32 i=0; i < CONTIG_SECTION_SIZE; ++i )
		m_startEndArray[i] = pair<uint64,uint64>(0,0);
	m_repeatSeqGap = false;
	m_repeatSeqGapSize = pair<uint64,uint64>(0,0);
}

boolean gnFileContig::SetRepeatSeqSize( const uint64 seqSize )
{
	if( !m_repeatSeqGap )
		return false;
	if( m_repeatSeqGapSize.first == seqSize )
		return true;
	if( m_repeatSeqGapSize.first == 0 )
	{
		m_repeatSeqGapSize.first = seqSize;
		return true;
	}
	m_repeatSeqGap = false;
	return false;
}
boolean gnFileContig::SetRepeatGapSize( const uint64 gapSize )
{
	if( !m_repeatSeqGap )
		return false;
	if( m_repeatSeqGapSize.second == gapSize )
		return true;
	if( m_repeatSeqGapSize.second == 0 )
	{
		m_repeatSeqGapSize.second = gapSize;
		return true;
	}
	m_repeatSeqGap = false;
	return false;
}

}	// end namespace genome

