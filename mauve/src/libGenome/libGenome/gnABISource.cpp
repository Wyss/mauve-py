/////////////////////////////////////////////////////////////////////////////
// File:            gnABISource.h
// Purpose:         Implements gnBaseSource for .ABI files
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

#include "libGenome/gnABISource.h"

using namespace std;
namespace genome {


gnABISource::gnABISource()
{
	m_openString = "";
	m_pFilter = gnFilter::fullDNASeqFilter();
}
gnABISource::gnABISource( const gnABISource& s ) : gnFileSource(s)
{
	vector< gnFileContig* >::const_iterator iter = s.m_contigList.begin();
	for( ; iter != s.m_contigList.end(); ++iter )
	{
		m_contigList.push_back( *iter );
	}
}
gnABISource::~gnABISource()
{
	m_ifstream.close();
	vector< gnFileContig* >::iterator iter = m_contigList.begin();
	for( ; iter != m_contigList.end(); ++iter )
	{
		gnFileContig* fg = *iter;
		*iter = 0;
		delete fg;
	}
}

// Contig Access methods	
boolean gnABISource::HasContig( const string& name ) const
{
	vector< gnFileContig* >::const_iterator iter = m_contigList.begin();
	for( ; iter != m_contigList.end(); ++iter )
	{
		if( name == (*iter)->GetName() )
			return true;
	}
	return false;
}
uint32 gnABISource::GetContigID( const string& name ) const
{
	vector< gnFileContig* >::const_iterator iter = m_contigList.begin();
	for( ; iter != m_contigList.end(); ++iter )
	{
		if( name == (*iter)->GetName() )
			return iter - m_contigList.begin();
	}
	return ALL_CONTIGS;
}

string gnABISource::GetContigName( uint32 i ) const{
	if( i < m_contigList.size() ){
		return m_contigList[i]->GetName();
	}
	return "";
}

gnSeqI gnABISource::GetContigSeqLength( uint32 i ) const{
	if( i < m_contigList.size() )
	{
		return m_contigList[i]->GetSeqLength();
	}else if( i == ALL_CONTIGS){
		gnSeqI seqlen = 0;
		for(int j=0; j < m_contigList.size(); j++)
			seqlen += m_contigList[j]->GetSeqLength();
		return seqlen;
	}
	return GNSEQI_ERROR;
}

boolean gnABISource::SeqRead( const gnSeqI start, char* buf, gnSeqI& bufLen, const uint32 contigI )
{
	return false;
}

boolean gnABISource::SeqSeek( const gnSeqI start, const uint32& contigI, uint64& startPos, uint64& readableBytes ){
	return SeqStartPos( start, *(m_contigList[contigI]), startPos, readableBytes );
}

boolean gnABISource::SeqStartPos( const gnSeqI start, gnFileContig& contig, uint64& startPos, uint64& readableBytes )
{
	return false;
}

gnFileContig* gnABISource::GetFileContig( const uint32 contigI ) const{
	if(m_contigList.size() > contigI)
		return m_contigList[contigI];
	return NULL;
}
 
boolean gnABISource::ParseStream( istream& fin )
{
	return false;
}


}	// end namespace genome

