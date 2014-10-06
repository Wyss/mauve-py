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

#ifndef _gnFileContig_h_
#define _gnFileContig_h_

#include "libGenome/gnDefs.h"

#include <string>
#include "libGenome/gnClone.h"
#include <utility>


namespace genome {

/**
 * gnFileContig is used by source classes to track the location of
 * sequence data on disk.  gnFileContig stores the start and end byte
 * offset, and the size of a repeated gap in the sequence data.
 * Also stores whether the sequence data is in the expected format or
 * if it is corrupted.
 */
class GNDLLEXPORT gnFileContig : public gnClone
{
public:
	gnFileContig();
	gnFileContig( std::string nameStr, const uint64 pos, const uint64 len );
	gnFileContig( const gnFileContig& fc );
	~gnFileContig();

	gnFileContig* Clone() const;
	void Clear();

	std::string GetName() const;
	gnSeqI GetSeqLength() const;
	std::pair<uint64,uint64> GetFileStartEnd() const;
	uint64 GetFileLength() const;
	std::pair<uint64,uint64> GetSectStartEnd( const gnContigSection i ) const;
	uint64 GetSectLength( gnContigSection i ) const;
	boolean HasRepeatSeqGap() const;
	std::pair<uint64,uint64> GetRepeatSeqGapSize() const;

	boolean SetName( std::string nameStr );
	boolean SetSeqLength( const gnSeqI len );
	boolean AddToSeqLength( const gnSeqI len );
	boolean SetFileStart( const uint64 s );
	boolean SetFileEnd( const uint64 e );
	boolean SetFileStartEnd( const std::pair<uint64,uint64> se );
	boolean SetSectStart( const gnContigSection i, const uint64 s );
	boolean SetSectEnd( const gnContigSection i, const uint64 e );
	boolean SetSectStartEnd( const gnContigSection i, const std::pair<uint64,uint64> se);
	boolean SetRepeatSeqGap( const boolean rsg );
	boolean SetRepeatSeqGapSize( const std::pair<uint64,uint64> rsgSize );
	boolean SetRepeatSeqSize( const uint64 seqSize );
	boolean SetRepeatGapSize( const uint64 gapSize );
private:
	std::string m_name;
	gnSeqI m_seqLength;
	std::pair<uint64,uint64> m_fileStartEnd;
	
	std::pair<uint64,uint64> m_startEndArray[CONTIG_SECTION_SIZE];
	// sequence access
	boolean m_repeatSeqGap;  // if true, use m_repeatSeqGapSize
	std::pair< uint64, uint64 > m_repeatSeqGapSize;
};// class gnFileContig

  // Clone
inline
gnFileContig* gnFileContig::Clone() const
{
	return new gnFileContig( *this );
}
  // GET
inline
std::string gnFileContig::GetName() const
{
	return m_name;
}
inline
gnSeqI gnFileContig::GetSeqLength() const
{
	return m_seqLength;
}
inline
std::pair<uint64,uint64> gnFileContig::GetFileStartEnd() const
{
	return m_fileStartEnd;
}
inline
uint64 gnFileContig::GetFileLength() const
{
	return m_fileStartEnd.second - m_fileStartEnd.first + 1;
}
inline
std::pair<uint64,uint64> gnFileContig::GetSectStartEnd( const gnContigSection i ) const
{
	if( (uint32)i < CONTIG_SECTION_SIZE )
		return m_startEndArray[(uint32)i];
	return std::pair<uint64,uint64>(0,0);
}
inline
uint64 gnFileContig::GetSectLength( gnContigSection i ) const
{
	if( (uint32)i < CONTIG_SECTION_SIZE )
		return m_startEndArray[(uint32)i].second - m_startEndArray[(uint32)i].first + 1;
	return 0;
}
inline
boolean gnFileContig::HasRepeatSeqGap() const
{
	return m_repeatSeqGap;
}
inline
std::pair<uint64,uint64> gnFileContig::GetRepeatSeqGapSize() const
{
	return m_repeatSeqGapSize;
}
  // SET
inline
boolean gnFileContig::SetName( std::string nameStr )
{
	m_name = nameStr;
	return true;
}
inline
boolean gnFileContig::SetSeqLength( const gnSeqI len )
{
	m_seqLength = len;
	return true;
}
inline
boolean gnFileContig::AddToSeqLength( const gnSeqI len )
{
	m_seqLength += len;
	return true;
}
inline
boolean gnFileContig::SetFileStart( const uint64 s )
{
	m_fileStartEnd.first = s;
	return true;
}
inline
boolean gnFileContig::SetFileEnd( const uint64 e )
{
	m_fileStartEnd.second = e;
	return true;
}
inline
boolean gnFileContig::SetFileStartEnd( const std::pair<uint64,uint64> se )
{
	m_fileStartEnd = se;
	return true;	
}
inline
boolean gnFileContig::SetSectStart( const gnContigSection i, const uint64 s )
{
	if( (uint32)i < CONTIG_SECTION_SIZE )
	{
		m_startEndArray[(uint32)i].first = s;
		return true;
	}
	return false;
}
inline
boolean gnFileContig::SetSectEnd( const gnContigSection i, const uint64 e )
{
	if( (uint32)i < CONTIG_SECTION_SIZE )
	{
		m_startEndArray[(uint32)i].second = e;
		return true;
	}
	return false;
}
inline
boolean gnFileContig::SetSectStartEnd( const gnContigSection i, const std::pair<uint64,uint64> se )
{
	if( (uint32)i < CONTIG_SECTION_SIZE )
	{
		m_startEndArray[(uint32)i] = se;
		return true;
	}
	return false;
}
inline
boolean gnFileContig::SetRepeatSeqGap( const boolean rsg )
{
	m_repeatSeqGap = rsg;
	return true;
}
inline
boolean gnFileContig::SetRepeatSeqGapSize( const std::pair<uint64,uint64> rsgSize )
{
	return  SetRepeatSeqSize( rsgSize.first ) && 
		SetRepeatGapSize( rsgSize.second );
}



}	// end namespace genome

#endif
	// _gnFileContig_h_
