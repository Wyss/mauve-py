/////////////////////////////////////////////////////////////////////////////
// File:            gnLocation.cpp
// Purpose:         Standard Location for Feature
// Description:     Feature location
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
#include "libGenome/gnLocation.h"
#include "libGenome/gnDebug.h"


using namespace std;
namespace genome {


gnLocation::gnLocation()
{
	Clear();
}
gnLocation::gnLocation( const gnLocation& s)
{
	SetBounds( s.m_start, s.m_startLength, s.m_end, s.m_endLength );
	m_type = s.m_type;
}
gnLocation::gnLocation( const gnSeqI start, const gnSeqI startLength, const gnSeqI end, const gnSeqI endLength, gnLocationType type, string contigName)
{
	SetBounds( start, startLength, end, endLength );
	m_type = type;
	m_name = contigName;
}
gnLocation::gnLocation( const gnSeqI start, const gnSeqI end, const gnLocationType type, string contigName )
{
	SetBounds( start, 0, end, 0 );
	m_type = type;
	m_name = contigName;
}

gnLocation* gnLocation::Clone() const
{
	return new gnLocation(*this);
}

void gnLocation::Clear()
{
	m_start = 0;
	m_end = 0;
	m_startLength = 0;
	m_endLength = 0;
	m_type = LT_Nothing;
}


void gnLocation::GetBounds( gnSeqI &s, gnSeqI &sl, gnSeqI &e, gnSeqI &el ) const
{
	s = m_start;
	sl = m_startLength;
	e = m_end;
	el = m_endLength;
}

void gnLocation::SetBounds( const gnSeqI start, const gnSeqI startLength, const gnSeqI end, const gnSeqI endLength)
{
	SetStart(start, startLength);
	SetEnd(end, endLength);
}
void gnLocation::SetBounds( const gnSeqI start, const gnSeqI end)
{
	m_start = start;
	m_end = end;
}

boolean gnLocation::CropTo( const gnLocation &l )
{
	gnSeqI tmp;
	gnSeqI start = l.GetStart();
	gnSeqI end = l.GetEnd();
	if(m_start < start){
		tmp = start < m_end ? start : m_end;
		m_startLength += tmp - m_start;
		m_start = tmp;
	}
	if(m_end < end){
		tmp = end > m_start ? end : m_start;
		m_endLength += m_end - tmp;
		m_end = tmp;
	}

	if( (l.GetFirst() > GetFirst()) )
	{
		if( l.GetFirst() <= m_end )
			m_startLength = m_start - l.GetFirst();
		else if( l.GetFirst() <= GetLast() )
		{
			m_end = l.GetFirst();
			m_start = l.GetFirst() + 1;
			m_startLength = 0;
		}	
		else
			Clear();
	}
	if( l.GetLast() < GetLast() )
	{
		if( l.GetLast() >= m_start )
			m_endLength = l.GetLast() - m_end;
		else if( l.GetLast() >= GetFirst() )
		{
			m_start = l.GetLast();
			m_end = l.GetLast() - 1;
			m_endLength = 0;
		}	
		else
			Clear();
	}
	if(m_start == m_end)
		return false;
	return true;
}

boolean gnLocation::CropStart( const gnSeqI start ){
	gnSeqI tmp;
	if(m_start < start){
		tmp = start < m_end ? start : m_end;
		m_startLength += tmp - m_start;
		m_start = tmp;
	}
	if(m_start == m_end)
		return false;
	return true;
}

boolean gnLocation::CropEnd( const gnSeqI end ){
	gnSeqI tmp;
	if(m_end > end){
		tmp = end > m_start ? end : m_start;
		m_endLength += m_end - tmp;
		m_end = tmp;
	}
	if(m_start == m_end)
		return false;
	return true;
}

// Intersects
boolean gnLocation::Intersects( const gnLocation &l, const intersectRestriction ir ) const{
	
	if( ir == determinedRegions )
	{
		if( (l.GetStart() <= m_end) && (l.GetEnd() >= m_start) )
			return true;
	}
	else if( ir == undeterminedRegions )
	{
		if( (l.GetFirst() <= m_start) && (l.GetLast() >= GetFirst()) )
			return true;
		if( (l.GetFirst() <= GetLast()) && (l.GetLast() >= m_end) )
			return true;
	}
	else if( ir == allRegions )
	{
		if( (l.GetFirst() <= GetLast()) && (l.GetLast() >= GetFirst()) )
			return true;
	}
	return false;
}

boolean gnLocation::Contains( const gnLocation &l, const intersectRestriction ir ) const{
	if(ir == determinedRegions)
		return m_start <= l.GetStart() && l.GetEnd() <= m_end;
	else if(ir == undeterminedRegions)
		return  (GetFirst() <= l.GetFirst() && l.GetLast() < m_start) ||
				(m_end < l.GetFirst() && l.GetLast() <= GetLast());
	return GetFirst() <= l.GetFirst() && l.GetLast() <= GetLast();
}


// Move
boolean gnLocation::MovePositive( const gnSeqI diff )
{
	if(m_start > GNSEQI_END - diff || m_end > GNSEQI_END - diff)
		return false;
	m_start += diff;
	m_end += diff;
	return true;
}
boolean gnLocation::MoveNegative( const gnSeqI diff )
{
	if(m_start < diff || m_end < diff)
		return false;
	m_start -= diff;
	m_end -= diff;
	return true;
}
boolean gnLocation::MoveTo( const int direction, const gnSeqI diff )
{
	if( direction > 0 )
		return MovePositive( diff );
	return MoveNegative( diff );
}

gnLocation gnLocation::GetUnion( const gnLocation &l ) const
{
	ErrorMsg("gnLocation::getUnion -- not implemented\n");
	return l;
}

// Intersects
gnLocation gnLocation::GetIntersection( const gnLocation &l, const intersectRestriction ir ) const{
	gnLocation inter_loc;
	if( ir == determinedRegions )
	{
		if( (l.GetStart() <= m_end) && (l.GetEnd() >= m_start) ){
			inter_loc.m_start = l.m_start > m_start ? l.m_start : m_start;
			inter_loc.m_end = l.m_end < m_end ? l.m_end : m_end;
		}
	}
	else if( ir == undeterminedRegions )
	{
		ErrorMsg("Not implemented!");
	}
	else if( ir == allRegions )
	{
		ErrorMsg("Not implemented!");
	}
	return inter_loc;
}

}	// end namespace genome

