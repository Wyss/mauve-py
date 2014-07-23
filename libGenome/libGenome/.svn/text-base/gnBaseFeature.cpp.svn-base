/////////////////////////////////////////////////////////////////////////////
// File:            gnBaseFeature.cpp
// Purpose:         implements the gnBaseFeature
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

#include "libGenome/gnBaseFeature.h"
#include "libGenome/gnLocation.h"
#include "libGenome/gnStringQualifier.h"
#include "libGenome/gnDebug.h"
#include <list>


using namespace std;
namespace genome {


gnBaseFeature::gnBaseFeature( )
{
	m_id = 0;
	m_name = "";
	m_locationType = gnLocation::LT_Nothing;
	m_broken = false;
	m_spec = NULL;
}
gnBaseFeature::gnBaseFeature( string& name, uint32 id, gnFragmentSpec* spec, gnLocation::gnLocationType lt, boolean broken )
{
	m_id = id;
	m_name = name;
	m_locationType = lt;
	m_broken = broken;
	m_spec = NULL;
}
gnBaseFeature::~gnBaseFeature()
{
	for( uint32 i=0; i < m_qualifierList.size(); i++)
		delete m_qualifierList[i];
}
// Clone
// Location Modification methods
boolean gnBaseFeature::MovePositive( const gnSeqI i )
{
	boolean still_valid = true;
	for(uint32 locationI=0; locationI< m_locationList.size(); locationI++){
		still_valid = still_valid && m_locationList[locationI].MovePositive(i);
	}
	return still_valid;
}
boolean gnBaseFeature::MoveNegative( const gnSeqI i )
{
	boolean still_valid = true;
	for(uint32 locationI=0; locationI< m_locationList.size(); locationI++){
		still_valid = still_valid && m_locationList[locationI].MoveNegative(i);
	}
	return still_valid;
}
//delete everything bigger than i
boolean gnBaseFeature::CropEnd( const gnSeqI i )
{
	boolean still_valid = true;
	for(uint32 locationI=0; locationI< m_locationList.size(); locationI++){
		still_valid = still_valid && m_locationList[locationI].CropEnd(i);
	}
	return still_valid;
}
//delete everything less than i
boolean gnBaseFeature::CropStart( const gnSeqI i ){
	boolean still_valid = true;
	for(uint32 locationI=0; locationI< m_locationList.size(); locationI++){
		still_valid = still_valid && m_locationList[locationI].CropStart(i);
	}
	return still_valid;
}
uint32 gnBaseFeature::FirstIndexOfQualifier( const string& name, uint32 listI ) const{
	if(listI >= m_qualifierList.size())
		return GNSEQI_END;
	
	uint32 i=listI;
	for(; i < m_qualifierList.size(); i++){
		string qname = m_qualifierList[i]->GetName();
		if(qname == name)
			break;
	}
	return i;
}
uint32 gnBaseFeature::LastIndexOfQualifier( const string& name, uint32 listI ) const{
	if(listI >= m_qualifierList.size())
		return GNSEQI_END;
	
	uint32 i=m_qualifierList.size()-1;
	for(; i >= listI; i--)
		if(m_qualifierList[i]->GetName() == name)
			break;
	return i;
}
boolean gnBaseFeature::RemoveQualifier( uint32 listI ){
	if(listI < m_qualifierList.size()){
		delete m_qualifierList[ listI ];
		m_qualifierList.erase(m_qualifierList.begin() + listI);
		return true;
	}
	return false;
}
boolean gnBaseFeature::SetQualifier( string& name, string& value, uint32 listI ){
	if(listI < m_qualifierList.size()){
		delete m_qualifierList[ listI ];
		m_qualifierList[listI] = new gnStringQualifier(name, value);
		return true;
	}
	return false;
}
boolean gnBaseFeature::SetQualifierName( string& name, uint32 listI ){
	if(listI < m_qualifierList.size()){
		gnStringQualifier* new_qual = new gnStringQualifier(name, m_qualifierList[listI]->GetValue());
		delete m_qualifierList[listI];
		m_qualifierList[listI] = new_qual;
		return true;
	}
	return false;
}
boolean gnBaseFeature::SetQualifierValue( string& value, uint32 listI ){
	if(listI < m_qualifierList.size()){
		gnStringQualifier* new_qual = new gnStringQualifier(m_qualifierList[listI]->GetName(), value);
		delete m_qualifierList[listI];
		m_qualifierList[listI] = new_qual;
		return true;
	}
	return false;
}
// Compare methods
boolean gnBaseFeature::Contains( gnSeqI i ) const{
	for(uint32 locationI=0; locationI< m_locationList.size(); locationI++){
		if((m_locationList[locationI].GetStart() <= i) &&
			(m_locationList[locationI].GetEnd() >= i))
			return true;
	}
	return false;
}
boolean gnBaseFeature::Contains( const gnLocation& l ) const{
	for(uint32 locationI=0; locationI< m_locationList.size(); locationI++){
		if(m_locationList[locationI].Contains(l))
			return true;
	}
	return false;
}
boolean gnBaseFeature::Contains( gnBaseFeature* feature ) const{
	for(uint32 locationI=0; locationI< m_locationList.size(); locationI++){
		uint32 i=0;
		for(; i < feature->GetLocationListLength(); i++){
			if(m_locationList[locationI].Contains(feature->GetLocation(i)))
				break;	//it's contained so we'll break the loop early
		}
		if(i == feature->GetLocationListLength())
			return false;	//the loop went full cycle.  we found an uncontained location
	}
	return true;
}
boolean gnBaseFeature::IsContainedBy( const gnLocation& l ) const{
	for(uint32 locationI=0; locationI< m_locationList.size(); locationI++){
		if(!l.Contains(m_locationList[locationI]))
			return false;
	}
	return true;
}

boolean gnBaseFeature::Intersects( const gnLocation& l ) const{
	for(uint32 locationI=0; locationI< m_locationList.size(); locationI++){
		if(!l.Intersects(m_locationList[locationI]))
			return false;
	}
	return true;
}
boolean gnBaseFeature::Intersects( gnBaseFeature* feature ) const{
	for(uint32 locationI=0; locationI< m_locationList.size(); locationI++){
		uint32 i=0;
		for(; i < feature->GetLocationListLength(); i++){
			if(m_locationList[locationI].Intersects(feature->GetLocation(i)))
				break;	//it's contained so we'll break the loop early
		}
		if(i == feature->GetLocationListLength())
			return false;	//the loop went full cycle.  we found an uncontained location
	}
	return true;
}
inline
boolean gnBaseFeature::AddLocation( const gnLocation& l, uint32 listI ){
	if(listI <= m_locationList.size()){
		m_locationList.insert(m_locationList.begin() + listI, l);
		return true;
	}
	return false;
}

gnLocation gnBaseFeature::GetLocation( uint32 listI ) const{
	if(listI < m_locationList.size())
		return m_locationList[listI];
	return gnLocation();
}

boolean gnBaseFeature::RemoveLocation( uint32 listI ){
	if(listI < m_locationList.size()){
		m_locationList.erase(m_locationList.begin() + listI);
		return true;
	}
	return false;
}

boolean gnBaseFeature::SetLocation( const gnLocation& l, uint32 listI ){
	if(listI < m_locationList.size()){
		m_locationList[listI] = l;
		return true;
	}
	return false;
}

boolean gnBaseFeature::AddQualifier( gnBaseQualifier* qualifier ){
	if(qualifier != NULL){
		m_qualifierList.push_back(qualifier);
		return true;
	}
	return false;
}

boolean gnBaseFeature::HasQualifier( const string& name ) const{
	for(uint32 i=0; i < m_qualifierList.size(); i++)
		if(m_qualifierList[i]->GetName() == name)
			return true;
	return false;
}

string gnBaseFeature::GetQualifierName( uint32 listI ) const{
	if(listI < m_qualifierList.size())
		return m_qualifierList[listI]->GetName();
	return string();
}

string gnBaseFeature::GetQualifierValue( uint32 listI ) const{
	if(listI < m_qualifierList.size())
		return m_qualifierList[listI]->GetValue();
	return string();
}

gnBaseQualifier* gnBaseFeature::GetQualifier( uint32 listI ){
	if(listI < m_qualifierList.size())
		return m_qualifierList[listI]->Clone();
	return NULL;
}


}	// end namespace genome

