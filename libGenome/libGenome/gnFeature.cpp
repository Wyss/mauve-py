/////////////////////////////////////////////////////////////////////////////
// File:            gnFeature.cpp
// Purpose:         implements the gnBaseFeature for generic features
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

#include "libGenome/gnFeature.h"
#include "libGenome/gnLocation.h"
#include "libGenome/gnStringQualifier.h"
#include "libGenome/gnDebug.h"
#include <list>


using namespace std;
namespace genome {


gnFeature::gnFeature()
{
/*	m_id = 0;
	m_name = "";
	m_locationType = gnLocation::LT_Nothing;
	m_broken = false;
*/}
gnFeature::gnFeature( string& name, uint32 id, gnLocation::gnLocationType lt, boolean broken )
 : gnBaseFeature(name, id, NULL, lt, broken)
{
}
gnFeature::gnFeature( const gnFeature& s )
{
	uint32 i;
	m_id = s.m_id;
	m_name = s.m_name;
	m_locationType = s.m_locationType;
	m_broken = s.m_broken;
	m_spec = s.m_spec;
	for( i=0; i < s.m_locationList.size(); i++)
		m_locationList.push_back(s.m_locationList[i]);
	for( i=0; i < s.m_qualifierList.size(); i++)
		m_qualifierList.push_back(s.m_qualifierList[i]->Clone());
}
gnFeature::~gnFeature()
{
}

}	// end namespace genome

