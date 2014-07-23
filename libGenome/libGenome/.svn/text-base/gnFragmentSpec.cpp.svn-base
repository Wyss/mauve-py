/////////////////////////////////////////////////////////////////////////////
// File:            gnFragmentSpec.cpp
// Purpose:         implements gnMultiSpec< gnContigSpec > for sequence fragments
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

#include "libGenome/gnFragmentSpec.h"
#include <string>


using namespace std;
namespace genome {


gnFragmentSpec::gnFragmentSpec()
{
	gnBaseSpec::Clear();
}

gnFragmentSpec::gnFragmentSpec( const gnFragmentSpec& s )
{
	m_sourceName = s.m_sourceName;
	m_name = s.m_name;
	m_reverseComplement = s.m_reverseComplement;
	m_circular = s.m_circular;
	//copy the header list.
	uint32 list_size =  s.m_headerList.size();
	m_headerList.reserve(list_size);
	for(uint32 i=0; i < list_size; i++)
		m_headerList.push_back(s.m_headerList[i]->Clone());
	//copy the contig list.
	list_size =  s.m_SpecList.size();
	m_SpecList.reserve(list_size);
	uint32 i=0;
	for(; i < list_size; i++)
		m_SpecList.push_back(s.m_SpecList[i]->Clone());
	//copy the feature list
	list_size =  s.m_featureList.size();
	m_featureList.reserve(list_size);
	for(i=0; i < list_size; i++)
		m_featureList.push_back(s.m_featureList[i]->Clone());
}
gnFragmentSpec::~gnFragmentSpec()
{
	Clear();
}
// Clone

void gnFragmentSpec::Clear()
{
	uint32 list_size = m_SpecList.size();
	for(uint32 i=0; i < list_size; i++)
		delete m_SpecList[i];
	m_SpecList.clear();
	list_size = m_featureList.size();
	for(uint32 i=0; i < list_size; i++)
		delete m_featureList[i];
	m_featureList.clear();
	gnMultiSpec< gnContigSpec >::Clear();
}
//Fragment
void gnFragmentSpec::GetContainedFeatures(const gnLocation& lt, vector<gnBaseFeature*>& feature_vector, vector<uint32>& index_vector) const{
	for(uint32 i=0; i < m_featureList.size(); i++){
		if(m_featureList[i]->IsContainedBy(lt)){
			feature_vector.push_back(m_featureList[i]->Clone());
			index_vector.push_back(i);
		}
	}
}
void gnFragmentSpec::GetIntersectingFeatures(const gnLocation& lt, vector<gnBaseFeature*>& feature_vector, vector<uint32>& index_vector) const{
	for(uint32 i=0; i < m_featureList.size(); i++){
		if(m_featureList[i]->Intersects(lt)){
			feature_vector.push_back(m_featureList[i]->Clone());
			index_vector.push_back(i);
		}
	}
}
void gnFragmentSpec::GetBrokenFeatures(const gnLocation& lt, vector<gnBaseFeature*>& feature_vector) const{
	for(uint32 i=0; i < m_featureList.size(); i++)
		if(m_featureList[i]->IsBroken() && m_featureList[i]->IsContainedBy(lt))
			feature_vector.push_back(m_featureList[i]->Clone());
}

//do feature cropping... done
void gnFragmentSpec::CropStart( gnSeqI cropLen ){
	uint32 flen = m_featureList.size();
	for(uint32 featureI = 0; featureI < flen; featureI++){
		m_featureList[featureI]->CropStart(cropLen);
	}
	gnMultiSpec< gnContigSpec >::CropStart(cropLen);
}

void gnFragmentSpec::CropEnd( gnSeqI cropLen ){
	uint32 flen = m_featureList.size();
	for(uint32 featureI = 0; featureI < flen; featureI++){
		m_featureList[featureI]->CropEnd(cropLen);
	}
	gnMultiSpec< gnContigSpec >::CropEnd(cropLen);
}

gnFragmentSpec* gnFragmentSpec::CloneRange( const gnSeqI startI, const gnSeqI len ) const{
	if(len == 0)
		return new gnFragmentSpec();

	//find the valid range of specs to copy
	uint32 firstSpec = GetSpecIndexByBase(startI);
	gnSeqI total_copylen = len;
	uint32 endSpec;
	if(len != GNSEQI_END){
		endSpec = GetSpecIndexByBase(startI + len - 1);
	}else{
		endSpec = GetSpecListLength() - 1;
		total_copylen = GetLength() - startI;
	}

	//find their starting and ending bases
	gnSeqI firstBase = startI - GetSpecStartBase(firstSpec);
	gnSeqI firstSpecLen = GetSpec(firstSpec)->GetLength();
	boolean spans_specs = true;
	gnSeqI firstCopyLen = firstSpecLen - firstBase;
	if(firstCopyLen >= total_copylen){
		spans_specs = false;
		firstCopyLen = total_copylen;
	}
	gnFragmentSpec* destSpec = new gnFragmentSpec();
	gnContigSpec* newSpec = m_SpecList[firstSpec]->CloneRange(firstBase, firstCopyLen);
	destSpec->AddSpec( newSpec );

	gnSeqI cur_copylen = firstCopyLen;
	//add all the completely covered specs in the middle
	for(uint32 specI = firstSpec + 2; specI <= endSpec; specI++){
		destSpec->AddSpec(GetSpec(specI-1)->Clone());
		cur_copylen += GetSpec(specI-1)->GetLength();
	}
	//add the last spec if necessary
	if(spans_specs){
		newSpec = m_SpecList[endSpec]->CloneRange( 0, total_copylen - cur_copylen);
		destSpec->AddSpec(newSpec);
	}
	
	//now clone all the appropriate features
	gnLocation lt;
	vector<gnBaseFeature*> feature_vector;
	vector<uint32> index_vector;
	lt.SetStart(startI);
	lt.SetEnd(startI + total_copylen);
	GetIntersectingFeatures(lt, destSpec->m_featureList, index_vector);
	
	return destSpec;
}

/*IMPLEMENT ME!!  ADD FEATURE Manipulation*/
void gnFragmentSpec::SetReverseComplement( const boolean value )
{
	if(value == m_reverseComplement)
		return;
	//reverse the spec list entries
	vector<gnContigSpec*> tmp_spec_list;
	for(uint32 i=0; i < GetSpecListLength(); i++){
		//transmit rev_comp down the tree
		GetSpec(i)->SetReverseComplement(!GetSpec(i)->IsReverseComplement());
		tmp_spec_list.insert(tmp_spec_list.begin(), GetSpec(i));
	}
	m_SpecList = tmp_spec_list;
	m_reverseComplement = value;
}

}	// end namespace genome

