/////////////////////////////////////////////////////////////////////////////
// File:            gnGenomeSpec.cpp
// Purpose:         implements gnMultiSpec< gnFragmentSpec > for genomes
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

#include "libGenome/gnGenomeSpec.h"
#include <string>


using namespace std;
namespace genome {


gnGenomeSpec::gnGenomeSpec()
{
	gnBaseSpec::Clear();
}

gnGenomeSpec::gnGenomeSpec( const gnGenomeSpec& s )
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
	//copy the fragment list.
	list_size =  s.m_SpecList.size();
	m_SpecList.reserve(list_size);
	for(uint32 i=0; i < list_size; i++)
		m_SpecList.push_back(s.m_SpecList[i]->Clone());
}
gnGenomeSpec::~gnGenomeSpec()
{
	Clear();
}
// Clone
void gnGenomeSpec::Clear()
{
	uint32 list_size = m_SpecList.size();
	for(uint32 i=0; i < list_size; i++)
		delete m_SpecList[i];
	m_SpecList.clear();
	gnMultiSpec< gnFragmentSpec >::Clear();
}

void gnGenomeSpec::SetReverseComplement( const boolean value )
{
	if(value == m_reverseComplement)
		return;
	//reverse the spec list entries
	vector<gnFragmentSpec*> tmp_spec_list;
	for(uint32 i=0; i < GetSpecListLength(); i++){
		//transmit rev_comp down to each fragment spec
		GetSpec(i)->SetReverseComplement(!GetSpec(i)->IsReverseComplement());
		tmp_spec_list.insert(tmp_spec_list.begin(), GetSpec(i));
	}
	m_SpecList = tmp_spec_list;
	m_reverseComplement = value;
}

gnGenomeSpec* gnGenomeSpec::CloneRange( const gnSeqI startI, const gnSeqI len ) const{
	if(len == 0)
		return new gnGenomeSpec();

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
	
	gnGenomeSpec* destSpec = new gnGenomeSpec();
	gnFragmentSpec* newSpec = m_SpecList[firstSpec]->CloneRange(firstBase, firstCopyLen);
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
	return destSpec;
}


//IMPLEMENT ME: What to do with headers?
void gnGenomeSpec::MergeFragments(const uint32 startF, const uint32 endF){
	if(startF > m_SpecList.size() || endF > m_SpecList.size())
		Throw_gnEx(FragmentIndexOutOfBounds());
	if(startF > endF)
		Throw_gnEx(FragmentIndexOutOfBounds());
	
	for(uint32 i = startF + 1; i < endF; i++){
		gnFragmentSpec* delspec = m_SpecList[startF + 1];
		m_SpecList.erase(m_SpecList.begin() + startF + 1);

		for(uint32 j = 0; j < delspec->GetSpecListLength(); j++)
			m_SpecList[startF]->AddSpec(delspec->GetSpec(j));
		delete delspec;
	}
}

uint32 gnGenomeSpec::AddFeature( gnBaseFeature* feat ){
	uint32 count = 0;
	uint32 len = 0;
	uint32 featureI = 0;
	uint32 specListLen = GetSpecListLength();

	for(uint32 specI = 0; specI < specListLen; specI++){
		len = GetSpec(specI)->GetLength();
		gnLocation lt(count, count+len);
		if(feat->IsContainedBy(lt)){
			return featureI + GetSpec(specI)->AddFeature(feat);
		}
		count += len;
		featureI += GetSpec(specI)->GetFeatureListLength();
	}
	//if we get this far then the feature has invalid coordinates
	Throw_gnEx(SeqIndexOutOfBounds());
}

uint32 gnGenomeSpec::GetFeatureListLength() const
{
	uint32 len = 0;
	for(uint32 i=0; i < GetSpecListLength(); i++)
		len += GetSpec(i)->GetFeatureListLength();
	return len;
}

gnBaseFeature* gnGenomeSpec::GetFeature(const uint32 i ) const
{
	uint32 count = 0;
	uint32 len = 0;
	for(uint32 specI=0; specI < GetSpecListLength(); specI++){
		len = GetSpec(specI)->GetFeatureListLength();
		if(count <= i && i < count + len){
			gnBaseFeature* feat = GetSpec(specI)->GetFeature(i - count);
			feat->MovePositive(GetSpecStartBase(specI));
			return feat;
		}
		count += len;
	}
	Throw_gnEx(FeatureIndexOutOfBounds());
}

void gnGenomeSpec::GetContainedFeatures(const gnLocation& lt, vector<gnBaseFeature*>& feature_vector, vector<uint32>& index_vector) const{
	uint32 ss_size = GetSpecListLength();
	uint32 fl_size = 0;
	gnSeqI start_base = 0;
	for(uint32 i=0; i < ss_size; i++){
		gnLocation sub_lt = lt;
		gnSeqI sub_len = GetSpec(i)->GetLength();
		sub_lt.MoveNegative(start_base);
		sub_lt.CropEnd(sub_len);
		GetSpec(i)->GetContainedFeatures(sub_lt, feature_vector, index_vector);
		uint32 fvs = feature_vector.size();
		for(uint32 j = 0; j < fvs; j++){
			feature_vector[j]->MovePositive(start_base);
			index_vector[j]+= fl_size;
		}
		if(fvs > 0)
			return;
		start_base += sub_len;
		fl_size += GetSpec(i)->GetFeatureListLength();
	}
}

void gnGenomeSpec::GetIntersectingFeatures(const gnLocation& lt, vector<gnBaseFeature*>& feature_vector, vector<uint32>& index_vector) const{
	uint32 ss_size = GetSpecListLength();
	uint32 fl_size = 0;
	gnSeqI start_base = 0;
	for(uint32 i=0; i < ss_size; i++){
		gnLocation sub_lt = lt;
		gnSeqI sub_len = GetSpec(i)->GetLength();
		sub_lt.MoveNegative(start_base);
		sub_lt.CropEnd(sub_len);
		GetSpec(i)->GetIntersectingFeatures(sub_lt, feature_vector, index_vector);
		uint32 fvs = feature_vector.size();
		for(uint32 j = 0; j < fvs; j++){
			feature_vector[j]->MovePositive(start_base);
			index_vector[j]+= fl_size;
		}
		if(fvs > 0)
			return;
		start_base += sub_len;
		fl_size += GetSpec(i)->GetFeatureListLength();
	}
}
void gnGenomeSpec::GetBrokenFeatures(const gnLocation& lt, vector<gnBaseFeature*>& feature_vector) const{
	uint32 ss_size = GetSpecListLength();
	uint32 fl_size = 0;
	gnSeqI start_base = 0;
	for(uint32 i=0; i < ss_size; i++){
		gnLocation sub_lt = lt;
		gnSeqI sub_len = GetSpec(i)->GetLength();
		sub_lt.MoveNegative(start_base);
		sub_lt.CropEnd(sub_len);
		GetSpec(i)->GetBrokenFeatures(sub_lt, feature_vector);
		uint32 fvs = feature_vector.size();
		for(uint32 j = 0; j < fvs; j++)
			feature_vector[j]->MovePositive(start_base);
		if(fvs > 0)
			return;
		start_base += sub_len;
		fl_size += GetSpec(i)->GetFeatureListLength();
	}
}

void gnGenomeSpec::RemoveFeature( const uint32 i ){
	uint32 count = 0;
	uint32 len = 0;
	uint32 specI=0;
	for(; specI < GetSpecListLength(); specI++){
		len = GetSpec(specI)->GetFeatureListLength();
		if(count <= i && i < count + len){
			GetSpec(specI)->RemoveFeature( i - count);
		}
		count += len;
	}
	Throw_gnEx(FeatureIndexOutOfBounds());
}

}	// end namespace genome

