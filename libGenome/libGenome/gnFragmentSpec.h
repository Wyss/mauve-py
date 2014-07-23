/////////////////////////////////////////////////////////////////////////////
// File:            gnFragmentSpec.h
// Purpose:         abstract Spec class
// Description:     Genome level spec class
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

#ifndef _gnFragmentSpec_h_
#define _gnFragmentSpec_h_

#include "libGenome/gnDefs.h"

#include <vector>
#include <string>

#include "libGenome/gnClone.h"
#include "libGenome/gnBaseFeature.h"
#include "libGenome/gnBaseHeader.h"
#include "libGenome/gnContigSpec.h"
#include "libGenome/gnMultiSpec.h"
#include "libGenome/gnException.h"


namespace genome {

/**
 * gnFragmentSpec contains a list of specs which make up a sequence
 * fragment.  It also contains a list of features which relate to
 * the sequence fragment.  Finally it contains a list of sequence
 * related header data.  This class is usually created and filled by
 * a file reader class like gnGBKSource.
 */
class GNDLLEXPORT gnFragmentSpec : public gnMultiSpec< gnContigSpec >
{
public:
	gnFragmentSpec();
	/**
	 * Destructor, frees memory.
	 */
	virtual ~gnFragmentSpec();
	/**
	 * Copy constructor.
	 * @param s the gnGenomeSpec to copy.
	 */
	gnFragmentSpec( const gnFragmentSpec& s);
	virtual gnFragmentSpec* Clone() const;
	virtual void Clear();
// Base Spec stuff
	virtual void SetReverseComplement( const boolean value );

//Multispec stuff
/*	virtual uint32 GetSpecListLength() const;
	virtual gnContigSpec* GetSpec( const uint32 i ) const;
	virtual gnContigSpec* GetSpecByBase( const gnSeqI baseI ) const;
	virtual void AddSpec( gnBaseSpec* spec, const uint32 i = UINT32_MAX );
	virtual void RemoveSpec( uint32 i );
*/
	virtual void CropStart( gnSeqI cropLen );
	virtual void CropEnd( gnSeqI cropLen );
	virtual uint32 AddFeature( gnBaseFeature* feat );
	virtual uint32 GetFeatureListLength() const;
	virtual gnBaseFeature* GetFeature( const uint32 i ) const;
	virtual void GetContainedFeatures(const gnLocation& lt, std::vector<gnBaseFeature*>& feature_vector, std::vector<uint32>& index_vector) const;
	virtual void GetIntersectingFeatures(const gnLocation& lt, std::vector<gnBaseFeature*>& feature_vector, std::vector<uint32>& index_vector) const;
	virtual void GetBrokenFeatures(const gnLocation& lt, std::vector<gnBaseFeature*>& feature_vector) const;
	virtual void RemoveFeature( const uint32 i );

	/**
	 * Copies a specified range of bases and returns a pointer to
	 * the resulting gnFragmentSpec.  You must delete the copy when you
	 * are finished with it.
	 * @param startI The first base pair to copy
	 * @param len The length of the piece to copy
	 * @return A copy of the gnFragmentSpec containing only the specified bases
	 */
	virtual gnFragmentSpec* CloneRange( const gnSeqI startI, const gnSeqI len ) const;

protected:
//	std::vector <gnContigSpec*> m_SpecList;
	//Feature stuff...
	std::vector <gnBaseFeature*> m_featureList;
	
}; // class gnFragmentSpec

inline
gnFragmentSpec* gnFragmentSpec::Clone() const{
	return new gnFragmentSpec(*this);
}
/*
inline
uint32 gnFragmentSpec::GetSpecListLength() const{
	return m_SpecList.size();
}

inline
gnContigSpec *gnFragmentSpec::GetSpec( const uint32 i ) const{
	if(i < m_SpecList.size())
		return m_SpecList[i];
	Throw_gnEx(ContigIndexOutOfBounds());
}

inline
gnContigSpec* gnFragmentSpec::GetSpecByBase( const gnSeqI baseI ) const{
	return m_SpecList[GetSpecIndexByBase(baseI)];
}

inline
void gnFragmentSpec::AddSpec( gnBaseSpec* spec, const uint32 i ){
	uint32 index = i == UINT32_MAX ? m_SpecList.size() : i;
	if(index <= m_SpecList.size()){
		m_SpecList.insert(m_SpecList.begin() + index, (gnContigSpec*)spec);
	}
}

inline
void gnFragmentSpec::RemoveSpec( uint32 i ){
	if(i < GetSpecListLength()){
		m_SpecList.erase(m_SpecList.begin() + i);
	}else
		Throw_gnEx(ContigIndexOutOfBounds());
}
*/
inline
uint32 gnFragmentSpec::AddFeature( gnBaseFeature* feat ) 
{
	m_featureList.push_back(feat);
	feat->SetSpec(this);
	return m_featureList.size()-1;
}
inline
uint32 gnFragmentSpec::GetFeatureListLength() const
{
	return m_featureList.size();
}
inline
gnBaseFeature* gnFragmentSpec::GetFeature(const uint32 i) const
{
	return m_featureList[i]->Clone();
}
inline
void gnFragmentSpec::RemoveFeature( const uint32 i) 
{
	if(i >= m_featureList.size())
		Throw_gnEx(FeatureIndexOutOfBounds());
	delete m_featureList[i];
	m_featureList.erase(m_featureList.begin() + i);
}


}	// end namespace genome

#endif
	// _gnFragmentSpec_h_
