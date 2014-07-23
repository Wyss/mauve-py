/////////////////////////////////////////////////////////////////////////////
// File:            gnGenomeSpec.h
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

#ifndef _gnGenomeSpec_h_
#define _gnGenomeSpec_h_

#include "libGenome/gnDefs.h"

#include <vector>
#include <string>

#include "libGenome/gnClone.h"
#include "libGenome/gnBaseFeature.h"
#include "libGenome/gnBaseHeader.h"
#include "libGenome/gnMultiSpec.h"
#include "libGenome/gnFragmentSpec.h"
#include "libGenome/gnException.h"


namespace genome {

/**
 * This class contains references to the sequence data, annotation,
 * and other related header data for an organism's genome.  The genome
 * is organized into a list of sequence fragments (specs), a list of
 * features, and a list of headers.
 */
class GNDLLEXPORT gnGenomeSpec : public gnMultiSpec< gnFragmentSpec >
{
public:
	gnGenomeSpec();
	/**
	 * Copy constructor.
	 * @param s the gnGenomeSpec to copy.
	 */
	gnGenomeSpec( const gnGenomeSpec& s);
	/**
	 * Destructor, frees memory.
	 */
	virtual ~gnGenomeSpec();
// Clone
	virtual gnGenomeSpec* Clone() const;
	virtual void Clear();
// Base Spec stuff
	virtual void SetReverseComplement( const boolean value );

//Multispec stuff
/*	virtual uint32 GetSpecListLength() const;
	virtual gnFragmentSpec* GetSpec( const uint32 i ) const;
	virtual gnFragmentSpec* GetSpecByBase( const gnSeqI baseI ) const;
	virtual void AddSpec( gnBaseSpec* spec, const uint32 i = UINT32_MAX );
	virtual void RemoveSpec( uint32 i );
*/
	virtual void MergeFragments( const uint32 startC, const uint32 endC);

	virtual uint32 AddFeature( gnBaseFeature* feat );
	virtual uint32 GetFeatureListLength() const;
	virtual gnBaseFeature* GetFeature( const uint32 i ) const;
	virtual void GetContainedFeatures(const gnLocation& lt, std::vector<gnBaseFeature*>& feature_vector, std::vector<uint32>& index_vector) const;
	virtual void GetIntersectingFeatures(const gnLocation& lt, std::vector<gnBaseFeature*>& feature_vector, std::vector<uint32>& index_vector) const;
	virtual void GetBrokenFeatures(const gnLocation& lt, std::vector<gnBaseFeature*>& feature_vector) const;
	virtual void RemoveFeature( const uint32 i );
	/**
	 * Copies a specified range of bases and returns a pointer to
	 * the resulting gnGenomeSpec.  You must delete the copy when you
	 * are finished with it.
	 * @param startI The first base pair to copy
	 * @param len The length of the piece to copy
	 * @return A copy of the gnGenomeSpec containing only the specified bases
	 */
	virtual gnGenomeSpec* CloneRange( const gnSeqI startI, const gnSeqI len ) const;

protected:
//	std::vector <gnFragmentSpec*> m_SpecList;
	
}; // class gnGenomeSpec

inline
gnGenomeSpec* gnGenomeSpec::Clone() const{
	return new gnGenomeSpec(*this);
}
/*
inline
uint32 gnGenomeSpec::GetSpecListLength() const{
	return m_SpecList.size();
}

inline
gnFragmentSpec* gnGenomeSpec::GetSpec( const uint32 i ) const{
	if(i < m_SpecList.size())
		return m_SpecList[i];
	Throw_gnEx(FragmentIndexOutOfBounds());
}

inline
gnFragmentSpec* gnGenomeSpec::GetSpecByBase( const gnSeqI baseI ) const{
	return m_SpecList[GetSpecIndexByBase(baseI)];
}

inline
void gnGenomeSpec::AddSpec( gnBaseSpec* spec, const uint32 i ){
	uint32 index = i == UINT32_MAX ? m_SpecList.size() : i;
	if(index <= m_SpecList.size()){
		m_SpecList.insert(m_SpecList.begin() + index, (gnFragmentSpec*)spec);
	}
}

inline
void gnGenomeSpec::RemoveSpec( uint32 i ){
	if(i < GetSpecListLength()){
		m_SpecList.erase(m_SpecList.begin() + i);
	}
}
*/

}	// end namespace genome

#endif
	// _gnGenomeSpec_h_
