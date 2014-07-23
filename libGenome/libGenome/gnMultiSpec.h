/////////////////////////////////////////////////////////////////////////////
// File:            gnMultiSpec.h
// Purpose:         multiple spec template class
// Description:     Template for specs which can contain multiple data sources
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

#ifndef _gnMultiSpec_h_
#define _gnMultiSpec_h_

#include "libGenome/gnDefs.h"
#include "libGenome/gnException.h"

#include <vector>
#include <string>

#include "libGenome/gnClone.h"
#include "libGenome/gnBaseFeature.h"
#include "libGenome/gnBaseHeader.h"
#include "libGenome/gnBaseSpec.h"


namespace genome {

/**
 * This class defines an interface for specs which have multiple subspecs
 * containing sequence data.
 */
template< class SubSpec >
class GNDLLEXPORT gnMultiSpec : public gnBaseSpec
{
public:
	gnMultiSpec(){}
	/**
	 * Destructor, frees memory.
	 */
	virtual ~gnMultiSpec(){}
	virtual gnMultiSpec* Clone() const = 0;

// Base Spec stuff
	virtual gnSeqI GetLength() const;
	virtual void CropStart( gnSeqI cropLen );
	virtual void CropEnd( gnSeqI cropLen );

	virtual std::string GetSourceName() const;
	virtual void SetSourceName(const std::string& sourceName);
	virtual void Clear();

	virtual boolean SeqRead(const gnSeqI start, gnSeqC* buf, gnSeqI& bufLen, const uint32 contigI ) const;

//Multispec stuff
	/**
	 * Returns the number of subspecs this spec contains.
	 * @return The number of subspecs this spec contains.
	 */
	virtual uint32 GetSpecListLength() const;
	/**
	 * Get the spec at index i in the list of subspecs.
	 * @param i The index of the spec
	 * @return The spec.
	 */
	virtual SubSpec* GetSpec( const uint32 i ) const;
	/**
	 * Get the spec from the list of subspecs which contains the specified base pair.
	 * @param baseI The base pair index.
	 * @return The spec.
	 */
	virtual SubSpec* GetSpecByBase( const gnSeqI baseI ) const;
	/**
	 * Get the index of the subspec which contains the specified base pair.
	 * @param baseI The base pair index.
	 * @return The subspec index.
	 */
	virtual uint32 GetSpecIndexByBase( const gnSeqI baseI ) const;
	/**
	 * Get the index of the spec which has the given name.
	 * @param name The name of the spec to find.
	 * @return The spec index.
	 */
	virtual uint32 GetSpecIndexByName( const std::string& name ) const;
	/**
	 * Get the starting base pair, in this spec's sequence, of the given subspec.
	 * @param specI The subspec index.
	 * @return The subspec's starting base pair.
	 */
	virtual gnSeqI GetSpecStartBase( const uint32 specI ) const;
	/**
	 * Get the ending base pair, in this spec's sequence, of the given subspec.
	 * @param specI The subspec index.
	 * @return The subspec's ending base pair.
	 */
	virtual gnSeqI GetSpecEndBase( const uint32 specI ) const;
	/**
	 * Add a subspec to this spec.
	 * Throws an exception if the insertion index i is out of range
	 * @param spec The subspec to add.
	 * @param i The subspec to insert before
	 * @return True if successful
	 */
	virtual void AddSpec( SubSpec* spec, const uint32 i = UINT32_MAX );
	/**
	 * Remove a subspec from this spec.
	 * @param i The index of the subspec to remove.
	 */
	virtual void RemoveSpec( uint32 i );

	/**
	 * Returns the number of headers this spec contains.
	 * @return The number of headers this spec contains.
	 */
	virtual uint32 GetHeaderListLength() const;
	/**
	 * Add a header to this spec, adds to the end of the header list by default.
	 * @param head The header to add.
	 * @param i The header to insert before.  If i is larger than the size of the header list
	 * the header will be added to the end of the list.
	 */
	virtual void AddHeader( gnBaseHeader *head, const uint32 i = UINT32_MAX);
	/**
	 * Get the headers at index i in the list of headers.
	 * Throws a HeaderIndexOutOfBounds exception if a nonexistant header is referenced
	 * @param i The index of the header
	 * @return The header.
	 */
	virtual gnBaseHeader* GetHeader( const uint32 i ) const;
	/**
	 * Find the first header with the specified name, starting at index i.
	 * @param name The name of the header to find.
	 * @param i The index to start looking at.
	 * @return The header or NULL if none was found.
	 */
	virtual gnBaseHeader* GetHeader( const std::string& name, uint32& i) const;
	/**
	 * Remove a header from this spec.
	 * Throws a HeaderIndexOutOfBounds exception if a nonexistant header is referenced
	 * @param i The index of the header to remove.
	 */
	virtual void RemoveHeader( uint32 i );

	/**
	 * Add a feature to this spec.
	 * @param feat The feature to add.
	 * @return True if successful
	 */
	uint32 AddFeature( gnBaseFeature* feat );
	/**
	 * Returns the number of features this spec contains.
	 * @return The number of features this spec contains.
	 */
	virtual uint32 GetFeatureListLength() const = 0;
	/**
	 * Get the feature at index i in the list of features.
	 * @param i The index of the feature
	 * @return The feature.
	 */
	virtual gnBaseFeature* GetFeature( const uint32 i ) const = 0;
	/**
	 * Creates a list of all features which are contained by coordinates specified.
	 * @param lt The coordinates containing the features to return.
	 * @param feature_vector The std::vector to store features in.
	 * @param index_vector A std::vector of indices which correspond to the features.
	 */
	virtual void GetContainedFeatures(const gnLocation& lt, std::vector<gnBaseFeature*>& feature_vector, std::vector<uint32>& index_vector) const = 0;
	/**
	 * Creates a list of all features which intersect the coordinates specified.
	 * @param lt The coordinates intersecting the features to return.
	 * @param feature_vector The std::vector to store features in.
	 * @param index_vector A std::vector of indices which correspond to the features.
	 */
	virtual void GetIntersectingFeatures(const gnLocation& lt, std::vector<gnBaseFeature*>& feature_vector, std::vector<uint32>& index_vector) const = 0;
	/**
	 * Creates a list of features which may have been broken by an edit.
	 * @param lt The coordinates containing the features to return.
	 * @param feature_vector The std::vector to store features in.
	 */
	virtual void GetBrokenFeatures(const gnLocation& lt, std::vector<gnBaseFeature*>& feature_vector) const = 0;
	/**
	 * Remove a feature from this spec.
	 * @param i The index of the feature to remove.
	 * @return True if successful
	 */
	virtual void RemoveFeature( const uint32 i ) = 0;
	
protected:
	std::string m_sourceName;
	
	std::vector < SubSpec* > m_SpecList;
	std::vector < gnBaseHeader* > m_headerList;
private:

}; // class gnMultiSpec

template< class SubSpec >
inline
uint32 gnMultiSpec< SubSpec >::GetHeaderListLength() const
{
	return m_headerList.size();
}

template< class SubSpec >
inline
gnBaseHeader* gnMultiSpec< SubSpec >::GetHeader(const uint32 i) const
{
	if(i < m_headerList.size()){
		return m_headerList[i];
	}
	return 0;
}

template< class SubSpec >
inline
std::string gnMultiSpec< SubSpec >::GetSourceName() const{
	return m_sourceName;
}

template< class SubSpec >
inline
void gnMultiSpec< SubSpec >::SetSourceName(const std::string& sourceName){
	m_sourceName = sourceName;
}

template< class SubSpec >
uint32 gnMultiSpec< SubSpec >::GetSpecListLength() const{
	return m_SpecList.size();
}

template< class SubSpec >
SubSpec* gnMultiSpec< SubSpec >::GetSpec( const uint32 i ) const{
	if(i < m_SpecList.size())
		return m_SpecList[i];
	Throw_gnEx(FragmentIndexOutOfBounds());
}

template< class SubSpec >
SubSpec* gnMultiSpec< SubSpec >::GetSpecByBase( const gnSeqI baseI ) const{
	return m_SpecList[GetSpecIndexByBase(baseI)];
}

template< class SubSpec >
void gnMultiSpec< SubSpec >::AddSpec( SubSpec* spec, const uint32 i ){
	uint32 index = i == UINT32_MAX ? m_SpecList.size() : i;
	if(index <= m_SpecList.size()){
		m_SpecList.insert(m_SpecList.begin() + index, spec);
	}
}

template< class SubSpec >
void gnMultiSpec< SubSpec >::RemoveSpec( uint32 i ){
	if(i < GetSpecListLength()){
		m_SpecList.erase(m_SpecList.begin() + i);
	}
}

template< class SubSpec >
void gnMultiSpec< SubSpec >::Clear()
{
	gnBaseSpec::Clear();
	uint32 list_size = m_headerList.size();
	for(uint32 i=0; i < list_size; i++)
		delete m_headerList[i];
	m_headerList.clear();
}

template< class SubSpec >
gnSeqI gnMultiSpec< SubSpec >::GetLength() const
{
	gnSeqI subLen = 0;	//aggregate the lower levels.
	for(uint32 i=0; i < GetSpecListLength(); i++)
		subLen += GetSpec(i)->GetLength();
	return subLen;
}

//Return the index of the subspec which contains the base in question
template< class SubSpec >
uint32 gnMultiSpec< SubSpec >::GetSpecIndexByBase( const gnSeqI baseI ) const{
	gnSeqI cur_length = 0;
	for(uint32 i=0; i < GetSpecListLength(); i++){
		cur_length += GetSpec(i)->GetLength();
		if(baseI < cur_length)
			return i;
	}
	//if we made it here then the base was out of bounds
	Throw_gnEx(SeqIndexOutOfBounds());
}

template< class SubSpec >
uint32 gnMultiSpec< SubSpec >::GetSpecIndexByName( const std::string& name ) const{
	for(uint32 i=0; i < GetSpecListLength(); i++){
		if(name == GetSpec(i)->GetName())
			return i;
	}
	Throw_gnEx(SpecIndexOutOfBounds());
}

//returns the starting base pair index of the given subspec.
template< class SubSpec >
gnSeqI gnMultiSpec< SubSpec >::GetSpecStartBase( const uint32 specI ) const{
	uint32 i;
	if(specI >= GetSpecListLength())
		Throw_gnEx(SpecIndexOutOfBounds());
	
	gnSeqI start_base = 0;
	for(i=0; i < specI; i++){
		start_base += GetSpec(i)->GetLength();
	}
	return start_base;
}

template< class SubSpec >
gnSeqI gnMultiSpec< SubSpec >::GetSpecEndBase( const uint32 specI ) const{
	uint32 i;
	if(specI >= GetSpecListLength())
		Throw_gnEx(SpecIndexOutOfBounds());
	
	gnSeqI end_base = 0;
	for(i=0; i <= specI; i++){
		end_base += GetSpec(i)->GetLength();
	}
	return end_base;
}

template< class SubSpec >
void gnMultiSpec< SubSpec >::CropStart( gnSeqI cropLen ){
	gnSeqI curbase = 0;
	for(uint32 specI = 0; specI < GetSpecListLength(); specI++){
		curbase += GetSpec(specI)->GetLength();
		if(curbase <= cropLen){
			//delete the spec completely
			gnBaseSpec *tmp_spec = GetSpec(specI);
			RemoveSpec(specI);
			delete tmp_spec;
			specI--;
		}else{
			//recurse and break
			gnSeqI sub_len = cropLen - (curbase - GetSpec(specI)->GetLength());
			GetSpec(specI)->CropStart(sub_len);
			break;
		}
	}
}

template< class SubSpec >
void gnMultiSpec< SubSpec >::CropEnd( gnSeqI cropLen ){
	gnSeqI curbase = 0;
	gnSeqI cropbase = GetLength() - cropLen;  //last base to keep
	boolean trash_the_rest = false;
	for(uint32 specI = 0; specI < GetSpecListLength(); specI++){
		curbase += GetSpec(specI)->GetLength();		
		if(trash_the_rest){
			//delete the spec entirely
			gnBaseSpec *tmp_spec = GetSpec(specI);
			RemoveSpec(specI);
			delete tmp_spec;
			specI--;
			continue;
		}else if(curbase > cropbase){
			GetSpec(specI)->CropEnd(curbase - cropbase);
			trash_the_rest = true;
		}else if(curbase == cropbase)
			trash_the_rest = true;
	}
}

template< class SubSpec >
void gnMultiSpec< SubSpec >::AddHeader( gnBaseHeader* head, const uint32 i) 
{
	uint32 index = i == UINT32_MAX ? m_headerList.size() : i;
	m_headerList.insert(m_headerList.begin() + index, head);
}

template< class SubSpec >
gnBaseHeader* gnMultiSpec< SubSpec >::GetHeader( const std::string& name, uint32& i) const{
	for(; i < m_headerList.size(); i++){
		if( m_headerList[i]->GetHeaderName() == name)
			return m_headerList[i];
	}
	Throw_gnEx(HeaderIndexOutOfBounds());
}

template< class SubSpec >
void gnMultiSpec< SubSpec >::RemoveHeader( uint32 i) 
{
	if(i <= m_headerList.size()){
		m_headerList.erase(m_headerList.begin() + i);
	}
	Throw_gnEx(HeaderIndexOutOfBounds());
}

template< class SubSpec >
boolean gnMultiSpec< SubSpec >::SeqRead(const gnSeqI start, gnSeqC* buf, gnSeqI& bufLen, const uint32 contigI ) const{
	if( bufLen == 0 )
		return true;
	if(contigI == ALL_CONTIGS){
		gnSeqI curpos = 0;
		gnSeqI readBytes = 0;
		gnSeqI remainingBytes = bufLen;
		uint32 curSpecI = 0;
		//seek to start spec.
		for(curSpecI=0; curSpecI < GetSpecListLength(); curSpecI++){
			curpos += GetSpec(curSpecI)->GetLength();
			if(curpos > start)
				break;
		}
		if(curpos <= start)
			Throw_gnEx(SeqIndexOutOfBounds());
		//read until we're done;
		while((remainingBytes > 0) && (curSpecI < GetSpecListLength())) {
			gnSeqI readable = GetSpec(curSpecI)->GetLength();
			//check for circular
			gnSeqI start_pos = readBytes == 0 ? start - (curpos - readable) : 0;
			gnSeqI to_read = readable - start_pos >= remainingBytes ? remainingBytes : readable - start_pos;
			boolean success = GetSpec(curSpecI)->SeqRead(start_pos, buf+readBytes, to_read, contigI);

			readBytes += to_read;
			remainingBytes -= to_read;
			if(!success)
				break;
			curSpecI++;
		}
		bufLen = readBytes;
		return true;
	}else{	//read from the specified contig.
		if(contigI < GetSpecListLength())
			return GetSpec(contigI)->SeqRead(start, buf, bufLen, ALL_CONTIGS);
		else
			Throw_gnEx(SpecIndexOutOfBounds());
	}
	return false;
}




}	// end namespace genome

#endif
	// _gnMultiSpec_h_
