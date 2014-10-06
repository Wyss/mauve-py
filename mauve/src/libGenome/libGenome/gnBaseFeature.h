/////////////////////////////////////////////////////////////////////////////
// File:            gnBaseFeature.h
// Purpose:         abstract Feature class
// Description:     Provides an interface for Features in memory and on disk.
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

#ifndef _gnBaseFeature_h_
#define _gnBaseFeature_h_

#include "libGenome/gnDefs.h"
#include <string>
#include <vector>
#include "libGenome/gnClone.h"
#include "libGenome/gnLocation.h"
#include "libGenome/gnBaseQualifier.h"

namespace genome {


class gnFragmentSpec;

/**
 * Features of DNA sequences can be accessed and manipulated using classes
 * derived from gnBaseFeature.  gnBaseFeature outlines a basic interface
 * and gives functionality its derived sequence classes, such as
 * gnCDSFeature.
 */
class GNDLLEXPORT gnBaseFeature : public gnClone
{
public:
	gnBaseFeature();
	gnBaseFeature( std::string& name, uint32 id = 0, gnFragmentSpec* spec = NULL, gnLocation::gnLocationType lt = gnLocation::LT_Nothing, boolean broken = false );
	/**
	 * Destructor, frees memory.
	 */
	~gnBaseFeature();

	virtual gnBaseFeature* Clone() const = 0;
	/**
	 * Gets the feature name (e.g. CDS, trna...).
	 * @return The feature name.
	 */
	virtual std::string GetName() const;
	/**
	 * Sets the feature name.
	 * @param name The feature name.
	 */
	virtual void SetName( const std::string& name );
	/**
	 * Gets this feature's ID.
	 * @return The feature's ID.
	 */
	virtual uint32 GetID() const;
	/**
	 * Sets this feature's ID.
	 * @param id The feature's ID.
	 */
	virtual void SetID(uint32 id);
	/**
	 * Gets the gnFragmentSpec to which this feature refers.
	 * @return The feature's ID.
	 * @see gnFragmentSpec
	 */
	virtual gnFragmentSpec* GetSpec() const;
	/**
	 * Sets the gnFragmentSpec to which this feature refers.
	 * @param spec A pointer to the fragment spec.
	 */
	virtual void SetSpec(gnFragmentSpec* spec);
	/**
	 * Gets this feature's location type.
	 * @return The feature's location type.
	 * @see gnLocationType
	 */
	virtual gnLocation::gnLocationType GetLocationType() const;
	/**
	 * Sets this feature's location type.
	 * LT_BetweenBases is not a valid feature location type.
	 * @param lType The feature's location type.
	 * @see gnLocationType
	 */
	virtual void SetLocationType( gnLocation::gnLocationType lType );
	/**
	 * Returns the number of locations this feature describes.
	 * @return The number of locations this feature describes.
	 */
	virtual uint32 GetLocationListLength() const;
	/**
	 * Adds a new location before the location at listI.
	 * @param l The gnLocation to add.
	 * @param listI The index into the location list to insert before.
	 * @return True if successful, false otherwise.
	 * @see gnLocation
	 */
	virtual boolean AddLocation( const gnLocation& l, uint32 listI = 0);
	/**
	 * Gets the location at listI.
	 * @param listI The index into the location list to get.
	 * @return The location.
	 */
	virtual gnLocation GetLocation( uint32 listI ) const;
	/**
	 * Deletes the location at listI.
	 * @param listI The index into the location list to delete.
	 * @return True if successful, false otherwise.
	 */
	virtual boolean RemoveLocation( uint32 listI );
	/**
	 * Sets the location at listI to a new location.
	 * @param l The new location
	 * @param listI The index into the location list to set.
	 * @return True if successful, false otherwise.
	 */
	virtual boolean SetLocation( const gnLocation& l, uint32 listI );
	/**
	 * Increases this feature's coordinates by a specific number of bases.
	 * @param i The number of bases to increase by.
	 * @return True if successful, false if the feature was broken by the change.
	 */
	virtual boolean MovePositive( const gnSeqI i );
	/**
	 * Decreases this feature's coordinates by a specific number of bases.
	 * @param i The number of bases to decrease by.
	 * @return True if successful, false if the feature was broken by the change.
	 */
	virtual boolean MoveNegative( const gnSeqI i );
	/**
	 * Crops the start locations of this feature by the specified amount.
	 * @param i The amount to crop.
	 * @return True if successful, false if the feature was broken by the change.
	 */
	virtual boolean CropStart( const gnSeqI i );
	/**
	 * Crops the end locations of this feature by the specified amount.
	 * @param i The amount to crop.
	 * @return True if successful, false if the feature was broken by the change.
	 */
	virtual boolean CropEnd( const gnSeqI i );
	/**
	 * Crops the locations of this feature to fit within the given location.
	 * @param l The location to crop reduce to.
	 * @return True if successful, false if the feature was broken by the change.
	 */
	virtual boolean Crop( const gnLocation& l );
	/**
	 * Returns true if the feature is broken
	 * @return True if the feature is broken
	 */
	virtual boolean IsBroken() const;
	/**
	 * Sets whether the feature is broken or not.
	 * @param broke True if the feature should be broken, false otherwise
	 */
	virtual void SetBroken(boolean broke);

	/**
	 * Checks if the given coordinate is contained by this feature
	 * @param i The coordinate to check
	 * @return True if "i" is contained by this feature.  False otherwise
	 */
	virtual boolean Contains( gnSeqI i ) const;
	/**
	 * Checks if the given location is contained by this feature
	 * @param l The location to check
	 * @return True if "l" is contained by this feature.  False otherwise
	 */
	virtual boolean Contains( const gnLocation& l ) const;
	/**
	 * Checks if the given feature is entirely contained by this feature
	 * @param feature The feature to check
	 * @return True if "feature" is contained by this feature.  False otherwise
	 */
	virtual boolean Contains( gnBaseFeature* feature ) const;
	/**
	 * Checks if this feature is entirely contained by the given location
	 * @param l The location to check
	 * @return True if "l" contains by this feature.  False otherwise
	 */
	virtual boolean IsContainedBy( const gnLocation& l ) const;
	/**
	 * Checks if the given location intersects this feature
	 * @param l The location to check
	 * @return True if "l" intersects this feature.  False otherwise
	 */
	virtual boolean Intersects( const gnLocation& l ) const;
	/**
	 * Checks if the given feature intersects this feature
	 * @param feature The location to check
	 * @return True if "feature" intersects this feature.  False otherwise
	 */
	virtual boolean Intersects( gnBaseFeature* feature ) const;

	/**
	 * Returns the number of qualifiers in this feature.
	 * @return The number of qualifiers in this feature.
	 */
	virtual uint32 GetQualifierListLength() const;
	/**
	 * Adds a new qualifier.
	 * @param qualifier The qualifier to add.
	 * @return True if successful, false otherwise.
	 */
	virtual boolean AddQualifier( gnBaseQualifier* qualifier );
	/**
	 * Looks for a qualifier by name.
	 * @param name The name of the qualifier to look for.
	 * @return True if a qualifier was found, false otherwise.
	 */
	virtual boolean HasQualifier( const std::string& name ) const;
	/**
	 * Searches for a qualifier by name, starting at the index listI.
	 * @param name The name of the qualifier to look for.
	 * @param listI The index into the qualifier list to start the search at.
	 * @return The index of the qualifier.
	 */
	virtual uint32 FirstIndexOfQualifier( const std::string& name, uint32 listI ) const;
	/**
	 * Searches for a qualifier by name, ending at the index listI.
	 * @param name The name of the qualifier to look for.
	 * @param listI The index into the qualifier list to end the search at.
	 * @return The index of the qualifier.
	 */
	virtual uint32 LastIndexOfQualifier( const std::string& name, uint32 listI ) const;
	/**
	 * Gets the name of the qualifier at the list index listI.
	 * @param listI The index into the qualifier list.
	 * @return The name of the qualifier.
	 */
	virtual std::string GetQualifierName( uint32 listI ) const;
	/**
	 * Gets the value of the qualifier at the list index listI.
	 * @param listI The index into the qualifier list.
	 * @return The value of the qualifier.
	 */
	virtual std::string GetQualifierValue( uint32 listI ) const;
	/**
	 * Gets the qualifier at the list index listI.
	 * @param listI The index into the qualifier list.
	 * @return The qualifier.
	 */
	virtual gnBaseQualifier* GetQualifier( uint32 listI );
	/**
	 * Deletes the qualifier at the list index listI.
	 * @param listI The index into the qualifier list.
	 * @return True if successful, false otherwise.
	 */
	virtual boolean RemoveQualifier( uint32 listI );
	/**
	 * Set the name and value of the qualifier at the list index listI.
	 * @param name The new name of the qualifier.
	 * @param value The new value of the qualifier.
	 * @param listI The index into the qualifier list.
	 * @return True if successful, false otherwise.
	 */
	virtual boolean SetQualifier( std::string& name, std::string& value, uint32 listI );
	/**
	 * Set the name of the qualifier at the list index listI.
	 * @param name The new name of the qualifier.
	 * @param listI The index into the qualifier list.
	 * @return True if successful, false otherwise.
	 */
	virtual boolean SetQualifierName( std::string& name, uint32 listI );
	/**
	 * Set the value of the qualifier at the list index listI.
	 * @param value The new value of the qualifier.
	 * @param listI The index into the qualifier list.
	 * @return True if successful, false otherwise.
	 */
	virtual boolean SetQualifierValue( std::string& value, uint32 listI );
protected:
	uint32 m_id;
	std::string m_name;
	boolean m_broken;
	gnLocation::gnLocationType m_locationType;
	std::vector< gnLocation > m_locationList;
	std::vector< gnBaseQualifier* > m_qualifierList;
	gnFragmentSpec* m_spec;
};// class gnBaseFeature

inline
std::string gnBaseFeature::GetName() const{
	return m_name;
}
inline
void gnBaseFeature::SetName( const std::string& name ){
	m_name = name;
}
inline
uint32 gnBaseFeature::GetID() const{
	return m_id;
}
inline
void gnBaseFeature::SetID(uint32 id){
	m_id = id;
}
inline
gnFragmentSpec* gnBaseFeature::GetSpec() const{
	return m_spec;
}
inline
void gnBaseFeature::SetSpec(gnFragmentSpec* spec){
	m_spec = spec;
}
inline
boolean gnBaseFeature::IsBroken() const{
	return m_broken;
}
inline
void gnBaseFeature::SetBroken(boolean broke){
	m_broken = broke;
}

inline
gnLocation::gnLocationType gnBaseFeature::GetLocationType() const{
	return m_locationType;
}
inline
void gnBaseFeature::SetLocationType( gnLocation::gnLocationType lType ){
	m_locationType = lType;
}
inline
uint32 gnBaseFeature::GetLocationListLength() const{
	return m_locationList.size();
}
inline
boolean gnBaseFeature::Crop( const gnLocation& l ){
	return CropStart(l.GetStart()) && CropEnd(l.GetEnd());
}
inline
uint32 gnBaseFeature::GetQualifierListLength() const{
	return m_qualifierList.size();
}


}	// end namespace genome

#endif
	// _gnBaseFeature_h_
