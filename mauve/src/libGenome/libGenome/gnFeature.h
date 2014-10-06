/////////////////////////////////////////////////////////////////////////////
// File:            gnFeature.h
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

#ifndef _gnFeature_h_
#define _gnFeature_h_

#include "libGenome/gnDefs.h"

#include <string>
#include <vector>
#include "libGenome/gnBaseFeature.h"
#include "libGenome/gnBaseQualifier.h"


namespace genome {

/**
 * gnFeature stores sequence features in memory.
 * It contains a list of locations and qualifiers which are used to describe this feature.
 * It can be referred to by ID.
 */
class GNDLLEXPORT gnFeature : public gnBaseFeature
{
public:
	/**
	 * Empty constructor.
	 */
	gnFeature( );
	/**
	 * Creates a memory feature with the specified name.
	 * @param name The name of the feature.
	 * @param id The id number of the feature
	 * @param lt The type of sequence location covered by this feature
	 * @param broken True if the feature was broken by some sequence manipulation, false otherwise.
	 */
	gnFeature( std::string& name, uint32 id = 0, gnLocation::gnLocationType lt = gnLocation::LT_Nothing, boolean broken = false );
	/**
	 * Copy constructor.
	 * @param s The gnFeature to copy.
	 */
	gnFeature( const gnFeature& s );
	/**
	 * Destructor, frees memory.
	 */
	~gnFeature();
// Clone
	gnFeature* Clone() const;

private:
}; // class gnFeature

inline
gnFeature* gnFeature::Clone() const
{
	return new gnFeature(*this);
}


}	// end namespace genome

#endif
	// _gnFeature_h_
