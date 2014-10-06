/////////////////////////////////////////////////////////////////////////////
// File:            gnLocation.h
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

#ifndef _gnLocation_h_
#define _gnLocation_h_

#include "libGenome/gnDefs.h"
#include <string>
#include <iostream>
#include "libGenome/gnClone.h"


namespace genome {

/**
 * This class is used to store sequence locations.  gnBaseFeature uses it
 * to track feature locations.  gnLocation is capable of representing any 
 * GenBank style location qualifier.  gnLocation stores the start coordinate
 * and a startLength value which represents the length of an undetermined
 * region of sequence prior to the starting coordinate.  A startLength of
 * GNSEQI_END represents an unknown start, whereas a startLength of 0
 * implies that the start coordinate is unambiguous.  Any value other value
 * indicates a range of coordinates for the starting position.  The range is
 * bounded by the given start coordinate and extends startLength characters
 * upstream.  Likewise, an endLength of GNSEQI_END represents an unknown 
 * ambiguous end coordinate, and an endLength of 0 implies that the end 
 * coordinate is unambiguous.  Any other value for endLength denotes a range
 * of ending coordinates starting at end and continuing for endLength characters.
 */

class GNDLLEXPORT gnLocation : public gnClone
{
public:
	enum intersectRestriction{
		determinedRegions,
		undeterminedRegions,
		allRegions
	};
	
	enum gnLocationType{
		LT_Standard,	//standard == join multiple locations.
		LT_BetweenBases,
		LT_Complement,
		LT_Order,
		LT_Group,
		LT_OneOf,
		LT_Nothing
	};

	static const gnSeqI Defined = 0;
	static const gnSeqI Unknown = GNSEQI_END;

public:
	/**
	 *	Empty constructor, does nothing.
	 */
	gnLocation();
	/**
	 *	Copy constructor, copies a gnLocation.
	 * @param s the location to copy.
	 */
	gnLocation( const gnLocation& s);
	/**
	 *	Constructor, creates a gnLocation.
	 * @param start The start position within the sequence.
	 * @param end The end position within the sequence.
	 * @param type The type of the location, LT_Standard by default.
	 * @param contigName The name of the contig associated with this location, empty by default.
	 * @see gnLocationType
	 */
	gnLocation( const gnSeqI start, const gnSeqI end, const gnLocationType type = LT_Standard, std::string contigName = "");
	/**
	 *	Constructor, creates a gnLocation.
	 * @param start The start position within the sequence.
	 * @param startLength GNSEQI_END if the starting location is unbounded to the left, 0 if defined, otherwise it specifies a range
	 * @param end The end position within the sequence.
	 * @param endLength 0 if defined, GNSEQI_END if unbounded to the right, otherwise specifies a range
	 * @param type The type of the location, LT_Standard by default.
	 * @param contigName The name of the contig associated with this location, empty by default.
	 * @see gnLocationType
	 */
	gnLocation( const gnSeqI start, const gnSeqI startLength, const gnSeqI end, gnSeqI endLength, const gnLocationType type = LT_Standard, std::string contigName = "");

	gnLocation* Clone() const;

	/**
	 * Resets this location to default values
	 */
	void Clear();
	/**
	 * Crops this location to fit within the specified location boundaries.
	 * @param l The location boundaries.
	 * @return True if the location still exists, false if the crop amount is larger than the location.
	 */
	boolean CropTo( const gnLocation &l );
	/**
	 * Crops the start location by the specified amount.
	 * @param start The amount to crop.
	 * @return True if the location still exists, false if the crop amount is larger than the location.
	 */
	boolean CropStart( const gnSeqI start );
	/**
	 * Crops the end location by the specified amount.
	 * @param end The amount to crop.
	 * @return True if the location still exists, false if the crop amount is larger than the location.
	 */
	boolean CropEnd( const gnSeqI end );

	/**
	 * Returns true if the given location intersects this location in the region specified by the parameter "ir".
	 * @param l The location to check for intersection.
	 * @param ir The region to check for intersection.  This can be either determinedRegions, undeterminedRegions, or allRegions
	 * @see intersectRestriction 
	 * @return True if the location still exists, false if the crop amount is larger than the location.
	 */
	boolean Intersects( const gnLocation &l, const intersectRestriction ir = determinedRegions ) const;
	/**
	 * Checks wether another gnLocation is contained by this one.
	 * @param l The gnLocation to be checked.
	 * @param cr A specification of the type of intersection to be checked.
	 * @see intersectRestriction 
	 * @return True if the location still exists, false if the crop amount is larger than the location.
	 */
	boolean Contains( const gnLocation& l, const intersectRestriction cr = determinedRegions ) const;
	/**
	 * Increases the location start and end.
	 * @param diff The amount to increase the location start and end.
	 */
	boolean MovePositive( const gnSeqI diff );
	/**
	 * Decreases the location start and end.
	 * @param diff The amount to decrease the location start and end.
	 */
	boolean MoveNegative( const gnSeqI diff );
	/**
	 * Moves the location start and end.
	 * @param direction Negative integer will decrease start and end, positive will increase them.
	 * @param diff The amount to move the location start and end.
	 */
	boolean MoveTo( const int direction, const gnSeqI diff );

	/**
	 * Gets the ending base pair.
	 * @return The ending base pair.
	 */
	gnSeqI GetEnd() const;
	/**
	 * Gets the ending base pair's definition value.
	 * @return -1 if the ending location is unbounded to the left, 0 if defined, 1 if unbounded to the right
	 */
	gnSeqI GetEndLength() const;
	gnSeqI GetLast() const;
	/**
	 * Gets the starting base pair.
	 * @return The starting base pair.
	 */
	gnSeqI GetStart() const;
	/**
	 * Gets the starting base pair's definition value.
	 * @return -1 if the starting location is unbounded to the left, 0 if defined, 1 if unbounded to the right
	 */
	gnSeqI GetStartLength() const;
	gnSeqI GetFirst() const;

	/**
	 * Gets the location type.
	 * @return the location type.
	 * @see gnLocationType
	 */
	gnLocationType GetType() const;

	/**
	 * Gets the location's boundaries.
	 * @param s The start index.
	 * @param sl The start length.
	 * @param e The end index.
	 * @param el The end length.
	 */
	void GetBounds( gnSeqI &s, gnSeqI &sl, gnSeqI &e, gnSeqI &el ) const;

	/**
	 * Returns true if the end is unbounded to the right.
	 * @return True if the end is unbounded to the right.
	 */
	bool IsEndBoundLonger() const;
	/**
	 * Returns true if the start is unbounded to the right.
	 * @return True if the start is unbounded to the right.
	 */
	bool IsStartBoundLonger() const;
	/**
	 * Returns true if the end is unbounded to the left.
	 * @return True if the end is unbounded to the left.
	 */
	bool IsEndBoundShorter() const;
	/**
	 * Returns true if the start is unbounded to the left.
	 * @return True if the start is unbounded to the left.
	 */
	bool IsStartBoundShorter() const;

	/**
	 * Sets the end location.
	 * @param end The end location.
	 */
	void SetEnd( const gnSeqI end );
	/**
	 * Sets the end location and definition.
	 * @param end The end location.
	 * @param endLength -1 if the ending location is unbounded to the left, 0 if defined, 1 if unbounded to the right
	 */
	void SetEnd( const gnSeqI end, const gnSeqI endLength );
	/**
	 * Sets the end definition.
	 * @param endLength -1 if the ending location is unbounded to the left, 0 if defined, 1 if unbounded to the right
	 */
	void SetEndLength( const gnSeqI endLength );

	/**
	 * Sets the start location/
	 * @param start The start location.
	 */
	void SetStart( const gnSeqI start );
	/**
	 * Sets the start location and definition.
	 * @param start The end location.
	 * @param startLength -1 if the ending location is unbounded to the left, 0 if defined, 1 if unbounded to the right
	 */
	void SetStart( const gnSeqI start, const gnSeqI startLength );
	/**
	 * Sets the start definition.
	 * @param startLength -1 if the ending location is unbounded to the left, 0 if defined, 1 if unbounded to the right
	 */
	void SetStartLength( const gnSeqI startLength );
	/**
	 * Sets the location type.
	 * @param lt The location type.
	 */
	void SetType( const gnLocationType lt );
	/**
	 * Sets the location's boundaries.
	 * @param start The start index.
	 * @param startLength -1 if the ending location is unbounded to the left, 0 if defined, 1 if unbounded to the right
	 * @param end The end index.
	 * @param endLength -1 if the ending location is unbounded to the left, 0 if defined, 1 if unbounded to the right
	 */
	void SetBounds( const gnSeqI start, const gnSeqI startLength, const gnSeqI end, const gnSeqI endLength);
	/**
	 * Sets the location's boundaries.
	 * @param start The start index.
	 * @param end The end index.
	 */
	void SetBounds( const gnSeqI start, const gnSeqI end);

	gnLocation GetUnion( const gnLocation &l ) const;
	
	gnLocation GetIntersection( const gnLocation &l, const intersectRestriction ir ) const;
	
private:
	std::string m_name;
	gnSeqI m_start;
	gnSeqI m_startLength;
	gnSeqI m_end;
	gnSeqI m_endLength;
	
	gnLocationType m_type;

}; // class gnLocation

// GET END
inline
gnSeqI gnLocation::GetEnd() const
{
	return m_end;
}
inline
gnSeqI gnLocation::GetEndLength() const
{
	return m_endLength;
}
inline
gnSeqI gnLocation::GetLast() const
{
	return m_end + m_endLength;
}
// GET START
inline
gnSeqI gnLocation::GetStart() const
{
	return m_start;
}
inline
gnSeqI gnLocation::GetStartLength() const
{
	return m_startLength;
}
inline
gnSeqI gnLocation::GetFirst() const
{
	return m_start > m_startLength ? m_start - m_startLength : 0;
}

inline
gnLocation::gnLocationType gnLocation::GetType() const
{
	return m_type;
}

inline
bool gnLocation::IsEndBoundLonger() const
{
	return m_endLength > 0;
}
inline
bool gnLocation::IsStartBoundLonger() const
{
	return m_startLength > 0;
}
inline
bool gnLocation::IsEndBoundShorter() const
{
	return m_endLength == GNSEQI_END;
}
inline
bool gnLocation::IsStartBoundShorter() const
{
	return m_startLength == GNSEQI_END;
}

inline
void gnLocation::SetEnd( const gnSeqI end )
{
	m_end = end;
}
inline
void gnLocation::SetEnd( const gnSeqI end, const gnSeqI endLength )
{
	m_end = end;
	m_endLength = endLength;
}
inline
void gnLocation::SetEndLength( const gnSeqI endLength )
{
	m_endLength = endLength;
}
inline
void gnLocation::SetStart( const gnSeqI start )
{
	m_start = start;
}
inline
void gnLocation::SetStart( const gnSeqI start, const gnSeqI startLength )
{
	m_start = start;
	m_startLength = startLength;
}
inline
void gnLocation::SetStartLength( const gnSeqI startLength )
{
	m_startLength = startLength;
}

inline
void gnLocation::SetType( const gnLocationType lt )
{
	m_type = lt;
}



}	// end namespace genome

#endif
	// _gnLocation_h_
