/////////////////////////////////////////////////////////////////////////////
// File:            gnCompare.h
// Purpose:         Compares all sequences
// Description:     Compares sequences
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

#ifndef _gnCompare_h_
#define _gnCompare_h_

#include "libGenome/gnDefs.h"

#include <string>
#include <cstring>
#include <cctype>
#include "libGenome/gnClone.h"


namespace genome {

class GNDLLEXPORT gnCompare : public gnClone
{
public:
	
	/**
	 * Returns a pointer to a ready to use Protein comparison object
	 * <b>Do not delete</b> the object when finished with it, it is a static object.
	 */
	static const gnCompare *ProteinSeqCompare();
	/**
	 * Returns a pointer to a ready to use DNA comparison object
	 * This comparator can handle comparisons between ambiguous DNA sequences.
	 * <b>Do not delete</b> the object when finished with it, it is a static object.
	 */
	static const gnCompare *DNASeqCompare();
	/**
	 * Returns a pointer to a ready to use RNA comparison object
	 * This comparator can handle comparisons between ambiguous RNA sequences.
	 * <b>Do not delete</b> the object when finished with it, it is a static object.
	 */
	static const gnCompare *RNASeqCompare();
	
	enum gnCompareType{
		ProteinSeqCompareType,
		DNASeqCompareType,
		RNASeqCompareType,
	};

	gnCompare();
	/**
	 * Creates a gnCompare for a predefined compare type.
	 * Used by the static compare constructors to avoid the "static initialization order fiasco"
	 * For general use of a predefined comaprison type, you should use one of the static
	 * member functions to get a pointer to that comparison object.
	 * @param c_type The type of comparator to create.
	 * @see gnCompareType
	 */
	gnCompare( const gnCompareType c_type );
	~gnCompare();
	
	gnCompare* Clone() const;

	std::string GetName() const;
	void SetName( std::string name );

	// Use less than comparisons since all other operators can be derived
	// gnSeqC less than operations
	boolean LessThan( gnSeqC ch, gnSeqC ch2, boolean case_sensitive = false) const;
	/** True if ch2 is equal to or contained within the scope of ch
	 *
	 */
	boolean Contains( gnSeqC ch, gnSeqC ch2, boolean case_sensitive = false ) const;
	// True if ch2 is general enough to equal ch
	// gnSeqC[] less than operations
	boolean LessThan( const gnSeqC* seq, const gnSeqC* seq2, const uint32 len, boolean case_sensitive = false ) const;
	boolean Contains( const gnSeqC* seq, const gnSeqC* seq2, const uint32 len, boolean case_sensitive = false ) const;
	// std::string
	boolean LessThan( const std::string &seq, const std::string &seq2, boolean case_sensitive = false) const;
	boolean Contains( const std::string &seq, const std::string &seq2, boolean case_sensitive = false) const;
	
	// fill map

	//adds a character which is equivalent to itself
	void SetSingle( const gnSeqC ch );
	//adds a pair of equivalent characters
	void SetPair( const gnSeqC ch, const gnSeqC ch2 );
	//adds ch as being wholly contained by ch2
	void SetContained( const gnSeqC ch, const gnSeqC ch2 );

	void RemoveSingle( const gnSeqC ch );
	void RemovePair( const gnSeqC ch, const gnSeqC ch2 );
	void RemoveContained( const gnSeqC ch, const gnSeqC ch2 );
	
private:
	gnCompare( const gnCompare& sf ){};
	void CreateProteinComparator();
	void CreateDNAComparator();
	void CreateRNAComparator();

	void AddArrayEntry( gnSeqC *array[GNSEQC_MAX], const gnSeqC ch, const gnSeqC ch2);
	void DelArrayEntry( gnSeqC *array[GNSEQC_MAX], const gnSeqC ch, const gnSeqC ch2);

	std::string m_name;
	boolean m_ignoreCase;
	
	gnSeqC* m_pairArray[GNSEQC_MAX];
	gnSeqC* m_containArray[GNSEQC_MAX];
	
};//class gnCompare

inline
gnCompare* gnCompare::Clone() const{
	return new gnCompare(*this);
}

inline
std::string gnCompare::GetName() const{
	return m_name;
}
inline
void gnCompare::SetName( std::string name ){
	m_name = name;
}

// gnSeqC
inline
boolean gnCompare::LessThan( gnSeqC ch, gnSeqC ch2, boolean case_sensitive) const
{
	if(!case_sensitive){
		ch = toupper(ch);
		ch2 = toupper(ch2);
	}
		
	if(strchr(m_pairArray[ch], ch2) == 0)
		return ch < ch2 ? true : false;
	return false;
}


// gnSeqC[]
inline
boolean gnCompare::LessThan( const gnSeqC* seq, const gnSeqC* seq2, const uint32 len, boolean case_sensitive ) const{
	for( uint32 i=0; i < len ; ++i )
		if(LessThan(seq[i], seq2[i], case_sensitive))
			return true;
	return false;
}

// std::string
inline
boolean gnCompare::LessThan( const std::string &seq, const std::string &seq2, boolean case_sensitive) const
{
	gnSeqI shorter_len = seq.length() < seq2.length() ? seq.length() : seq2.length();
	return LessThan( (gnSeqC*)seq.data(), (gnSeqC*)seq2.data(), shorter_len, case_sensitive );
}

// fill map
inline
void gnCompare::SetSingle( const gnSeqC ch ){
	AddArrayEntry(m_pairArray, ch, ch);
	AddArrayEntry(m_containArray, ch, ch);
}
inline
void gnCompare::SetPair( const gnSeqC ch, const gnSeqC ch2 ){
	AddArrayEntry(m_pairArray, ch, ch2);
	AddArrayEntry(m_pairArray, ch2, ch);
}
inline
void gnCompare::SetContained( const gnSeqC ch, const gnSeqC ch2 ){
	AddArrayEntry(m_containArray, ch2, ch);
}

inline
void gnCompare::RemoveSingle( const gnSeqC ch )
{
	DelArrayEntry(m_pairArray, ch, ch);
	DelArrayEntry(m_containArray, ch, ch);
}
inline
void gnCompare::RemovePair( const gnSeqC ch, const gnSeqC ch2 )
{
	DelArrayEntry(m_pairArray, ch, ch2);
	DelArrayEntry(m_pairArray, ch2, ch);
}
inline
void gnCompare::RemoveContained( const gnSeqC ch, const gnSeqC ch2 )
{
	DelArrayEntry(m_containArray, ch2, ch);
}


}	// end namespace genome

#endif
	// _gnSeqCompare_h_
