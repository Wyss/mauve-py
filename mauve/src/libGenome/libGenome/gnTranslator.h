/////////////////////////////////////////////////////////////////////////////
// File:            gnTranslator.h
// Purpose:         Translator for all Sequences
// Description:     Translates DNA and protein sequences
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

#ifndef _gnTranslator_h_
#define _gnTranslator_h_


#include "libGenome/gnDefs.h"

#include <string>
#include <vector>
#include "libGenome/gnClone.h"
#include "libGenome/gnBaseFilter.h"
#include "libGenome/gnCompare.h"

namespace genome {

/**
 * The gnTranslator class facilitates the translation of ambiguous DNA and RNA 
 * sequences to protein sequence.  It can also perform reverse translation from
 * Protein to RNA or DNA though the utility of this functionality is somewhat
 * questionable since it does not assign nucleotides probabilistically.
 * Do not attempt to instantiate this class unless you are defining a new
 * translation type.  Instead use the static member access functions ProteinDNATranslator(), 
 * ProteinRNATranslator(), DNAProteinTranslator(), and RNAProteinTranslator()
 * In general, the gnFastTranslator class should be used instead of this one.
 * @see gnFastTranslator
 */
class GNDLLEXPORT gnTranslator : public gnBaseFilter
{
public:

	static const gnTranslator *ProteinDNATranslator();
	static const gnTranslator *ProteinRNATranslator();
	static const gnTranslator *DNAProteinTranslator();
	static const gnTranslator *RNAProteinTranslator();

	enum gnTranslatorType{
		ProteinDNATranslatorType,
		ProteinRNATranslatorType,
		DNAProteinTranslatorType,
		RNAProteinTranslatorType,
	};
	
	gnTranslator();
	gnTranslator( gnTranslatorType t_type );
	gnTranslator( const gnTranslator& sf );

	gnTranslator* Clone() const;

	// gnSeqC 
	virtual gnSeqC Filter( const gnSeqC ch ) const;

	virtual void Filter( gnSeqC** seq, gnSeqI& len ) const;
	// std::string
	virtual void Filter( std::string &seq ) const;

	// Default gnSeqC
	/**
	 * Sets the default character to insert when no valid translation exists.
	 * @param ch1 The default character
	 */
	void SetDefaultChar( const gnSeqC ch1 );
	/**
	 * Gets the default character which is inserted into the destination
	 * sequence when no valid translation exists.
	 * @return The default character
	 */
	gnSeqC GetDefaultChar() const;
	/**
	 * Set whether the default character is inserted upon translation failure.
	 * @param use True if the default character should be used.
	 */
	void UseDefaultChar( const boolean use = true);
	/**
	 * Set the expected number of characters in each unit of translation.  For DNA to
	 * Protein, for instance, this is 3 because each codon is 3 characters of input.
	 * @param defaultInputWidth The expected input width
	 */
	void SetDefaultInputWidth( const uint32 defaultInputWidth);
	/**
	 * Get the expected number of characters in each unit of translation.  For DNA to
	 * Protein, for instance, this is 3 because each codon is 3 characters of input.
	 * @return The expected input width
	 */
	uint32 GetDefaultInputWidth() const;

	// fill map
	/**
	 * Set a translation mapping.  The first value is the input, the second value is the
	 * output.  During translation, any occurrences of the first std::string will be transformed
	 * to the second std::string.  A gnCompare object is used to facilitate ambiguity matching.
	 * @param input The input std::string
	 * @param output The std::string to output every time the input std::string is encountered
	 */
	void SetPair( const std::string& input, const std::string& output );
	/**
	 * Removes a translation mapping.  RemovePair removes the translation mapping
	 * corresponding to the given input std::string.
	 * @param input The input std::string to delete
	 */
	void RemovePair( const std::string& input );
	/**
	 * Set the comparator to use.  Ambiguous base matching is facilitated using a gnCompare
	 * object.  This must be set to allow sequences to be compared.
	 * @param comp A pointer to the gnCompare object to use.
	 */
	void SetCompare( const gnCompare* comp );
private:
	void CreateProteinDNATranslator();
	void CreateProteinRNATranslator();
	void CreateDNAProteinTranslator();
	void CreateRNAProteinTranslator();

	//for each entry in the input table there is a corresponding
	//entry in the output table.
	std::vector<std::string> m_inputTable, m_outputTable;

	const gnCompare* compare;
	
	boolean use_default;
	gnSeqC m_defaultChar;
	uint32 m_defaultInputWidth;
};//class gnTranslator

inline
gnTranslator* gnTranslator::Clone() const
{
	return new gnTranslator(*this);
}

inline
void gnTranslator::SetDefaultChar( const gnSeqC ch1 )
{
	m_defaultChar = ch1;
	use_default = true;
}
inline
gnSeqC gnTranslator::GetDefaultChar() const
{
	return m_defaultChar;
}

inline
void gnTranslator::UseDefaultChar(const boolean use)
{
	use_default = use;
}

inline
void gnTranslator::SetDefaultInputWidth( const uint32 defaultInputWidth){
	m_defaultInputWidth = defaultInputWidth;
}

inline
uint32 gnTranslator::GetDefaultInputWidth() const{
	return m_defaultInputWidth;
}

inline
void gnTranslator::SetCompare( const gnCompare* comp ){
	compare = comp;
}


}	// end namespace genome


#endif // _gnTranslator_h_
