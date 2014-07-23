/////////////////////////////////////////////////////////////////////////////
// File:            gnFastTranslator.h
// Purpose:         Fast translator for all Sequences
// Description:     Caches translations of each possible sequence in a tree
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

#ifndef _gnFastTranslator_h_
#define _gnFastTranslator_h_


#include "libGenome/gnDefs.h"

#include <string>
#include <vector>
#include <map>
#include "libGenome/gnClone.h"
#include "libGenome/gnBaseFilter.h"
#include "libGenome/gnTranslator.h"


namespace genome {

class GNDLLEXPORT gnFastTranslator : public gnBaseFilter
{
public:

	static const gnFastTranslator *ProteinDNATranslator();
	static const gnFastTranslator *DNAProteinTranslator();

	enum gnTranslatorType{
		ProteinDNATranslatorType,
		DNAProteinTranslatorType,
	};

	gnFastTranslator();
	gnFastTranslator( gnTranslatorType t_type );
	gnFastTranslator( const gnFastTranslator& sf );
	gnFastTranslator& operator= (const gnFastTranslator& sf);
	gnFastTranslator* Clone() const;
	
	/**
	 * Queries the specified gnTranslator for every possible combination translation
	 * of the characters specified in the inputs std::string.  An input_width may be specified
	 * so that every possible combination of "input_width" characters in "inputs" will be
	 * cached.  This is useful for DNA to protein translations, for example.
	 * @param tranny The gnTranslator to cache.
	 * @param inputs The characters to cache from tranny.
	 * @param input_width The number of characters in each query to make to tranny.
	 */
	virtual void CacheTranslator(const gnTranslator* tranny, std::string inputs, const gnSeqI input_width);

	// gnSeqC 
	virtual gnSeqC Filter( const gnSeqC ch ) const;

	virtual void Filter( gnSeqC** seq, gnSeqI& len ) const;
	// std::string
	virtual void Filter( std::string &seq ) const;

	// Default gnSeqC
	void SetDefaultChar( const gnSeqC ch1 );
	gnSeqC GetDefaultChar() const;
	void UseDefaultChar( const boolean use = true);
	// fill map
	void SetPair( const std::string& ch1, const std::string& ch2 );
	void RemovePair( const std::string& ch );

private:

	void CreateProteinDNATranslator();
	void CreateDNAProteinTranslator();

	//map an input std::string to an output std::string
	std::map<std::string, std::string> m_transCache;
	const gnTranslator * m_translator;
	
	boolean use_default;
	gnSeqC m_defaultChar;
};//class gnFastTranslator

inline
gnFastTranslator* gnFastTranslator::Clone() const
{
	return new gnFastTranslator(*this);
}

inline
void gnFastTranslator::SetDefaultChar( const gnSeqC ch1 )
{
	m_defaultChar = ch1;
	use_default = true;
}
inline
gnSeqC gnFastTranslator::GetDefaultChar() const
{
	return m_defaultChar;
}

inline
void gnFastTranslator::UseDefaultChar(const boolean use)
{
	use_default = use;
}


}	// end namespace genome

#endif // _gnFastTranslator_h_
