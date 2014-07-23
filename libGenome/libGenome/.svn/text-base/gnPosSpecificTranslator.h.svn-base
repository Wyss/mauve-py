/////////////////////////////////////////////////////////////////////////////
// File:            gnPosSpecificTranslator.h
// Purpose:         Special case ORF translation
// Description:     Used to translate sequences differently based on the position of
//					input characters.  Useful for tranlating genes because the first
//					codon is translated differently
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

#ifndef _gnPosSpecificTranslator_h_
#define _gnPosSpecificTranslator_h_


#include "libGenome/gnDefs.h"

#include <string>
#include <vector>
#include <map>
#include "libGenome/gnClone.h"
#include "libGenome/gnBaseFilter.h"
#include "libGenome/gnFastTranslator.h"


namespace genome {

/** 
 * Used to translate sequences differently based on the position of
 * input characters.  Useful for tranlating genes because the first
 * codon is translated differently.
 */
class GNDLLEXPORT gnPosSpecificTranslator : public gnBaseFilter
{
public:

	static const gnPosSpecificTranslator *ProteinDNATranslator();
	static const gnPosSpecificTranslator *DNAProteinTranslator();

	enum gnTranslatorType{
		ProteinDNATranslatorType,
		DNAProteinTranslatorType,
	};

	gnPosSpecificTranslator();
	gnPosSpecificTranslator( gnTranslatorType t_type );
	gnPosSpecificTranslator( const gnPosSpecificTranslator& sf );
	gnPosSpecificTranslator& operator= (const gnPosSpecificTranslator& sf);
	gnPosSpecificTranslator* Clone() const;
	

	// gnSeqC 
	virtual gnSeqC Filter( const gnSeqC ch ) const;

	virtual void Filter( gnSeqC** seq, gnSeqI& len ) const;
	// std::string
	virtual void Filter( std::string &seq ) const;

private:
	gnTranslatorType m_type;	// defines the type of translator this is
	const gnBaseFilter* filter;		// this is the filter used to do the translation
	void CreateProteinDNATranslator();
	void CreateDNAProteinTranslator();

};//class gnPosSpecificTranslator

inline
gnPosSpecificTranslator* gnPosSpecificTranslator::Clone() const
{
	return new gnPosSpecificTranslator(*this);
}


}	// end namespace genome

#endif // _gnPosSpecificTranslator_h_
