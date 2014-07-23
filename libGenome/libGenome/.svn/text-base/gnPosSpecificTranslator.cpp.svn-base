/////////////////////////////////////////////////////////////////////////////
// File:            gnPosSpecificTranslator.cpp
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

#include "libGenome/gnPosSpecificTranslator.h"


using namespace std;
namespace genome {


//	static data access, avoids static initialization order fiasco
const gnPosSpecificTranslator *gnPosSpecificTranslator::ProteinDNATranslator(){
	const static gnPosSpecificTranslator* t_trans = new gnPosSpecificTranslator(ProteinDNATranslatorType);
	return t_trans;
}
const gnPosSpecificTranslator *gnPosSpecificTranslator::DNAProteinTranslator(){
	const static gnPosSpecificTranslator* t_trans = new gnPosSpecificTranslator(DNAProteinTranslatorType);
	return t_trans;
}

//	public:
gnPosSpecificTranslator::gnPosSpecificTranslator()
{
}

gnPosSpecificTranslator::gnPosSpecificTranslator( const gnPosSpecificTranslator &sf )
{
	m_name = sf.m_name;
}
gnPosSpecificTranslator::gnPosSpecificTranslator( gnTranslatorType t_type )
{
	m_type = t_type;
	switch(t_type){
		case ProteinDNATranslatorType:
			// not special done for this case
			filter = gnFastTranslator::ProteinDNATranslator();
			break;
		case DNAProteinTranslatorType:
			filter = gnFastTranslator::DNAProteinTranslator();
			break;
	}
}

// gnSeqC
gnSeqC gnPosSpecificTranslator::Filter( const gnSeqC ch ) const{
	return filter->Filter(ch);
}

void gnPosSpecificTranslator::Filter( gnSeqC** seq, gnSeqI& len ) const{
	return filter->Filter(seq, len);
}
	// string
void gnPosSpecificTranslator::Filter( string &seq ) const{
	switch(m_type){
		case ProteinDNATranslatorType:
			filter->Filter( seq );
			break;
		case DNAProteinTranslatorType:
			string first_codon = seq.substr(0,3);
			string original_sequence = seq;
			filter->Filter( seq );
			for(int charI = 0; charI < first_codon.length(); charI++ )
				first_codon[charI] = tolower(first_codon[charI]);
			if( first_codon == "ttg" || first_codon == "gtg" )
				seq[0] = 'M';
			// replace . with selenocysteine U where appropriate
			string::size_type dot_pos = seq.find( "." );
			while( dot_pos != string::npos ){
				if( dot_pos == seq.size() - 1 )
					break;  // dont convert the final residue!
				string dot_codon = original_sequence.substr( dot_pos * 3, 3 );
				dot_codon[0] = tolower( dot_codon[0] );
				dot_codon[1] = tolower( dot_codon[1] );
				dot_codon[2] = tolower( dot_codon[2] );
				if( dot_codon == "tga" )
					seq[ dot_pos ] = 'U';
				dot_pos = seq.find( ".", dot_pos + 1 );
			}
			break;
	}
}

}	// end namespace genome

