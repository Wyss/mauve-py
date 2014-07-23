/////////////////////////////////////////////////////////////////////////////
// File:            gnTranslator.h
// Purpose:         Filter for all Sequences
// Description:     translates, converts sequence
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

#include "libGenome/gnTranslator.h"
#include "libGenome/gnCompare.h"
#include <string>

using namespace std;
using namespace std;
namespace genome {


//	static data access, avoids static initialization order fiasco
const gnTranslator *gnTranslator::ProteinDNATranslator(){
	const static gnTranslator* t_trans = new gnTranslator(ProteinDNATranslatorType);
	return t_trans;
}
const gnTranslator *gnTranslator::ProteinRNATranslator(){
	const static gnTranslator* t_trans = new gnTranslator(ProteinRNATranslatorType);
	return t_trans;
}
const gnTranslator *gnTranslator::DNAProteinTranslator(){
	const static gnTranslator* t_trans = new gnTranslator(DNAProteinTranslatorType);
	return t_trans;
}
const gnTranslator *gnTranslator::RNAProteinTranslator(){
	const static gnTranslator* t_trans = new gnTranslator(RNAProteinTranslatorType);
	return t_trans;
}

//	public:
gnTranslator::gnTranslator()
{
	use_default = false;
	m_defaultChar = 0;
	m_defaultInputWidth = 1;
}

gnTranslator::gnTranslator( const gnTranslator &sf )
{
	m_name = sf.m_name;
	use_default = sf.use_default;
	m_defaultChar = sf.m_defaultChar;
	compare = sf.compare;
	m_inputTable = sf.m_inputTable;
	m_outputTable = sf.m_outputTable;
	m_defaultInputWidth = sf.m_defaultInputWidth;
}
gnTranslator::gnTranslator( gnTranslatorType t_type )
{
	use_default = false;
	m_defaultChar = 0;
	switch(t_type){
		case ProteinDNATranslatorType:
			CreateProteinDNATranslator();
			break;
		case ProteinRNATranslatorType:
			CreateProteinRNATranslator();
			break;
		case DNAProteinTranslatorType:
			CreateDNAProteinTranslator();
			break;
		case RNAProteinTranslatorType:
			CreateRNAProteinTranslator();
			break;
	}
}

	// gnSeqC 
gnSeqC gnTranslator::Filter( const gnSeqC ch ) const{
	for(uint32 i=0; i < m_inputTable.size(); i++){
		if(m_inputTable[i].length() == 1)
			if(compare->Contains(m_inputTable[i][0], ch))
				return m_outputTable[i][0];
	}
	return m_defaultChar;
}

void gnTranslator::Filter( gnSeqC** seq, gnSeqI& len ) const{
	uint32 curpos = 0;
	string output;
	while(curpos < len){
		uint32 i=0;
		for(; i < m_inputTable.size(); i++){
			//don't compare if there aren't enough chars
			uint32 curlen = m_inputTable[i].length();
			if(len - curpos < curlen)
				continue;
			if(compare->Contains(m_inputTable[i].data(), *seq + curpos, curlen)){
				output += m_outputTable[i];
				curpos += curlen;
				break;
			}
		}
		if(i == m_inputTable.size()){
			//no match was found.  
			if(use_default)  //fill with the default char?
				output += m_defaultChar;
			curpos += m_defaultInputWidth;
		}
	}
	if(output.length() > len){
		delete[] *seq;
		*seq = new gnSeqC[output.length()];
	}
	len = output.length();
	memcpy(*seq, output.data(), len);
}
	// string
void gnTranslator::Filter( string &seq ) const{
	uint32 curpos = 0;
	uint32 len = seq.length();
	string output;
	while(curpos < len){
		uint32 i=0;
		for(; i < m_inputTable.size(); i++){
			//don't compare if there aren't enough chars
			uint32 curlen = m_inputTable[i].length();
			if(len - curpos < curlen)
				continue;
			if(compare->Contains(m_inputTable[i], seq.substr(curpos, curlen))){
				output += m_outputTable[i];
				curpos += curlen;
				break;
			}
		}
		if(i == m_inputTable.size()){
			//no match was found.  
			if(use_default)  //fill with the default char?
				output += m_defaultChar;
			curpos += m_defaultInputWidth;
		}
	}
	seq = output;
}

// fill map
void  gnTranslator::SetPair( const string& ch1, const string& ch2 )
{
	if(ch1.length() == 0)
		return;	//cant have an empty input, empty output is ok

	m_inputTable.push_back(ch1);
	m_outputTable.push_back(ch2);
}

void gnTranslator::RemovePair( const string& ch )
{
	for(uint32 i=0; i < m_inputTable.size(); i++){
		if(m_inputTable[i] == ch){
			m_inputTable.erase(m_inputTable.begin()+i);
			m_outputTable.erase(m_outputTable.begin()+i);
		}
	}
}

// standard comparators
void gnTranslator::CreateProteinDNATranslator(){
	SetName( "Protein to DNA Translator" );
	
	SetDefaultChar('X');
	SetCompare(gnCompare::ProteinSeqCompare());
	m_defaultInputWidth = 1;
	SetPair( "F", "TTY" );
	SetPair( "L", "YTX" );	//fix this somehow.  how?
	SetPair( "I", "ATH" );
	SetPair( "M", "ATG" );
	SetPair( "V", "GTX" );
	SetPair( "P", "CCX" );
	SetPair( "T", "ACX" );
	SetPair( "A", "GCX" );
	SetPair( "Y", "TAY" );
	SetPair( ".", "TRR" );//fix this somehow.  how?
	SetPair( "H", "CAY" );
	SetPair( "Q", "CAR" );
	SetPair( "N", "AAY" );
	SetPair( "K", "AAR" );
	SetPair( "D", "GAY" );
	SetPair( "E", "GAR" );
	SetPair( "C", "TGY" );
	SetPair( "W", "TGG" );
	SetPair( "G", "GGX" );

	SetPair( "S", "TCX" );
	SetPair( "S", "AGY");
	SetPair( "R", "CGX");
	SetPair( "R", "AGR");
}

void gnTranslator::CreateProteinRNATranslator(){
	SetName( "Protein to RNA Translator" );
	SetDefaultChar('X');
	SetCompare(gnCompare::ProteinSeqCompare());
	m_defaultInputWidth = 1;

	SetPair( "F", "UUY" );
	SetPair( "L", "YUX" );	//fix this somehow.  how?
	SetPair( "I", "AUH" );
	SetPair( "M", "AUG" );
	SetPair( "V", "GUX" );
	SetPair( "P", "CCX" );
	SetPair( "U", "ACX" );
	SetPair( "A", "GCX" );
	SetPair( "Y", "UAY" );
	SetPair( ".", "URR" );//fix this somehow.  how?
	SetPair( "H", "CAY" );
	SetPair( "Q", "CAR" );
	SetPair( "N", "AAY" );
	SetPair( "K", "AAR" );
	SetPair( "D", "GAY" );
	SetPair( "E", "GAR" );
	SetPair( "C", "UGY" );
	SetPair( "W", "UGG" );
	SetPair( "G", "GGX" );

	SetPair( "S", "UCX" );
	SetPair( "S", "AGY");
	SetPair( "R", "CGX");
	SetPair( "R", "AGR");
}

void gnTranslator::CreateDNAProteinTranslator(){
	SetName( "DNA to Protein Translator" );
	SetDefaultChar('X');
	SetCompare(gnCompare::DNASeqCompare());
	m_defaultInputWidth = 3;
	use_default = true;
	
	SetPair( "TTY", "F" );
	SetPair( "CTX", "L" );
	SetPair( "TTR", "L" );
	SetPair( "ATH", "I" );
	SetPair( "ATG", "M" );
	SetPair( "GTX", "V" );
	SetPair( "CCX", "P" );
	SetPair( "ACX", "T" );
	SetPair( "GCX", "A" );
	SetPair( "TAY", "Y" );
	SetPair( "TGG", "W" );
	SetPair( "TGA", "." );
	SetPair( "TAR", "." );
	SetPair( "CAY", "H" );
	SetPair( "CAR", "Q" );
	SetPair( "AAY", "N" );
	SetPair( "AAR", "K" );
	SetPair( "GAY", "D" );
	SetPair( "GAR", "E" );
	SetPair( "TGY", "C" );
	SetPair( "GGX", "G" );

	SetPair( "TCX", "S" );
	SetPair( "AGY", "S" );
	SetPair( "CGX", "R" );
	SetPair( "AGR", "R" );
	
	SetPair( "tty", "F" );
	SetPair( "ctx", "L" );
	SetPair( "ttr", "L" );
	SetPair( "ath", "I" );
	SetPair( "atg", "M" );
	SetPair( "gtx", "V" );
	SetPair( "ccx", "P" );
	SetPair( "acx", "T" );
	SetPair( "gcx", "A" );
	SetPair( "tay", "Y" );
	SetPair( "tgg", "W" );
	SetPair( "tga", "." );
	SetPair( "tar", "." );
	SetPair( "cay", "H" );
	SetPair( "car", "Q" );
	SetPair( "aay", "N" );
	SetPair( "aar", "K" );
	SetPair( "gay", "D" );
	SetPair( "gar", "E" );
	SetPair( "tgy", "C" );
	SetPair( "ggx", "G" );

	SetPair( "tcx", "S" );
	SetPair( "agy", "S" );
	SetPair( "cgx", "R" );
	SetPair( "agr", "R" );

}

void gnTranslator::CreateRNAProteinTranslator(){
	SetName( "RNA to Protein Translator" );
	SetDefaultChar('X');
	SetCompare(gnCompare::RNASeqCompare());
	m_defaultInputWidth = 3;
	use_default = true;
	
	SetPair( "UUY", "F" );
	SetPair( "CUX", "L" );
	SetPair( "UUR", "L" );
	SetPair( "AUH", "I" );
	SetPair( "AUG", "M" );
	SetPair( "GUX", "V" );
	SetPair( "CCX", "P" );
	SetPair( "ACX", "T" );
	SetPair( "GCX", "A" );
	SetPair( "UAY", "Y" );
	SetPair( "UGG", "W" );
	SetPair( "UGA", "." );
	SetPair( "UAR", "." );
	SetPair( "CAY", "H" );
	SetPair( "CAR", "Q" );
	SetPair( "AAY", "N" );
	SetPair( "AAR", "K" );
	SetPair( "GAY", "D" );
	SetPair( "GAR", "E" );
	SetPair( "UGY", "C" );
	SetPair( "GGX", "G" );

	SetPair( "UCX", "S" );
	SetPair( "AGY", "S" );
	SetPair( "CGX", "R" );
	SetPair( "AGR", "R" );


	SetPair( "uuy", "F" );
	SetPair( "cux", "L" );
	SetPair( "uur", "L" );
	SetPair( "auh", "I" );
	SetPair( "aug", "M" );
	SetPair( "gux", "V" );
	SetPair( "ccx", "P" );
	SetPair( "acx", "T" );
	SetPair( "gcx", "A" );
	SetPair( "uay", "Y" );
	SetPair( "ugg", "W" );
	SetPair( "uga", "." );
	SetPair( "uar", "." );
	SetPair( "cay", "H" );
	SetPair( "car", "Q" );
	SetPair( "aay", "N" );
	SetPair( "aar", "K" );
	SetPair( "gay", "D" );
	SetPair( "gar", "E" );
	SetPair( "ugy", "C" );
	SetPair( "ggx", "G" );

	SetPair( "ucx", "S" );
	SetPair( "agy", "S" );
	SetPair( "cgx", "R" );
	SetPair( "agr", "R" );
}

}	// end namespace genome

