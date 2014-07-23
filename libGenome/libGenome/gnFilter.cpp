/////////////////////////////////////////////////////////////////////////////
// File:            gnFilter.h
// Purpose:         Filter for all Sequences
// Description:     Filters sequences, translates, reverse complement, converts
//                   additions, etc.
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

#include "libGenome/gnFilter.h"
#include "libGenome/gnDebug.h"
#include <cstring>


using namespace std;
namespace genome {


//	public:
const gnFilter *gnFilter::alphabetCharacterFilter(){
	const static gnFilter* t_filt = new gnFilter(alphabetCharacterFilterType);
	return t_filt;
}

const gnFilter *gnFilter::numberCharacterFilter(){
	const static gnFilter* t_filt = new gnFilter(numberCharacterFilterType);
	return t_filt;
}


const gnFilter *gnFilter::proteinSeqFilter(){
	const static gnFilter* t_filt = new gnFilter(proteinSeqFilterType);
	return t_filt;
}

const gnFilter *gnFilter::basicDNASeqFilter(){
	const static gnFilter* t_filt = new gnFilter(basicDNASeqFilterType);
	return t_filt;
}

const gnFilter *gnFilter::fullDNASeqFilter(){
	const static gnFilter* t_filt = new gnFilter(fullDNASeqFilterType);
	return t_filt;
}

const gnFilter *gnFilter::basicRNASeqFilter(){
	const static gnFilter* t_filt = new gnFilter(basicRNASeqFilterType);
	return t_filt;
}

const gnFilter *gnFilter::fullRNASeqFilter(){
	const static gnFilter* t_filt = new gnFilter(fullRNASeqFilterType);
	return t_filt;
}

const gnFilter *gnFilter::DNAtoRNAFilter(){
	const static gnFilter* t_filt = new gnFilter(DNAtoRNAFilterType);
	return t_filt;
}

const gnFilter *gnFilter::RNAtoDNAFilter(){
	const static gnFilter* t_filt = new gnFilter(RNAtoDNAFilterType);
	return t_filt;
}

const gnFilter *gnFilter::DNAComplementFilter(){
	const static gnFilter* t_filt = new gnFilter(DNAComplementFilterType);
	return t_filt;
}

const gnFilter *gnFilter::RNAComplementFilter(){
	const static gnFilter* t_filt = new gnFilter(RNAComplementFilterType);
	return t_filt;
}


//	public:
gnFilter::gnFilter()
{
	m_defaultChar = 'n';
	m_rDefaultChar = 'n';
	for( gnSeqC i = 0; i < GNSEQC_MAX; ++i )
		m_pairArray[i] = NO_REVCOMP_CHAR;
}
gnFilter::gnFilter( const gnSeqC defaultChar, const gnSeqC rdefaultChar )
{
	m_defaultChar = defaultChar; 
	m_rDefaultChar = rdefaultChar;
	for( gnSeqC i = 0; i < GNSEQC_MAX; ++i )
		m_pairArray[i] = NO_REVCOMP_CHAR;
}

gnFilter::gnFilter( const gnFilter &sf )
{
	m_name = sf.m_name;
	for( gnSeqC i = 0; i < GNSEQC_MAX; ++i )
		m_pairArray[i] = sf.m_pairArray[i];
	m_defaultChar = sf.m_defaultChar;
	m_rDefaultChar = sf.m_rDefaultChar;
}

gnFilter::gnFilter( const gnFilterType f_type ){
	for( gnSeqC i = 0; i < GNSEQC_MAX; ++i )
		m_pairArray[i] = NO_REVCOMP_CHAR;
	switch(f_type){
		case alphabetCharacterFilterType:
			CreateAlphabetCharacterFilter();
			break;
		case numberCharacterFilterType:
			CreateNumberCharacterFilter();
			break;
		case proteinSeqFilterType:
			CreateProteinFilter();
			break;
		case basicDNASeqFilterType:
			CreateBasicDNAFilter();
			break;
		case fullDNASeqFilterType:
			CreateFullDNAFilter();
			break;
		case basicRNASeqFilterType:
			CreateBasicRNAFilter();
			break;
		case fullRNASeqFilterType:
			CreateFullRNAFilter();
			break;
		case DNAtoRNAFilterType:
			CreateDNAtoRNAFilter();
			break;
		case RNAtoDNAFilterType:
			CreateRNAtoDNAFilter();
			break;
		case DNAComplementFilterType:
			CreateDNAComplementFilter();
			break;		
		case RNAComplementFilterType:
			CreateRNAComplementFilter();
			break;
	}
}


gnFilter::~gnFilter()
{
}

inline
void gnFilter::Filter( gnSeqC** seq, gnSeqI& len ) const
{
	Array<gnSeqC> array_buf( len );
	gnSeqC* tmp = array_buf.data;
	gnSeqI c=0;
	for(uint32 i=0; i < len; i++)
		if(IsValid((*seq)[i]))
			tmp[c++] = m_pairArray[(*seq)[i]];
	len = c;
	memcpy(*seq, tmp, len);
}

void gnFilter::ReverseFilter( gnSeqC** seq, gnSeqI& len ) const
{
	gnSeqC tmp, dum;
	uint32 halfLen = len/2;
	uint32 end = len - 1;
	uint32 curB = 0;
	uint32 curE = end;
	for( uint32 i=0; i < halfLen ; ++i )
	{
		tmp = m_pairArray[(*seq)[i]];
		dum = m_pairArray[(*seq)[ end - i ]];
		if(dum != NO_REVCOMP_CHAR)
			(*seq)[ curB++ ] = dum;
		if(tmp != NO_REVCOMP_CHAR)
			(*seq)[ curE-- ] = tmp;
	}
	if(len&0x1){
		tmp = m_pairArray[(*seq)[halfLen]];
		if(tmp != NO_REVCOMP_CHAR)
			(*seq)[curB++] = tmp;
	}
	// now for the memmove
	if(curE >= curB){
		memmove(*seq+curB, *seq+curE+1, end - curE);
		len = end - curE + curB;
	}

}

void gnFilter::Filter( string &seq ) const
{
	gnSeqI c=0;
	for(uint32 i=0; i < seq.length(); i++)
		if(IsValid(seq[i]))
			seq[c++] = m_pairArray[seq[i]];
}

void gnFilter::ReverseFilter( string &seq ) const
{
	gnSeqC tmp, dum;
	uint32 halfLen = seq.length()/2;
	uint32 end = seq.length() - 1;
	uint32 curB = 0;
	uint32 curE = end;
	for( uint32 i=0; i < halfLen ; ++i )
	{
		tmp = m_pairArray[seq[i]];
		dum = m_pairArray[seq[ end - i ]];
		if(dum != NO_REVCOMP_CHAR)
			seq[ curB++ ] = dum;
		if(tmp != NO_REVCOMP_CHAR)
			seq[ curE-- ] = tmp;
	}
	if(seq.length()&0x1){
		tmp = m_pairArray[seq[halfLen]];
		if(tmp != NO_REVCOMP_CHAR)
			seq[curB++] = tmp;
	}
	// now for the memmove
	if(curE >= curB){
		seq.erase(curB, curE-curB);
	}
}

// standard filters
void gnFilter::CreateAlphabetCharacterFilter()
{
	SetDefaultChar( 0, 0 );
	SetName( "Alphabet Character Filter" );
	SetPair( 'A', 'a' );
	SetPair( 'B', 'b' );
	SetPair( 'C', 'c' );
	SetPair( 'D', 'd' );
	SetPair( 'E', 'e' );
	SetPair( 'F', 'f' );
	SetPair( 'G', 'g' );
	SetPair( 'H', 'h' );
	SetPair( 'I', 'i' );
	SetPair( 'J', 'j' );
	SetPair( 'K', 'k' );
	SetPair( 'L', 'l' );
	SetPair( 'M', 'm' );
	SetPair( 'N', 'n' );
	SetPair( 'O', 'o' );
	SetPair( 'P', 'p' );
	SetPair( 'Q', 'q' );
	SetPair( 'R', 'r' );
	SetPair( 'S', 's' );
	SetPair( 'T', 't' );
	SetPair( 'U', 'u' );
	SetPair( 'V', 'v' );
	SetPair( 'W', 'w' );
	SetPair( 'X', 'x' );
	SetPair( 'Y', 'y' );
	SetPair( 'Z', 'z' );
}

void gnFilter::CreateNumberCharacterFilter()
{
	SetDefaultChar( 0, 0 );
	SetName( "Number Character Filter" );
	SetSingle( '0' );
	SetSingle( '1' );
	SetSingle( '2' );
	SetSingle( '3' );
	SetSingle( '4' );
	SetSingle( '5' );
	SetSingle( '6' );
	SetSingle( '7' );
	SetSingle( '8' );
	SetSingle( '9' );
}

void gnFilter::CreateProteinFilter()
{
	SetDefaultChar( 'u', 'u' );
	SetName( "Protein Filter" );
	SetSingle( 'A' );
	SetSingle( 'R' );
	SetSingle( 'N' );
	SetSingle( 'D' );
	SetSingle( 'C' );
	SetSingle( 'Q' );
	SetSingle( 'E' );
	SetSingle( 'G' );
	SetSingle( 'H' );
	SetSingle( 'I' );
	SetSingle( 'L' );
	SetSingle( 'K' );
	SetSingle( 'M' );
	SetSingle( 'F' );
	SetSingle( 'P' );
	SetSingle( 'S' );
	SetSingle( 'T' );
	SetSingle( 'W' );
	SetSingle( 'Y' );
	SetSingle( 'V' );
	
	SetSingle( 'a' );
	SetSingle( 'r' );
	SetSingle( 'n' );
	SetSingle( 'd' );
	SetSingle( 'c' );
	SetSingle( 'q' );
	SetSingle( 'e' );
	SetSingle( 'g' );
	SetSingle( 'h' );
	SetSingle( 'i' );
	SetSingle( 'l' );
	SetSingle( 'k' );
	SetSingle( 'm' );
	SetSingle( 'f' );
	SetSingle( 'p' );
	SetSingle( 's' );
	SetSingle( 't' );
	SetSingle( 'w' );
	SetSingle( 'y' );
	SetSingle( 'v' );
}

void gnFilter::CreateBasicDNAFilter()
{
	SetDefaultChar( 'n', 'n' );
	SetName( "Basic DNA Filter" );
	SetSingle( 'a' );
	SetSingle( 'c' );
	SetSingle( 'g' );
	SetSingle( 't' );
	SetSingle( 'A' );
	SetSingle( 'C' );
	SetSingle( 'G' );
	SetSingle( 'T' );
	SetSingle( 'n' );
	SetSingle( 'N' );
	SetSingle( 'x' );
	SetSingle( 'X' );
	SetSingle( '-' );
}
void gnFilter::CreateFullDNAFilter()
{	
	SetDefaultChar( 'n', 'n' );
	SetName( "Full DNA Filter" );
	SetSingle( 'a' );
	SetSingle( 'c' );
	SetSingle( 'g' );
	SetSingle( 't' );
	SetSingle( 'A' );
	SetSingle( 'C' );
	SetSingle( 'G' );
	SetSingle( 'T' );
	SetSingle( 'r' );
	SetSingle( 'y' );
	SetSingle( 'k' );
	SetSingle( 'm' );
	SetSingle( 'b' );
	SetSingle( 'v' );
	SetSingle( 'd' );
	SetSingle( 'h' );
	SetSingle( 'R' );
	SetSingle( 'Y' );
	SetSingle( 'K' );
	SetSingle( 'M' );
	SetSingle( 'B' );
	SetSingle( 'V' );
	SetSingle( 'D' );
	SetSingle( 'H' );
	SetSingle( 's' );
	SetSingle( 'S' );
	SetSingle( 'w' );
	SetSingle( 'W' );
	SetSingle( 'n' );
	SetSingle( 'N' );
	SetSingle( 'x' );
	SetSingle( 'X' );
	SetSingle( '-' );
}
void gnFilter::CreateBasicRNAFilter()
{
	SetDefaultChar( 'n', 'n' );
	SetName( "Basic RNA Filter" );
	SetSingle( 'a' );
	SetSingle( 'c' );
	SetSingle( 'g' );
	SetSingle( 'u' );
	SetSingle( 'A' );
	SetSingle( 'C' );
	SetSingle( 'G' );
	SetSingle( 'U' );
	SetSingle( 'n' );
	SetSingle( 'N' );
	SetSingle( '-' );
}
void gnFilter::CreateFullRNAFilter()
{
	SetDefaultChar( 'n', 'n' );
	SetName( "Full RNA Filter" );
	SetSingle( 'a' );
	SetSingle( 'c' );
	SetSingle( 'g' );
	SetSingle( 'u' );
	SetSingle( 'A' );
	SetSingle( 'C' );
	SetSingle( 'G' );
	SetSingle( 'U' );
	SetSingle( 'r' );
	SetSingle( 'y' );
	SetSingle( 'k' );
	SetSingle( 'm' );
	SetSingle( 'b' );
	SetSingle( 'v' );
	SetSingle( 'd' );
	SetSingle( 'h' );
	SetSingle( 'R' );
	SetSingle( 'Y' );
	SetSingle( 'K' );
	SetSingle( 'M' );
	SetSingle( 'B' );
	SetSingle( 'V' );
	SetSingle( 'D' );
	SetSingle( 'H' );
	SetSingle( 's' );
	SetSingle( 'S' );
	SetSingle( 'w' );
	SetSingle( 'W' );
	SetSingle( 'n' );
	SetSingle( 'N' );
	SetSingle( '-' );
}


void gnFilter::CreateDNAtoRNAFilter(){
	SetDefaultChar( 'n', 'n' );
	SetName( "Full DNA to RNA Filter" );
	SetSingle( 'a' );
	SetSingle( 'c' );
	SetSingle( 'g' );
	SetPair( 't', 'u' );
	SetSingle( 'A' );
	SetSingle( 'C' );
	SetSingle( 'G' );
	SetPair( 'T', 'U' );
	SetSingle( 'r' );
	SetSingle( 'y' );
	SetSingle( 'k' );
	SetSingle( 'm' );
	SetSingle( 'b' );
	SetSingle( 'v' );
	SetSingle( 'd' );
	SetSingle( 'h' );
	SetSingle( 'R' );
	SetSingle( 'Y' );
	SetSingle( 'K' );
	SetSingle( 'M' );
	SetSingle( 'B' );
	SetSingle( 'V' );
	SetSingle( 'D' );
	SetSingle( 'H' );
	SetSingle( 's' );
	SetSingle( 'S' );
	SetSingle( 'w' );
	SetSingle( 'W' );
	SetSingle( 'n' );
	SetSingle( 'N' );
	SetSingle( '-' );
}

void gnFilter::CreateRNAtoDNAFilter(){
	SetDefaultChar( 'n', 'n' );
	SetName( "Full RNA to DNA Filter" );
	SetSingle( 'a' );
	SetSingle( 'c' );
	SetSingle( 'g' );
	SetPair( 'u', 't' );
	SetSingle( 'A' );
	SetSingle( 'C' );
	SetSingle( 'G' );
	SetPair( 'U', 'T' );
	SetSingle( 'r' );
	SetSingle( 'y' );
	SetSingle( 'k' );
	SetSingle( 'm' );
	SetSingle( 'b' );
	SetSingle( 'v' );
	SetSingle( 'd' );
	SetSingle( 'h' );
	SetSingle( 'R' );
	SetSingle( 'Y' );
	SetSingle( 'K' );
	SetSingle( 'M' );
	SetSingle( 'B' );
	SetSingle( 'V' );
	SetSingle( 'D' );
	SetSingle( 'H' );
	SetSingle( 's' );
	SetSingle( 'S' );
	SetSingle( 'w' );
	SetSingle( 'W' );
	SetSingle( 'n' );
	SetSingle( 'N' );
	SetSingle( '-' );
}

void gnFilter::CreateDNAComplementFilter(){
	SetDefaultChar( 'n', 'n' );
	SetName( "Full DNA Complement Filter" );
	SetPair( 'a', 't' );
	SetPair( 'A', 'T' );
	SetPair( 't', 'a' );
	SetPair( 'T', 'A' );
	SetPair( 'c', 'g' );
	SetPair( 'C', 'G' );
	SetPair( 'g', 'c' );
	SetPair( 'G', 'C' );
	SetPair( 'r', 'y' );
	SetPair( 'R', 'Y' );
	SetPair( 'y', 'r' );
	SetPair( 'Y', 'R' );
	SetPair( 'k', 'm' );
	SetPair( 'K', 'M' );
	SetPair( 'm', 'k' );
	SetPair( 'M', 'K' );
	SetSingle( 's' );
	SetSingle( 'S' );
	SetSingle( 'w' );
	SetSingle( 'W' );
	SetPair( 'b', 'v' );
	SetPair( 'B', 'V' );
	SetPair( 'v', 'b' );
	SetPair( 'V', 'B' );
	SetPair( 'd', 'h' );
	SetPair( 'D', 'H' );
	SetPair( 'h', 'd' );
	SetPair( 'H', 'D' );
	SetSingle( 'n' );
	SetSingle( 'N' );
	SetSingle( 'x' );
	SetSingle( 'X' );
	SetSingle( '-' );
}

void gnFilter::CreateRNAComplementFilter(){
	SetDefaultChar( 'n', 'n' );
	SetName( "Full RNA Complement Filter" );
	SetPair( 'a', 'u' );
	SetPair( 'A', 'U' );
	SetPair( 'u', 'a' );
	SetPair( 'U', 'A' );
	SetPair( 'c', 'g' );
	SetPair( 'C', 'G' );
	SetPair( 'g', 'c' );
	SetPair( 'G', 'C' );
	SetPair( 'r', 'y' );
	SetPair( 'R', 'Y' );
	SetPair( 'y', 'r' );
	SetPair( 'Y', 'R' );
	SetPair( 'k', 'm' );
	SetPair( 'K', 'M' );
	SetPair( 'm', 'k' );
	SetPair( 'M', 'K' );
	SetSingle( 's' );
	SetSingle( 'S' );
	SetSingle( 'w' );
	SetSingle( 'W' );
	SetPair( 'b', 'v' );
	SetPair( 'B', 'V' );
	SetPair( 'v', 'b' );
	SetPair( 'V', 'B' );
	SetPair( 'd', 'h' );
	SetPair( 'D', 'H' );
	SetPair( 'h', 'd' );
	SetPair( 'H', 'D' );
	SetSingle( 'n' );
	SetSingle( 'N' );
	SetSingle( '-' );
}

}	// end namespace genome

