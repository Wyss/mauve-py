/////////////////////////////////////////////////////////////////////////////
// File:            gnCompare.cpp
// Purpose:         Coparator for all Sequences
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

#include "libGenome/gnCompare.h"
#include <cstring>

//	public:


using namespace std;
namespace genome {


const gnCompare *gnCompare::ProteinSeqCompare(){
	const static gnCompare* t_comp = new gnCompare(ProteinSeqCompareType);
	return t_comp;
}
const gnCompare *gnCompare::DNASeqCompare(){
	const static gnCompare* t_comp = new gnCompare(DNASeqCompareType);
	return t_comp;
}
const gnCompare *gnCompare::RNASeqCompare(){
	const static gnCompare* t_comp = new gnCompare(RNASeqCompareType);
	return t_comp;
}

gnCompare::gnCompare( const gnCompareType c_type ){
	for( gnSeqC i = 0; i < GNSEQC_MAX; ++i ){
		m_pairArray[i] = new gnSeqC[1];
		m_pairArray[i][0] = 0;
		m_containArray[i] = new gnSeqC[1];
		m_containArray[i][0] = 0;
	}
	switch(c_type){
		case ProteinSeqCompareType:
			CreateProteinComparator();
			break;
		case DNASeqCompareType:
			CreateDNAComparator();
			break;
		case RNASeqCompareType:
			CreateRNAComparator();
			break;
	}
}

boolean gnCompare::Contains( gnSeqC ch, gnSeqC ch2, boolean case_sensitive) const
{
	if(!case_sensitive){
		ch = toupper(ch);
		ch2 = toupper(ch2);
	}
	if(strchr(m_containArray[ch], ch2) == 0)
		return false;
	return true;
}

boolean gnCompare::Contains( const gnSeqC* seq, const gnSeqC* seq2, const uint32 len, boolean case_sensitive ) const{
	for( uint32 i=0; i < len ; ++i )
		if(!Contains(seq[i], seq2[i], case_sensitive))
			return false;
	return true;
}

boolean gnCompare::Contains( const string &seq, const string &seq2, boolean case_sensitive) const
{
	gnSeqI shorter_len = seq.length() < seq2.length() ? seq.length() : seq2.length();
	return Contains( (gnSeqC*)seq.data(), (gnSeqC*)seq2.data(), shorter_len, case_sensitive );
}
//	public:
gnCompare::gnCompare()
{
	for( gnSeqC i = 0; i < GNSEQC_MAX; ++i ){
		m_pairArray[i] = new gnSeqC[1];
		m_pairArray[i][0] = 0;
		m_containArray[i] = new gnSeqC[1];
		m_containArray[i][0] = 0;
	}
}
/*
gnCompare::gnCompare( const gnCompare &sf )
{
	m_name = sf.m_name;
	for( gnSeqC i = 0; i < GNSEQC_MAX; ++i ){
		m_pairArray[i] = new gnSeqC[strlen(sf.m_pairArray[i])+1];
		strcpy(m_pairArray[i], sf.m_pairArray[i]);
		m_containArray[i] = new gnSeqC[strlen(sf.m_containArray[i])+1];
		strcpy(m_containArray[i], sf.m_containArray[i]);
	}
}
*/
gnCompare::~gnCompare()
{
	for( gnSeqC i = 0; i < GNSEQC_MAX; ++i ){
		delete[] m_pairArray[i];
		delete[] m_containArray[i];
	}
}


void gnCompare::AddArrayEntry( gnSeqC *array[GNSEQC_MAX], const gnSeqC ch, const gnSeqC ch2){
	uint32 curlen = strlen(array[ch]);
	gnSeqC* tmp = new gnSeqC[curlen + 2];
	strcpy(tmp, array[ch]);
	tmp[curlen] = ch2;
	tmp[curlen+1] = 0;
	delete[] array[ch];
	array[ch] = tmp;
}

void gnCompare::DelArrayEntry( gnSeqC *array[GNSEQC_MAX], const gnSeqC ch, const gnSeqC ch2){
	//check that the pair exists
	gnSeqC* loc = strchr(m_containArray[ch], ch2);
	uint32 count = 0;
	while(loc != NULL){
		count++;
		loc = strchr(loc+1, ch2);
	}
	if(count == 0)
		return;

	uint32 curlen = strlen(array[ch]);
	gnSeqC* tmp = new gnSeqC[curlen - count];
	uint32 tmppos = 0;
	for(uint32 i=0; i < curlen; i++)
		if(m_containArray[ch][i] != ch2)
			tmp[tmppos++] = m_containArray[ch][i];
	tmp[tmppos] = 0;
	delete[] array[ch];
	array[ch] = tmp;
}


void gnCompare::CreateProteinComparator()
{
	SetName( "Protein Comparator" );
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
	SetSingle( '.' );
	
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

void gnCompare::CreateDNAComparator()
{
	SetName( "Full DNA Comparator" );

	SetSingle( 'a' );
	SetSingle( 'c' );
	SetSingle( 'g' );
	SetSingle( 't' );
	SetSingle( 'r' );
	SetSingle( 'k' );
	SetSingle( 's' );
	SetSingle( 'm' );
	SetSingle( 'y' );
	SetSingle( 'w' );
	SetSingle( 'b' );
	SetSingle( 'v' );
	SetSingle( 'd' );
	SetSingle( 'h' );
	SetSingle( 'n' );
	SetSingle( 'x' );

	SetSingle( 'A' );
	SetSingle( 'C' );
	SetSingle( 'G' );
	SetSingle( 'T' );
	SetSingle( 'R' );
	SetSingle( 'K' );
	SetSingle( 'S' );
	SetSingle( 'M' );
	SetSingle( 'Y' );
	SetSingle( 'W' );
	SetSingle( 'B' );
	SetSingle( 'V' );
	SetSingle( 'D' );
	SetSingle( 'H' );
	SetSingle( 'N' );
	SetSingle( 'X' );

	SetPair( 'g', 'r' );
	SetPair( 'g', 'k' );
	SetPair( 'g', 's' );
	SetPair( 'g', 'd' );
	SetPair( 'g', 'v' );
	SetPair( 'g', 'b' );
	SetPair( 'g', 'x' );
	SetPair( 'g', 'n' );
	SetPair( 'G', 'R' );
	SetPair( 'G', 'K' );
	SetPair( 'G', 'S' );
	SetPair( 'G', 'D' );
	SetPair( 'G', 'V' );
	SetPair( 'G', 'B' );
	SetPair( 'G', 'X' );
	SetPair( 'G', 'N' );
	SetPair( 'a', 'r' );
	SetPair( 'a', 'w' );
	SetPair( 'a', 'm' );
	SetPair( 'a', 'd' );
	SetPair( 'a', 'v' );
	SetPair( 'a', 'h' );
	SetPair( 'a', 'x' );
	SetPair( 'a', 'n' );
	SetPair( 'A', 'R' );
	SetPair( 'A', 'W' );
	SetPair( 'A', 'M' );
	SetPair( 'A', 'D' );
	SetPair( 'A', 'V' );
	SetPair( 'A', 'H' );
	SetPair( 'A', 'B' );
	SetPair( 'A', 'X' );
	SetPair( 'A', 'N' );
	SetPair( 'c', 's' );
	SetPair( 'c', 'm' );
	SetPair( 'c', 'y' );
	SetPair( 'c', 'v' );
	SetPair( 'c', 'b' );
	SetPair( 'c', 'h' );
	SetPair( 'c', 'x' );
	SetPair( 'c', 'n' );
	SetPair( 'C', 'S' );
	SetPair( 'C', 'M' );
	SetPair( 'C', 'Y' );
	SetPair( 'C', 'V' );
	SetPair( 'C', 'B' );
	SetPair( 'C', 'H' );
	SetPair( 'C', 'X' );
	SetPair( 'C', 'N' );
	SetPair( 't', 'k' );
	SetPair( 't', 'w' );
	SetPair( 't', 'y' );
	SetPair( 't', 'd' );
	SetPair( 't', 'b' );
	SetPair( 't', 'h' );
	SetPair( 't', 'x' );
	SetPair( 't', 'n' );
	SetPair( 'T', 'K' );
	SetPair( 'T', 'W' );
	SetPair( 'T', 'Y' );
	SetPair( 'T', 'D' );
	SetPair( 'T', 'B' );
	SetPair( 'T', 'H' );
	SetPair( 'T', 'X' );
	SetPair( 'T', 'N' );

	SetPair( 'r', 'x' );
	SetPair( 'y', 'x' );
	SetPair( 'w', 'x' );
	SetPair( 's', 'x' );
	SetPair( 'k', 'x' );
	SetPair( 'm', 'x' );
	SetPair( 'b', 'x' );
	SetPair( 'd', 'x' );
	SetPair( 'h', 'x' );
	SetPair( 'v', 'x' );
	SetPair( 'n', 'x' );

	SetPair( 'R', 'X' );
	SetPair( 'Y', 'X' );
	SetPair( 'W', 'X' );
	SetPair( 'S', 'X' );
	SetPair( 'K', 'X' );
	SetPair( 'M', 'X' );
	SetPair( 'B', 'X' );
	SetPair( 'D', 'X' );
	SetPair( 'H', 'X' );
	SetPair( 'V', 'X' );
	SetPair( 'N', 'X' );

	SetPair( 'r', 'n' );
	SetPair( 'y', 'n' );
	SetPair( 'w', 'n' );
	SetPair( 's', 'n' );
	SetPair( 'k', 'n' );
	SetPair( 'm', 'n' );
	SetPair( 'b', 'n' );
	SetPair( 'd', 'n' );
	SetPair( 'h', 'n' );
	SetPair( 'v', 'n' );
	SetPair( 'x', 'n' );

	SetPair( 'R', 'N' );
	SetPair( 'Y', 'N' );
	SetPair( 'W', 'N' );
	SetPair( 'S', 'N' );
	SetPair( 'K', 'N' );
	SetPair( 'M', 'N' );
	SetPair( 'B', 'N' );
	SetPair( 'D', 'N' );
	SetPair( 'H', 'N' );
	SetPair( 'V', 'N' );
	SetPair( 'X', 'N' );

	SetPair( 'k', 'd' );
	SetPair( 'k', 'b' );
	SetPair( 'r', 'd' );
	SetPair( 'r', 'v' );
	SetPair( 's', 'v' );
	SetPair( 's', 'b' );
	SetPair( 'm', 'v' );
	SetPair( 'm', 'h' );
	SetPair( 'w', 'd' );
	SetPair( 'w', 'h' );
	SetPair( 'y', 'b' );
	SetPair( 'y', 'h' );

	SetPair( 'K', 'D' );
	SetPair( 'K', 'B' );
	SetPair( 'R', 'D' );
	SetPair( 'R', 'V' );
	SetPair( 'S', 'V' );
	SetPair( 'S', 'B' );
	SetPair( 'M', 'V' );
	SetPair( 'M', 'H' );
	SetPair( 'W', 'D' );
	SetPair( 'W', 'H' );
	SetPair( 'Y', 'B' );
	SetPair( 'Y', 'H' );
//set the containment properties
	SetContained( 'g', 'r' );
	SetContained( 'g', 'k' );
	SetContained( 'g', 's' );
	SetContained( 'g', 'd' );
	SetContained( 'g', 'v' );
	SetContained( 'g', 'b' );
	SetContained( 'g', 'x' );
	SetContained( 'g', 'n' );
	SetContained( 'G', 'R' );
	SetContained( 'G', 'K' );
	SetContained( 'G', 'S' );
	SetContained( 'G', 'D' );
	SetContained( 'G', 'V' );
	SetContained( 'G', 'B' );
	SetContained( 'G', 'X' );
	SetContained( 'G', 'N' );
	SetContained( 'a', 'r' );
	SetContained( 'a', 'w' );
	SetContained( 'a', 'm' );
	SetContained( 'a', 'd' );
	SetContained( 'a', 'v' );
	SetContained( 'a', 'h' );
	SetContained( 'a', 'x' );
	SetContained( 'a', 'n' );
	SetContained( 'A', 'R' );
	SetContained( 'A', 'W' );
	SetContained( 'A', 'M' );
	SetContained( 'A', 'D' );
	SetContained( 'A', 'V' );
	SetContained( 'A', 'H' );
	SetContained( 'A', 'B' );
	SetContained( 'A', 'X' );
	SetContained( 'A', 'N' );
	SetContained( 'c', 's' );
	SetContained( 'c', 'm' );
	SetContained( 'c', 'y' );
	SetContained( 'c', 'v' );
	SetContained( 'c', 'b' );
	SetContained( 'c', 'h' );
	SetContained( 'c', 'x' );
	SetContained( 'c', 'n' );
	SetContained( 'C', 'S' );
	SetContained( 'C', 'M' );
	SetContained( 'C', 'Y' );
	SetContained( 'C', 'V' );
	SetContained( 'C', 'B' );
	SetContained( 'C', 'H' );
	SetContained( 'C', 'X' );
	SetContained( 'C', 'N' );
	SetContained( 't', 'k' );
	SetContained( 't', 'w' );
	SetContained( 't', 'y' );
	SetContained( 't', 'd' );
	SetContained( 't', 'b' );
	SetContained( 't', 'h' );
	SetContained( 't', 'x' );
	SetContained( 't', 'n' );
	SetContained( 'T', 'K' );
	SetContained( 'T', 'W' );
	SetContained( 'T', 'Y' );
	SetContained( 'T', 'D' );
	SetContained( 'T', 'B' );
	SetContained( 'T', 'H' );
	SetContained( 'T', 'X' );
	SetContained( 'T', 'N' );

	SetContained( 'r', 'x' );
	SetContained( 'y', 'x' );
	SetContained( 'w', 'x' );
	SetContained( 's', 'x' );
	SetContained( 'k', 'x' );
	SetContained( 'm', 'x' );
	SetContained( 'b', 'x' );
	SetContained( 'd', 'x' );
	SetContained( 'h', 'x' );
	SetContained( 'v', 'x' );
	SetContained( 'n', 'x' );

	SetContained( 'R', 'X' );
	SetContained( 'Y', 'X' );
	SetContained( 'W', 'X' );
	SetContained( 'S', 'X' );
	SetContained( 'K', 'X' );
	SetContained( 'M', 'X' );
	SetContained( 'B', 'X' );
	SetContained( 'D', 'X' );
	SetContained( 'H', 'X' );
	SetContained( 'V', 'X' );
	SetContained( 'N', 'X' );

	SetContained( 'r', 'n' );
	SetContained( 'y', 'n' );
	SetContained( 'w', 'n' );
	SetContained( 's', 'n' );
	SetContained( 'k', 'n' );
	SetContained( 'm', 'n' );
	SetContained( 'b', 'n' );
	SetContained( 'd', 'n' );
	SetContained( 'h', 'n' );
	SetContained( 'v', 'n' );
	SetContained( 'x', 'n' );

	SetContained( 'R', 'N' );
	SetContained( 'Y', 'N' );
	SetContained( 'W', 'N' );
	SetContained( 'S', 'N' );
	SetContained( 'K', 'N' );
	SetContained( 'M', 'N' );
	SetContained( 'B', 'N' );
	SetContained( 'D', 'N' );
	SetContained( 'H', 'N' );
	SetContained( 'V', 'N' );
	SetContained( 'X', 'N' );

	SetContained( 'k', 'd' );
	SetContained( 'k', 'b' );
	SetContained( 'r', 'd' );
	SetContained( 'r', 'v' );
	SetContained( 's', 'v' );
	SetContained( 's', 'b' );
	SetContained( 'm', 'v' );
	SetContained( 'm', 'h' );
	SetContained( 'w', 'd' );
	SetContained( 'w', 'h' );
	SetContained( 'y', 'b' );
	SetContained( 'y', 'h' );

	SetContained( 'K', 'D' );
	SetContained( 'K', 'B' );
	SetContained( 'R', 'D' );
	SetContained( 'R', 'V' );
	SetContained( 'S', 'V' );
	SetContained( 'S', 'B' );
	SetContained( 'M', 'V' );
	SetContained( 'M', 'H' );
	SetContained( 'W', 'D' );
	SetContained( 'W', 'H' );
	SetContained( 'Y', 'B' );
	SetContained( 'Y', 'H' );

}

void gnCompare::CreateRNAComparator()
{
	SetName( "Full RNA Comparator" );

	SetSingle( 'a' );
	SetSingle( 'c' );
	SetSingle( 'g' );
	SetSingle( 'u' );
	SetSingle( 'r' );
	SetSingle( 'k' );
	SetSingle( 's' );
	SetSingle( 'm' );
	SetSingle( 'y' );
	SetSingle( 'w' );
	SetSingle( 'b' );
	SetSingle( 'v' );
	SetSingle( 'd' );
	SetSingle( 'h' );
	SetSingle( 'n' );
	SetSingle( 'x' );

	SetSingle( 'A' );
	SetSingle( 'C' );
	SetSingle( 'G' );
	SetSingle( 'U' );
	SetSingle( 'R' );
	SetSingle( 'K' );
	SetSingle( 'S' );
	SetSingle( 'M' );
	SetSingle( 'Y' );
	SetSingle( 'W' );
	SetSingle( 'B' );
	SetSingle( 'V' );
	SetSingle( 'D' );
	SetSingle( 'H' );
	SetSingle( 'N' );
	SetSingle( 'X' );

	SetPair( 'g', 'r' );
	SetPair( 'g', 'k' );
	SetPair( 'g', 's' );
	SetPair( 'g', 'd' );
	SetPair( 'g', 'v' );
	SetPair( 'g', 'b' );
	SetPair( 'g', 'x' );
	SetPair( 'g', 'n' );
	SetPair( 'G', 'R' );
	SetPair( 'G', 'K' );
	SetPair( 'G', 'S' );
	SetPair( 'G', 'D' );
	SetPair( 'G', 'V' );
	SetPair( 'G', 'B' );
	SetPair( 'G', 'X' );
	SetPair( 'G', 'N' );
	SetPair( 'a', 'r' );
	SetPair( 'a', 'w' );
	SetPair( 'a', 'm' );
	SetPair( 'a', 'd' );
	SetPair( 'a', 'v' );
	SetPair( 'a', 'h' );
	SetPair( 'a', 'x' );
	SetPair( 'a', 'n' );
	SetPair( 'A', 'R' );
	SetPair( 'A', 'W' );
	SetPair( 'A', 'M' );
	SetPair( 'A', 'D' );
	SetPair( 'A', 'V' );
	SetPair( 'A', 'H' );
	SetPair( 'A', 'B' );
	SetPair( 'A', 'X' );
	SetPair( 'A', 'N' );
	SetPair( 'c', 's' );
	SetPair( 'c', 'm' );
	SetPair( 'c', 'y' );
	SetPair( 'c', 'v' );
	SetPair( 'c', 'b' );
	SetPair( 'c', 'h' );
	SetPair( 'c', 'x' );
	SetPair( 'c', 'n' );
	SetPair( 'C', 'S' );
	SetPair( 'C', 'M' );
	SetPair( 'C', 'Y' );
	SetPair( 'C', 'V' );
	SetPair( 'C', 'B' );
	SetPair( 'C', 'H' );
	SetPair( 'C', 'X' );
	SetPair( 'C', 'N' );
	SetPair( 'u', 'k' );
	SetPair( 'u', 'w' );
	SetPair( 'u', 'y' );
	SetPair( 'u', 'd' );
	SetPair( 'u', 'b' );
	SetPair( 'u', 'h' );
	SetPair( 'u', 'x' );
	SetPair( 'u', 'n' );
	SetPair( 'U', 'K' );
	SetPair( 'U', 'W' );
	SetPair( 'U', 'Y' );
	SetPair( 'U', 'D' );
	SetPair( 'U', 'B' );
	SetPair( 'U', 'H' );
	SetPair( 'U', 'X' );
	SetPair( 'U', 'N' );

	SetPair( 'r', 'x' );
	SetPair( 'y', 'x' );
	SetPair( 'w', 'x' );
	SetPair( 's', 'x' );
	SetPair( 'k', 'x' );
	SetPair( 'm', 'x' );
	SetPair( 'b', 'x' );
	SetPair( 'd', 'x' );
	SetPair( 'h', 'x' );
	SetPair( 'v', 'x' );
	SetPair( 'n', 'x' );

	SetPair( 'R', 'X' );
	SetPair( 'Y', 'X' );
	SetPair( 'W', 'X' );
	SetPair( 'S', 'X' );
	SetPair( 'K', 'X' );
	SetPair( 'M', 'X' );
	SetPair( 'B', 'X' );
	SetPair( 'D', 'X' );
	SetPair( 'H', 'X' );
	SetPair( 'V', 'X' );
	SetPair( 'N', 'X' );

	SetPair( 'r', 'n' );
	SetPair( 'y', 'n' );
	SetPair( 'w', 'n' );
	SetPair( 's', 'n' );
	SetPair( 'k', 'n' );
	SetPair( 'm', 'n' );
	SetPair( 'b', 'n' );
	SetPair( 'd', 'n' );
	SetPair( 'h', 'n' );
	SetPair( 'v', 'n' );
	SetPair( 'x', 'n' );

	SetPair( 'R', 'N' );
	SetPair( 'Y', 'N' );
	SetPair( 'W', 'N' );
	SetPair( 'S', 'N' );
	SetPair( 'K', 'N' );
	SetPair( 'M', 'N' );
	SetPair( 'B', 'N' );
	SetPair( 'D', 'N' );
	SetPair( 'H', 'N' );
	SetPair( 'V', 'N' );
	SetPair( 'X', 'N' );

	SetPair( 'k', 'd' );
	SetPair( 'k', 'b' );
	SetPair( 'r', 'd' );
	SetPair( 'r', 'v' );
	SetPair( 's', 'v' );
	SetPair( 's', 'b' );
	SetPair( 'm', 'v' );
	SetPair( 'm', 'h' );
	SetPair( 'w', 'd' );
	SetPair( 'w', 'h' );
	SetPair( 'y', 'b' );
	SetPair( 'y', 'h' );

	SetPair( 'K', 'D' );
	SetPair( 'K', 'B' );
	SetPair( 'R', 'D' );
	SetPair( 'R', 'V' );
	SetPair( 'S', 'V' );
	SetPair( 'S', 'B' );
	SetPair( 'M', 'V' );
	SetPair( 'M', 'H' );
	SetPair( 'W', 'D' );
	SetPair( 'W', 'H' );
	SetPair( 'Y', 'B' );
	SetPair( 'Y', 'H' );
//set the containment properties
	SetContained( 'g', 'r' );
	SetContained( 'g', 'k' );
	SetContained( 'g', 's' );
	SetContained( 'g', 'd' );
	SetContained( 'g', 'v' );
	SetContained( 'g', 'b' );
	SetContained( 'g', 'x' );
	SetContained( 'g', 'n' );
	SetContained( 'G', 'R' );
	SetContained( 'G', 'K' );
	SetContained( 'G', 'S' );
	SetContained( 'G', 'D' );
	SetContained( 'G', 'V' );
	SetContained( 'G', 'B' );
	SetContained( 'G', 'X' );
	SetContained( 'G', 'N' );
	SetContained( 'a', 'r' );
	SetContained( 'a', 'w' );
	SetContained( 'a', 'm' );
	SetContained( 'a', 'd' );
	SetContained( 'a', 'v' );
	SetContained( 'a', 'h' );
	SetContained( 'a', 'x' );
	SetContained( 'a', 'n' );
	SetContained( 'A', 'R' );
	SetContained( 'A', 'W' );
	SetContained( 'A', 'M' );
	SetContained( 'A', 'D' );
	SetContained( 'A', 'V' );
	SetContained( 'A', 'H' );
	SetContained( 'A', 'B' );
	SetContained( 'A', 'X' );
	SetContained( 'A', 'N' );
	SetContained( 'c', 's' );
	SetContained( 'c', 'm' );
	SetContained( 'c', 'y' );
	SetContained( 'c', 'v' );
	SetContained( 'c', 'b' );
	SetContained( 'c', 'h' );
	SetContained( 'c', 'x' );
	SetContained( 'c', 'n' );
	SetContained( 'C', 'S' );
	SetContained( 'C', 'M' );
	SetContained( 'C', 'Y' );
	SetContained( 'C', 'V' );
	SetContained( 'C', 'B' );
	SetContained( 'C', 'H' );
	SetContained( 'C', 'X' );
	SetContained( 'C', 'N' );
	SetContained( 'u', 'k' );
	SetContained( 'u', 'w' );
	SetContained( 'u', 'y' );
	SetContained( 'u', 'd' );
	SetContained( 'u', 'b' );
	SetContained( 'u', 'h' );
	SetContained( 'u', 'x' );
	SetContained( 'u', 'n' );
	SetContained( 'U', 'K' );
	SetContained( 'U', 'W' );
	SetContained( 'U', 'Y' );
	SetContained( 'U', 'D' );
	SetContained( 'U', 'B' );
	SetContained( 'U', 'H' );
	SetContained( 'U', 'X' );
	SetContained( 'U', 'N' );

	SetContained( 'r', 'x' );
	SetContained( 'y', 'x' );
	SetContained( 'w', 'x' );
	SetContained( 's', 'x' );
	SetContained( 'k', 'x' );
	SetContained( 'm', 'x' );
	SetContained( 'b', 'x' );
	SetContained( 'd', 'x' );
	SetContained( 'h', 'x' );
	SetContained( 'v', 'x' );
	SetContained( 'n', 'x' );

	SetContained( 'R', 'X' );
	SetContained( 'Y', 'X' );
	SetContained( 'W', 'X' );
	SetContained( 'S', 'X' );
	SetContained( 'K', 'X' );
	SetContained( 'M', 'X' );
	SetContained( 'B', 'X' );
	SetContained( 'D', 'X' );
	SetContained( 'H', 'X' );
	SetContained( 'V', 'X' );
	SetContained( 'N', 'X' );

	SetContained( 'r', 'n' );
	SetContained( 'y', 'n' );
	SetContained( 'w', 'n' );
	SetContained( 's', 'n' );
	SetContained( 'k', 'n' );
	SetContained( 'm', 'n' );
	SetContained( 'b', 'n' );
	SetContained( 'd', 'n' );
	SetContained( 'h', 'n' );
	SetContained( 'v', 'n' );
	SetContained( 'x', 'n' );

	SetContained( 'R', 'N' );
	SetContained( 'Y', 'N' );
	SetContained( 'W', 'N' );
	SetContained( 'S', 'N' );
	SetContained( 'K', 'N' );
	SetContained( 'M', 'N' );
	SetContained( 'B', 'N' );
	SetContained( 'D', 'N' );
	SetContained( 'H', 'N' );
	SetContained( 'V', 'N' );
	SetContained( 'X', 'N' );

	SetContained( 'k', 'd' );
	SetContained( 'k', 'b' );
	SetContained( 'r', 'd' );
	SetContained( 'r', 'v' );
	SetContained( 's', 'v' );
	SetContained( 's', 'b' );
	SetContained( 'm', 'v' );
	SetContained( 'm', 'h' );
	SetContained( 'w', 'd' );
	SetContained( 'w', 'h' );
	SetContained( 'y', 'b' );
	SetContained( 'y', 'h' );

	SetContained( 'K', 'D' );
	SetContained( 'K', 'B' );
	SetContained( 'R', 'D' );
	SetContained( 'R', 'V' );
	SetContained( 'S', 'V' );
	SetContained( 'S', 'B' );
	SetContained( 'M', 'V' );
	SetContained( 'M', 'H' );
	SetContained( 'W', 'D' );
	SetContained( 'W', 'H' );
	SetContained( 'Y', 'B' );
	SetContained( 'Y', 'H' );
}


}	// end namespace genome

