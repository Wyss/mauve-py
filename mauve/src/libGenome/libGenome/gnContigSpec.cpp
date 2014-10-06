/////////////////////////////////////////////////////////////////////////////
// File:            gnContigSpec.cpp
// Purpose:         Abstract Contig Spec class
// Description:     Defines an interface for contig specs
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

#include "libGenome/gnContigSpec.h"
#include "libGenome/gnFilter.h"


using namespace std;
namespace genome {


void gnContigSpec::CropStart( gnSeqI cropLen ){
	m_start = m_reverseComplement ? (GetSourceLength() + m_start - cropLen) % GetSourceLength() : (m_start + cropLen) % GetSourceLength();
	m_length -= cropLen;
	return;
}

void gnContigSpec::CropEnd( gnSeqI cropLen ){
	m_length -= cropLen;
}

void gnContigSpec::SetReverseComplement( const boolean value )
{
	//translate coordinates if revComp has changed.
	if((m_reverseComplement != value) && (m_length > 0))
		m_start = (m_start + m_length) % GetSourceLength();
	m_reverseComplement = value;
}

boolean gnContigSpec::SeqRead(const gnSeqI start, gnSeqC* buf, gnSeqI& bufLen, const uint32 contigI ) const{
	boolean success;
	bufLen = bufLen < m_length - start? bufLen : m_length - start;	//trim the read down to size.
	gnSeqI readable = bufLen;
	gnSeqI read_start = start;	// coordinate translation

	if(contigI == ALL_CONTIGS)
		read_start = m_reverseComplement ? (m_start - start - readable + GetSourceLength()) % GetSourceLength() : start + m_start;

	success = Read(read_start, buf, readable);

	if(m_circular){
		gnSeqI beginread = bufLen - readable;  //read whatever couldn't be read before
		read_start = m_reverseComplement ? (m_start - readable + GetSourceLength()) % GetSourceLength() : m_start;
		success = Read(read_start , buf + readable, beginread);
		readable += beginread;
	}

	bufLen = readable;

	if(m_reverseComplement)
		gnFilter::DNAComplementFilter()->ReverseFilter(&buf, bufLen);

	return success;
}

void gnContigSpec::Clear(){
	gnBaseSpec::Clear();
	m_start = 0;
	m_length = 0;
	m_SourceContigIndex = ALL_CONTIGS;
}


}	// end namespace genome

