/*******************************************************************************
 * $Id: DNAFileSML.cpp,v 1.4 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing, copying, and modification
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnSequence.h"
#include "libGenome/gnFilter.h"
#include "libMems/DNAFileSML.h"

using namespace std;
using namespace genome;
namespace mems {

DNAFileSML::DNAFileSML() : FileSML(){
	FileSML::header.version = FormatVersion();
}

DNAFileSML::DNAFileSML(const string& fname, const uint8* table, const uint32 alpha_bits){
	header.alphabet_bits = alpha_bits;
	memcpy(header.translation_table, table, UINT8_MAX);
	filename = fname;
	header.version = FormatVersion();
}

DNAFileSML& DNAFileSML::operator=(const DNAFileSML& msa ){
	FileSML::operator=(msa);
	return *this;
}

DNAFileSML* DNAFileSML::Clone() const{
	DNAFileSML *bdsa = new DNAFileSML();
	(*bdsa) = *this;
	return bdsa;
}

uint64 DNAFileSML::GetNeededMemory(gnSeqI len){
	uint64 neededmem = (len * FileSML::header.alphabet_bits) / 8;
	//forward and reverse copies of the sequence
	neededmem += len * 2;
	neededmem += sizeof(bmer) * len;
	return neededmem;
}

uint32 DNAFileSML::CalculateMaxMerSize() const{
	return 62 / header.alphabet_bits;
}

uint64 DNAFileSML::GetMer(gnSeqI position) const{
	return GetDnaMer( position );
}

uint64 DNAFileSML::GetSeedMer( gnSeqI offset ) const{
	return GetDnaSeedMer( offset );
}

void DNAFileSML::FillSML(const gnSequence& seq, vector<bmer>& sml_array)
{
	FillDnaSML(seq, sml_array);
}

} // namespace mems
