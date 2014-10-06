/*******************************************************************************
 * $Id: DNAMemorySML.cpp,v 1.3 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing, copying, and modification
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnSequence.h"
#include "libMems/DNAMemorySML.h"

using namespace std;
using namespace genome;
namespace mems {

DNAMemorySML::DNAMemorySML(const uint8* table, const uint32 alpha_bits) : 
MemorySML( table, alpha_bits )
{}

DNAMemorySML& DNAMemorySML::operator=(const DNAMemorySML& msa ){
	MemorySML::operator=(msa);
	return *this;
}

DNAMemorySML* DNAMemorySML::Clone() const{
	DNAMemorySML *bdsa = new DNAMemorySML();
	(*bdsa) = *this;
	return bdsa;
}

uint64 DNAMemorySML::GetMer(gnSeqI position) const{
	return GetDnaMer( position );
}

uint64 DNAMemorySML::GetSeedMer( gnSeqI offset ) const{
	return GetDnaSeedMer( offset );
}

void DNAMemorySML::FillSML(const gnSequence& seq, vector<bmer>& sml_array)
{
	FillDnaSML(seq, sml_array);
}

} // namespace mems
