/*******************************************************************************
 * $Id: MemorySML.cpp,v 1.8 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing, copying, and modification
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/MemorySML.h"
#include <algorithm>

using namespace std;
using namespace genome;
namespace mems {

MemorySML::MemorySML(const uint8* table, const uint32 alpha_bits ){
	header.alphabet_bits = alpha_bits;
	memcpy(header.translation_table, table, UINT8_MAX);
	header.version = 0;
}

MemorySML::MemorySML(const MemorySML& msa){
	*this = msa;
}

MemorySML& MemorySML::operator=(const MemorySML& msa ) {
 	SortedMerList::operator=( msa );
 	positions = msa.positions;
 	return *this;
}

MemorySML* MemorySML::Clone() const{
	return new MemorySML(*this);
}

void MemorySML::Clear(){
	SortedMerList::Clear();
	positions.clear();
}

void MemorySML::Create(const gnSequence& seq, const uint64 seed ){
	SortedMerList::Create(seq, seed);

	vector<bmer> sml_array;
	boolean is_spaced_seed = header.seed_length != header.seed_weight;
	if( is_spaced_seed )
		FillDnaSeedSML( seq, sml_array );
	else
		FillSML( seq, sml_array );
	sort( sml_array.begin(), sml_array.end(), &bmer_lessthan );
	positions.reserve( sml_array.size() );
	for(gnSeqI merI = 0; merI < sml_array.size(); merI++ ){
		positions.push_back( sml_array[merI].position );
	}

}

boolean MemorySML::Read(vector<bmer>& readVector, gnSeqI size, gnSeqI offset )
{
	readVector.clear();
	if( offset > positions.size() )
		return false;

	gnSeqI last_mer = offset + size;
	boolean success = true;
	if( last_mer > positions.size() ){
		last_mer = positions.size();
		success = false;
	}

	bmer cur_mer;
	for(gnSeqI merI = offset; merI < last_mer; merI++){
		cur_mer.position = positions[merI];
		cur_mer.mer = GetSeedMer( cur_mer.position );
		readVector.push_back( cur_mer );
	}
	return success;
}

void MemorySML::Merge(SortedMerList& sa, SortedMerList& sa2){

}

bmer MemorySML::operator[](gnSeqI index)
{
	bmer cur_mer;
	cur_mer.position = positions[index];
	cur_mer.mer = GetSeedMer( cur_mer.position );
	return cur_mer;
}

} // namespace mems
