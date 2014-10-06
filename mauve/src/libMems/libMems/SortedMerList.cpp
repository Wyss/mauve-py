/*******************************************************************************
 * $Id: SortedMerList.cpp,v 1.23 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing, copying, and modification
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/SortedMerList.h"

using namespace std;
using namespace genome;
namespace mems {

const uint8* SortedMerList::BasicDNATable(){
	static const uint8* const bdt = SortedMerList::CreateBasicDNATable();
	return bdt;
}

const uint8* SortedMerList::ProteinTable(){
	static const uint8* const bdt = SortedMerList::CreateProteinTable();
	return bdt;
}

const uint8* SortedMerList::CreateBasicDNATable(){
	uint8* bdt = new uint8[UINT8_MAX];
	memset(bdt, 0, UINT8_MAX);
	bdt['c'] = 1;
	bdt['C'] = 1;
	bdt['b'] = 1;
	bdt['B'] = 1;
	bdt['y'] = 1;
	bdt['Y'] = 1;
	bdt['g'] = 2;
	bdt['G'] = 2;
	bdt['s'] = 2;
	bdt['S'] = 2;
	bdt['k'] = 2;
	bdt['K'] = 2;
	bdt['t'] = 3;
	bdt['T'] = 3;
	return bdt;
}

const uint8* SortedMerList::CreateProteinTable(){
	uint8* pt = new uint8[UINT8_MAX];
	memset(pt, 0, UINT8_MAX);
	pt['A'] = 0;
	pt['R'] = 1;
	pt['N'] = 2;
	pt['D'] = 3;
	pt['C'] = 4;
	pt['Q'] = 5;
	pt['E'] = 6;
	pt['G'] = 7;
	pt['H'] = 8;
	pt['I'] = 9;
	pt['L'] = 10;
	pt['K'] = 11;
	pt['M'] = 12;
	pt['F'] = 13;
	pt['P'] = 14;
	pt['S'] = 15;
	pt['T'] = 16;
	pt['W'] = 17;
	pt['Y'] = 18;
	pt['V'] = 19;
	
	pt['a'] = 0;
	pt['r'] = 1;
	pt['n'] = 2;
	pt['d'] = 3;
	pt['c'] = 4;
	pt['q'] = 5;
	pt['e'] = 6;
	pt['g'] = 7;
	pt['h'] = 8;
	pt['i'] = 9;
	pt['l'] = 10;
	pt['k'] = 11;
	pt['m'] = 12;
	pt['f'] = 13;
	pt['p'] = 14;
	pt['s'] = 15;
	pt['t'] = 16;
	pt['w'] = 17;
	pt['y'] = 18;
	pt['v'] = 19;
	return pt;
}

SortedMerList::SortedMerList(){
	//default to BasicDNA settings
	header.length = 0;
	header.alphabet_bits = 2;
	header.unique_mers = NO_UNIQUE_COUNT;
	memcpy(header.translation_table, BasicDNATable(), UINT8_MAX);
	header.description[0] = 0;
	header.seed_length = DNA_MER_SIZE;
	header.id = 0;
	header.circular = false;
	mask_size = DNA_MER_SIZE;
	mer_mask = 0;
	seed_mask = 0;
	// init sequence data to null
	sequence = NULL;
	binary_seq_len = 0;
}

SortedMerList::SortedMerList( const SortedMerList& sa ){
	sequence = NULL;
	*this = sa;
}

SortedMerList& SortedMerList::operator=(const SortedMerList& sa)
{
	header = sa.header;
	mer_mask = sa.mer_mask;
	seed_mask = sa.seed_mask;
	mask_size = sa.mask_size;
	binary_seq_len = sa.binary_seq_len;

	// copy binary sequence data
	if( sa.sequence != NULL ){
		if( sequence != NULL )
			delete[] sequence;
		sequence = new uint32[binary_seq_len];
		memcpy(sequence, sa.sequence, sizeof(uint32) * binary_seq_len);
	}else
		sequence = NULL;

	return *this;
}

SortedMerList::~SortedMerList(){
	if( sequence != NULL )
		delete[] sequence;
}

void SortedMerList::Clear(){
	//default to BasicDNA settings
	header.length = 0;
	header.alphabet_bits = 2;
	header.unique_mers = NO_UNIQUE_COUNT;
	memcpy(header.translation_table, BasicDNATable(), UINT8_MAX);
	header.description[0] = 0;
	header.seed_length = DNA_MER_SIZE;
	header.id = 0;
	header.circular = false;
	mask_size = DNA_MER_SIZE;
	mer_mask = 0;
	seed_mask = 0;
	// delete sequence data
	if( sequence != NULL ){
		delete[] sequence;
		sequence = NULL;
	}
	binary_seq_len = 0;
}

uint32 SortedMerList::CalculateMaxMerSize() const{
	bmer tmp;
	return (sizeof(tmp.mer) * 8) / header.alphabet_bits;
}

boolean SortedMerList::FindMer(const uint64 query_mer, gnSeqI& result){
	bmer merle;
	merle.mer = query_mer;
	gnSeqI last_pos = Length();
	if( last_pos == 0 || (last_pos < header.seed_length && !header.circular) )
		return false;
	last_pos -= header.circular ? 1 : header.seed_length;
	result = bsearch(merle, 0, last_pos );
	return ((*this)[result].mer == merle.mer);
}

boolean SortedMerList::Find(const string& query_seq, gnSeqI& result) {
	struct bmer merle;
	merle.mer = 0;

	//check the length to make sure it is small enough
	gnSeqI len = query_seq.length() * header.alphabet_bits < 64 ? 
		query_seq.length() : 64 / header.alphabet_bits;
		
	translate((uint8*)&merle.mer, query_seq.c_str(), len);
	return FindMer( merle.mer, result );
}

void SortedMerList::FindAll(const string& query_seq, vector<gnSeqI> result) {
	struct bmer merle;
	merle.mer = 0;

	//check the length to make sure it is small enough
	gnSeqI len = query_seq.length() * header.alphabet_bits < 64 ? 
		query_seq.length() : 64 / header.alphabet_bits;
		
	translate((uint8*)&merle.mer, query_seq.c_str(), len);
	
	//find the first match then start filling forward.
	gnSeqI matchI = 0;
	gnSeqI last_pos = Length();
	last_pos -= header.circular ? 1 : header.seed_length;
	bmer matchmer;
	matchI = bsearch(merle, 0, last_pos);

	//first seek backwards
	int64 cur_matchI = matchI;
	matchmer = (*this)[matchI];
	while(cur_matchI >= 0 && matchmer.mer == merle.mer){
		cur_matchI--;
		matchmer = (*this)[cur_matchI];
	}
	int64 first_matchI = cur_matchI+1;

	//now seek forwards
	cur_matchI = matchI+1;
	matchmer = (*this)[cur_matchI];
	while(cur_matchI < GNSEQI_END && matchmer.mer == merle.mer){
		cur_matchI++;
		matchmer = (*this)[cur_matchI];
	}
	//fill the result array
	for(matchI = first_matchI; matchI < cur_matchI; matchI++)
		result.push_back(matchI);
}

string SortedMerList::Description() const{
	return header.description;
}

void SortedMerList::SetDescription(const string& d){
	strncpy(header.description, d.c_str(), DESCRIPTION_SIZE-1);
}

uint SortedMerList::SeedLength() const{
	return header.seed_length;
}
/**
 * Returns the weight of the seed that this SML was sorted on.
 */
uint SortedMerList::SeedWeight() const{
	return header.seed_weight;
}
/**
 * Returns the seed pattern that this SML was sorted on.
 */
uint64 SortedMerList::Seed() const{
	return header.seed;
}

boolean SortedMerList::IsCircular() const{
	return header.circular;
}

uint64 SortedMerList::GetMerMask() const{
	return mer_mask;
}

uint64 SortedMerList::GetSeedMask() const{
	return seed_mask;
}

uint32 SortedMerList::GetMerMaskSize() const{
	return mask_size;
}

void SortedMerList::SetMerMaskSize(uint32 mer_size){
	if(mer_size > header.seed_length)
		mask_size = header.seed_length;
	else
		mask_size = mer_size;

	// calculate the mer mask
	mer_mask = UINT32_MAX;
	mer_mask <<= 32;
	mer_mask |= UINT32_MAX;
	mer_mask <<= (64 - header.alphabet_bits * mer_size);
}

gnSeqI SortedMerList::Length() const{
	return header.length;
}

gnSeqI SortedMerList::SMLLength() const{
	// make sure there was at least one seed
	if( header.length < header.seed_length )
		return 0;
	if( !header.circular )
		return header.length - header.seed_length + 1;
	return header.length;
}

sarID_t SortedMerList::GetID() const{
	return header.id;
}
void SortedMerList::SetID(const sarID_t d){
	header.id = d;
}

#define OPT_HEADER_ALPHABET_BITS DNA_ALPHA_BITS

void SortedMerList::SetSequence(gnSeqC* seq_buf, gnSeqI seq_len){
	binary_seq_len = (seq_len * header.alphabet_bits) / 32;
	if((seq_len * header.alphabet_bits) % 32 != 0)
		binary_seq_len++;

	binary_seq_len+=2;	// zero-pad the end for extra working room

	if( sequence != NULL )
		delete[] sequence;
	sequence = new uint32[binary_seq_len];
	translate32(sequence, seq_buf, seq_len);
}

// this should return a mer containing all characters covered by the
// spaced seed
uint64 SortedMerList::GetMer(gnSeqI position) const
{
	//check this for access violations.
	uint64 mer_a;
	gnSeqI mer_word, mer_bit;
	uint32 merle;
	//get mer_a
	mer_a = 0;
	mer_word = (position * (gnSeqI)OPT_HEADER_ALPHABET_BITS) / (gnSeqI)32;
	mer_bit = (position * (gnSeqI)OPT_HEADER_ALPHABET_BITS) % (gnSeqI)32;
	mer_a |= sequence[mer_word++];
	mer_a <<= 32;
	mer_a |= sequence[mer_word++];
	if(mer_bit > 0){
		merle = sequence[mer_word];
		merle >>= 32 - mer_bit;
		mer_a <<= mer_bit;
		mer_a |= merle;
	}
	mer_a &= mer_mask;
	return mer_a;
}

//potential buffer overflows here.  make dest extra big.
void SortedMerList::GetBSequence(uint32* dest, const gnSeqI len, const gnSeqI offset){
	//first determine the byte offset of the sequence within the file.
	if(offset >= header.length){
		Throw_gnEx( IndexOutOfBounds() );
	}
	uint64 startpos = (offset * OPT_HEADER_ALPHABET_BITS) / 32;
	int begin_remainder = (offset * OPT_HEADER_ALPHABET_BITS) % 32;
	uint64 readlen = offset + len < header.length ? len : header.length - offset;

	gnSeqI word_read_len = (readlen * OPT_HEADER_ALPHABET_BITS) / 32;
	int end_remainder = (readlen * OPT_HEADER_ALPHABET_BITS) % 32;
	if(begin_remainder + (readlen * OPT_HEADER_ALPHABET_BITS) > 32
	   && end_remainder > 0)
		word_read_len++;
	if(begin_remainder > 0)
		word_read_len++;
	
	//now do the actual read
	memcpy((char*)dest, (char*)sequence + (startpos * 4), word_read_len * 4);
	
	//now shift if needed
	ShiftWords(dest, word_read_len, -begin_remainder);
	
	//now mask if needed
	if(end_remainder > begin_remainder){
		uint32 mask = 0xFFFFFFFF;
		mask <<= 32 - (end_remainder - begin_remainder);
		dest[word_read_len-1] &= mask;
	}else if(end_remainder < begin_remainder){
		uint32 mask = 0xFFFFFFFF;
		mask <<= (begin_remainder - end_remainder);
		dest[word_read_len-2] &= mask;
	}
}

gnSeqI SortedMerList::bsearch(const struct bmer& query_mer, const gnSeqI start, const gnSeqI end) {

	gnSeqI middle = (start + end) / 2;
	struct bmer midmer = (*this)[middle];
	if(midmer.mer == query_mer.mer)
		return middle;
	else if((midmer.mer < query_mer.mer) && (middle < end))
		return bsearch(query_mer, middle + 1, end);
	else if((midmer.mer > query_mer.mer) && (start < middle))
		return bsearch(query_mer, start, middle - 1);
	
	//if we get here then the mer was not found.
	//return where it would be if it existed.
	return middle;
}

//translate the character sequence to binary form based on the
//translation table.
void SortedMerList::translate(uint8* dest, const gnSeqC* src, const gnSeqI len) const{
	uint8 start_bit = 0;
	gnSeqI cur_byte = 0;
	const uint32 alpha_bits = OPT_HEADER_ALPHABET_BITS;
	dest[cur_byte] = 0;
	for(uint32 i=0; i < len; i++){
		uint8 tmp = header.translation_table[src[i]];
		if(start_bit + alpha_bits <= 8){
			tmp <<= 8 - start_bit - alpha_bits;
			dest[cur_byte] |= tmp;
		}else{
			uint8 over_bits = (start_bit + alpha_bits) % 8;
			uint8 tmp2 = tmp;
			tmp2 <<= 8 - over_bits;
			tmp >>= over_bits;
			dest[cur_byte] |= tmp;
			dest[cur_byte+1] |= tmp2;
		}
		start_bit += alpha_bits;
		if(start_bit >= 8){
			start_bit %= 8;
			cur_byte++;
			dest[cur_byte] = 0;
		}
	}
}

void SortedMerList::translate32(uint32* dest, const gnSeqC* src, const gnSeqI len) const{
	if( len == 0 )
		return;
	uint8 start_bit = 0;
	gnSeqI cur_word = 0;
	const uint32 alpha_bits = OPT_HEADER_ALPHABET_BITS;
	dest[cur_word] = 0;
	for(uint32 i=0; i < len; i++){
		if(src[i]=='-'){
			cerr << "ERROR! gap character encountered at genome sequence position " << i << std::endl;
			cerr << "Input sequences must be unaligned and ungapped!\n";
			throw "Gap in genome sequence\n";
		}
		uint32 tmp = header.translation_table[src[i]];
		if(start_bit + alpha_bits <= 32){
			tmp <<= 32 - start_bit - alpha_bits;
			dest[cur_word] |= tmp;
			start_bit += alpha_bits;
			if(start_bit >= 32 && i < len - 1){
				start_bit %= 32;
				cur_word++;
				dest[cur_word] = 0;
			}
		}else{
			uint8 over_bits = (start_bit + alpha_bits) % 32;
			uint32 tmp2 = tmp;
			tmp2 <<= 32 - over_bits;
			tmp >>= over_bits;
			dest[cur_word] |= tmp;
			cur_word++;
			dest[cur_word] = 0;
			dest[cur_word] |= tmp2;
			start_bit = over_bits;
		}
	}
}
SMLHeader SortedMerList::GetHeader() const{
	return header;
}

gnSeqI SortedMerList::UniqueMerCount(){
	if(header.unique_mers != NO_UNIQUE_COUNT)
		return header.unique_mers;

	uint32 MER_BUFFER_SIZE = 16384;  //not quite arbitrary (2^14)
	gnSeqI cur_pos = 0;
	vector<bmer> mer_vector;
	bmer prev_mer;
	gnSeqI m_unique = 0;
	gnSeqI report_interval = MER_BUFFER_SIZE * 212;
	while(cur_pos < header.length){
		if(!Read(mer_vector, MER_BUFFER_SIZE, cur_pos)){
			break;
//			DebugMsg("SortedMerList::UniqueMerCount: Error reading bmer vector.");
//			return NO_UNIQUE_COUNT;
		}
		uint32 mer_count = mer_vector.size();
		if(mer_count == 0)
			break;
		if(cur_pos > 0 && prev_mer.mer != mer_vector[0].mer)
			m_unique++;
		
		//count them up.
		uint32 i = 0;
		for(uint32 j = 1; j < mer_count; j++){
			if((mer_vector[i].mer & mer_mask) != (mer_vector[j].mer & mer_mask) )
				m_unique++;
			i++;
		}
		prev_mer = mer_vector[i];
		cur_pos += mer_count;
		if( cur_pos % report_interval == 0 ){
//			cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
			cout << m_unique << "/" << cur_pos << endl;
		}
	}
	cout << endl;
	m_unique++;
	header.unique_mers = m_unique;
	return header.unique_mers;
}

//will not handle more than 8GB sequence on 32-bit systems
void SortedMerList::ShiftWords(unsigned int* data, uint32 length, int32 bits)
{
	int32 word_bits = 8 * sizeof(unsigned int);
	if(bits > 0 && bits < word_bits){
		//shift everything right starting at the end
		data[length - 1] >>= bits;
		for(int i=length-2; i >= 0; i--){
			uint32 tmp = data[i];
			tmp <<= word_bits - bits;
			data[i+1] |= tmp;
			data[i] >>= bits;
		}
	}else if(bits < 0 && bits > (-1)*word_bits){
		bits *= -1;
		//shift everything left
		data[0] <<= bits;
		for(uint32 i=0; i < length; i++){
			uint32 tmp = data[i+1];
			tmp >>= word_bits - bits;
			data[i] |= tmp;
			data[i+1] <<= bits;
		}
	}
}

void SortedMerList::FillSML(gnSeqC* seq_buf, gnSeqI seq_len, boolean circular, vector<bmer>& sml_array){
	const uint32 alpha_bits = OPT_HEADER_ALPHABET_BITS;
	const uint32 mer_size = header.seed_length;
	gnSeqI sar_len = seq_len;
	if(!circular)
		sar_len -= header.seed_length - 1;
	sml_array.reserve(sar_len);

	bmer cur_suffix;
	cur_suffix.mer = 0;
	cur_suffix.position = 0;

	/* now fill in the suffix array with the forward sequence*/
	for(gnSeqI i=0; i < mer_size; i++){
		cur_suffix.mer <<= alpha_bits;
		cur_suffix.mer |= header.translation_table[seq_buf[i]];
	}
	uint8 dead_bits = 64 - (mer_size * alpha_bits);
	cur_suffix.mer <<= dead_bits;

	sml_array.push_back(cur_suffix);

	//fill sml_array with mers
	for(gnSeqI seqI = 1; seqI < sar_len; seqI++){//already added the
													//first one
		cur_suffix.position++;
		cur_suffix.mer <<= alpha_bits;
		uint64 new_mer = header.translation_table[seq_buf[seqI+(mer_size-1)]];
		new_mer <<= dead_bits;
		cur_suffix.mer |= new_mer;
		sml_array.push_back(cur_suffix);
	}
}

void SortedMerList::FillSML(const gnSequence& seq, vector<bmer>& sml_array){
	gnSeqI seq_len = seq.length();
	Array<gnSeqC> seq_buf( seq_len );
	seq.ToArray(seq_buf.data, seq_len);
	FillSML(seq_buf.data, seq_len, seq.isCircular(), sml_array);
}

void SortedMerList::FillSML(gnSeqI seq_len, vector<gnSeqI>& pos_array){
	pos_array.clear();
	pos_array.reserve( seq_len );
	for(gnSeqI seqI = 0; seqI < seq_len; seqI++ )
		pos_array.push_back(seqI);
}

uint64 SortedMerList::GetDnaMer(gnSeqI offset) const
{
	// get the forward orientation mer
	uint64 mer_a = SortedMerList::GetMer( offset );
	//find the reverse complement of mer_a and return it if it's
	//smaller
	uint64 mer_c = RevCompMer( mer_a, header.seed_length );	//mer_c will be the reverse complement
	
	// for debugging
//	if( mer_c < mer_a )
//		return mer_c;
	return mer_a < mer_c ? mer_a : mer_c;
}

#define OPT_ALPHA_MASQ 0x00000003

uint64 SortedMerList::RevCompMer( uint64 mer_a, int mer_length ) const
{
	//find the reverse complement of mer_a and return it if it's
	//smaller
	uint64 mer_b, mer_c = 0;	//mer_c will be the reverse complement
	mer_b = ~mer_a;
//	uint32 masq = 0xffffffff;
//	masq >>= 32 - header.alphabet_bits;
	for(uint32 i = 0; i < 64; i += OPT_HEADER_ALPHABET_BITS){
		mer_c |= mer_b & OPT_ALPHA_MASQ;
//		mer_c |= mer_b & masq;
		mer_b >>= OPT_HEADER_ALPHABET_BITS;
		mer_c <<= OPT_HEADER_ALPHABET_BITS;
	}
	mer_c <<= 64 - (OPT_HEADER_ALPHABET_BITS * (mer_length+1));
	mer_c |= 1;
	return mer_c;
}


void SortedMerList::FillDnaSML(const gnSequence& seq, vector<bmer>& sml_array){
	/* now fill in the suffix array with the forward sequence*/
	const uint32 alpha_bits = OPT_HEADER_ALPHABET_BITS;
	const uint32 mer_size = header.seed_length;
	gnSeqI sar_len = seq.length();
	if( sar_len < header.seed_length )
		return;	// can't have an sml if there ain't enough sequence
	if( !seq.isCircular() )
		sar_len -= ( header.seed_length - 1);
	sml_array.reserve(sar_len);

	uint32 dead_bits = 64 - (mer_size * alpha_bits);
	uint64 create_mask = UINT32_MAX;
	create_mask <<= 32;
	create_mask |= UINT32_MAX;
	create_mask <<= dead_bits;

	bmer cur_suffix, rcur_suffix;
	cur_suffix.mer = sequence[0];
	cur_suffix.mer <<= 32;
	cur_suffix.mer |= sequence[1];
	cur_suffix.mer &= create_mask;
	cur_suffix.position = 0;
	rcur_suffix.mer = 0;
	rcur_suffix.position = 0;
	
	//find the reverse complement of cur_suffix.mer and return it if it's
	//smaller
	uint64 mer_b = 0;
	mer_b = ~cur_suffix.mer;
//	uint32 masq = 0xffffffff;
//	masq >>= 32 - alpha_bits;
	for(uint32 i = 0; i < 64; i += alpha_bits){
//		rcur_suffix.mer |= mer_b & masq;
		rcur_suffix.mer |= mer_b & OPT_ALPHA_MASQ;
		mer_b >>= alpha_bits;
		rcur_suffix.mer <<= alpha_bits;
	}
	rcur_suffix.mer <<= dead_bits - alpha_bits;
	rcur_suffix.mer |= 1;

	//add the first mer
	if(cur_suffix.mer < rcur_suffix.mer)
		sml_array.push_back(cur_suffix);
	else
		sml_array.push_back(rcur_suffix);

	//fill sml_array with mers
	gnSeqI 	endI = sar_len + mer_size;
	if(seq.isCircular())
		endI += mer_size;

	uint32 rdead_bits = 64 - alpha_bits - dead_bits;
	uint64 tmp_rseq = 0;
	uint32 seqI = (mer_size * alpha_bits) / 32;
	int32 cur_bit = 32 - alpha_bits - ((mer_size * alpha_bits) % 32);
	uint32 cur_seq = sequence[seqI];
	uint64 tmp_seq;
//	uint32 alpha_mask = 0xFFFFFFFF;
//	alpha_mask >>= 32 - alpha_bits;
	uint64 revalpha_mask = OPT_ALPHA_MASQ;
	revalpha_mask <<= dead_bits;

	//which is slower? a memory operation or a conditional?
	//probably a memory operation.
	for(gnSeqI cur_pos = mer_size + 1; cur_pos < endI; cur_pos++){//already added the
													//first one
		//increment positions
		cur_suffix.position++;
		rcur_suffix.position++;
		//extract the next character
		tmp_seq = cur_seq;
		tmp_seq >>= cur_bit;
		tmp_seq &= OPT_ALPHA_MASQ;
		tmp_seq <<= dead_bits;
		
		//add it to the forward mer
		cur_suffix.mer <<= alpha_bits;
		cur_suffix.mer |= tmp_seq;

		//do the reverse complement mer
		tmp_seq = ~tmp_seq;
		tmp_seq &= revalpha_mask;
		tmp_rseq = tmp_seq;
		tmp_rseq <<= rdead_bits;
		rcur_suffix.mer >>= alpha_bits;
		rcur_suffix.mer |= tmp_rseq;
		rcur_suffix.mer &= create_mask;
		rcur_suffix.mer |= 1;
		if(cur_suffix.mer < rcur_suffix.mer)
			sml_array.push_back(cur_suffix);
		else
			sml_array.push_back(rcur_suffix);

		cur_bit -= alpha_bits;
		if(cur_bit < 0){
			cur_bit += alpha_bits;
			cur_seq <<= 16;		//trade bitwise ops for conditional
			cur_seq <<= 16 - (cur_bit);
			seqI++;
			tmp_seq = sequence[seqI];
			tmp_seq >>= cur_bit;
			cur_seq |= tmp_seq;
			cur_bit += 32 - alpha_bits;
		}
	}
}


uint64 SortedMerList::GetSeedMer( gnSeqI offset ) const
{
	//check this for access violations.
	uint64 mer_a = SortedMerList::GetMer( offset );
	uint64 mer_b = SortedMerList::GetMer( offset + 1 );
	uint64 seed_mer = 0;
	uint64 alpha_mask = 1;
	alpha_mask <<= OPT_HEADER_ALPHABET_BITS;
	alpha_mask--;
	alpha_mask <<= 62;
	uint64 cur_alpha_mask = alpha_mask;
	uint64 char_mask = 1;
	char_mask <<= header.seed_length - 1;
	uint64 cur_mer = mer_a;
	const int mer_transition = 64 / OPT_HEADER_ALPHABET_BITS;
	int patternI = 0;
	int rshift_amt = 64 - OPT_HEADER_ALPHABET_BITS;
	for( ; patternI < header.seed_length; patternI++ ){
		if( patternI == mer_transition ){
			cur_mer = mer_b;
			cur_alpha_mask = alpha_mask;
			rshift_amt = 64 - OPT_HEADER_ALPHABET_BITS;
		}
		if( (header.seed & char_mask) != 0 ){
			uint64 char_tmp = cur_mer & cur_alpha_mask;
			char_tmp >>= rshift_amt;
			seed_mer <<= OPT_HEADER_ALPHABET_BITS;
			seed_mer |= char_tmp;
		}
		cur_alpha_mask >>= OPT_HEADER_ALPHABET_BITS;
		char_mask >>= 1;
		rshift_amt -= OPT_HEADER_ALPHABET_BITS;
	}

	seed_mer <<= 64 - (OPT_HEADER_ALPHABET_BITS * header.seed_weight);
	return seed_mer;
}

uint64 SortedMerList::GetDnaSeedMer( gnSeqI offset ) const
{
	uint64 seed_mer = SortedMerList::GetSeedMer( offset );
	uint64 rev_mer = RevCompMer( seed_mer, header.seed_weight );
	return seed_mer < rev_mer ? seed_mer : rev_mer;
}

void SortedMerList::FillDnaSeedSML(const gnSequence& seq, vector<bmer>& sml_array){
	// first get the length of the sequence
	gnSeqI sar_len = SMLLength();
	if( sar_len == 0 )
		return;	// can't have an sml if there ain't enough sequence
	sml_array.resize(sar_len);
	
	/* now fill in the sml_array with the forward sequence */
	for( gnSeqI seedI = 0; seedI < sar_len; seedI++ ){
		sml_array[seedI].mer = GetDnaSeedMer( seedI );
		sml_array[seedI].position = seedI;
	}
}


void SortedMerList::Create(const gnSequence& seq, const uint64 seed){
	
	if(CalculateMaxMerSize() == 0)
		Throw_gnExMsg( SMLCreateError(), "Alphabet size is too large" );

	int seed_length = getSeedLength( seed );
	int seed_weight = getSeedWeight( seed );
	
	if(seed_length > CalculateMaxMerSize())
		Throw_gnExMsg( SMLCreateError(), "Mer size is too large" );

	if(seed_length == 0)
		Throw_gnExMsg( SMLCreateError(), "Can't have 0 seed length" );

	//determine sequence and sar length and read in sequence
	gnSeqI seq_len = seq.length();
	if(!seq.isCircular()){
		header.circular = false;
	}else
		header.circular = true;
	// use the nifty Array class as a wrapper for the buffer to ensure correct deallocation
	gnSeqI buf_len = seq.isCircular() ? seq_len + seed_length : seq_len;
	Array<gnSeqC> seq_buf( buf_len );
	seq.ToArray(seq_buf.data, seq_len);
	if( seq.isCircular() )
		seq.ToArray(seq_buf.data + seq_len, seed_length-1);

	// set header information
	header.length = seq_len;
	header.seed_length = seed_length;
	header.seed_weight = seed_weight;
	header.seed = seed;

	SetMerMaskSize( seed_weight );
	seed_mask = mer_mask;
	SetMerMaskSize( seed_length );

	SetSequence( seq_buf.data, buf_len );
}

} // namespace mems
