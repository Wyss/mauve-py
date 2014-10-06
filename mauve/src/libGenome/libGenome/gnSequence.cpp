/////////////////////////////////////////////////////////////////////////////
// File:            gnSequence.h
// Purpose:         Source Seq, where sequence is stored in source file
// Description:     implements the baseSeq interface for dna source files.
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

#include "libGenome/gnSequence.h"
#include "libGenome/gnSourceFactory.h"
#include "libGenome/gnBaseSource.h"
#include "libGenome/gnGenomeSpec.h"
#include "libGenome/gnFragmentSpec.h"
#include "libGenome/gnSourceSpec.h"
#include "libGenome/gnStringSpec.h"
#include "libGenome/gnStringHeader.h"
#include <cstring>


using namespace std;
namespace genome {


// Constructor
gnSequence::gnSequence( )
{
	spec = new gnGenomeSpec();
	comparator = gnCompare::DNASeqCompare();
}
gnSequence::gnSequence( const gnSeqC* seq){
	spec = new gnGenomeSpec();
	if (seq[0] != 0){
		gnFragmentSpec *fragSpec = new gnFragmentSpec();
		spec->AddSpec(fragSpec);
		fragSpec->AddSpec(new gnStringSpec(seq));
	}
	comparator = gnCompare::DNASeqCompare();
}
gnSequence::gnSequence( const string& str){
	spec = new gnGenomeSpec();
	if (str.length() != 0){
		gnFragmentSpec *fragSpec = new gnFragmentSpec();
		spec->AddSpec(fragSpec);
		fragSpec->AddSpec(new gnStringSpec(str));
	}
	comparator = gnCompare::DNASeqCompare();
}
gnSequence::gnSequence( const gnGenomeSpec& gngs){
	spec = gngs.Clone();
	comparator = gnCompare::DNASeqCompare();
}
gnSequence::gnSequence( const gnFragmentSpec& gnfs){
	spec = new gnGenomeSpec();
	spec->AddSpec(gnfs.Clone());
	comparator = gnCompare::DNASeqCompare();
}
gnSequence::gnSequence( const gnContigSpec& gcs){
	spec = new gnGenomeSpec();
	gnFragmentSpec *fragSpec = new gnFragmentSpec();
	fragSpec->AddSpec(gcs.Clone());
	comparator = gnCompare::DNASeqCompare();
}
gnSequence::gnSequence( gnSeqC *bases, gnSeqI length){
	spec = new gnGenomeSpec();
	if (length > 0){
		gnFragmentSpec *fragSpec = new gnFragmentSpec();
		spec->AddSpec(fragSpec);
		fragSpec->AddSpec(new gnStringSpec(string(bases, length)));
	}
	comparator = gnCompare::DNASeqCompare();
}	
gnSequence::gnSequence( const gnSequence& seq )
{
	*this = seq;
}

// Destructor
gnSequence::~gnSequence()
{
	if (spec != NULL) 
		delete spec;
}

gnSequence& gnSequence::operator=( const gnSequence& seq){
	spec = seq.spec->Clone();
	filter_list = seq.filter_list;
	comparator = seq.comparator;
	return *this;
}

// Clone
gnSequence* gnSequence::Clone() const
{
	return new gnSequence( *this );
}
//IMPLEMENT ME!! Consider ambiguities!
int gnSequence::compare(const string& str) const{
	STACK_TRACE_START
		gnSeqI len = length();
		gnSeqI seq_len = str.length();
		gnSeqI comparelen = len < seq_len ? len : seq_len;
		gnSeqI compared = 0;
		while(comparelen > 0){
			gnSeqI curlen = comparelen < BUFFER_SIZE ? comparelen : BUFFER_SIZE;
			string bases = ToString(compared, curlen);
			string seq_bases = str.substr(compared, curlen);
			if(comparator->LessThan(bases, seq_bases))
				return -1;
			else if(comparator->LessThan(seq_bases, bases))
				return 1;

			compared += curlen;
			comparelen -= curlen;
		}
		if(len < seq_len)
			return -1;
		else if(len > seq_len)
			return 1;
		return 0;
	STACK_TRACE_END
}

int gnSequence::compare(const gnSequence& seq) const{
	STACK_TRACE_START
		gnSeqI len = length();
		gnSeqI seq_len = seq.length();
		gnSeqI comparelen = len < seq_len ? len : seq_len;
		gnSeqI compared = 0;
		while(comparelen > 0){
			gnSeqI curlen = comparelen < BUFFER_SIZE ? comparelen : BUFFER_SIZE;
			string bases = ToString(curlen, compared+1);
			string seq_bases = seq.ToString(curlen, compared+1);
			if(comparator->LessThan(bases, seq_bases))
				return -1;
			else if(comparator->LessThan(seq_bases, bases))
				return 1;

			compared+= curlen;
			comparelen -= curlen;
		}
		if(len < seq_len)
			return -1;
		else if(len > seq_len)
			return 1;
		return 0;
	STACK_TRACE_END
}

void gnSequence::insert( const gnSeqI offset, const gnSeqC *bases, const gnSeqI len){
	STACK_TRACE_START
		string str(bases, len);
		gnStringSpec gpbs(str);
		insert(offset, gpbs);
	STACK_TRACE_END
}
//Offset is the gene index, not a computer index, to insert before.
//This is the default insert.  It inserts entire contigs
void gnSequence::insert( const gnSeqI offset, const gnGenomeSpec& gpbs){
	STACK_TRACE_START
		if(offset == 0)
			Throw_gnEx(SeqIndexOutOfBounds());
		if(offset == GNSEQI_END || spec->GetLength() < offset){ //simple append.
			for(uint32 i=0; i < gpbs.GetSpecListLength(); i++)
				spec->AddSpec(gpbs.GetSpec(i)->Clone());
			return;
		}
		//clone this sequence
		gnGenomeSpec* tmpSpec = spec->Clone();
		//crop off the end of this sequence past the insert point
		//crop off the beginning of the cloned sequence up to the insert
		gnSeqI real_offset = offset - 1;
		gnSeqI endCrop = spec->GetLength() - real_offset;
		spec->CropEnd(endCrop);
		tmpSpec->CropStart(real_offset);
		//append the inserted sequence and the end of this sequence.
		insert(GNSEQI_END, gpbs);
		insert(GNSEQI_END, *tmpSpec);
		delete tmpSpec;
	STACK_TRACE_END
}
// Concatenation operators
gnSequence const gnSequence::operator+(const gnSequence& seq) const{
	STACK_TRACE_START
		gnSequence rval(*this);
		rval.append(seq);
		return rval;
	STACK_TRACE_END
}
// Subsequences and base deletion
gnSequence gnSequence::subseq(const gnSeqI offset, const gnSeqI length) const{
	STACK_TRACE_START
		if(length == 0)
			return gnSequence();
		if(offset == 0)
			Throw_gnEx(SeqIndexOutOfBounds());
		
		gnSequence tmpSeq;
		delete tmpSeq.spec;
		tmpSeq.spec = spec->CloneRange(offset - 1, length);
		return tmpSeq;
	STACK_TRACE_END
}
void gnSequence::erase( const gnSeqI offset, const gnSeqI len ){
	STACK_TRACE_START
		if(offset == 0)
			Throw_gnEx(SeqIndexOutOfBounds());
		//range checking, etc.
		gnSeqI current_length = length();
		if(offset > current_length)
			Throw_gnEx(SeqIndexOutOfBounds());
		gnSeqI endBase = offset + len - 1 < current_length ? offset + len - 1 : current_length;
		gnSeqI startBase = offset - 1;

		gnGenomeSpec* tmpSpec = spec->CloneRange(endBase, GNSEQI_END);
		uint32 lennard = tmpSpec->GetLength();
		spec->CropEnd(current_length - startBase);
		lennard = length();
		insert(GNSEQI_END, *tmpSpec);
		delete tmpSpec;
	STACK_TRACE_END
}
void gnSequence::splitContig(const gnSeqI splitI, const uint32 contigI){
	STACK_TRACE_START
		gnSeqI splitBase = splitI;
		gnSeqI current_length = length();
		if(splitI == 0)
			Throw_gnEx(SeqIndexOutOfBounds());
		if(contigI == ALL_CONTIGS && splitI > current_length)
				Throw_gnEx(SeqIndexOutOfBounds());
		else
			localToGlobal(contigI, splitBase);
		gnGenomeSpec* tmpSpec = spec->Clone();
		tmpSpec->CropStart(splitBase);
		spec->CropEnd(current_length - splitBase);
		
		insert(GNSEQI_END, *tmpSpec);
		delete tmpSpec;
	STACK_TRACE_END
}

// IO operators
istream& operator>>(istream& is, gnSequence& gps){	//read from source.
	string bases;
	is >> bases;
	gps.append(bases);
	return is;
}
ostream& operator<<(ostream& os, const gnSequence& gps){ //write to source.
	os << gps.ToString();
	return os;
}
string gnSequence::ToString( const gnSeqI len, const gnSeqI offset ) const
{
	STACK_TRACE_START
		string str;
		ToString(str, len, offset);
		return str;
	STACK_TRACE_END
}
boolean gnSequence::ToString( string& str, const gnSeqI len, const gnSeqI offset ) const
{
	STACK_TRACE_START
		gnSeqI real_offset = offset - 1;
		gnSeqI m_length = length();
		gnSeqI readSize = len > m_length - real_offset ? m_length - real_offset : len;
		Array<char> array_buf( readSize+1 );
		char *buf = array_buf.data;
		boolean success = spec->SeqRead( real_offset, buf, readSize, ALL_CONTIGS);
		buf[readSize] = '\0';
		str = buf;

		//now filter the string
		list<const gnBaseFilter*>::const_iterator iter = filter_list.begin();
		for(; iter != filter_list.end(); iter++)
			(*iter)->Filter(str);

		if( success )
			return true;
		return false;
	STACK_TRACE_END
}
boolean gnSequence::ToArray( gnSeqC* pSeqC, gnSeqI length, const gnSeqI offset ) const
{
	STACK_TRACE_START
		gnSeqI real_offset = offset - 1;
		if(offset == GNSEQI_END)
			return false;
		gnSeqC* tmp = new gnSeqC[length];
		boolean success = spec->SeqRead( real_offset, tmp, length, ALL_CONTIGS);
		
		//now filter the array
		list<const gnBaseFilter*>::const_iterator iter = filter_list.begin();
		gnSeqC** tomp = &tmp;
		for(; iter != filter_list.end(); iter++)
			(*iter)->Filter(tomp, length);
		memcpy(pSeqC, *tomp, length);
		delete[] *tomp;
		
		if( success )
			return true;

		return false;
	STACK_TRACE_END
}
gnSeqC gnSequence::GetSeqC( const gnSeqI offset ) const
{
	STACK_TRACE_START
		char block;
		gnSeqI readLen = 1;
		gnSeqI real_offset = offset - 1;
		boolean success = spec->SeqRead( real_offset, &block, readLen, ALL_CONTIGS);

		//now filter the char
		list<const gnBaseFilter*>::const_iterator iter = filter_list.begin();
		for(; iter != filter_list.end(); iter++)
			block = (*iter)->Filter(block);

		if( success )
			return block;

		return GNSEQC_NULL;
	STACK_TRACE_END
}
gnSeqC gnSequence::operator[]( const gnSeqI offset ) const
{
	STACK_TRACE_START
		return GetSeqC( offset );
	STACK_TRACE_END
}

gnSeqI gnSequence::contigListSize() const{
	STACK_TRACE_START
		return spec->GetSpecListLength();
	STACK_TRACE_END
}
gnSeqI gnSequence::contigListLength() const{
	STACK_TRACE_START
		return spec->GetSpecListLength();
	STACK_TRACE_END
}

//find the contig associated with base
uint32 gnSequence::contigIndexByBase( const gnSeqI baseI) const{
	STACK_TRACE_START
		return spec->GetSpecIndexByBase(baseI-1);
	STACK_TRACE_END
}
gnSequence gnSequence::contig( const uint32 contigI) const{
	STACK_TRACE_START
		if(contigI == ALL_CONTIGS)
			return *this;
		return gnSequence(*spec->GetSpec(contigI));
	STACK_TRACE_END
}
//returns a gnSequence pointer containing the specified contig.
gnSequence gnSequence::contigByBase( const gnSeqI baseI) const{
	STACK_TRACE_START
		return gnSequence(*spec->GetSpecByBase(baseI-1));
	STACK_TRACE_END
}
uint32 gnSequence::contigIndexByName( string& contigName) const{
	STACK_TRACE_START
		return spec->GetSpecIndexByName(contigName);
	STACK_TRACE_END
}

gnSeqI gnSequence::contigStart( const uint32 contigI) const{
	STACK_TRACE_START
		int64 leafI = contigI;
		// add one for the stupid geneticists whose indices start at 1
		return spec->GetSpecStartBase( leafI ) + 1;
	STACK_TRACE_END
}

gnSeqI gnSequence::contigLength( const uint32 contigI) const{
	STACK_TRACE_START
		return spec->GetSpec( contigI )->GetLength();
	STACK_TRACE_END
}

string gnSequence::contigName( const uint32 contigI) const{
	STACK_TRACE_START
		return spec->GetSpec(contigI)->GetName();
	STACK_TRACE_END
}

void gnSequence::globalToLocal(uint32& contigI, gnSeqI& baseI) const{
	STACK_TRACE_START
		contigI = contigIndexByBase(baseI);
		baseI -= (contigStart(contigI) - 1);
	STACK_TRACE_END
}

void gnSequence::localToGlobal(const uint32 contigI, gnSeqI& baseI) const{
	STACK_TRACE_START
		if(baseI > contigLength(contigI))
			Throw_gnEx(SeqIndexOutOfBounds());
		baseI += contigStart(contigI) - 1;
	STACK_TRACE_END
}

void gnSequence::globalToSource(uint32& contigI, gnSeqI& baseI) const{
	STACK_TRACE_START
		baseI--;	//convert from geneticist coordinates
		gnSeqI seq_contig = spec->GetSpecIndexByBase(baseI);
		gnSeqI new_base = baseI - spec->GetSpecStartBase(seq_contig);
		gnFragmentSpec *fragment_spec = spec->GetSpec(seq_contig);
		seq_contig = fragment_spec->GetSpecIndexByBase(new_base);
		new_base = new_base - fragment_spec->GetSpecStartBase(seq_contig);
		gnContigSpec* contig_spec = fragment_spec->GetSpec(seq_contig);
		contigI = contig_spec->GetSourceContigIndex();
		gnSeqI contig_start = contig_spec->GetStart();
		if(contig_spec->IsReverseComplement()){
			gnSeqI source_len = contig_spec->GetSourceLength();
			//it may seem counter intuitive, but we can add the source_len and mod by source_len
			//to deal with circularity
	//		contig_start = (contig_start - contig_spec->GetLength() + source_len) % source_len;
			//if we are reverse complement then the contig_start will be the ending base pair
			//we want to subtract new_base from it. and we can add/mod by source_len to deal
			//with circles.
			baseI = (contig_start - new_base - 1 + source_len) % source_len;
		}else
			baseI =  contig_start + new_base + 1;
	STACK_TRACE_END
}

void gnSequence::localToSource(uint32& contigI, gnSeqI& baseI) const{
	STACK_TRACE_START
		localToGlobal(contigI, baseI);
		globalToSource(contigI, baseI);
	STACK_TRACE_END
}

gnSequence gnSequence::contigByName( string& contigName) const{
	STACK_TRACE_START
		uint32 contigIndex = spec->GetSpecIndexByName(contigName);
		return gnSequence(*spec->GetSpec(contigIndex));
	STACK_TRACE_END
}

void gnSequence::setContigName( const uint32 contigI, const string& contig_name) {
	STACK_TRACE_START
		if(contigI == ALL_CONTIGS)
			spec->SetName(contig_name);
		else
			spec->GetSpec(contigI)->SetName(contig_name);
	STACK_TRACE_END
}

void gnSequence::setReverseComplement( const boolean revComp, const uint32 contigI){
	STACK_TRACE_START
		if(contigI == ALL_CONTIGS)
			spec->SetReverseComplement(revComp);
		else{
			gnFragmentSpec* contig = spec->GetSpec(contigI);
			contig->SetReverseComplement(revComp);
		}
	STACK_TRACE_END
}

boolean gnSequence::isReverseComplement( const uint32 contigI ){
	STACK_TRACE_START
		if(contigI == ALL_CONTIGS)
			return spec->IsReverseComplement();
		gnFragmentSpec* contig = spec->GetSpec(contigI);
		return contig->IsReverseComplement();
	STACK_TRACE_END
}

uint32 gnSequence::getHeaderListLength(const uint32 contigI) const{
	STACK_TRACE_START
		if(contigI == ALL_CONTIGS)
			return spec->GetHeaderListLength();
		else{
			return spec->GetSpec(contigI)->GetHeaderListLength();
		}
	STACK_TRACE_END
}

gnBaseHeader* gnSequence::getHeader(const uint32 contigI, const uint32 headerI) const{
	STACK_TRACE_START
		if(contigI == ALL_CONTIGS)
			return spec->GetHeader(headerI);
		else{
			return spec->GetSpec(contigI)->GetHeader(headerI);
		}
	STACK_TRACE_END
}

void gnSequence::addHeader(const uint32 contigI, gnBaseHeader* header, const uint32 headerI){
	STACK_TRACE_START
		if(contigI == ALL_CONTIGS)
			spec->AddHeader(header, headerI);
		else{
			spec->GetSpec(contigI)->AddHeader(header, headerI);
		}
	STACK_TRACE_END
}

void gnSequence::removeHeader(const uint32 contigI, const uint32 headerI){
	STACK_TRACE_START
		if(contigI == ALL_CONTIGS)
			spec->RemoveHeader(headerI);
		else
			spec->GetSpec(contigI)->RemoveHeader(headerI);
	STACK_TRACE_END
}

void gnSequence::merge(const gnSeqI startBase, const gnSeqI endBase){
	STACK_TRACE_START
	STACK_TRACE_END
}

void gnSequence::mergeContigs(const uint32 startC, const uint32 endC){
	STACK_TRACE_START
		spec->MergeFragments(startC, endC);
	STACK_TRACE_END
}

bool gnSequence::LoadSource(const string sourcename){
	STACK_TRACE_START
		gnSourceFactory* m_sourceFactory = gnSourceFactory::GetSourceFactory();
		gnBaseSource *m_pSource = m_sourceFactory->AddSource(sourcename, true);
		if (m_pSource!=NULL)
		{
			if(spec != NULL)
				delete spec;
			spec = m_pSource->GetSpec();
			return true;
		}
		return false;
	STACK_TRACE_END
}

gnSeqI gnSequence::find(const gnSequence& search, const gnSeqI offset)const{
	STACK_TRACE_START
		//this is really adHoc... should probably be switched
		string searchIn=ToString();
		string find= search.ToString();
		string::size_type pos = searchIn.find(find, offset);
		if (pos == string::npos)
			return GNSEQI_ERROR;
		else{
			return pos++; 
		}
	STACK_TRACE_END
}

}	// end namespace genome

