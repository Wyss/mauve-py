/////////////////////////////////////////////////////////////////////////////
// File:            gnGBKSource.h
// Purpose:         Implements gnBaseSource for GenBank sequences
// Description:     
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
#include "libGenome/gnFeature.h"
#include "libGenome/gnGBKSource.h"
#include "libGenome/gnSourceSpec.h"
#include "libGenome/gnSourceHeader.h"
#include "libGenome/gnSourceQualifier.h"
#include "libGenome/gnLocation.h"
#include "libGenome/gnStringTools.h"
#include "libGenome/gnDebug.h"
#include "libGenome/gnStringQualifier.h"
#include <string>
#include <cstring>


using namespace std;
namespace genome {


gnGBKSource::gnGBKSource()
{
	m_openString = "";
	m_pFilter = gnFilter::proteinSeqFilter();
	if(m_pFilter == NULL){
		DebugMsg("Error using static sequence filters.");
	}
}
gnGBKSource::gnGBKSource( const gnGBKSource& s ) : gnFileSource(s)
{
	vector< gnFileContig* >::const_iterator iter = s.m_contigList.begin();
	for( ; iter != s.m_contigList.end(); ++iter )
	{
		m_contigList.push_back( (*iter)->Clone() );
	}
}
gnGBKSource::~gnGBKSource()
{
	m_ifstream.close();
	vector< gnFileContig* >::iterator iter = m_contigList.begin();
	for( ; iter != m_contigList.end(); ++iter )
	{
		gnFileContig* fg = *iter;
		*iter = 0;
		delete fg;
	}
}
boolean gnGBKSource::HasContig( const string& name ) const
{
	for(uint32 i = 0 ; i <= m_contigList.size(); i++ )
	{
		if( name == m_contigList[i]->GetName() )
			return true;
	}
	return false;
}
uint32 gnGBKSource::GetContigID( const string& name ) const
{
	for(uint32 i = 0 ; i <= m_contigList.size(); i++ )
	{
		if( name == m_contigList[i]->GetName() )
			return i;
	}
	return ALL_CONTIGS;
}
string gnGBKSource::GetContigName( const uint32 i ) const
{
	if( i < m_contigList.size() )
	{
		return m_contigList[i]->GetName();
	}
	return "";
}
gnSeqI gnGBKSource::GetContigSeqLength( const uint32 i ) const
{
	if( i == ALL_CONTIGS)
		return m_spec->GetLength();
	if( i < m_contigList.size() )
	{
		return m_contigList[i]->GetSeqLength();
	}
	return GNSEQI_ERROR;
}

boolean gnGBKSource::SeqRead( const gnSeqI start, char* buf, gnSeqI& bufLen, const uint32 contigI )
{
	boolean result = false;
#pragma omp critical
{
	result = SeqReadImpl( start, buf, bufLen, contigI );
}
	return result;
}

boolean gnGBKSource::SeqReadImpl( const gnSeqI start, char* buf, gnSeqI& bufLen, const uint32 contigI ){

	uint64 startPos = 0;
	uint64 readableBytes = 0;
	if( !SeqSeek( start, contigI, startPos, readableBytes ) )
	{
		bufLen = 0;
		return false;
	}
	
	if( contigI == ALL_CONTIGS )
	{
		uint32 curLen = 0;
		uint64 bytesRead = 0;
		while (curLen < bufLen)
		{
//SeqSeek to start, Figure out how much can be read before SeqSeeking again.
			if(readableBytes <= 0)	//Look out for zero length contigs!  IMPLEMENT ME
				if( !SeqSeek( start + curLen, contigI, startPos, readableBytes ) ){
					bufLen = curLen;
					return true;
				}
			//readLen is the amount to read on this pass
			uint64 readLen = (bufLen - curLen) < readableBytes ? (bufLen - curLen) : readableBytes;	
			Array<gnSeqC> array_buf( readLen );
			gnSeqC* tmpBuf = array_buf.data;

			// read chars and filter
			m_ifstream.read(tmpBuf, readLen);
			uint64 gc = m_ifstream.gcount();
			bytesRead += gc;
			readableBytes -= gc;
			for(uint32 i=0; i < gc; i++){
				if( m_pFilter->IsValid(tmpBuf[i]) ){
					buf[curLen] = tmpBuf[i];
					curLen++;
				}
			}
			if(m_ifstream.eof()){	//we hit the end of the file.  bail out.
				m_ifstream.clear();
				bufLen = curLen;
				return true;
			}
		}
		bufLen = curLen;
	}
	else if( contigI < m_contigList.size() )
	{
		uint32 curLen = 0;
		//check to see if the buffer is bigger than the contig.  if so truncate it.
		gnSeqI contigSize = m_contigList[contigI]->GetSeqLength();
		bufLen = bufLen < contigSize ? bufLen : contigSize;
		while (curLen < bufLen)
		{
			uint64 readLen = bufLen - curLen;	//the amount to read on this pass
			Array<gnSeqC> array_buf( readLen );
			gnSeqC* tmpBuf = array_buf.data;

			// read chars and filter
			m_ifstream.read(tmpBuf, readLen);
			uint64 gc = m_ifstream.gcount();
//			cout << "Read " << gc << " chars from " << m_openString << "\n";
//			cout << "Checking character validity on: " << tmpBuf << "\n";
			for(uint32 i=0; i < gc; i++){
				if( m_pFilter->IsValid(tmpBuf[i]) ){
					buf[curLen] = tmpBuf[i];
					curLen++;
				}
			}
			if(m_ifstream.eof()){	//we hit the end of the file.  bail out.
				m_ifstream.clear();
				bufLen = curLen;
				return true;
			}
		}
		bufLen = curLen;
	}
	return true;

}
// private:
// figures out which contig the sequence starts at then calls SeqStartPos to get the offset within that contig
// returns startPos, the file offset where the sequence starts
// returns true if successful, false otherwise
boolean gnGBKSource::SeqSeek( const gnSeqI start, const uint32& contigI, uint64& startPos, uint64& readableBytes )
{
	if( contigI == ALL_CONTIGS )
	{
		// find first contig
		gnSeqI curIndex = 0;
		vector< gnFileContig* >::iterator iter = m_contigList.begin();
		for( ; iter != m_contigList.end(); ++iter )
		{
			uint64 len = (*iter)->GetSeqLength();
			if( (curIndex + len) > start )
				break;
			curIndex += len;
		}
		if( iter == m_contigList.end() )
			return false;
		// seek to start
		gnSeqI startIndex = start - curIndex;  //startIndex is starting pos. within the contig
		return SeqStartPos( startIndex, *(*iter), startPos, readableBytes );
	}
	else if( contigI < m_contigList.size() )
	{
		return SeqStartPos( start, *(m_contigList[contigI]), startPos, readableBytes );
	}
	return false;
}
//Returns startPos, the file offset where the sequence starts.
boolean gnGBKSource::SeqStartPos( const gnSeqI start, gnFileContig& contig, uint64& startPos, uint64& readableBytes )
{
	readableBytes = 0;
	uint32 curLen = 0;
	//seek to the file offset where the contig starts
	startPos = contig.GetSectStartEnd(gnContigSequence).first;	//set startPos to start where the contig starts
	m_ifstream.seekg( startPos, ios::beg );
	if( m_ifstream.eof() ){
		ErrorMsg("ERROR in gnGBKSource::Incorrect contig start position, End of file reached!\n");
		return false;
	}
	while( true )
	{
		  // READ the rest of the contig skipping over invalid characters until we get to the starting base pair.
		  // startPos will contain the file offset with the starting base pair
		uint32 tmpbufsize = contig.GetSectStartEnd(gnContigSequence).second - startPos;
		if(tmpbufsize == 0){
			ErrorMsg("ERROR in gnGBKSource: stored contig size is incorrect.");
			return false;
		}
		uint64 startOffset = start;
		if(contig.HasRepeatSeqGap()){	//check for sequence integrity
			startOffset += (9 + m_newlineSize) * (start / 60 + 1) + start / 10 + 1;
			if( m_newlineSize == 2 )	// test compensation for strange bug
				startOffset--;
			startPos+=startOffset;
			m_ifstream.seekg(startPos , ios::beg);
			readableBytes = contig.GetSectStartEnd(gnContigSequence).second - startPos;
			return true;
		}

		//sequence is corrupt.  read in base by base
		tmpbufsize = tmpbufsize < BUFFER_SIZE ? tmpbufsize : BUFFER_SIZE;  //read in the smaller of the two.
		Array<char> array_buf( tmpbufsize );
		char* tmpbuf = array_buf.data;

		m_ifstream.read( tmpbuf, tmpbufsize );
		if( m_ifstream.eof() ){
			ErrorMsg("ERROR in gnGBKSource::Read End of file reached!\n");
			return false;
		}
		for( uint32 i=0; i < tmpbufsize; ++i ){
			if( m_pFilter->IsValid(tmpbuf[i]) ){
				if( curLen >= start ){ //stop when we reach the starting offset within this contig
					startPos += i;
					m_ifstream.seekg( startPos, ios::beg );  //seek to startPos
					readableBytes = contig.GetSectStartEnd(gnContigSequence).second - startPos;
					return true;
				}
				++curLen;  //each time we read a valid b.p., increment the sequence length
			}
		}
		startPos += tmpbufsize;
	}
	return true;
}

void gnGBKSource::FormatString(string& data, uint32 offset, uint32 width){
	//first remove newlines and corresponding whitespace
	string::size_type newline_loc = data.find_first_of('\n', 0);
	while(newline_loc != string::npos){
		if(data[newline_loc-1] == '\r')
			newline_loc--;
		string::size_type text_loc = newline_loc;
		while((data[text_loc] == ' ') ||(data[text_loc] == '	')||(data[text_loc] == '\n')||(data[text_loc] == '\r')){
			text_loc++;
			if(text_loc+1 == data.length())
				break;
		}
		data = (data.substr(0, newline_loc) + " " + data.substr(text_loc));
		newline_loc = data.find_first_of('\n', 0);
	}
	//now reformat with newlines and whitespace, observing word boundaries...
	string output_string = "";
	for(uint32 charI = 0; charI < data.length();){
		//get the substring to append and increment charI
		string::size_type base_loc = charI;
		string append_string;
		while(base_loc - charI <= width){
			string::size_type space_loc = data.find_first_of(' ', base_loc+1);
			if(space_loc - charI < width)
				base_loc = space_loc;
			else if(base_loc == charI){
				//word is too big for one line.  split it.
				append_string = data.substr(charI, width);
				charI+=width;
			}else{
				append_string = data.substr(charI, base_loc - charI);
				charI = base_loc;
			}
		}
		output_string += string(offset, ' ') + append_string;
		if(charI + width < data.length())
			output_string += "\r\n";
	}
	data = output_string;
}

template< class SubSpec >
void WriteHeader(gnMultiSpec< SubSpec >* spec, const string& hdr, ofstream& m_ofstream) {
	gnBaseHeader* gpbh = NULL;
	uint32 header_index = 0;
	try{
		do{
			gpbh = spec->GetHeader(hdr, header_index);
			m_ofstream << gpbh->GetHeader();
			header_index++;
		}while(gpbh != NULL);
	}catch(gnException& gne){}
}

boolean gnGBKSource::Write(gnSequence& seq, const string& filename){
	ofstream m_ofstream(filename.c_str(), ios::out | ios::binary);
	if(!m_ofstream.is_open())
		return false;

	string newline = "\r\n";
	gnGenomeSpec* spec = seq.GetSpec();

	// output general file header first if one exists.
	if(spec->GetHeaderListLength() == 1){
		gnBaseHeader *gpbh = spec->GetHeader(0);
		string name = gpbh->GetHeaderName();
		//IMPLEMENT ME!  Is platform specific newline substitution necessary?
		if(string::npos != name.find(".SEQ")){
			string header = gpbh->GetHeader();
			m_ofstream << header;
		}
	}
	// TODO:  Figure out where the buffer overflow is and reduce
	// this back to BUFFER_SIZE -- overflow appears not to corrupt
	// sequence data
	Array<gnSeqC> array_buf( 2 * BUFFER_SIZE );
	gnSeqC *bases = array_buf.data;

	for(uint32 specI = 0; specI < spec->GetSpecListLength(); specI++){
		gnFragmentSpec* subSpec = spec->GetSpec(specI);
		
		//write out contig headers.  start with LOCUS...
		m_ofstream << "LOCUS       ";
		//write Locus Name
		string contigName = subSpec->GetName();
		if(contigName.length() > SEQ_LOCUS_NAME_LENGTH)
			contigName = contigName.substr(0, SEQ_LOCUS_NAME_LENGTH);
		uint32 filler_size = SEQ_LOCUS_NAME_LENGTH - contigName.length();
		m_ofstream << contigName << string(filler_size, ' ');
		//write Locus Length
		string length_string = uintToString(subSpec->GetLength());
		filler_size = SEQ_LOCUS_SIZE_LENGTH - length_string.size();
		m_ofstream << string(filler_size, ' ') << length_string << " bp ";
		//write dnatype
		string dnatype = string(SEQ_LOCUS_DNATYPE_LENGTH, ' ');
		uint32 head_look_i = 0;
		gnBaseHeader* gpbh = NULL;
		try{
			gpbh = subSpec->GetHeader("LOCUS", head_look_i);
		}catch(gnException& gne){}
		if( gpbh != NULL )
			dnatype = gpbh->GetHeader().substr(SEQ_LOCUS_DNATYPE_OFFSET, SEQ_LOCUS_DNATYPE_LENGTH);
		m_ofstream << dnatype << string(2, ' ');
		//write circularity
		string circular = subSpec->IsCircular() ? string("circular  ") : string(10, ' ');
		m_ofstream << circular;
		//write division code
		string division = string(SEQ_LOCUS_DIVCODE_LENGTH, ' ');
		if(gpbh != NULL)
			division = gpbh->GetHeader().substr(SEQ_LOCUS_DIVCODE_OFFSET, SEQ_LOCUS_DIVCODE_LENGTH);
		m_ofstream << division;
		//write date -- IMPLEMENT ME!  format the real date to spec! dd-mmm-yyyy
		string date = string(SEQ_LOCUS_DATE_LENGTH, ' ');
		if(gpbh != NULL)
			date = gpbh->GetHeader().substr(SEQ_LOCUS_DATE_OFFSET, SEQ_LOCUS_DATE_LENGTH);
		m_ofstream << string(7, ' ') << date << "\r\n";
		
		//write out the rest of the headers if they were supplied!
		WriteHeader(subSpec, "DEFINITION", m_ofstream);
		WriteHeader(subSpec, "ACCESSION", m_ofstream);
		WriteHeader(subSpec, "VERSION", m_ofstream);
		WriteHeader(subSpec, "KEYWORDS", m_ofstream);
		WriteHeader(subSpec, "SEGMENT", m_ofstream);
		WriteHeader(subSpec, "SOURCE", m_ofstream);
		WriteHeader(subSpec, "REFERENCE", m_ofstream);
		WriteHeader(subSpec, "COMMENT", m_ofstream);

		//write out feature table!
		m_ofstream << "FEATURES             Location/Qualifiers" << "\r\n";
		for(uint32 featureI = 0; featureI < subSpec->GetFeatureListLength(); featureI++){
			//write a feature tag
			gnBaseFeature *gpmf = subSpec->GetFeature(featureI);
			string featureName = gpmf->GetName();
			m_ofstream << string(SEQ_SUBTAG_COLUMN, ' ') << featureName;
			m_ofstream << string(SEQ_FEATURE_LOC_OFFSET - featureName.length() - SEQ_SUBTAG_COLUMN, ' ');
			//ready to output location b.s.
			uint32 location_count = gpmf->GetLocationListLength();
			uint32 line_pos = SEQ_FEATURE_LOC_OFFSET;
			uint32 parenthesis_count = 0;
			if(location_count > 1){
				m_ofstream << "join(";
				line_pos += 5;
				parenthesis_count++;
			}
			gnLocation::gnLocationType loc_type = gpmf->GetLocationType();
			switch(loc_type){
				case gnLocation::LT_Standard:
					break;
				case gnLocation::LT_Complement:
					m_ofstream << "complement(";
					line_pos += 11;
					parenthesis_count++;
					break;
				case gnLocation::LT_Order:
					m_ofstream << "order(";
					line_pos += 6;
					parenthesis_count++;
					break;
				case gnLocation::LT_Group:
					m_ofstream << "group(";
					parenthesis_count++;
					line_pos += 6;
					break;
				case gnLocation::LT_OneOf:
					m_ofstream << "one-of(";
					parenthesis_count++;
					line_pos += 7;
					break;				
				default:
					break;
			}
			//create the location string, then see if it will fit on the line
			string location;
			for(uint32 locationI = 0; locationI < location_count; locationI++){
				gnLocation gpl = gpmf->GetLocation(locationI);
				if(gpl.IsStartBoundLonger())
					location += ">";
				if(gpl.IsStartBoundShorter())
					location += "<";
				location += uintToString(gpl.GetStart());
				gnSeqI end_loc = gpl.GetEnd();
				if(end_loc != 0){
					switch(gpl.GetType()){
						case gnLocation::LT_BetweenBases:
							location += "^";
							break;
						case gnLocation::LT_OneOf:
							location += ".";
							break;
						default:
							location += "..";
							break;
					}
					if(gpl.IsEndBoundShorter())
						location += "<";
					if(gpl.IsEndBoundLonger())
						location += ">";
					location+= uintToString(end_loc);
				}
				if(locationI +1 < location_count)
					location += ",";
				else{	//append necessary parenthesis
					for(;parenthesis_count > 0; parenthesis_count--)
						location += ")";
				}
				//put it on this line if it fits.  otherwise make a new line.
				if(line_pos + location.length() < SEQ_COLUMN_WIDTH){
					m_ofstream << location;
					line_pos += location.length();
				}else{
					m_ofstream << "\r\n" << string(SEQ_FEATURE_LOC_OFFSET, ' ') << location;
					line_pos = SEQ_FEATURE_LOC_OFFSET + location.length();
				}
				location = "";
			}
			m_ofstream << "\r\n";
			//now output qualifiers!  yaay!
			//god damn this is a big ugly piece of code.
			
			uint32 qualifier_count = gpmf->GetQualifierListLength();
			for(uint32 qualifierI = 0; qualifierI < qualifier_count; qualifierI++){
				m_ofstream << string(SEQ_FEATURE_LOC_OFFSET, ' ');
				gnBaseQualifier* qualifier = gpmf->GetQualifier(qualifierI);
				m_ofstream << "/" << qualifier->GetName() << "=";
				//IMPLEMENT ME! do a better word wrap on this bitch.
				string qually = string(qualifier->GetValue());
//				FormatString(qually, SEQ_FEATURE_LOC_OFFSET, 80 - SEQ_FEATURE_LOC_OFFSET);
//				qually = qually.substr(SEQ_FEATURE_LOC_OFFSET);
				m_ofstream << qually << "\r\n";
			}
			if(gpmf != NULL)
				delete gpmf;
		}
		
		//get information about the sequence we're writing out.
		gnSeqI readOffset = seq.contigStart(specI);
		gnSeqI readLength = seq.contigLength(specI);

		//finally - output base count and origin
		m_ofstream << "BASE COUNT ";
		gnSeqI a_count=0, c_count=0, g_count=0, t_count=0, other_count=0;
		gnSeqI countLen = readLength + readOffset;
		for(gnSeqI countI = readOffset; countI < countLen;){
			gnSeqI writeLen = countLen - countI < BUFFER_SIZE ? countLen - countI : BUFFER_SIZE;
			if(!seq.ToArray(bases, writeLen, countI))
				return false;
			gnSeqI a, c, g, t, other;
			BaseCount(string(bases, writeLen), a, c, g, t, other);
			a_count += a;
			c_count += c;
			g_count += g;
			t_count += t;
			other_count += other;
			countI += writeLen;
		}
		m_ofstream << uintToString(a_count) << " a ";
		m_ofstream << uintToString(c_count) << " c ";
		m_ofstream << uintToString(g_count) << " g ";
		m_ofstream << uintToString(t_count) << " t ";
		m_ofstream << uintToString(other_count) << " others" << "\r\n";

		string origin = "ORIGIN\r\n";
		head_look_i = 0;
		try{
			gpbh = subSpec->GetHeader("ORIGIN", head_look_i);
			origin = gpbh->GetHeader();
			m_ofstream << origin;
		}catch(gnException& gne){
			m_ofstream << "ORIGIN" << endl;
		}
		//write out the sequence
		gnSeqI contig_bases = 0;
		while(readLength > 0){	//buffer the read/writes
			gnSeqI writeLen = readLength < BUFFER_SIZE + 20 ? readLength : BUFFER_SIZE + 20;
			boolean success = seq.ToArray(bases, writeLen, readOffset);
			if(!success)
				return false;
			//print each 60 on their own lines...
			for(gnSeqI curbaseI = 0; curbaseI < writeLen; curbaseI += 60){
				string baseIndexStr = uintToString(contig_bases + curbaseI +1);
				m_ofstream << string(SEQ_BASES_INDEX_END - baseIndexStr.length(), ' ');
				m_ofstream << baseIndexStr;
				for(gnSeqI base_offset = 0; base_offset <= 50; base_offset+=10){
					if(writeLen <= curbaseI + base_offset)
						break;
					int64 print_length = writeLen - (curbaseI + base_offset);
					print_length = print_length > 10 ? 10 : print_length;
					m_ofstream << ' ' << string(bases + curbaseI + base_offset, print_length);
				}
				m_ofstream << "\r\n";
			}
			readLength -= writeLen;
			readOffset += writeLen;
			contig_bases += writeLen;
		}
		m_ofstream << "//\r\n";
	}
	
	m_ofstream.close();
	return true;
}

gnFileContig* gnGBKSource::GetFileContig( const uint32 contigI ) const{
	if(m_contigList.size() > contigI)
		return m_contigList[contigI];
	return NULL;
}

//File parsing access routine
boolean gnGBKSource::ParseStream( istream& fin )
{
	// INIT temp varables
	uint32 readState = 0;
	uint32 lineStart = 0;
	// INIT buffer
	uint32 sectionStart = 0;
	uint64 streamPos = 0;
	uint64 bufReadLen = 0;
	uint64 remainingBuffer = 0;
	Array<char> array_buf( BUFFER_SIZE );
	char* buf = array_buf.data;
	gnFragmentSpec* curFrag = 0;
	gnSourceSpec* curSpec = 0;
	gnSourceHeader *curHeader;
	gnFeature* curFeature;
	gnFileContig* curContig = 0;
	gnLocation::gnLocationType curBaseLocationType;
	gnSeqI curLocationStart;
	int32 curStartLength = 0;
	int32 curEndLength = 0;
	string curLocContig = "";
	string curQualifierName;
	uint64 curQualifierStart;
	string curContigName = "";
	gnSeqI seqLength = 0;
	gnSeqI seqChunk, seqChunkCount, gapChunk;
	boolean corruptWarning = false;
	
	//decide what type of newlines we have
	DetermineNewlineType();

	m_spec = new gnGenomeSpec();
	while( !fin.eof() )
	{
		size_t newstart = sectionStart < lineStart ? sectionStart : lineStart;
		if(sectionStart > 0){
			if(readState == 14)
				sectionStart = lineStart;
			else if( sectionStart >= lineStart )
				sectionStart -= lineStart;
			remainingBuffer = bufReadLen - newstart;
			memmove(buf, buf+newstart, remainingBuffer);
		}
		  // read chars
		fin.read( buf + remainingBuffer, BUFFER_SIZE - remainingBuffer);
		streamPos -= remainingBuffer;
		lineStart -= newstart;
		bufReadLen = fin.gcount();
		bufReadLen += remainingBuffer;
		
		for( uint32 i=remainingBuffer ; i < bufReadLen ; i++ )
		{
			char ch = buf[i];
			switch( readState )
			{
				case 0: 	//Assume we are in header at the start of a new line.  
							//Look for keywords starting in column 1
					if((ch == '\n')&&(buf[lineStart] != ' ')&&(buf[lineStart] != '\t') && buf[lineStart] != '\r' && buf[lineStart] != '\n'){  //not equal to whitespace
						if(curSpec == NULL){
							curSpec = new gnSourceSpec(this, m_spec->GetSpecListLength());
							curFrag = new gnFragmentSpec();
							curFrag->AddSpec(curSpec);
							curSpec->SetSourceName(m_openString);
							m_spec->AddSpec(curFrag);
						}
						if(lineStart != sectionStart){	//Add the previous header to our list
							uint32 j = SEQ_HEADER_NAME_LENGTH-1;
							for(; j > 0; j--)	
								if((buf[sectionStart+j] != ' ')&&(buf[sectionStart+j] != '	'))
									break;
							string header_name = string(buf+sectionStart, j+1);
							curHeader = new gnSourceHeader(this, header_name, sectionStart + streamPos, lineStart - sectionStart);
							//if this is header info _before_ a locus statement then its a general file header.
							if(strncmp(&buf[lineStart], "LOCUS", 5) == 0)
								m_spec->AddHeader(curHeader);
							else	//otherwise its a fragment header.
								curFrag->AddHeader(curHeader);
							sectionStart = lineStart;
						}
						
						if(strncmp(&buf[lineStart], "FEATURES", 8) == 0){
							sectionStart = i + 1;
							readState = 1;  //read in features
						}else if(strncmp(&buf[lineStart], "ORIGIN", 6) == 0){
							curHeader = new gnSourceHeader(this, string("ORIGIN"), sectionStart + streamPos, i - sectionStart + 1);
							curFrag->AddHeader(curHeader);
							curContig = new gnFileContig();
							curContig->SetName(curContigName);
							curContigName = "";
							readState = 13;  //read in base pairs
						}else if(strncmp(&buf[lineStart], "LOCUS", 5) == 0){
							if(strncmp(&buf[lineStart+SEQ_LOCUS_CIRCULAR_COLUMN-1], "circular", 8) == 0)
								curFrag->SetCircular(true);
							string locus_line = string(buf+lineStart, 80);
							if( locus_line.find(" DNA ") == string::npos &&
								locus_line.find(" dna ") == string::npos )
								m_pFilter = gnFilter::proteinSeqFilter();
							else
								m_pFilter = gnFilter::fullDNASeqFilter();

							uint32 j = SEQ_LOCUS_NAME_LENGTH+1;
							for(; j > 0; j--)	
								if((buf[lineStart+SEQ_LOCUS_NAME_COLUMN+j-1] != ' ')&&(buf[sectionStart+SEQ_LOCUS_NAME_COLUMN+j-1] != '	'))
									break;
							curContigName = string(buf+lineStart+SEQ_LOCUS_NAME_COLUMN-1, j+1);
							curFrag->SetName(curContigName);
						}
					}
					if(ch == '\n'){
						lineStart = i + 1;
					}
					break;
				case 1:	//look for feature tag in column six.  ignore whitespace before feature.
					if((ch == ' ')||(ch == '	')){
						break;
					}else if(ch == '\n'){
						lineStart = i + 1;
						sectionStart = i + 1;
						break;
					}else if(sectionStart == i){ //there was no whitespace, we hit a TAG instead
						i--;
						readState = 0; //Deal with a Header TAG
						sectionStart = i + 1;
						break;
					}else if((i - lineStart == SEQ_SUBTAG_COLUMN)||((buf[lineStart]=='	')&&(i==lineStart+1))){
						sectionStart = i;
						readState = 2;
					} //
				case 2:  //Get the feature name.  stop on whitespace
					if((ch == ' ')||(ch == '	')){
						string featureName(buf+sectionStart, i - sectionStart);
						curFeature = new gnFeature(featureName);
						curFrag->AddFeature(curFeature);
						sectionStart = i + 1;
						readState = 3;
					}
					break;
				case 3:   //Ignore whitespace before feature location
					if((ch == ' ')||(ch == '	')){
						break;
					}else if((ch == '\r')||(ch == '\n')){
						lineStart = i+1;
						break;
					}
					sectionStart = i;
					readState = 4;
				// KNOWN BUG HERE!
				// if JOIN is outside a group of complemented coordinates the feature type
				// is incorrectly set to LT_COMPLEMENT.  Instead each location's type should be set to LT_COMPLEMENT 
				case 4:		//Read a location start.  stop on (<.:^ and whitespace
					if((ch == ' ')||(ch == '	')||(ch == '(')||(ch == '.')||(ch=='^')||(ch==':')){
						string starter(buf+sectionStart, i - sectionStart);
						if(ch == '('){
							if(starter == "complement")
								curFeature->SetLocationType(gnLocation::LT_Complement);
							else if(starter == "order")
								curFeature->SetLocationType(gnLocation::LT_Order);
							else if(starter == "group")
								curFeature->SetLocationType(gnLocation::LT_Group);
							else if(starter == "one-of")
								curFeature->SetLocationType(gnLocation::LT_OneOf);
							sectionStart = i + 1;	//ignore join since it is default.
							break;
						}else if(ch == ':'){
							curLocContig = starter;
							sectionStart = i + 1;
							break;
						}
						curLocationStart = atoi(starter.c_str());
						readState = 6;	//read in end base by default.
						if(ch == '.'){
							//go to special state to look for another one.
							readState = 5;
							sectionStart = i + 1;
							break;
						}else if(ch == '^'){
							curBaseLocationType = gnLocation::LT_BetweenBases;
						}else if((ch == ' ')||(ch == '	')){
							//no end location go to qualifier
							gnLocation curLocation(curLocationStart, curLocationStart);
							curFeature->AddLocation(curLocation, curFeature->GetLocationListLength());
							readState = 7;
						}
						sectionStart = i + 1;

					}else if(ch == '<'){
						curStartLength = -1;
						sectionStart = i + 1;
					}else if(ch == '>'){
						curStartLength = 1;
						sectionStart = i + 1;
					}
					break;
				case 5: //look for another period or location start.
					if(ch == '.'){
						curBaseLocationType = gnLocation::LT_Standard;
						readState = 6;
						sectionStart = i + 1;
						break;
					}
					curBaseLocationType = gnLocation::LT_OneOf;
				case 6:	//see if there's a second location value.  stop on >, and whitespace
					if(ch == '>'){
						curEndLength = 1;
						sectionStart = i + 1;
					}else if(ch == '<'){
						curEndLength = -1;
						sectionStart = i + 1;
					}else if((ch == ' ')||(ch == '	')||(ch == ',')||(ch == '\r')||(ch == '\n')){
						//read end location
						string ender(buf+sectionStart, i - sectionStart);
						gnSeqI curLocationEnd = atoi(ender.c_str());
						gnLocation curLocation(curLocationStart, curStartLength, curLocationEnd, curEndLength, curBaseLocationType);
						curEndLength = 0;
						curStartLength = 0;
						curFeature->AddLocation(curLocation, curFeature->GetLocationListLength());
						readState = ch == ',' ? 3 : 7;  //read another loc if we need to.
						sectionStart = i+1;
						if( ch == '\n' )
							lineStart = i+1;
					}
					break;
				case 7:  //skip to start of qualifier
					if((ch != ' ')&&(ch != '	')&&(lineStart == i)){
						sectionStart = i;	// Hit a tag.  go deal with it.
						readState = 0;
						i--;
					}else if((ch != ' ')&&(ch != '	')&&((lineStart == i - SEQ_SUBTAG_COLUMN)||((buf[lineStart]=='	')&&(i==lineStart+1)))){
						sectionStart = i;	// Hit a feature.  go deal with it.
						readState = 2;
						i--;
					}else if(ch == ','){  //oops!  another location to read!
						sectionStart = i+1;
						readState = 3;
					}else if(ch == '/'){  //finally, a qualifier.
						sectionStart = i+1;
						readState = 8;
					}else if(ch == '\n')
						lineStart = i + 1;
					break;
				case 8:		//get a qualifier, stop on =
					if(ch == '='){
						curQualifierName = string(buf+sectionStart, i - sectionStart);
						readState = 9;
						sectionStart = i+1;
					}else if( ch == '\r' || ch == '\n' ){
						// this is a value-less qualifier
						curQualifierName = string(buf+sectionStart, i - sectionStart);
						curFeature->AddQualifier( new gnStringQualifier( curQualifierName, "" ));
						readState = 7;
						sectionStart = i+1;
						if( ch == '\n' )
							lineStart = i + 1;
					}
					break;
				case 9:		//are we getting a string? look for " or [
					if(ch == '"'){
						readState = 10;
						sectionStart = i;
						curQualifierStart = i + streamPos;
					}else if(ch == '['){
						readState = 12;
						sectionStart = i;
					}else if((ch == '\r')||(ch == '\n')){
						curFeature->AddQualifier(new gnSourceQualifier(this, curQualifierName, sectionStart + streamPos, i - sectionStart));
						sectionStart = i+1;
						readState = 7; //look for another qualifier
						if( ch == '\n' )
							lineStart = i + 1;
					}
					break;
				case 10:		//read until the end of the quotation. look out for escaped quotes
					if(ch == '"')
						readState = 11;
					if(ch == '\n'){
						lineStart = i + 1;
					}
					break;
				case 11:
					if(ch != '"'){
						gnSourceQualifier* gnsq = new gnSourceQualifier(this, curQualifierName, curQualifierStart, i - sectionStart);
						curFeature->AddQualifier(gnsq);
						sectionStart = i+1;
						readState = 7;	//look for another qualifier.
						if(ch == '\n')
							lineStart = i + 1;
					}else
						readState = 10;  //quote was escaped.  look for another.
					break;
				case 12:
					if(ch == ']'){
						curFeature->AddQualifier(new gnSourceQualifier(this, curQualifierName, sectionStart + streamPos, i - sectionStart));
						sectionStart = i+1;
						readState = 7;	//look for another qualifier.
					}
					break;
				case 13:	//start the sequence read.
					curContig->SetSectStart(gnContigSequence, i - 1 + streamPos);
					curContig->SetRepeatSeqGap(true);
					seqChunk = 0;
					seqChunkCount = 0;
					gapChunk = m_newlineSize + 1;
					readState = 14;
					break;
				case 14:
					while(i < bufReadLen){
						ch = buf[i];
						if((ch == '/')&&(i==lineStart)){
							readState = 15;	//end of this sequence
							break;
						}else if(m_pFilter->IsValid(ch)){
							if(gapChunk > 0){
								if((gapChunk > 1 && seqChunkCount > 0) ||
								  (gapChunk != 10 + m_newlineSize && seqChunkCount == 0)){
								  	if( !corruptWarning ){
										ErrorMsg("File is corrupt.  Proceed with caution.");
										corruptWarning = true;
									}
									curContig->SetRepeatSeqGap(false);
								}
								gapChunk = 0;
							}
							seqChunk++;
							seqLength++;
						}else{
							gapChunk++;
							if(seqChunk == 10){
								seqChunk = 0;
								seqChunkCount++;
								if(seqChunkCount == 6){
									//got a complete line.  start over
									seqChunkCount = 0;
								}
							}
							if(ch == '\n')
								lineStart = i + 1;
						}
						i++;
					}
					break;
				case 15:
					if((ch == '\n')&&(buf[lineStart+1] == '/')){
						curContig->SetSectEnd(gnContigSequence, lineStart - m_newlineSize + streamPos);
						curContig->SetSeqLength(seqLength);
						m_contigList.push_back(curContig);
						curContig = 0;
						curSpec->SetLength(seqLength);
						curSpec = 0;
						seqLength = 0;
						lineStart = i + 1;
						sectionStart = i + 1;
						readState = 0;
					}
					break;
			}
		}
		streamPos += bufReadLen;
	}
	if(curContig != 0){
		curContig->SetSectEnd(gnContigSequence, streamPos - 1);
		curContig->SetSeqLength(seqLength);
		m_contigList.push_back(curContig);
		curSpec->SetLength(seqLength);
	}
	if(curSpec != 0)
		if((curFrag->GetFeatureListLength() == 0) && (curFrag->GetHeaderListLength() == 0)
			&&(curSpec->GetLength() == 0)){
			m_spec->RemoveSpec(m_spec->GetSpecListLength() - 1);
			delete curFrag;
		}
	m_ifstream.clear();
	return true;
}

}	// end namespace genome

