/////////////////////////////////////////////////////////////////////////////
// File:            gnSEQSource.h
// Purpose:         Implements gnBaseSource for .SEQ files
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
#include "libGenome/gnSEQSource.h"
#include "libGenome/gnGBKSource.h"
#include "libGenome/gnSourceSpec.h"
#include "libGenome/gnSourceHeader.h"
#include "libGenome/gnSourceQualifier.h"
#include "libGenome/gnLocation.h"
#include "libGenome/gnStringTools.h"
#include "libGenome/gnDebug.h"
#include <cstring>


using namespace std;
namespace genome {


gnSEQSource::gnSEQSource()
{
	m_openString = "";
	m_pFilter = gnFilter::fullDNASeqFilter();
	if(m_pFilter == NULL){
		DebugMsg("Error using static sequence filters.\n");
	}
}
gnSEQSource::gnSEQSource( const gnSEQSource& s ) : gnFileSource(s)
{
	vector< gnFileContig* >::const_iterator iter = s.m_contigList.begin();
	for( ; iter != s.m_contigList.end(); ++iter )
	{
		m_contigList.push_back( (*iter)->Clone() );
	}
}
gnSEQSource::~gnSEQSource()
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
boolean gnSEQSource::HasContig( const string& name ) const
{
	for(uint32 i = 0 ; i <= m_contigList.size(); i++ )
	{
		if( name == m_contigList[i]->GetName() )
			return true;
	}
	return false;
}
uint32 gnSEQSource::GetContigID( const string& name ) const
{
	for(uint32 i = 0 ; i <= m_contigList.size(); i++ )
	{
		if( name == m_contigList[i]->GetName() )
			return i;
	}
	return ALL_CONTIGS;
}
string gnSEQSource::GetContigName( const uint32 i ) const
{
	if( i < m_contigList.size() )
	{
		return m_contigList[i]->GetName();
	}
	return "";
}
gnSeqI gnSEQSource::GetContigSeqLength( const uint32 i ) const
{
	if( i == ALL_CONTIGS)
		return m_spec->GetLength();
	if( i < m_contigList.size() )
	{
		return m_contigList[i]->GetSeqLength();
	}
	return GNSEQI_ERROR;
}

boolean gnSEQSource::SeqRead( const gnSeqI start, char* buf, gnSeqI& bufLen, const uint32 contigI ){
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
			gnSeqC* tmpBuf = new gnSeqC[readLen];	//read into tmpBuf, then filter tmpBuf into curBuf

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
			delete[] tmpBuf;
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
			gnSeqC* tmpBuf = new gnSeqC[readLen];	//read into tmpBuf, then filter tmpBuf into curBuf

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
			delete[] tmpBuf;
		}
		bufLen = curLen;
	}
	return true;

}
// private:
// figures out which contig the sequence starts at then calls SeqStartPos to get the offset within that contig
// returns startPos, the file offset where the sequence starts
// returns true if successful, false otherwise
boolean gnSEQSource::SeqSeek( const gnSeqI start, const uint32& contigI, uint64& startPos, uint64& readableBytes )
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
boolean gnSEQSource::SeqStartPos( const gnSeqI start, gnFileContig& contig, uint64& startPos, uint64& readableBytes )
{
	readableBytes = 0;
	uint32 curLen = 0;
	//seek to the file offset where the contig starts
	startPos = contig.GetSectStartEnd(gnContigSequence).first;	//set startPos to start where the contig starts
	m_ifstream.seekg( startPos, ios::beg );
	if( m_ifstream.eof() ){
		DebugMsg("ERROR in gnSEQSource::Incorrect contig start position, End of file reached!\n");
		return false;
	}
	while( true )
	{
		  // READ the rest of the contig skipping over invalid characters until we get to the starting base pair.
		  // startPos will contain the file offset with the starting base pair
		uint32 tmpbufsize = contig.GetSectStartEnd(gnContigSequence).second - startPos;
		if(tmpbufsize == 0){
			DebugMsg("ERROR in gnSEQSource: stored contig size is incorrect.");
			return false;
		}
		uint64 startOffset = start;
		if(contig.HasRepeatSeqGap()){
			if(contig.GetRepeatSeqGapSize().first > 0){
				if(contig.GetRepeatSeqGapSize().second > 0){
					startOffset += (start*contig.GetRepeatSeqGapSize().second)/contig.GetRepeatSeqGapSize().first;
					startPos+=startOffset;
					m_ifstream.seekg(startPos , ios::beg);
					readableBytes = contig.GetSectStartEnd(gnContigSequence).second - startPos;
					return true;
				}
			}else{
				startPos+=start;
				m_ifstream.seekg(startPos , ios::beg);
				readableBytes = contig.GetSectStartEnd(gnContigSequence).second - startPos;
				return true;
			}
		}
		tmpbufsize = tmpbufsize < BUFFER_SIZE ? tmpbufsize : BUFFER_SIZE;  //read in the smaller of the two.
		char *tmpbuf = new char[tmpbufsize];
		m_ifstream.read( tmpbuf, tmpbufsize );
		if( m_ifstream.eof() ){
			ErrorMsg("ERROR in gnSEQSource::Read End of file reached!\n");
			delete[] tmpbuf;
			return false;
		}
		for( uint32 i=0; i < tmpbufsize; ++i ){
			if( m_pFilter->IsValid(tmpbuf[i]) ){
				if( curLen >= start ){ //stop when we reach the starting offset within this contig
					startPos += i;
					m_ifstream.seekg( startPos, ios::beg );  //seek to startPos
					readableBytes = contig.GetSectStartEnd(gnContigSequence).second - startPos;
					delete[] tmpbuf;
					return true;
				}
				++curLen;  //each time we read a valid b.p., increment the sequence length
			}
		}
		startPos += tmpbufsize;
		delete[] tmpbuf;
	}
	return true;
}


//IMPLEMENT ME!  move these static methods somewhere else!  especially basecount!
void gnSEQSource::BaseCount(const string& bases, gnSeqI& a_count, gnSeqI& c_count, gnSeqI& g_count, gnSeqI& t_count, gnSeqI& other_count){
	a_count = 0;
	c_count = 0;
	g_count = 0;
	t_count = 0;
	other_count = 0;
	for(uint32 i=0; i < bases.length(); i++){
		if((bases[i] == 'a')||(bases[i] == 'A'))
			a_count++;
		else if((bases[i] == 'c')||(bases[i] == 'C'))
			c_count++;
		else if((bases[i] == 'g')||(bases[i] == 'G'))
			g_count++;
		else if((bases[i] == 't')||(bases[i] == 'T'))
			t_count++;
		else
			other_count++;
	}
}

void gnSEQSource::FormatString(string& data, uint32 offset, uint32 width){
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

boolean gnSEQSource::Write(gnGenomeSpec *spec, const string& filename){
	ErrorMsg("Writing DNAStar SEQ files is not supported at this time.  Try again next week.\n");
	return false;
}

gnFileContig* gnSEQSource::GetFileContig( const uint32 contigI ) const{
	if(m_contigList.size() > contigI)
		return m_contigList[contigI];
	return NULL;
}

//File parsing access routine
boolean gnSEQSource::ParseStream( istream& fin )
{
	// INIT temp varables
	uint32 readState = 0;
	uint32 lineStart = 0;
	int64 gapstart = -1;
	// INIT buffer
	uint32 sectionStart = 0;
	uint64 streamPos = 0;
	uint64 bufReadLen = 0;
	uint64 remainingBuffer = 0;
	char* buf = new char[BUFFER_SIZE];
	gnFragmentSpec* curFrag = 0;
	gnSourceSpec* curSpec = 0;
	gnSourceHeader *curHeader;
	gnBaseFeature* curFeature;
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
	gnSeqI lineSeqSize = 0;
	
	m_spec = new gnGenomeSpec();
	while( !fin.eof() )
	{
		if(sectionStart > 0){
			if(readState == 15)
				sectionStart = bufReadLen;
			else if(readState == 16)
				sectionStart = lineStart;
			remainingBuffer = bufReadLen - sectionStart;
			memmove(buf, buf+sectionStart, remainingBuffer);
		}
		  // read chars
		fin.read( buf + remainingBuffer, BUFFER_SIZE - remainingBuffer);
		streamPos -= remainingBuffer;
		lineStart -= sectionStart;
		if(gapstart > 0)
			gapstart -= sectionStart;
		sectionStart = 0;
		bufReadLen = fin.gcount();
		bufReadLen += remainingBuffer;
		
		for( uint32 i=remainingBuffer ; i < bufReadLen ; i++ )
		{
			char ch = buf[i];
			switch( readState )
			{
				case 0: 	//Assume we are in header at the start of a new line.  
							//Look for keywords starting in column 1
					if((ch == '\n')&&(buf[lineStart] != ' ')&&(buf[lineStart] != '	')){  //not equal to space or tab
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
							uint32 j = SEQ_LOCUS_NAME_LENGTH + 1;
							for(; j > 0; j--)	
								if((buf[lineStart+SEQ_LOCUS_NAME_COLUMN+j-1] != ' ')&&(buf[sectionStart+SEQ_LOCUS_NAME_COLUMN+j-1] != '	'))
									break;
							curContigName = string(buf+lineStart+SEQ_LOCUS_NAME_COLUMN-1, j+1);
							curFrag->SetName(curContigName);
						}else if(strncmp(&buf[lineStart], "^^", 2) == 0){
							//start the sequence.
							if(curContig == NULL){
								curContig = new gnFileContig();
								curContig->SetName(curContigName);
								curContigName = "";
							}
							i--;
							readState = 14;
							break;
						}
					}
					if(ch == '\n')
						lineStart = i + 1;
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
					}else if((ch == ' ')||(ch == '	')||(ch == ',')){
						//read end location
						string ender(buf+sectionStart, i - sectionStart);
						gnSeqI curLocationEnd = atoi(ender.c_str());
						gnLocation curLocation(curLocationStart, curStartLength, curLocationEnd, curEndLength, curBaseLocationType);
						curEndLength = 0;
						curStartLength = 0;
						curFeature->AddLocation(curLocation, curFeature->GetLocationListLength());
						readState = ch == ',' ? 3 : 7;  //read another loc if we need to.
						sectionStart = i+1;
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
					}
					break;
				case 9:		//are we getting a string? look for " or [
					if(ch == '"'){
						readState = 10;
						sectionStart = i;
						curQualifierStart = i + streamPos;
					}else if(ch == '['){
						readState = 11;
						sectionStart = i;
					}else if((ch == '\r')||(ch == '\n')){
						curFeature->AddQualifier(new gnSourceQualifier(this, curQualifierName, sectionStart + streamPos, i - sectionStart));
						sectionStart = i+1;
						readState = 7; //look for another qualifier
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
						curFeature->AddQualifier(new gnSourceQualifier(this, curQualifierName, curQualifierStart, i - sectionStart));
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
					if(ch == '^')	//stupid blattlab file format.
						readState = 14;
					else{
						curContig->SetSectStart(gnContigSequence, i + 1 + streamPos);
						readState = 16;
					}
					break;
				case 14:	//wait for newline before sequence starts.
					if(ch == '\n'){
						curContig->SetRepeatSeqGap(true);
						lineStart = i + 1;
						sectionStart = i + 1;
						curContig->SetSectStart(gnContigSequence, i + 1 + streamPos);
						readState = 15;
					}
					break;
				case 15:
					if(m_pFilter->IsValid(ch))
						seqLength++;
					else
						curContig->SetRepeatSeqGap(false);
					break;
				case 16:
					if((ch == '/')&&(i==lineStart)){
						readState = 17;
					}else if(m_pFilter->IsValid(ch)){
						seqLength++;
						lineSeqSize++;
						if(gapstart >= 0){
							curContig->SetRepeatGapSize(i - gapstart);
							gapstart = -1;
						}
					}else if(ch == '\n'){	//IMPLEMENT ME! Needs consistent gap size checking
						if(sectionStart == lineStart){
							curContig->SetRepeatSeqGap(true);
							curContig->SetRepeatSeqSize(seqLength);
							gapstart = i;
							for(; gapstart >= lineStart; gapstart--)
								if(m_pFilter->IsValid(buf[gapstart]))
									break;
							gapstart++;
						}else if(lineSeqSize != curContig->GetRepeatSeqGapSize().first)
							curContig->SetRepeatSeqGap(false);
						lineSeqSize = 0;
						lineStart = i + 1;
					}
					break;
				case 17:
					if((ch == '\n')&&(buf[lineStart+1] == '/')){
						curContig->SetSectEnd(gnContigSequence, lineStart - 2 + streamPos);
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
	if(curSpec != NULL)
		if((curFrag->GetFeatureListLength() == 0) && (curFrag->GetHeaderListLength() == 0)
			&&(curSpec->GetLength() == 0)){
			m_spec->RemoveSpec(m_spec->GetSpecListLength() - 1);
			delete curFrag;
		}
	m_ifstream.clear();
	delete[] buf;
	return true;
}

}	// end namespace genome

