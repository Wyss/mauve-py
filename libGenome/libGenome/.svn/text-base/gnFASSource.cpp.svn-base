/////////////////////////////////////////////////////////////////////////////
// File:            gnFASSource.h
// Purpose:         Implements gnBaseSource for .FAS files
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

#include "libGenome/gnFASSource.h"
#include "libGenome/gnBaseSpec.h"
#include "libGenome/gnFilter.h"
#include "libGenome/gnStringTools.h"
#include "libGenome/gnSourceSpec.h"
#include "libGenome/gnSourceHeader.h"
#include "libGenome/gnDebug.h"


using namespace std;
namespace genome {


gnFASSource::gnFASSource()
{
	m_openString = "";
	m_pFilter = gnFilter::fullDNASeqFilter();
	if(m_pFilter == NULL){
		DebugMsg("Error using static sequence filters.");
	}
}
gnFASSource::gnFASSource( const gnFASSource& s ) : gnFileSource(s)
{
	vector< gnFileContig* >::const_iterator iter = s.m_contigList.begin();
	for( ; iter != s.m_contigList.end(); ++iter )
	{
		m_contigList.push_back( (*iter)->Clone() );
	}
}
gnFASSource::~gnFASSource()
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
boolean gnFASSource::HasContig( const string& nameStr ) const
{
	for(uint32 i = 0 ; i <= m_contigList.size(); i++ )
	{
		if( nameStr == m_contigList[i]->GetName() )
			return true;
	}
	return false;
}
uint32 gnFASSource::GetContigID( const string& name ) const
{
	for(uint32 i = 0 ; i <= m_contigList.size(); i++ )
	{
		if( name == m_contigList[i]->GetName() )
			return i;
	}
	return ALL_CONTIGS;
}
string gnFASSource::GetContigName( const uint32 i ) const
{
	if( i < m_contigList.size() ){
		return m_contigList[i]->GetName();
	}
	return "";
}
gnSeqI gnFASSource::GetContigSeqLength( const uint32 i ) const
{
	if( i < m_contigList.size() ){
		return m_contigList[i]->GetSeqLength();
	}else if( i == ALL_CONTIGS){
		gnSeqI seqlen = 0;
		for(uint32 j=0; j < m_contigList.size(); j++)
			seqlen += m_contigList[j]->GetSeqLength();
		return seqlen;
	}
	return GNSEQI_ERROR;
}
gnFileContig* gnFASSource::GetContig( const uint32 i ) const
{
	if( i < m_contigList.size() ){
		return m_contigList[i];
	}
	return 0;
}

boolean gnFASSource::SeqRead( const gnSeqI start, char* buf, gnSeqI& bufLen, const uint32 contigI ) 
{
	boolean result = true;
#pragma omp critical
{
	result = SeqReadImpl( start, buf, bufLen, contigI );
}
	return result;
}

boolean gnFASSource::SeqReadImpl( const gnSeqI start, char* buf, gnSeqI& bufLen, const uint32 contigI ) 
{
	m_ifstream.clear();
	uint32 contigIndex = contigI;
	uint64 startPos = 0;
	uint64 readableBytes = 0;
	if( !SeqSeek( start, contigIndex, startPos, readableBytes ) ){
		bufLen = 0;
		return false;
	}
	
	if( contigI == ALL_CONTIGS ){
		uint32 curLen = 0;
		uint64 bytesRead = 0;
		while (curLen < bufLen){
//SeqSeek to start, Figure out how much can be read before SeqSeeking again.
			if(readableBytes <= 0)	//Look out for zero length contigs!  IMPLEMENT ME
				if( !SeqSeek( start + curLen, contigIndex, startPos, readableBytes ) ){
					bufLen = curLen;
					return true;
				}
			//readLen is the amount to read on this pass
			uint64 readLen = (bufLen - curLen) < readableBytes ? (bufLen - curLen) : readableBytes;
			Array<gnSeqC> array_buf( readLen );
			gnSeqC* tmpBuf = array_buf.data;	//read into tmpBuf, then filter tmpBuf into curBuf

			// read chars and filter
			m_ifstream.read(tmpBuf, readLen);
			uint64 gc = m_ifstream.gcount();
			bytesRead += gc;
			readableBytes -= gc;
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
	else if( contigI < m_contigList.size() )
	{
		uint32 curLen = 0;
		//check to see if the buffer is bigger than the contig.  if so truncate it.
		gnSeqI contigSize = m_contigList[contigI]->GetSeqLength();
		bufLen = bufLen < contigSize ? bufLen : contigSize;
		while (curLen < bufLen)
		{
			uint64 readLen = bufLen - curLen;	//the amount to read on this pass
			Array<gnSeqC> array_buf( readLen );	// use Array template to ensure proper deallocation
			gnSeqC* tmpBuf = array_buf.data;	//read into tmpBuf, then filter tmpBuf into curBuf

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
boolean gnFASSource::SeqSeek( const gnSeqI start, const uint32 contigI, uint64& startPos, uint64& readableBytes )
{
	if( contigI == ALL_CONTIGS ){
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
boolean gnFASSource::SeqStartPos( const gnSeqI start, gnFileContig& contig, uint64& startPos, uint64& readableBytes ){
	readableBytes = 0;
	uint32 curLen = 0;
	//seek to the file offset where the contig starts
	startPos = contig.GetSectStartEnd(gnContigSequence).first;	//set startPos to start where the contig starts

	if(contig.HasRepeatSeqGap())
		if(contig.GetRepeatSeqGapSize().first > 0)
			if(contig.GetRepeatSeqGapSize().second > 0){
//check this	
				startPos += start + (start/contig.GetRepeatSeqGapSize().first)*contig.GetRepeatSeqGapSize().second;
				readableBytes = contig.GetSectStartEnd(gnContigSequence).second - startPos;
				m_ifstream.seekg( startPos, ios::beg );
				return true;
			}
	m_ifstream.seekg( startPos, ios::beg );
	if( m_ifstream.eof() ){
		ErrorMsg("ERROR in gnFASSource::Incorrect contig start position, End of file reached!\n");
		return false;
	}
	while( true ){
		  // READ the rest of the contig skipping over invalid characters until we get to the starting base pair.
		  // startPos will contain the file offset with the starting base pair
		uint32 tmpbufsize = contig.GetSectStartEnd(gnContigSequence).second - startPos;
		if(tmpbufsize == 0){
			ErrorMsg("ERROR in gnFASSource: stored contig size is incorrect.\n");
			return false;
		}
		tmpbufsize = tmpbufsize < BUFFER_SIZE ? tmpbufsize : BUFFER_SIZE;  //read in the smaller of the two.
		Array<char> array_buf( tmpbufsize );	// use Array template to ensure proper deallocation
		char *tmpbuf = array_buf.data;
		m_ifstream.read( tmpbuf, tmpbufsize );
		if( m_ifstream.eof() ){
			ErrorMsg("ERROR in gnFASSource::Read End of file reached!\n");
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

//write out the given source in this file format.
boolean gnFASSource::Write(gnBaseSource *source, const string& filename){
	ofstream m_ofstream(filename.c_str(), ios::out | ios::binary);
	if(!m_ofstream.is_open())
		return false;
	uint32 contigCount = source->GetContigListLength();
	for(uint32 contigI = 0; contigI < contigCount; contigI++){
		string contigName = source->GetContigName(contigI);
		m_ofstream << ">" << contigName << ";\n";
		gnSeqI seqLength = source->GetContigSeqLength(contigI);
		while(seqLength > 0){
			gnSeqI writeLen = seqLength < BUFFER_SIZE ? seqLength : BUFFER_SIZE;		
			Array<gnSeqC> array_buf( writeLen+1 );	// use Array template to ensure proper deallocation
			gnSeqC *bases = array_buf.data;
			boolean success = source->SeqRead(0, bases, writeLen, contigI);
			if(!success)
				return false;
			bases[writeLen] = 0;
			m_ofstream << bases << "\n";
			seqLength -= writeLen;
		}
	}
	m_ofstream.close();
	return true;
}

//write out the given sequence in this file format.
void gnFASSource::Write(gnSequence& seq, const string& filename, boolean write_coords, boolean enforce_unique_names){
	ofstream m_ofstream(filename.c_str(), ios::out | ios::binary);
	if(!m_ofstream.is_open())
		Throw_gnEx(FileNotOpened());
	Write(seq, m_ofstream, write_coords, enforce_unique_names);
	m_ofstream.close();
}

void gnFASSource::Write(gnSequence& seq, ostream& m_ostream, boolean write_coords, boolean enforce_unique_names){
	vector<string> contigNameList;
	Array<gnSeqC> array_buf( BUFFER_SIZE );	// use Array template to ensure proper deallocation
	gnSeqC *bases = array_buf.data;
	gnGenomeSpec* spec = seq.GetSpec();
	
	//set the newline type.
	string newline = "\r\n";
	gnSeqI readOffset = 1;

	for(uint32 fragI = 0; fragI < seq.contigListLength(); fragI++){
		//write out contig name and header
		string contigName = seq.contigName(fragI);

		if(enforce_unique_names){
			uint32 name_count = 0;
			for(uint32 i=0; i < contigNameList.size(); i++)
				if(contigNameList[i] == contigName)
					name_count++;
			contigNameList.push_back(contigName);
			if(name_count > 0)
				contigName += "_" + uintToString(name_count);
		}

		gnFragmentSpec* subSpec = spec->GetSpec(fragI);
		gnBaseHeader *gpbh = subSpec->GetHeader(0);
		string header = "";
		if(gpbh != NULL){
			header = gpbh->GetHeader();
			//delete everything after the first newline.
			string::size_type newlinepos = header.find_first_of('\n', 0);
			if(newlinepos != string::npos){
				if(header[newlinepos-1] == '\r')
					newlinepos--;
				header = header.substr(0, newlinepos);
			}
		}	//IMPLEMENT ME!! Does header line need to be < 80 chars?

		//write out the sequence
		gnSeqI readLength = seq.contigLength(fragI);
		m_ostream << ">" << contigName;
		if(write_coords)
			m_ostream << " " << "(" << readOffset << ", " << readOffset + readLength - 1 << ")";
		m_ostream << "; " << header << endl;

		gnSeqI linePos = 0;
		while(readLength > 0){	//buffer the read/writes
			gnSeqI read_chunk_size = (BUFFER_SIZE / FAS_LINE_WIDTH) * FAS_LINE_WIDTH;
			gnSeqI writeLen = readLength < read_chunk_size ? readLength : read_chunk_size;
			
			if(!seq.ToArray(bases, writeLen, readOffset))
				return;
			for(gnSeqI curbase = 0; curbase < writeLen; curbase += FAS_LINE_WIDTH){
				gnSeqI writeout_size = writeLen - curbase < FAS_LINE_WIDTH ? writeLen - curbase : FAS_LINE_WIDTH;
				m_ostream << string(bases + curbase, writeout_size) << endl;
			}
			readLength -= writeLen;
			readOffset += writeLen;
		}
		if(linePos != 0)
			m_ostream << endl;
	}

}

gnGenomeSpec *gnFASSource::GetSpec() const{
	gnGenomeSpec *spec = new gnGenomeSpec();
	for(uint32 i=0; i < m_contigList.size(); i++){
		//create specs for the fragment and contig levels
		gnFragmentSpec *fragmentSpec = new gnFragmentSpec();
		gnSourceSpec *contigSpec = new gnSourceSpec((gnBaseSource*)this, i);
		//set up the spec tree structure
		spec->AddSpec(fragmentSpec, i);
		fragmentSpec->AddSpec(contigSpec);

		fragmentSpec->SetName(m_contigList[i]->GetName());
		fragmentSpec->SetSourceName(m_openString);
		contigSpec->SetName(m_contigList[i]->GetName());
		contigSpec->SetSourceName(m_openString);
		
		pair<uint64,uint64> headsect = m_contigList[i]->GetSectStartEnd(gnContigHeader);
		if(headsect.first != headsect.second){
			gnSourceHeader *gpsh = new gnSourceHeader((gnBaseSource *)this, string(""), headsect.first, headsect.second - headsect.first);
			fragmentSpec->AddHeader(gpsh, 0);
		}
	}
	return spec;
}

gnFileContig* gnFASSource::GetFileContig( const uint32 contigI ) const{
	if(m_contigList.size() > contigI)
		return m_contigList[contigI];
	return NULL;
}


boolean gnFASSource::ParseStream( istream& fin )
{
	// INIT temp varables
	uint32 readState = 0;
	gnFileContig* currentContig = 0;
	string nameFStr;
	uint64 seqLength = 0, gapLength = 0;
	// INIT buffer
	uint64 streamPos = 0;
	uint64 bufReadLen = 0;
	Array<gnSeqC> array_buf( BUFFER_SIZE );	// use Array template to ensure proper deallocation
	char* buf = array_buf.data;
	boolean paren_hit = false;
	uint32 repeatSeqSize = 0;
	boolean corrupt_msg = false;

	//decide what type of newlines we have
	DetermineNewlineType();

	
	while( !fin.eof() )
	{
		  // read chars
		fin.read( buf, BUFFER_SIZE);
		streamPos += bufReadLen;
		bufReadLen = fin.gcount();
		for( uint32 i=0 ; i < bufReadLen ; i++ )
		{
			char ch = buf[i];
			switch( readState )
			{
				case 0: // CHECK file validity
					if( !((buf[0] == '>') || !m_pFilter->IsValid(buf[0])) )
					{
						return false; // NOT a .fas file
					}
					readState = 1;
				case 1: // BRANCH
					if( ch == '>' )
					{
						seqLength = 0; gapLength = 0; // reset count
						currentContig = new gnFileContig();
						currentContig->SetFileStart( streamPos + i );
						currentContig->SetRepeatSeqGap(true);
						currentContig->SetRepeatSeqSize( repeatSeqSize );
						currentContig->SetRepeatGapSize( m_newlineSize );
						readState = 2;
						paren_hit = false;
						nameFStr = "";
					}
					else
						++gapLength;
					break;
				case 2: // > CONTIG NAME
					if( isNewLine(ch) || ch == ';')
					{
						currentContig->SetName( nameFStr );
						currentContig->SetSectStart( gnContigHeader, streamPos + i + 1 );
						if( ch == ';' )
							readState = 3;
						else{
							readState = 4;
							if(ch == '\r')
								currentContig->SetSectStart( gnContigHeader, streamPos + i + 2 );
						}
					}
					else if( ch == '(' ){
						if(isSpace(buf[i-1])){
							//delete the last char in nameFStr
							nameFStr = nameFStr.substr(0, nameFStr.length()-1);
						}
						paren_hit = true;
					}else if((!isSpace(ch) || nameFStr.length()) && !paren_hit )
						nameFStr += ch;
					break;
				case 3: // >... ; REMARK
					if( isNewLine(ch) )
						readState = 4;
					break;
				case 4: // >... ; /n NEWLINE BRANCH
					if( ch == '>' )
					{
						readState = 3;
					}
					else if( m_pFilter->IsValid(ch) )
					{
						currentContig->SetSectEnd( gnContigHeader, streamPos + i );
						currentContig->SetSectStart( gnContigSequence, streamPos + i);
						seqLength = 1; gapLength = 0;
						readState = 5;
					}
					break;
				case 5: // SEQUENCE
					while(i < bufReadLen){
						ch = buf[i];
						if( m_pFilter->IsValid(ch) ){
							if(gapLength > 0){
								if(seqLength != repeatSeqSize){
									if( !corrupt_msg ){
										corrupt_msg = true;
										ErrorMsg( "Sequence file appears corrupt, proceeding with caution\n" );
									}
									currentContig->SetRepeatSeqGap(false);
								}if(gapLength != m_newlineSize){
									//file is corrupt, can't do jumps.
									if( !corrupt_msg ){
										corrupt_msg = true;
										ErrorMsg( "Sequence file appears corrupt, proceeding with caution\n" );
									}
									currentContig->SetRepeatSeqGap(false);
								}
								currentContig->AddToSeqLength( seqLength );
								seqLength = 0;
								gapLength = 0;
							}
							seqLength++;
						}else if( ch == '>'){
							currentContig->AddToSeqLength( seqLength );
							currentContig->SetSectEnd( gnContigSequence, streamPos + i - 1);
							currentContig->SetFileEnd( streamPos + i - 1 );
							m_contigList.push_back(currentContig);
							readState = 1;
							i--;
							break;
						}
						else if( isNewLine(ch) ){
							if( repeatSeqSize == 0 ){
								repeatSeqSize = seqLength;
								currentContig->SetRepeatSeqSize( repeatSeqSize );
							}
							gapLength++;
						}
						else{
							// Error cannot do nice jumps
							currentContig->SetRepeatSeqGap(false);
							if( !corrupt_msg ){
								corrupt_msg = true;
								ErrorMsg( "Sequence file appears corrupt, proceeding with caution\n" );
							}
						}
						i++;
					}
					break;
				default:
					ErrorMsg("ERROR");
					return false;
					break;
			}
		}// for all buf
	}// while !eof
	// CLEAN UP
	if( currentContig != 0 )
	{
		if( readState == 2 )
		{
			currentContig->SetName( nameFStr );
		}
		if( (readState >= 2) && (readState < 5) )
		{
			currentContig->SetSectEnd( gnContigHeader, streamPos + bufReadLen );
		}
		else if( readState ==  5 )
		{
			currentContig->AddToSeqLength( seqLength );
			currentContig->SetSectEnd( gnContigSequence, streamPos + bufReadLen );
		}			
		currentContig->SetFileEnd( streamPos + bufReadLen );
		m_contigList.push_back(currentContig);
	}
	m_ifstream.clear();
	return true;
}


}	// end namespace genome

