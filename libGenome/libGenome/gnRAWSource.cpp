/////////////////////////////////////////////////////////////////////////////
// File:            gnRAWSource.h
// Purpose:         Implements gnBaseSource for raw data files
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
#include "libGenome/gnRAWSource.h"
#include "libGenome/gnGenomeSpec.h"
#include "libGenome/gnFragmentSpec.h"
#include "libGenome/gnSourceSpec.h"
#include "libGenome/gnStringTools.h"
#include "libGenome/gnDebug.h"


using namespace std;
namespace genome {


gnRAWSource::gnRAWSource()
{
	m_openString = "";
	m_contig = NULL;
	m_pFilter = NULL;
}

gnRAWSource::gnRAWSource( const gnRAWSource& s ) : gnFileSource(s)
{
	m_contig = NULL;
	if(s.m_contig != NULL)
		m_contig = s.m_contig->Clone();
}

gnRAWSource::~gnRAWSource()
{
	m_ifstream.close();
	delete m_contig;
}

boolean gnRAWSource::HasContig( const string& name ) const
{
	if( name.length() == 0 )
		return true;
	return false;
}

uint32 gnRAWSource::GetContigID( const string& name ) const
{
	return ALL_CONTIGS;
}

string gnRAWSource::GetContigName( const uint32 i ) const
{
	return "";
}

gnSeqI gnRAWSource::GetContigSeqLength( const uint32 i ) const
{
	if( m_contig && (i == 0 || i == ALL_CONTIGS))
		return m_contig->GetSeqLength();
	return GNSEQI_ERROR;
}

boolean gnRAWSource::SeqRead( const gnSeqI start, char* buf, gnSeqI& bufLen, const uint32 contigI ){
	return Read( start, buf, bufLen );
}

gnGenomeSpec *gnRAWSource::GetSpec() const{
	return m_spec->Clone();
}

boolean gnRAWSource::Write(gnSequence& seq, const string& filename){
	ofstream m_ofstream(filename.c_str(), ios::out | ios::binary);
	if(!m_ofstream.is_open())
		return false;

	gnSeqC buf[BUFFER_SIZE + 1];
	buf[BUFFER_SIZE] = 0;
	gnSeqI readOffset = 0;
	gnSeqI readLength = seq.length();
	while(readLength > 0){	//buffer the read/writes
		gnSeqI writeLen = readLength < BUFFER_SIZE ? readLength : BUFFER_SIZE;
		if(!seq.ToArray(buf, writeLen, readOffset + 1))
			return false;
		m_ofstream.write( buf, writeLen );
		readLength -= writeLen;
		readOffset += writeLen;
	}
	m_ofstream.flush();
	m_ofstream.close();
	return true;
}

gnFileContig* gnRAWSource::GetFileContig( const uint32 contigI ) const{
	if(contigI > 0)
		return NULL;
	return m_contig;
}

//File parsing access routine
boolean gnRAWSource::ParseStream( istream& fin )
{
	// init variables
	uint64 streamPos = 0;
	uint64 bufReadLen = 0;
	Array<char> array_buf( BUFFER_SIZE );
	char* buf = array_buf.data;
	gnSeqI seqLength = 0;
	
	if( m_contig == NULL )
		m_contig = new gnFileContig();
	m_contig->SetName( "Raw Data" );
	m_contig->SetRepeatSeqGap(true);
	m_contig->SetSectStart(gnContigSequence, 0);
	
	uint64 offset = 0;
	if( !CheckRawData() ){
		fin.seekg( 0, ios::end );
		offset = fin.tellg();
	}
	else
	{
		while( !fin.eof() )
		{
			  // read chars
			fin.read( buf , BUFFER_SIZE );
			bufReadLen = fin.gcount();
			
			for( uint32 i=0 ; i < bufReadLen ; i++ )
			{				
				if(m_pFilter == NULL || m_pFilter->IsValid(buf[i]))
					seqLength++;
				else{
					m_contig->SetRepeatSeqGap(false);
				}
			}
			streamPos += bufReadLen;
		}
	}
	m_contig->SetSectEnd(gnContigSequence, streamPos + offset);
	m_contig->SetSeqLength(seqLength + offset );
	m_spec = new gnGenomeSpec();
	gnFragmentSpec* fragspec = new gnFragmentSpec();
	gnSourceSpec* sspec = new gnSourceSpec(this);
	sspec->SetSourceName(m_openString);
	m_spec->AddSpec(fragspec);
	fragspec->AddSpec(sspec);

	m_ifstream.clear();
	return true;
}

}	// end namespace genome

