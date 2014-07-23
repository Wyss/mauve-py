/////////////////////////////////////////////////////////////////////////////
// File:            gnFileSource.h
// Purpose:         Implements generic gnBaseSource methods
// Description:     Provides a general implementation for accessing DNA
//                  sequence contig data files.
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

#include "libGenome/gnFileSource.h"
#include <fstream>
#include "libGenome/gnFilter.h"


using namespace std;
namespace genome {

gnFileSource::gnFileSource() :
m_pFilter(gnFilter::fullDNASeqFilter())
{}

//copy constructor
gnFileSource::gnFileSource(const gnFileSource& gnfs){
	m_openString = gnfs.m_openString;
	m_pFilter = gnfs.m_pFilter;
	m_newlineType = gnfs.m_newlineType;
	m_newlineSize = gnfs.m_newlineSize;
#pragma omp critical
{
	m_ifstream.open( m_openString.c_str(), ios::in | ios::binary );
	if( !m_ifstream.is_open() )
		m_ifstream.clear();
}
}

// Open, Close	
void gnFileSource::Open( string openString )
{
	boolean opened = true;
#pragma omp critical
{
	m_ifstream.open(openString.c_str(), ios::in | ios::binary );
	if( m_ifstream.is_open() )
	{
		m_openString = openString;
		if( ParseStream(m_ifstream) )
		{
			;
		}
		else{
			m_ifstream.clear();
			m_ifstream.close();
		}
	}else{
		m_ifstream.clear();
		opened = false;
	}
}
	if(!opened)
		Throw_gnEx(FileNotOpened());
}
void gnFileSource::Open( )
{
	bool opened = true;
#pragma omp critical
{
	m_ifstream.open( m_openString.c_str(), ios::in | ios::binary );
	if( !m_ifstream.is_open() ){
		opened = false;
		m_ifstream.clear();
	}
}
	if(!opened)
		Throw_gnEx(FileNotOpened());
}
void gnFileSource::Close()
{
#pragma omp critical
	m_ifstream.close();
}

boolean gnFileSource::Read( const uint64 pos, char* buf, gnSeqI& bufLen) 
{
	bool failure = false;
#pragma omp critical
{
	m_ifstream.seekg(pos, ios::beg);
	m_ifstream.read(buf, bufLen);
	failure = m_ifstream.fail();
}
	if(failure){
#pragma omp critical
		m_ifstream.clear();
		return false;
	}
	return true;
}

void gnFileSource::DetermineNewlineType()
{
	// set default values
	m_newlineType = gnNewlineUnix;
	m_newlineSize = 1;

	//decide what type of newlines we have
	char buf[ BUFFER_SIZE ];
	m_ifstream.getline( buf, BUFFER_SIZE);
	m_ifstream.seekg(-2, ios::cur);
	m_ifstream.read( buf, 2);
	m_ifstream.seekg(0);
	if(buf[1] == '\n'){
		if(buf[0] == '\r'){
			m_newlineType = gnNewlineWindows;
			m_newlineSize = 2;
		}else{ 
			if(buf[1] == '\r')
				m_newlineType = gnNewlineMac;
			else
				m_newlineType = gnNewlineUnix;
		}
	}
}

}	// end namespace genome

