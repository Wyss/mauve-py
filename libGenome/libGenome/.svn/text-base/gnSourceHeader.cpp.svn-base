/////////////////////////////////////////////////////////////////////////////
// File:            gnSourceHeader.cpp
// Purpose:         Source Header class
// Description:     Provides an interface for Headers on disk.
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


#include <string>
#include "libGenome/gnBaseSource.h"
#include "libGenome/gnSourceHeader.h"


using namespace std;
namespace genome {


gnSourceHeader::gnSourceHeader(){
	m_source = NULL;
	m_start = 0;
	m_length = 0;
}
gnSourceHeader::gnSourceHeader( gnBaseSource* source, const string& name, const uint32 begin, const uint32 length ){
	m_source = source;
	m_name = string(name);
	m_start = begin;
	m_length = length;
}
gnSourceHeader::gnSourceHeader(const gnSourceHeader& s){
	m_source = s.m_source;
	m_start = s.m_start;
	m_length = s.m_length;
	m_name = string(s.m_name);
}
gnSourceHeader::~gnSourceHeader(){
};
string gnSourceHeader::GetHeader() const{
	Array<char> array_buf( m_length );
	char* buf = array_buf.data;
	gnSeqI readBytes = m_length;
	m_source->Read(m_start, buf, readBytes);
	string rval(buf, readBytes);
	return rval;
}

}	// end namespace genome

