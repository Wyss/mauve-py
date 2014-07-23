/////////////////////////////////////////////////////////////////////////////
// File:            gnSourceQualifier.cpp
// Purpose:         Source Qualifier class
// Description:     Provides an interface for Qualifier on disk.
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
#include "libGenome/gnSourceQualifier.h"


using namespace std;
namespace genome {


gnSourceQualifier::gnSourceQualifier(){
	m_source = NULL;
	m_name = "";
	m_start = 0;
	m_length = 0;
}
gnSourceQualifier::gnSourceQualifier( gnBaseSource* source, string& name, uint32 begin, uint32 length ){
	m_source = source;
	m_name = name;
	m_start = begin;
	m_length = length;
}
gnSourceQualifier::gnSourceQualifier(const gnSourceQualifier& s){
	m_source = s.m_source;
	m_start = s.m_start;
	m_length = s.m_length;
	m_name = string(s.m_name);
}
gnSourceQualifier::~gnSourceQualifier(){
};
string gnSourceQualifier::GetValue() const{
	Array<char> array_buf( m_length );
	char* buf = array_buf.data;
	gnSeqI readBytes = m_length;
	m_source->Read(m_start, buf, readBytes);
	string rval(buf, readBytes);
	return rval;
}


}	// end namespace genome

