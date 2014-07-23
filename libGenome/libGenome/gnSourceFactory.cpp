/////////////////////////////////////////////////////////////////////////////
// File:            gnSourceFactory.cpp
// Purpose:         Source Factory for all Sources
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

#include "libGenome/gnSourceFactory.h"
#include "libGenome/gnBaseSource.h"
#include "libGenome/gnStringTools.h"
#include "libGenome/gnFASSource.h"
#include "libGenome/gnGBKSource.h"
#include "libGenome/gnSEQSource.h"
#include "libGenome/gnRAWSource.h"
#include "libGenome/gnException.h"

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif


using namespace std;
namespace genome {


gnSourceFactory::gnSourceFactory()
{
	m_sourceClassList.insert( map< string, gnBaseSource* >::value_type(".fas", new gnFASSource()));
	m_sourceClassList.insert( map< string, gnBaseSource* >::value_type(".FAS", new gnFASSource()));
	m_sourceClassList.insert( map< string, gnBaseSource* >::value_type(".seq", new gnSEQSource()));
	m_sourceClassList.insert( map< string, gnBaseSource* >::value_type(".SEQ", new gnSEQSource()));
	m_sourceClassList.insert( map< string, gnBaseSource* >::value_type(".gbk", new gnGBKSource()));
	m_sourceClassList.insert( map< string, gnBaseSource* >::value_type(".GBK", new gnGBKSource()));
	m_sourceClassList.insert( map< string, gnBaseSource* >::value_type(".gb", new gnGBKSource()));
	m_sourceClassList.insert( map< string, gnBaseSource* >::value_type(".GB", new gnGBKSource()));
	m_sourceClassList.insert( map< string, gnBaseSource* >::value_type(".raw", new gnRAWSource()));
	m_sourceClassList.insert( map< string, gnBaseSource* >::value_type(".RAW", new gnRAWSource()));
	m_pDefaultSourceClass = new gnFASSource();
}

gnSourceFactory::~gnSourceFactory()
{
	vector< gnBaseSource* >::iterator iter = m_sourceList.begin();
	for(; iter != m_sourceList.end(); ++iter )
	{
		delete *iter;
	}

	map< string, gnBaseSource* >::iterator cIter = m_sourceClassList.begin();
	for(; cIter != m_sourceClassList.end(); ++cIter )
	{
		delete cIter->second;
	}
}

gnSourceFactory* gnSourceFactory::GetSourceFactory()
{
	//use construct on first use method to avoid the static constructor
	//initialization fiasco...
	static gnSourceFactory* m_sSourceFactory = new gnSourceFactory();
	return m_sSourceFactory;
}


boolean gnSourceFactory::DelSourceClass( const string& ext ){
	map< string, gnBaseSource* >::iterator iter = m_sourceClassList.find( ext );
	if( iter != m_sourceClassList.end() ){
		m_sourceClassList.erase( iter );
		return true;
	}
	return false;
}
gnBaseSource* gnSourceFactory::GetSourceClass( const string& ext ) const{
	map< string, gnBaseSource* >::const_iterator iter = m_sourceClassList.find( ext );
	if( iter != m_sourceClassList.end() ){
		return iter->second;
	}
	return m_pDefaultSourceClass;
}
gnBaseSource* gnSourceFactory::MatchSourceClass( const string& sourceStr ) const{
	string::size_type dot_loc = sourceStr.rfind('.');
	if(dot_loc != string::npos){
		string ext = sourceStr.substr(dot_loc, sourceStr.length() - dot_loc);
		return GetSourceClass(ext);
	}
	return m_pDefaultSourceClass;
}
boolean gnSourceFactory::HasSourceClass( const string& ext ) const{
	if( m_sourceClassList.find(ext) != m_sourceClassList.end() ){
		return true;
	}
	return false;
}
boolean gnSourceFactory::SetSourceClass( const string& ext, const gnBaseSource& source ){
	map< string, gnBaseSource* >::iterator iter = m_sourceClassList.find( ext );
	if( iter == m_sourceClassList.end() ){
		m_sourceClassList.insert( 
			map< string, gnBaseSource* >::value_type( ext, source.Clone() ) );
	}else{
		iter->second = source.Clone();
	}
	return true;
}
boolean gnSourceFactory::AddPath( const string& path ){
	if( PathExists( path ) && !HasPath(path) )
	{
		m_pathList.push_back( path );	
		return true;
	}
	return false;
}
boolean gnSourceFactory::DelPath( uint32 i ){
	if( i < m_pathList.size() ){
		vector<string>::iterator iter = m_pathList.begin() + i;
		m_pathList.erase( iter );
		return true;
	}
	return false;
}
boolean gnSourceFactory::InsPath( const string& path, uint32 i ){
	if( (i < m_pathList.size()) && PathExists( path ) )
	{
		vector<string>::iterator iter = m_pathList.begin() + i;
		m_pathList.insert( iter, path );
		return true;
	}
	return false;
}
string gnSourceFactory::GetPath( uint32 i ) const{
	if( i < m_pathList.size() ){
		return m_pathList[i];
	}
	return "";
}
boolean gnSourceFactory::HasPath( string path ) const{
	standardizePathString( path );
	for( uint32 i = 0; i < m_pathList.size(); ++i ){
		if( m_pathList[i] == path )
			return true;
	}
	return false;
}

boolean gnSourceFactory::GetURL( const string& urlStr, string& localFile ){
	return false;
}
  // Sources
gnBaseSource* gnSourceFactory::AddSource( const string& sourceStr, boolean searchPaths )
{
	string openStr = sourceStr;
	// Check if Exists
	gnBaseSource* source = HasSource( sourceStr, false );
	if( source != 0 )
		return source;
	// Check if File Present and valid
	gnBaseSource* newSource = MatchSourceClass( sourceStr )->Clone();
	if(newSource == NULL)
		return NULL;

	// Check the URL to see if we need to cache the file locally
	if(sourceStr.substr(0, 7) == "http://"){
		ErrorMsg("Sorry, no HTTP support.\n");
		return NULL;
	}else if(sourceStr.substr(0, 6) == "ftp://"){
		ErrorMsg("Sorry, no FTP support.\n");
		return NULL;
	}else if(sourceStr.substr(0, 8) == "file:///")
		openStr = sourceStr.substr(8);

	// Now try to open the local file
	try{
		newSource->Open( openStr );
		m_sourceList.push_back( newSource );
		return newSource;
	}catch(gnException& e){
		if(e.GetCode() != FileNotOpened()){
			delete newSource;
			e.AddCaller(__PRETTY_FUNCTION__);
			throw e;
		}
	}
	
	if( searchPaths )
	{
		// Check other paths for Exists or Presents/valid
		string file = getFileString( openStr );
		vector< string >::iterator iter = m_pathList.begin();
		for( ; iter != m_pathList.end(); ++iter )
		{
			string newSourceStr = *iter + file;
			// Check if Exists
			source = HasSource( newSourceStr, false );
			if( source != 0 )
			{
				delete newSource;
				return source;
			}
			// Check if File Present and valid
			try{
				newSource->Open( newSourceStr );
				m_sourceList.push_back( newSource );
				return newSource;
			}catch(gnException& e){
				if(e.GetCode() != FileNotOpened()){
					delete newSource;
					e.AddCaller(__PRETTY_FUNCTION__);
					throw e;
				}
			}
		}
	}
	// Cannot find a valid source
	delete newSource;
	Throw_gnEx(FileNotOpened());
}
gnBaseSource* gnSourceFactory::GetSource( uint32 i ) const{
	if( i < m_sourceList.size() ){
		return *(m_sourceList.begin() + i);
	}
	return 0;
}
void gnSourceFactory::DelSource( uint32 i ){
	if( i >= m_sourceList.size() )
		Throw_gnEx(IndexOutOfBounds());
	vector< gnBaseSource* >::iterator iter = m_sourceList.begin() + i;
	gnBaseSource* source = *iter;
	try{
		source->Close();
		m_sourceList.erase(iter);
		delete source;
	}catch(gnException& e){
		e.AddCaller(__PRETTY_FUNCTION__);
		throw e;
	}
}

boolean gnSourceFactory::DelSource( const gnBaseSource* source ){
	vector< gnBaseSource* >::iterator iter = m_sourceList.begin();
	for( ; iter != m_sourceList.end(); ++iter ){
		if( *iter == source ){
			gnBaseSource* source = *iter;
			try{
				source->Close();
				m_sourceList.erase(iter);
				delete source;
				return true;
			}catch(gnException& e){
				e.AddCaller(__PRETTY_FUNCTION__);
				throw e;
			}
		}
	}
	return false;
}
gnBaseSource* gnSourceFactory::HasSource( string sourceStr, boolean searchPaths ) const{
	standardizePathString( sourceStr );
	vector< gnBaseSource* >::const_iterator iter1 = m_sourceList.begin();
	for( ; iter1 != m_sourceList.end(); ++iter1 )
	{
		if( (*iter1)->GetOpenString() == sourceStr )
			return *iter1;
	}
	
	if( searchPaths )
	{
		string file = getFileString( sourceStr );
		vector< string >::const_iterator iter2 = m_pathList.begin();
		for( ; iter2 != m_pathList.end(); ++iter2 )
		{
			iter1 = m_sourceList.begin();
			for( ; iter1 != m_sourceList.end(); ++iter1 )
			{
				if( (*iter1)->GetOpenString() == (*iter2 + file) )
					return *iter1;
			}
		}
	}
	return 0;
}
	
// private:
boolean gnSourceFactory::PathExists( string path ) const{
#ifdef HAVE_UNISTD_H
	standardizePathString( path );
	char folder[FILENAME_MAX], *f2;
	f2 = getcwd( folder, FILENAME_MAX );
	
	if( chdir( path.c_str() ) ){
		return false;
	}
	int err = chdir( folder );
#endif
	return true;
}


}	// end namespace genome

