/////////////////////////////////////////////////////////////////////////////
// DataBase:        gnDataBaseSource.h
// Purpose:         Implements gnBaseSource for .DataBase files
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

#ifndef _gnDataBaseSource_h_
#define _gnDataBaseSource_h_

#include "libGenome/gnDefs.h"

#include <string>
#include "gnBaseSource.h"


namespace genome {

/**
 * Not yet implemented.
 */
class gnDataBaseSource : public gnBaseSource
{
public:
	virtual ~gnDataBaseSourcee(){}

	virtual gnDataBaseSource* Clone() const = 0;

	virtual gnBaseSource* Clone() const = 0;

	virtual boolean Open( std::string openString ) = 0;
	virtual boolean Open() = 0;
	virtual boolean Close() = 0;
	virtual std::string GetOpenString() const = 0;
	virtual uint32 GetContigListLength() const = 0;
	virtual boolean HasContig( const std::string& name ) const = 0;
	virtual uint32 GetContigID( const std::string& name ) const = 0;
	virtual std::string GetContigName( const uint32 i ) const = 0;
	virtual gnSeqI GetContigSeqLength( const uint32 i ) const = 0;
	virtual const gnFilter* GetFilter() const = 0;
	virtual boolean SetFilter( gnFilter* filter ) = 0;
	virtual boolean Read( const uint64 pos, char* buf, uint32& bufLen) = 0;
	virtual boolean SeqRead( const gnSeqI start, char* buf, uint32& bufLen, const uint32 contigI=ALL_CONTIGS ) = 0;
	virtual gnGenomeSpec *GetSpec() const = 0;
private:
	gnDataBaseSource(){}
};// class gnDataBaseSource


}	// end namespace genome

#endif
	// _gnDataBaseSource_h_
