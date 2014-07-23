/////////////////////////////////////////////////////////////////////////////
// File:            gnSourceFactory.h
// Purpose:         Manage data sources
// Description:     Manages data sources by tracking open files/databases/URLs
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

#ifndef _gnSourceFactory_h_
#define _gnSourceFactory_h_

#include "libGenome/gnDefs.h"

#include <string>
#include <vector>
#include <map>
#include "libGenome/gnBaseSource.h"

namespace genome {


/**
 * gnSourceFactory is the middle man when acessing a sequence data source
 * It tracks all data sources currently in use, ensuring that a particular
 * data source is only opened and parsed once.  When opening a data source
 * it first tries to interpret the source location as a URL, opening the
 * specified file, or downloading it if necessary.  If that fails, it will
 * attempt to open the source as a file on the local disk.  
 * gnSourceFactory uses the file extension to determine file format, so a
 * file which ends with .fas will be opened by gnFASSource.  Finally,
 * gnSourceFactory can be given directory paths to search when opening a
 * file whose path is not specified. 
 * IMPORTANT: Do not try to instantiate this class.  
 * To use this class do the following:  
 * gnSourceFactory* mySourceFactory = gnSourceFactory::GetSourceFactory();
 */
class GNDLLEXPORT gnSourceFactory{
public:
	~gnSourceFactory();
	
	/**
	 * Returns the current source factory.
	 * @return The current source factory.
	 */
	static gnSourceFactory* GetSourceFactory();
	  // Plugin Sources
	/**
	 * Returns the number of file extension to class mappings.
	 * @return The list size.
	 */
	uint32 GetSourceClassListSize() const;
	/**
	 * Deletes a file extension to class mapping.
	 * @param ext The extension to delete.
	 * @return True if successful.
	 */
	boolean DelSourceClass( const std::string& ext );
	/**
	 * Gets the source class which is mapped to the specified file extension.
	 * @param ext The extension to delete.
	 * @return The class associated with the extension.
	 */
	gnBaseSource* GetSourceClass( const std::string& ext ) const;
	/**
	 * Gets the source class which would be mapped to the std::string.
	 * @param sourceStr The std::string to check, usually a filename.
	 * @return The class associated with the extension.
	 */
	gnBaseSource* MatchSourceClass( const std::string& sourceStr ) const;
	/**
	 * Checks if the specified file extension is recognized.
	 * @param ext The extension to check.
	 * @return True if the extension has a class.
	 */
	boolean HasSourceClass( const std::string& ext ) const;
	/**
	 * Maps the specified file extension to the given source class.
	 * e.g. ".fas" to gnFASSource
	 * @param ext The extension to map.
	 * @param source The class to map
	 * @return True if successful.
	 */
	boolean SetSourceClass( const std::string& ext, const gnBaseSource& source );
	/**
	 * Sets a source class to be the default class for unknown file extensions.
	 * @param source The default class to map
	 * @return True if successful.
	 */
	boolean SetDefaultSourceClass( const gnBaseSource* source );
	/**
	 * Gets the source class which is the default class for unknown file extensions.
	 * @return The default class
	 */
	gnBaseSource* GetDefaultSourceClass() const;
	  // Directory paths to search for sources
	/**
	 * Returns the number of directory paths to search for files.
	 * @return The list size.
	 */
	uint32 GetPathListSize() const;
	/**
	 * Adds the directory to the search path.
	 * @param path The path to add.
	 * @return True if successful.
	 */
	boolean AddPath( const std::string& path );
	/**
	 * Deletes the directory path at index i from the search path list.
	 * @param i The index of the path to delete.
	 * @return True if successful, false if i is out of range.
	 */
	boolean DelPath( uint32 i );
	/**
	 * Inserts the directory path at index i in the search path list.
	 * @param path The path to insert.
	 * @param i The index of the path to insert before.
	 * @return True if successful, false if i is out of range.
	 */
	boolean InsPath( const std::string& path, uint32 i );
	/**
	 * Gets the directory path at index i in the path list.
	 * @param i The index of the path to get.
	 * @return The path or an empty std::string if i is out of range.
	 */
	std::string GetPath( uint32 i ) const;
	/**
	 * Checks the path list for the given path.
	 * @param path The path to look for.
	 * @return True if the path is in the path list.
	 */
	boolean HasPath( std::string path ) const;
	  // Sources
	/**
	 * Returns the number of open data sources.
	 * @return The list size.
	 */
	uint32 GetSourceListSize() const;
	/**
	 * Opens and returns a pointer to a source of genetic sequence data.
	 * If the source has already been opened, AddSource() returns a copy of the existing source class.
	 * @param sourceStr The file name or URL where the source is located.
	 * @param searchPaths Should the path list be searched if the file can't be found.
	 * @return A pointer to the source.
	 */
	gnBaseSource* AddSource( const std::string& sourceStr, boolean searchPaths = true );
	/**
	 * Gets the source at index i in the source list.
	 * @param i The index of the source to get.
	 * @return The source.
	 */
	gnBaseSource* GetSource( uint32 i ) const;
	/**
	 * Deletes the source at index i in the source list.
	 * This will close the associated file, network, or database connection.
	 * @param i The index of the source to delete.
	 * @throws IndexOutOfBounds if i is too large
	 */
	void DelSource( uint32 i );
	/**
	 * Deletes the given source from the source list.
	 * This will close the associated file, network, or database connection.
	 * @param source The source to close.
	 * @return True if successful.
	 */
	boolean DelSource( const gnBaseSource* source );
	/**
	 * Gets the source if it has already been opened.
	 * @param sourceStr The file name or URL where the source is located.
	 * @param searchPaths Should the path list be searched if the file can't be found.
	 * @return A pointer to the source.
	 */
	gnBaseSource* HasSource( std::string sourceStr, boolean searchPaths = true ) const;
		
private:
	gnSourceFactory();
	gnSourceFactory(gnSourceFactory& gnsf);
	gnSourceFactory& operator=(gnSourceFactory& gnsf);

	boolean PathExists( std::string path ) const;
	static boolean GetURL( const std::string& urlStr, std::string& localFile );
	
	std::vector< std::string > m_pathList;
	std::vector< gnBaseSource* > m_sourceList;
	std::map< std::string, gnBaseSource* > m_sourceClassList;
	gnBaseSource* m_pDefaultSourceClass;
};//class gnSourceFactory

// Plugin Sources
inline
uint32 gnSourceFactory::GetSourceClassListSize() const
{
	return m_sourceClassList.size();
}
inline
boolean gnSourceFactory::SetDefaultSourceClass( const gnBaseSource* source )
{
	if(m_pDefaultSourceClass != NULL){
		delete m_pDefaultSourceClass;
	}
	m_pDefaultSourceClass = source->Clone();
	return true;
}
inline
gnBaseSource* gnSourceFactory::GetDefaultSourceClass() const
{
	return m_pDefaultSourceClass;
}

// Directory paths to search for sources
inline
uint32 gnSourceFactory::GetPathListSize() const
{
	return m_pathList.size();
}
inline
uint32 gnSourceFactory::GetSourceListSize() const
{
	return m_sourceList.size();
}



}	// end namespace genome

#endif
	// _gnSourceFactory_h_
