/*******************************************************************************
 * $Id: Files.h,v 1.23 2004/04/19 23:10:13 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef __libMems_Files_h__
#define __libMems_Files_h__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// for CreateTempFilename
#ifdef WIN32
#include "windows.h"
#else
#include "unistd.h"
#endif

#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/exception.hpp"
#include "boost/algorithm/string.hpp"
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>


/**
 * Register a file name to be deleted before the process exits
 * When passed an empty string, it does not add to the list of files to delete
 * @param fname	The name of a file to delete, empty strings are ignored
 * @return A vector of file names registered for deletion
 */
std::vector< std::string >& registerFileToDelete( std::string fname = "" );

inline
std::vector< std::string >& registerFileToDelete( std::string fname ) {
	// since this vector is needed when atexit() is called we allocate it
	// on the heap so its destructor won't get called
	static std::vector< std::string >* files = new std::vector< std::string >();
#pragma omp critical
{
	if( fname != "" )
		files->push_back( fname );
}
	return *files;
}

void deleteRegisteredFiles();
inline
void deleteRegisteredFiles() {
	// don't be a slob, clean up after yourself:
	// delete any files that are laying around
	std::vector< std::string >& del_files = registerFileToDelete();
	for( int fileI = 0; fileI < del_files.size(); fileI++ )
		boost::filesystem::remove( del_files[ fileI ] );
	del_files.clear();	// clear the deleted files from the list
}


/**
 * Create a temporary file
 */
std::string CreateTempFileName(const std::string& prefix);


/* shamelessly ripped from wxWidgets and boostified*/
inline
std::string CreateTempFileName(const std::string& prefix)
{
   std::string dir, name, ret_path;
#ifdef WIN32
        char buf[MAX_PATH + 1];
#else
        char buf[PATH_MAX + 1];
#endif
        boost::filesystem::path path( prefix );
        dir = path.branch_path().string();
        name = path.leaf();
        if( name == "/" )
        {
                dir += name;
                name.clear();
        }
#if defined(WIN32)

    if ( dir.size() == 0 )
    {
                strncpy(buf, dir.c_str(), MAX_PATH);
        if ( !::GetTempPath(MAX_PATH, buf) )
			std::cerr << "GetTempPath\n";

                dir = buf;
        if ( dir.size()==0 )
            dir = ".";  // GetTempFileName() fails if we pass it an emptystd::string
    }
    else // we have a dir to create the file in
    {
        // ensure we use only the back slashes as GetTempFileName(), unlike all
        // the other APIs, is picky and doesn't accept the forward ones
                boost::algorithm::replace_all( dir, "/", "\\" );
    }

        strncpy(buf, path.string().c_str(), MAX_PATH);
        if ( !::GetTempFileName(dir.c_str(), name.c_str(), 0, buf) )
    {
        std::cerr << "GetTempFileName\n";
                path = boost::filesystem::path();
    }
        ret_path = buf;

#else // !Windows
    if ( dir.empty() )
    {
        char* env_val = getenv("TMP");
        dir = env_val != NULL ? env_val : "";

        if ( dir.size() == 0 ){
            env_val = getenv("TMPDIR");
            dir = env_val != NULL ? env_val : "";
        }

        if ( dir.size() == 0 ){
            env_val = getenv("TEMP");
            dir = env_val != NULL ? env_val : "";
        }

        if ( dir.size()==0 )
        {
            // default
            #ifdef __DOS__
                dir = ".";
            #else
                dir = "/tmp";
            #endif
        }
    }

    path = dir;
    path /= name;

    // we need to copy the path to the buffer in which mkstemp() can modify it
       std::string path_str = path.string();
        path_str += "XXXXXX";  // scratch space for mkstemp()
        strncpy( buf, path_str.c_str(), path_str.size()+1 );

#if defined(HAVE_MKSTEMP)
    // cast is safe because thestd::string length doesn't change
    int fdTemp = mkstemp( buf );
    if ( fdTemp == -1 )
    {
        // this might be not necessary as mkstemp() on most systems should have
        // already done it but it doesn't hurt neither...
//        path.clear();
    }
    else // mkstemp() succeeded
    {
                ret_path = buf;
        close(fdTemp);
    }
#else // !HAVE_MKSTEMP

#ifdef HAVE_MKTEMP
    // same as above
    if ( mktemp( buf ) )
                ret_path = buf;

#else // !HAVE_MKTEMP (includes __DOS__)
    // generate the unique file name ourselves
    unsigned my_pid = 0;
    #ifndef __DOS__
    my_pid = getpid();
    #endif

        std::ostringstream oss;

    std::string oss_str;
    static const size_t numTries = 1000;
    for ( size_t n = 0; n < numTries; n++ )
    {
	std::ostringstream oss;
                oss << path.string() << my_pid << "." << std::setfill('0') << std::setw(3) << n;
        // 3 hex digits is enough for numTries == 1000 < 4096
                boost::filesystem::path pathTry( oss.str() );
        oss_str = oss.str();
        if ( !boost::filesystem::exists(pathTry) )
            break;

    }

    ret_path = oss_str;
#endif // HAVE_MKTEMP/!HAVE_MKTEMP

#endif // HAVE_MKSTEMP/!HAVE_MKSTEMP

#endif // Windows/!Windows

        return ret_path;
}


#endif // __libMems_Files_h__

