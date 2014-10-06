/////////////////////////////////////////////////////////////////////////////
// File:            gnStringTools.h
// Purpose:         Random std::string manipulation tools
// Description:     Random std::string manipulation tools
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

#ifndef _gnStringTools_h_
#define _gnStringTools_h_

#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include "libGenome/gnDefs.h"

namespace genome {


unsigned int removeSpace(std::string &str);
void removeEndSpace(std::string &str);

bool isNewLine(char ch);
bool isWhiteSpace(char ch);
bool isSpace(char ch);

std::string uintToString(unsigned int value);

std::string ulongToString(unsigned long value);

std::string charArrayToString( char *charArray, unsigned int length);

bool isBase(char base);
void BaseCount(const std::string& bases, gnSeqI& a_count, gnSeqI& c_count, gnSeqI& g_count, gnSeqI& t_count, gnSeqI& other_count);

unsigned int parseValue(std::string &valueString);
int parseUintValue(std::string &valueString);
int parseIntValue(std::string &valueString);

std::vector< std::string > tokenizeString( const std::string &str, char delimiter  = '\t' );
std::vector< std::string > tokenizeString( const char* str, unsigned int len, char delimiter = '\t' );

// **** MAKE SURE YOU STANDARDIZE ALL PATHS ****
void standardizePathString( std::string &oFileName );

// **** ONLY '/' IS USED for file seperation on all below functions ****

// make more robust to handle dir and files
std::string getPathString( std::string oFileName );
// look for '/' at end
std::string getFileString( std::string oFileName );
// look for '/' before '.'
std::string getExtString( std::string oFileName );
// look for '/' before '.'
std::string getFileNoExtString( std::string oFileName );


}	// end namespace genome

#endif  // _gnStringTools_h_
