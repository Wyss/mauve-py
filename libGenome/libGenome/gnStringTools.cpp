/////////////////////////////////////////////////////////////////////////////
// File:            gnStringTools.cpp
// Purpose:         Random string manipulation tools
// Description:     Random string manipulation tools
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

#include "libGenome/gnStringTools.h"

using namespace std;
namespace genome {


void BaseCount(const string& bases, gnSeqI& a_count, gnSeqI& c_count, gnSeqI& g_count, gnSeqI& t_count, gnSeqI& other_count){
	a_count = 0;
	c_count = 0;
	g_count = 0;
	t_count = 0;
	other_count = 0;
	for(uint32 i=0; i < bases.length(); i++){
		if((bases[i] == 'a')||(bases[i] == 'A'))
			a_count++;
		else if((bases[i] == 'c')||(bases[i] == 'C'))
			c_count++;
		else if((bases[i] == 'g')||(bases[i] == 'G'))
			g_count++;
		else if((bases[i] == 't')||(bases[i] == 'T'))
			t_count++;
		else
			other_count++;
	}
}

// removes white space, keeps only one return for multiple returns
unsigned int removeSpace(string &str)
{
	bool onSpace = true;
	unsigned int nbrSpace = 0;
	bool containsReturn = false;
	unsigned int i;
	for( i = str.length(); i > 0 ; i--)
	{
		if( isspace(str[i-1]) )
		{
			nbrSpace++;
			if( (str[i-1] == '\n') || (str[i-1] == '\r') )
				containsReturn = true;
			onSpace = true;
		}
		else
		{
			onSpace = false;
			if( nbrSpace > 0 )
			{
				str.erase( i, nbrSpace-1);
				str[i] = (containsReturn?'\n':' ');
			}
			containsReturn = false;
			nbrSpace = 0;
		}
	}
	if( nbrSpace > 0 )
	{
		str.erase( i, nbrSpace);
	}
	if( str.length() > 0 )
	{
		if( isspace(str[str.length()-1]) ) str.erase(str.length()-1, 1);
	}
	return nbrSpace;
}

void removeEndSpace(string &str)
{
	unsigned int nbrSpace = 0;
	unsigned int i;
	for( i = str.length()-1; i > 0 ; i--)
	{
		if( !isSpace(str[i]) )
			break;
		nbrSpace++;
	}
	if( i != str.length() )
	{
		str.erase(i+1, nbrSpace);
	}
}


bool isNewLine(char ch)
{
	if( (ch == '\n') || (ch == '\r') )
	{
		return true;
	}
	return false;
}

bool isWhiteSpace(char ch)
{
	if( (ch == ' ') || (ch == '\t') )
	{
		return true;
	}
	return false;
}

bool isSpace(char ch)
{
	if( isWhiteSpace(ch) || isNewLine(ch) )
	{
		return true;
	}
	return false;
}

string uintToString(unsigned int value)
{
	string str = "";
	char ch = '\0';
	unsigned int b = 0;
	if( value == 0 )
		str = "0";
	while( value != 0 )
	{
		b = value % 10;
		value /= 10;
		ch = b + 48;
		str = ch + str;
	}
	return str;
}
string ulongToString(unsigned long value)
{
	string str = "";
	char ch = '\0';
	unsigned long b = 0;
	if( value == 0 )
		str = "0";
	while( value != 0 )
	{
		b = value % 10;
		value /= 10;
		ch = b + 48;
		str = ch + str;
	}
	return str;
}

unsigned int parseValue(string &valueString)
{
	unsigned int retValue = 0;
	unsigned int length = valueString.length();
	for( unsigned int i=0; i < length;  i++)
	{
		retValue = (retValue * 10) + (valueString[i] - '0');
	}
	return retValue;
}

int parseUintValue(string &valueString)
{
	int retValue = 0;
	unsigned int length = valueString.length();
	for( unsigned int i=0; i < length;  i++)
	{
		if( isdigit( valueString[i] ) )
		{
			retValue = (retValue * 10) + (valueString[i] - '0');
		}
		else
			break;
	}
	return retValue;
}

int parseIntValue(string &valueString)
{
	int sign = 1;
	int retValue = 0;
	unsigned int length = valueString.length();
	unsigned int i=0;
	for( ; i < length;  i++)
	{
		if( valueString[i] == '-' )
		{
			sign = -1;
			break;
		}
		else if( isdigit( valueString[i] ) )
		{
			retValue = (retValue * 10) + sign * (valueString[i] - '0');
			break;
		}
	}
	i++;
	for( ; i < length;  i++)
	{
		if( isdigit( valueString[i] ) )
		{
			retValue = (retValue * 10) + sign * (valueString[i] - '0');
		}
		else
			break;
	}
	return retValue;
}

vector< string > tokenizeString( const string &str, char delimiter )
{
	return tokenizeString( str.c_str(), str.length(), delimiter );
}

vector< string > tokenizeString( const char* str, unsigned int len, char delimiter )
{
	unsigned int lastIndex = 0 ;
	vector< string > tokenizeVector;
	unsigned int i=0;
	for( i = 0; i < len ; ++i )
	{
		if( str[i] == delimiter )
		{
			if( i > (lastIndex + 1) )
			{
				tokenizeVector.push_back( string( str + lastIndex, i - lastIndex ) );
			}
			lastIndex = i + 1;
		}
	}
	if( i > (lastIndex + 1) )
	{
		tokenizeVector.push_back( string( str + lastIndex, i - lastIndex ) );
	}
	return tokenizeVector;
}


void standardizePathString( string &oFileName )
{
	unsigned int len = oFileName.size();
	for( unsigned int i=0; i < len ; ++i )
	{
		if( oFileName[i] == '\\' )
			oFileName[i] = '/';
	}
}
string getPathString( string oFileName )
{
	string::size_type i = oFileName.rfind('/');
	if( i != string::npos )
		oFileName.erase(i+1, oFileName.length() - (i+1));
//	else
//		oFileName.clear();
	return oFileName;
}
string getFileString( string oFileName )
{
	string::size_type i = oFileName.rfind('/');
	if( i != string::npos )
		oFileName.erase(0, i + 1);
	return oFileName;
}
string getExtString( string oFileName )
{
	string::size_type i = oFileName.rfind('.');
	if( i != string::npos )
		oFileName.erase( 0, i+1);
//	else
//		oFileName.clear();
	return oFileName;
}

string getFileNoExtString( string oFileName )
{
	string::size_type i = oFileName.rfind('/');
	if( i != string::npos )
		oFileName.erase(0, i + 1);
	i = oFileName.rfind('.');
	if( i != string::npos )
		oFileName.erase( i, string::npos);
	return oFileName;
}

}	// end namespace genome

