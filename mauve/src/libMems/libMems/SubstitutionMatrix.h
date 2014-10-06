/*******************************************************************************
 * $Id: SubstitutionMatrix.h,v 1.7 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef __SubstitutionMatrix_h__
#define __SubstitutionMatrix_h__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnSequence.h"
#include <iostream>
#include <sstream>

namespace mems {

typedef int score_t;
static const score_t hoxd_matrix[4][4] = 
{ 
	{91,	-114,	-31,	-123}, // A

	{-114,	100,	-125,	-31}, // C

	{-31,	-125,	100,	-114}, // G

	{-123,	-31,	-114,	91}, // T
};

static const score_t default_gap_open = -400;
static const score_t default_gap_extend = -30;

class PairwiseScoringScheme
{
public:
	score_t matrix[4][4];	/**< 4x4 nucleotide substitution matrix */
	score_t gap_open;	/**< gap open penalty */
	score_t gap_extend;	/**< gap extend penalty */

	PairwiseScoringScheme( const score_t matrix[4][4], score_t gap_open, score_t gap_extend )
	{
		setMatrix(matrix);
		this->gap_open = gap_open;
		this->gap_extend = gap_extend;
	}

	PairwiseScoringScheme(){ *this = PairwiseScoringScheme( hoxd_matrix, default_gap_open, default_gap_extend ); }
	PairwiseScoringScheme& operator=( const PairwiseScoringScheme& pss )
	{
		setMatrix(pss.matrix);
		this->gap_open = pss.gap_open;
		this->gap_extend = pss.gap_extend;
		return *this;
	}
	void setMatrix( const score_t matrix[4][4] )
	{
		for( int i = 0; i < 4; ++i )
			for( int j = 0; j < 4; ++j )
				this->matrix[i][j] = matrix[i][j];
	}
};

static PairwiseScoringScheme& getDefaultScoringScheme()
{
	static PairwiseScoringScheme pss( hoxd_matrix, default_gap_open, default_gap_extend );
	return pss;
}

void readSubstitutionMatrix( std::istream& is, score_t matrix[4][4] );

inline
void readSubstitutionMatrix( std::istream& is, score_t matrix[4][4] )
{
	std::string tmp;
	std::getline( is, tmp );	// first line contains header info
	std::getline( is, tmp );	// second line contains sub mat column labels
	std::stringstream ss( tmp );
	std::string letter;
	bool format_ok = true;
	ss >> letter;
	format_ok = format_ok && letter == "A";
	ss >> letter;
	format_ok = format_ok && letter == "C";
	ss >> letter;
	format_ok = format_ok && letter == "G";
	ss >> letter;
	format_ok = format_ok && letter == "T";
	ss >> letter;
	format_ok = format_ok && letter == "N";
	if( !format_ok )
	{
		std::cerr << "Invalid substitution matrix format\n";
		throw "Invalid substitution matrix format\n";
	}

	for( int i = 0; i < 4; i++ )
	{
		is >> letter;	// the first character on each line should be a letter
		for( int j = 0; j < 4; j++ )
			is >> matrix[i][j];
		is >> letter;	// this should be the N sub score (which gets ignored)
	}
}

}

#endif // __SubstitutionMatrix_h__
