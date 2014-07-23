/////////////////////////////////////////////////////////////////////////////
// File:            gnProteinSequence.h
// Purpose:         Sequence class
// Description:     Provides a high level sequence interface to all types of
//					sequence data.
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

#ifndef _gnProteinSequence_h_
#define _gnProteinSequence_h_

#include "libGenome/gnDefs.h"

#include <string>
#include <list>
#include "libGenome/gnSequence.h"
#include "libGenome/gnFilter.h"


namespace genome {

/**
 * gnProteinSequence is a special kind of gnSequence which can be used for RNA sequences
 * It sets the default filters and comparators to the RNA filters and comparators.
 */

class GNDLLEXPORT gnProteinSequence : public gnSequence
{
public:
	/**
	 * Empty Constructor, creates an empty gnProteinSequence.
	 */
	gnProteinSequence();
	/**
	 * Creates a gnProteinSequence with a single contig containing the bases in "seq".
	 * @param seq The null terminated array of base pairs to use.
	 */
	gnProteinSequence( const gnSeqC* seq );
	/**
	 * Creates a gnProteinSequence with a single contig containing the bases in "str".
	 * @param str The base pairs to use.
	 */
	gnProteinSequence( const std::string& str );
	/**
	 * Creates a gnProteinSequence with the contigs stored in "gngs".
	 * @param gngs the gnGenomeSpec to get contigs from.
	 */
	gnProteinSequence( const gnGenomeSpec& gngs );
	/**
	 * Creates a gnProteinSequence with the contigs stored in "gnfs".
	 * @param gnfs the gnFragmentSpec to get contigs from.
	 */
	gnProteinSequence( const gnFragmentSpec& gnfs );
	/**
	 * Creates a gnProteinSequence with the contigs stored in "gncs".
	 * @param gncs the gnContigSpec to get contigs from.
	 */
	gnProteinSequence( const gnContigSpec& gncs );
	/**
	 * Creates a gnProteinSequence with a single contig containing the bases in "bases".
	 * @param bases The base pairs to use
	 * @param length The length of the base pair array.
	 */
	gnProteinSequence( gnSeqC *bases, const gnSeqI length);
	/**
	 * Copies the gnProteinSequence "seq".
	 * @param seq The gnProteinSequence to copy.
	 */
	gnProteinSequence( const gnProteinSequence& seq);
private:
	gnGenomeSpec *spec;
	list<const gnBaseFilter*> filter_list;
	const gnCompare* comparator;
}; // class gnProteinSequence

inline
gnProteinSequence::gnProteinSequence() : gnSequence(){
	filter_list.push_back(gnFilter::proteinSeqFilter());
	comparator = gnCompare::ProteinSeqCompare();
}
inline
gnProteinSequence::gnProteinSequence( const gnSeqC* seq ) : gnSequence(seq){
	filter_list.push_back(gnFilter::proteinSeqFilter());
	comparator = gnCompare::ProteinSeqCompare();
}
inline
gnProteinSequence::gnProteinSequence( const std::string& str ) : gnSequence(str){
	filter_list.push_back(gnFilter::proteinSeqFilter());
	comparator = gnCompare::ProteinSeqCompare();
}
inline
gnProteinSequence::gnProteinSequence( const gnGenomeSpec& gngs ) : gnSequence(gngs){
	filter_list.push_back(gnFilter::proteinSeqFilter());
	comparator = gnCompare::ProteinSeqCompare();
}
inline
gnProteinSequence::gnProteinSequence( const gnFragmentSpec& gnfs ) : gnSequence(gnfs){
	filter_list.push_back(gnFilter::proteinSeqFilter());
	comparator = gnCompare::ProteinSeqCompare();
}
inline
gnProteinSequence::gnProteinSequence( const gnContigSpec& gncs ) : gnSequence(gncs){
	filter_list.push_back(gnFilter::proteinSeqFilter());
	comparator = gnCompare::ProteinSeqCompare();
}
inline
gnProteinSequence::gnProteinSequence( gnSeqC *bases, const gnSeqI length) : gnSequence(bases, length){
	filter_list.push_back(gnFilter::proteinSeqFilter());
	comparator = gnCompare::ProteinSeqCompare();
}
inline
gnProteinSequence::gnProteinSequence( const gnProteinSequence& seq) : gnSequence(seq){
	filter_list.push_back(gnFilter::proteinSeqFilter());
	comparator = gnCompare::ProteinSeqCompare();
}


}	// end namespace genome

#endif
	// _gnProteinSequence_h_
