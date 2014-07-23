/////////////////////////////////////////////////////////////////////////////
// File:            gnRNASequence.h
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

#ifndef _gnRNASequence_h_
#define _gnRNASequence_h_

#include "libGenome/gnDefs.h"

#include <string>
#include <list>
#include "libGenome/gnSequence.h"
#include "libGenome/gnFilter.h"

namespace genome {


/**
 * gnRNASequence is a special kind of gnSequence which can be used for RNA sequences
 * It sets the default filters and comparators to the RNA filters and comparators.
 */

class GNDLLEXPORT gnRNASequence : public gnSequence
{
public:
	/**
	 * Empty Constructor, creates an empty gnRNASequence.
	 */
	gnRNASequence();
	/**
	 * Creates a gnRNASequence with a single contig containing the bases in "seq".
	 * @param seq The null terminated array of base pairs to use.
	 */
	gnRNASequence( const gnSeqC* seq );
	/**
	 * Creates a gnRNASequence with a single contig containing the bases in "str".
	 * @param str The base pairs to use.
	 */
	gnRNASequence( const std::string& str );
	/**
	 * Creates a gnRNASequence with the contigs stored in "gngs".
	 * @param gngs the gnGenomeSpec to get contigs from.
	 */
	gnRNASequence( const gnGenomeSpec& gngs );
	/**
	 * Creates a gnRNASequence with the contigs stored in "gnfs".
	 * @param gnfs the gnFragmentSpec to get contigs from.
	 */
	gnRNASequence( const gnFragmentSpec& gnfs );
	/**
	 * Creates a gnRNASequence with the contigs stored in "gncs".
	 * @param gncs the gnContigSpec to get contigs from.
	 */
	gnRNASequence( const gnContigSpec& gncs );
	/**
	 * Creates a gnRNASequence with a single contig containing the bases in "bases".
	 * @param bases The base pairs to use
	 * @param length The length of the base pair array.
	 */
	gnRNASequence( gnSeqC *bases, const gnSeqI length);
	/**
	 * Copies the gnRNASequence "seq".
	 * @param seq The gnRNASequence to copy.
	 */
	gnRNASequence( const gnRNASequence& seq);
private:
	gnGenomeSpec *spec;
	list<const gnBaseFilter*> filter_list;
	const gnCompare* comparator;
}; // class gnRNASequence

inline
gnRNASequence::gnRNASequence() : gnSequence(){
	filter_list.push_back(gnFilter::fullRNASeqFilter());
	comparator = gnCompare::RNASeqCompare();
}
inline
gnRNASequence::gnRNASequence( const gnSeqC* seq ) : gnSequence(seq){
	filter_list.push_back(gnFilter::fullRNASeqFilter());
	comparator = gnCompare::RNASeqCompare();
}
inline
gnRNASequence::gnRNASequence( const std::string& str ) : gnSequence(str){
	filter_list.push_back(gnFilter::fullRNASeqFilter());
	comparator = gnCompare::RNASeqCompare();
}
inline
gnRNASequence::gnRNASequence( const gnGenomeSpec& gngs ) : gnSequence(gngs){
	filter_list.push_back(gnFilter::fullRNASeqFilter());
	comparator = gnCompare::RNASeqCompare();
}
inline
gnRNASequence::gnRNASequence( const gnFragmentSpec& gnfs ) : gnSequence(gnfs){
	filter_list.push_back(gnFilter::fullRNASeqFilter());
	comparator = gnCompare::RNASeqCompare();
}
inline
gnRNASequence::gnRNASequence( const gnContigSpec& gncs ) : gnSequence(gncs){
	filter_list.push_back(gnFilter::fullRNASeqFilter());
	comparator = gnCompare::RNASeqCompare();
}
inline
gnRNASequence::gnRNASequence( gnSeqC *bases, const gnSeqI length) : gnSequence(bases, length){
	filter_list.push_back(gnFilter::fullRNASeqFilter());
	comparator = gnCompare::RNASeqCompare();
}
inline
gnRNASequence::gnRNASequence( const gnRNASequence& seq) : gnSequence(seq){
	filter_list.push_back(gnFilter::fullRNASeqFilter());
	comparator = gnCompare::RNASeqCompare();
}


}	// end namespace genome

#endif
	// _gnRNASequence_h_
