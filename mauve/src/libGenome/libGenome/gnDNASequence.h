/////////////////////////////////////////////////////////////////////////////
// File:            gnDNASequence.h
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

#ifndef _gnDNASequence_h_
#define _gnDNASequence_h_

#include "libGenome/gnDefs.h"

#include <string>
#include <list>
#include "libGenome/gnSequence.h"
#include "libGenome/gnFilter.h"

namespace genome {

/**
 * gnDNASequence is a special kind of gnSequence which can be used for DNA sequences
 * It sets the default filters and comparators to the DNA filters and comparators.
 */

class GNDLLEXPORT gnDNASequence : public gnSequence
{
public:
	/**
	 * Empty Constructor, creates an empty gnDNASequence.
	 */
	gnDNASequence();
	/**
	 * Creates a gnDNASequence with a single contig containing the bases in "seq".
	 * @param seq The null terminated array of base pairs to use.
	 */
	gnDNASequence( const gnSeqC* seq );
	/**
	 * Creates a gnDNASequence with a single contig containing the bases in "str".
	 * @param str The base pairs to use.
	 */
	gnDNASequence( const std::string& str );
	/**
	 * Creates a gnDNASequence with the contigs stored in "gngs".
	 * @param gngs the gnGenomeSpec to get contigs from.
	 */
	gnDNASequence( const gnGenomeSpec& gngs );
	/**
	 * Creates a gnDNASequence with the contigs stored in "gnfs".
	 * @param gnfs the gnFragmentSpec to get contigs from.
	 */
	gnDNASequence( const gnFragmentSpec& gnfs );
	/**
	 * Creates a gnDNASequence with the contigs stored in "gncs".
	 * @param gncs the gnContigSpec to get contigs from.
	 */
	gnDNASequence( const gnContigSpec& gncs );
	/**
	 * Creates a gnDNASequence with a single contig containing the bases in "bases".
	 * @param bases The base pairs to use
	 * @param length The length of the base pair array.
	 */
	gnDNASequence( gnSeqC *bases, const gnSeqI length);
	/**
	 * Copies the gnDNASequence "seq".
	 * @param seq The gnDNASequence to copy.
	 */
	gnDNASequence( const gnDNASequence& seq);
private:
	gnGenomeSpec *spec;
	list<const gnBaseFilter*> filter_list;
	const gnCompare* comparator;
}; // class gnDNASequence

inline
gnDNASequence::gnDNASequence() : gnSequence(){
	filter_list.push_back(gnFilter::fullDNASeqFilter());
	comparator = gnCompare::DNASeqCompare();
}
inline
gnDNASequence::gnDNASequence( const gnSeqC* seq ) : gnSequence(seq){
	filter_list.push_back(gnFilter::fullDNASeqFilter());
	comparator = gnCompare::DNASeqCompare();
}
inline
gnDNASequence::gnDNASequence( const std::string& str ) : gnSequence(str){
	filter_list.push_back(gnFilter::fullDNASeqFilter());
	comparator = gnCompare::DNASeqCompare();
}
inline
gnDNASequence::gnDNASequence( const gnGenomeSpec& gngs ) : gnSequence(gngs){
	filter_list.push_back(gnFilter::fullDNASeqFilter());
	comparator = gnCompare::DNASeqCompare();
}
inline
gnDNASequence::gnDNASequence( const gnFragmentSpec& gnfs ) : gnSequence(gnfs){
	filter_list.push_back(gnFilter::fullDNASeqFilter());
	comparator = gnCompare::DNASeqCompare();
}
inline
gnDNASequence::gnDNASequence( const gnContigSpec& gncs ) : gnSequence(gncs){
	filter_list.push_back(gnFilter::fullDNASeqFilter());
	comparator = gnCompare::DNASeqCompare();
}
inline
gnDNASequence::gnDNASequence( gnSeqC *bases, const gnSeqI length) : gnSequence(bases, length){
	filter_list.push_back(gnFilter::fullDNASeqFilter());
	comparator = gnCompare::DNASeqCompare();
}
inline
gnDNASequence::gnDNASequence( const gnDNASequence& seq) : gnSequence(seq){
	filter_list.push_back(gnFilter::fullDNASeqFilter());
	comparator = gnCompare::DNASeqCompare();
}


}	// end namespace genome

#endif
	// _gnDNASequence_h_
