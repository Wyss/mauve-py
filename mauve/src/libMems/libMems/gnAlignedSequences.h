/*******************************************************************************
 * $Id: gnAlignedSequences.h,v 1.5 2004/02/27 23:08:55 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

/////////////////////////////////////////////////////////////////////////////
// File:            gnAlignedSequences.h
// Purpose:         Aligned Sequences class
// Discription:     Provides an alignment interface for any number of alignable
//					sequences (the data of each of which is contained in a 
//                  genome::gnSequence object).
//                  Currently only compatible with ClustalW alignment files.
// Revisions:       
// Version:         A
// Created:         August 3, 2000, 11:55am
// Author:          Brian Gettler
// Last Edited:     May 3, 2001, 4:25pm
// Modified by:     
// Copyright:       (c)
// Licences:         
/////////////////////////////////////////////////////////////////////////////
#ifndef __gnAlignedSequences_h__
#define __gnAlignedSequences_h__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnSequence.h"
#include "libGenome/gnFilter.h"
#include <list>
#include <fstream>
#include <vector>

namespace mems {

// the number of characters in each row of an alignment file
const int MEGA_ALIGN_COLUMNS = 60;

/**
 * gnAlignedSequences allows for the manipulation of aligned sequence
 * data. 
 */

class gnAlignedSequences// : blClone
{
public:
	/**
	 * Empty Constructor, creates a default gnAlignedSequences.
	 */
	gnAlignedSequences();
	/**
	 * Copy Constructor, creates a copy of toCopy.
	 */
	gnAlignedSequences(const gnAlignedSequences &toCopy);


	/**
	 * Returns a vector of supported format names
	 */
	static const std::vector< std::string >& getSupportedFormats();

	/**
	 * Checks whether a particular format name is supported
	 */
	static boolean isSupportedFormat( const std::string& format_name );

	/**
	 * Writes out this sequence alignment in the specified format, 
	 * assuming the format is supported
	 */
	void output( const std::string& format_name, std::ostream& os ) const;

// sequence alignment loading
	/**
	 * Loads the data held in file alignedFileName (in ClustalW format).
	 * @param alignedFileName name of a file containing an alignment.
	 */
	void constructFromClustalW(std::string alignedFileName);
	/**
	 * Loads the data held in file alignedFileName (in Phylip format).
	 * @param alignedFileName name of a file containing an alignment.
	 */
	void constructFromPhylip(std::string alignedFileName);
	/**
	 * Loads the data held in file alignedFileName (in MSF format).
	 * @param alignedFileName name of a file containing an alignment.
	 */
	void constructFromMSF(std::string alignedFileName);
	/**
	 * Loads the data held in file alignedFileName (in Nexus format).
	 * @param alignedFileName name of a file containing an alignment.
	 */
	void constructFromNexus(std::string alignedFileName);
	/**
	 * Loads the data held in file alignedFileName (in Mega format).
	 * @param alignedFileName name of a file containing an alignment.
	 */
	void constructFromMega(std::string alignedFileName);
	
	/**
	 * Reads a single sequence entry in relaxed NEXUS format.  Assumes that
	 * the #NEXUS has already been read off.
	 * @param align_stream	The stream to read data from
	 */
	void constructFromRelaxedNexus( std::istream& align_stream );

	/**
	 * Assigns a file name to the alignment data for purposes of output.
	 * @param name the name of the file.
	 */
	void assignFileName(std::string name);

// output
	/**
	 * Writes alignment using the given output stream (in Phylip format).
	 * @param os the output stream.
	 * @return true if successful.
	 */
	bool outputPhylip(std::ostream& os) const;
	/**
	 * Writes alignment using the given output stream (in ClustalW format).
	 * @param os the output stream.
	 * @return true if successful.
	 */
	bool outputClustalW(std::ostream& os) const;
	/**
	 * Writes alignment using the given output stream (in MSF format).
	 * @param os the output stream.
	 * @return true if successful.
	 */
	bool outputMSF(std::ostream& os) const;
	/**
	 * Writes alignment using the given output stream (in Nexus format).
	 * @param os the output stream.
	 * @return true if successful.
	 */
	bool outputNexus(std::ostream& os) const;
	/**
	 * Writes alignment using the given output stream (in Mega format).
	 * @param os the output stream.
	 * @return true if successful.
	 */
	bool outputMega(std::ostream& os) const;
	/**
	 * Writes alignment in 3-base, codon segments using the given output
	 * stream (in Phylip format).
	 * @param os the output stream.
	 * @return true if successful.
	 */
	bool outputCodon(std::ostream& os) const;
	/**
	 * Writes alignment with consensus using the given output stream
	 * (in Phylip format).
	 * @param os the output stream.
	 * @return true if successful.
	 */
	bool outputWithConsensus(std::ostream& os);

// alignment manipulators that create new gnAlignedSequences
	/**
	 * Create a new alignment that is comprised of all of the segments
	 * in the initial alignment from start to stop (inclusive)
	 * if stop == 0, the alignment ends at the end
	 * @param start the beginning point of the segment.
	 * @param stop the end point of the segment.
	 * @return the new gnAlignedSequences that is created
	 */
	gnAlignedSequences getAlignedSegment(unsigned start, unsigned stop);
	/**
	 * Extracts every codonMultiple-th codon in the reading
	 * frame readingFrame beginning with startCodon
	 * reading frames supported: 1, 2 & 3 (no reverse complementing)
	 * @param readingFrame the codon reading frame.
	 * @param startCodon the number codon in readingFrame with which to begin.
	 * @param codonMultiple the multiple with which codons in readingFrame
	 *        are selected.
	 * @return the new gnAlignedSequences that is created
	 */
	gnAlignedSequences getCodons(int readingFrame, int startCodon, int codonMultiple);
	
	/**
	 * Returns the name of the file associated with this gnAlignedSequences.
	 * @return the alignment file name.
	 */
	std::string getAlignedSequenceFileName();
	/**
	 * Returns the size of each sequence in the alignment (all are identical).
	 * @return the size of the aligned sequences.
	 */
	gnSeqI alignedSeqsSize() const;

	/**
	 * Removes a single sequence from the alignment.
	 * @param seqName the name of the sequence to remove.
	 * @return true if successful (a sequence called seqName exists).
	 */
	bool removeAlignedSeq(std::string seqName);
	/**
	 * Removes a single sequence from the alignment.
	 * @param index the position in the of the sequence to be removed.
	 * @return true if successful (a sequence at index exists).
	 */
	bool removeAlignedSeq(unsigned index);

	/**
	 * Concatenates 2 alignmnets.
	 * @param toConcat the sequence which is appended to *this.
	 */
	void concatenateAlignedSequences(gnAlignedSequences toConcat);
	
	/**
	 * Extracts the variable sites from *this.
	 * @param variableSites the alignment consisting of all variable sites.
	 * @param countGapsAsMismatches true if gaps are to be considered mismatches.
	 */
	void extractVariableSites(gnAlignedSequences &variableSites, bool countGapsAsMismatches);

	/**
	 * Collapses the alignment accross all sequences.
	 * Sequences are compared in terms of base content -
	 * if the sequences of base pairs of equal, the sequences are identical
	 * @return true if there exist identical sequences that are collapsed.
	 */
	bool collapseIdenticalSequences();
	/**
	 * Accesses the alignment and returns the bases at that position in all
	 * sequences.
	 * @param offset the position in the alignment to access.
	 * @return a vector of characters at position offset.
	 */
	std::vector <char> operator[]( const int offset ); //const;
	
	/**
	 * Adds a sequence to the alignment.
	 * @param seqToAdd the sequence data.
	 * @param seqName the sequence's name.
	 */
	void addSequence(std::string& seqToAdd, std::string& seqName);
	/**
	 * Adds a sequence to the alignment.
	 * @param seqToAdd the sequence data.
	 * @param seqName the sequence's name.
	 */
	void addSequence(genome::gnSequence& seqToAdd, std::string& seqName);

	std::list <std::pair <std::string*, std::string*> > alignedSequences;
	std::vector< std::string > sequences;
	std::vector< std::string > names;
	std::vector< int64 > positions;		/**< If this is part of a larger alignment this vector stores start positions within that alignment */
	void seq( uint seqI );

private:

	/**
	 * Reads a relaxed NEXUS format alignment.
	 * @return true if successful.
	 */
	bool readRelaxedNexusAlignment( std::istream& align_stream );
	/**
	 * Aids constructFromClustalW.
	 * @return true if successful.
	 */
	bool readClustalWAlignment();
	/**
	 * Aids constructFromPhylip.
	 * @return true if successful.
	 */
	bool readPhylipAlignment();
	/**
	 * Aids constructFromMSF.
	 * @return true if successful.
	 */
	bool readMSFAlignment();
	/**
	 * Aids constructFromNexus.
	 * @return true if successful.
	 */
	bool readNexusAlignment();
	/**
	 * Aids constructFromMega.
	 * @return true if successful.
	 */
	bool readMegaAlignment();

	/**
	 * Aids readClustalWAlignment.
	 * @param alignmentFile the file that contains the alignment.
	 * @return true if successful.
	 */
	bool constructClustalWAlignedSequenceList(std::ifstream& alignmentFile);
	/**
	 * Aids readPhylipAlignment.
	 * @param alignmentFile the file that contains the alignment.
	 * @return true if successful.
	 */
	bool constructPhylipAlignedSequenceList(std::ifstream& alignmentFile);
	/**
	 * Aids readMSFAlignment.
	 * @param alignmentFile the file that contains the alignment.
	 * @return true if successful.
	 */
	bool constructMSFAlignedSequenceList(std::ifstream& alignmentFile);
	/**
	 * Aids readNexusAlignment.
	 * @param alignmentFile the file that contains the alignment.
	 * @return true if successful.
	 */
	bool constructNexusAlignedSequenceList(std::ifstream& alignmentFile);
	/**
	 * Aids readMegaAlignment.
	 * @param alignmentFile the file that contains the alignment.
	 * @return true if successful.
	 */
	bool constructMegaAlignedSequenceList(std::ifstream& alignmentFile);

	/**
	 * Determines whether a sequence of the given name is present in the list..
	 * @param sequenceName the name to be found.
	 * @param sequenceItr the list iterator to be employed.
	 * @return true if sequenceName is present.
	 */
	bool sequenceNameInList(std::string sequenceName, std::list <std::pair <std::string*, std::string*> >::iterator &sequenceItr);

	/**
	 * Determines whether a sequence of the given name is present in the list.
	 * @param sequenceName the name to be found.
	 * @return the index in the list or -1 if not present
	 */
	int sequenceNameInList( std::string& sequenceName );

	/**
	 * Reads all sequences in the alignment and creates a consensus.
	 * @return true if successful.
	 */
	bool buildConsensus();

	/**
	 * Adds a sequence to the alignment.
	 * @param seqToAdd the sequence data.
	 * @param seqName the sequence's name.
	 * @param consensusStart position in consensus to add sequence.
	 * @param originalConsensus the alignment's consensus.
	 */
	void addSequence(genome::gnSequence seqToAdd, std::string seqName, int consensusStart, std::string originalConsensus);

	/**
	 * Adds all segments in *this to the given alignment.
	 * @param alignment sequences to add.
	 * @param start segment start point.
	 * @param stop segment stop point.
	 */
	void addAllSegments(gnAlignedSequences &alignment, unsigned start, unsigned stop);
	/**
	 * Adds all segments in *this to the given alignment -
	 * replaces gaps with consensus data.
	 * @param alignment sequences to add.
	 * @param start segment start point.
	 * @param stop segment stop point.
	 */
	void addAllSegmentsReplaceGaps(gnAlignedSequences &alignment, unsigned start, unsigned stop);
	/**
	 * Removes all segments across all sequences in *this.
	 * @param start segment start point.
	 * @param stop segment stop point.
	 */
	void removeAllSegments(unsigned start, unsigned stop);
	
	/**
	 * Computes an index for a given base (0-25: a-z).
	 * @param base base to be converted.
	 * @return an index.
	 */
	int determineBaseIndex(char base);
	
	/**
	 * Searches given line for coordinates.
	 * @param line line to search.
	 * @return true if coordinates.
	 */
	bool coordinates(std::string line);

	std::string alignedSequenceFileName;
//	list <pair <string*, genome::gnSequence*> > alignedSequences;
	std::string consensus;
	std::vector <int> indexPositions; // 1->n if a standard alignment, variable for varible sites
}; // gnAlignedSequences


inline
void gnAlignedSequences::assignFileName(std::string name) {alignedSequenceFileName=name;}

inline
std::string gnAlignedSequences::getAlignedSequenceFileName() {return alignedSequenceFileName;}

}

#endif	// __gnAlignedSequences_h__
