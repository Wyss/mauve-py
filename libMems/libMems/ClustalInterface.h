/*******************************************************************************
 * $Id: ClustalInterface.h,v 1.12 2004/04/19 23:10:50 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _ClustalInterface_h_
#define _ClustalInterface_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/NumericMatrix.h"
#include "libGenome/gnFilter.h"
#include "libGenome/gnSequence.h"
#include "libMems/GappedAlignment.h"
#include "libMems/GappedAligner.h"

// attempt to auto-link the ClustalW library on windows
#if defined(WIN64)&&defined(NDEBUG)&&!defined(FASTDEBUG)&&defined(_OPENMP)
#pragma comment(lib, "ClustalW64omp.lib")
#endif
#if defined(WIN64)&&defined(FASTDEBUG)&&defined(_OPENMP)
#pragma comment(lib, "ClustalW64fdomp.lib")
#endif
#if defined(WIN32)&&!defined(WIN64)&&defined(NDEBUG)&&!defined(FASTDEBUG)&&defined(_OPENMP)
#pragma comment(lib, "ClustalWomp.lib")
#endif
#if defined(WIN32)&&!defined(WIN64)&&defined(FASTDEBUG)&&defined(_OPENMP)
#pragma comment(lib, "ClustalWfdomp.lib")
#endif
#if defined(WIN64)&&defined(NDEBUG)&&!defined(FASTDEBUG)&&!defined(_OPENMP)
#pragma comment(lib, "ClustalW64.lib")
#endif
#if defined(WIN64)&&defined(FASTDEBUG)&&!defined(_OPENMP)
#pragma comment(lib, "ClustalW64fd.lib")
#endif
#if defined(WIN32)&&!defined(WIN64)&&defined(NDEBUG)&&!defined(FASTDEBUG)&&!defined(_OPENMP)
#pragma comment(lib, "ClustalW.lib")
#endif
#if defined(WIN32)&&!defined(WIN64)&&defined(FASTDEBUG)&&!defined(_OPENMP)
#pragma comment(lib, "ClustalWfd.lib")
#endif


namespace mems {

class ClustalInterface : public GappedAligner {
public:
	/**
	 * Returns a reference to a usable ClustalInterface
	 */
	static ClustalInterface& getClustalInterface();
	/**
	 * Attempts to perform a multiple alignment using ClustalW between
	 * <code>r_begin</code> and <code>r_end</code>
	 */
	boolean Align( GappedAlignment& cr, Match* r_begin, Match* r_end, std::vector< genome::gnSequence* >& seq_table );
	/**
	 * Set the distance matrix to use when computing alignments, writes the guide tree to the location
	 * specified in <code>tree_filename</code>
	 * @param distance_matrix An NxN distance matrix for the sequences
	 * @param tree_filename The output file name for the guide tree
	 */
	void SetDistanceMatrix( NumericMatrix< double >& distance_matrix, std::string& tree_filename );
	/**
	 * Set the minimum flank size used to anchor alignments on the sequences
	 */
	void SetMinFlankSize( gnSeqI min_flank ){ min_flank_size = min_flank; }
	
	/**
	 * Try using the guide tree in the file given by tree_filename.  Throws an
	 * exception if the tree file couldn't be loaded
	 * @param tree_filename		The path to the guide tree file
	 * @param dist_mat			The distance matrix relating sequences
	 * @param seq_count			The number of genomes in the guide tree file
	 */
	void setGuideTree( std::string& tree_filename, NumericMatrix< double >& dist_mat, uint seq_count );
	
	/** returns true if a guide tree has been loaded already */
	boolean guideTreeLoaded() const { return distance_matrix.cols() > 0; };
	
	void SetDistanceMatrix( NumericMatrix< double >& distance_matrix, std::string& tree_filename, boolean reread_tree );
protected:
	boolean CallClustal( std::vector< std::string >& seq_table );
	NumericMatrix< double > distance_matrix;
	gnSeqI min_flank_size;
	int clustal_score_cutoff;
	bool allocated_aln;
private:
	ClustalInterface( const ClustalInterface& ci ){ *this = ci; }
	ClustalInterface& operator=( const ClustalInterface& ci );
	ClustalInterface();
};

}

#endif // _ClustalInterface_h_
