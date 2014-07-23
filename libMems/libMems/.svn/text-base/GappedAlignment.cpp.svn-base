/*******************************************************************************
 * $Id: GappedAlignment.cpp,v 1.27 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * Please see the file called COPYING for licensing, copying, and modification
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/GappedAlignment.h"
#include <sstream>
#include "libGenome/gnFilter.h"

#include <fstream>

using namespace std;
using namespace genome;
namespace mems {

GappedAlignment::GappedAlignment() : 
AbstractGappedAlignment< SparseAbstractMatch<> >()
{}

GappedAlignment::GappedAlignment( uint seq_count, gnSeqI align_length ) : 
AbstractGappedAlignment< SparseAbstractMatch<> >( seq_count, align_length )
{
	align_matrix.resize(seq_count);
}

void GappedAlignment::SetAlignment( const vector< string >& seq_align ){
	align_matrix = seq_align;
	if( seq_align.size() > 0 )
		SetAlignmentLength(seq_align[0].size());
	else
		SetAlignmentLength(0);
}

std::ostream& operator<<( std::ostream& os, const GappedAlignment& ga ); //write to source.
std::ostream& operator<<( std::ostream& os, const GappedAlignment& ga ){
	os << "GappedAlignmentSeqs: " << ga.SeqCount() << endl;
	os << ga.AlignmentLength();
	for( uint seqI = 0; seqI < ga.SeqCount(); seqI++ )
		os << '\t' << ga.Start( seqI );
	os << endl;
	for( uint seqI = 0; seqI < ga.SeqCount(); seqI++ ){
		os << ga.align_matrix[ seqI ] << endl;
	}
	return os;
};

std::istream& operator>>( std::istream& is, GappedAlignment& ga ); // read from source
std::istream& operator>>( std::istream& is, GappedAlignment& ga ){
	uint seq_count;
	string nuffin;
	is >> nuffin;
	is >> seq_count;
	ga = GappedAlignment( seq_count, 0 );
	is >> nuffin;
	for( uint seqI = 0; seqI < seq_count; seqI++ ){
		int64 startI;
		is >> startI;
		ga.SetStart( seqI, startI );
	}
	for( uint seqI = 0; seqI < seq_count; seqI++ ){
		string seq;
		is >> seq;
		ga.align_matrix.push_back( seq );
	}
	if( ga.align_matrix.size() > 0 )
		ga.SetAlignmentLength( ga.align_matrix[ 0 ].length() );
	return is;
};

}
