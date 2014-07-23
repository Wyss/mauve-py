/*******************************************************************************
 * $Id: GappedAligner.h,v 1.12 2004/04/19 23:10:50 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _GappedAligner_h_
#define _GappedAligner_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnSequence.h"
#include "libMems/GappedAlignment.h"
#include "libMems/Match.h"

namespace mems {

class GappedAligner {
public:
	GappedAligner(){ max_alignment_length = 10000; }	// default to something
	GappedAligner& operator=( const GappedAligner& ga )
	{ 
		max_alignment_length = ga.max_alignment_length;
		return *this;
	}
	/**
	 * Set the maximum allowed length for a gapped alignment.  Sequences above this length
	 * threshold will be ignored.
	 * @param max_length The maximum length
	 */
	void SetMaxAlignmentLength( gnSeqI len ){max_alignment_length = len;}
	virtual boolean Align( GappedAlignment& cr, Match* r_begin, Match* r_end, std::vector< genome::gnSequence* >& seq_table ) = 0;
protected:
	gnSeqI max_alignment_length;
};





boolean getInterveningCoordinates( std::vector< genome::gnSequence* >& seq_table, Match* r_begin, Match* r_end, uint seqI, int64& gap_lend, int64& gap_rend );

inline
boolean getInterveningCoordinates( std::vector< genome::gnSequence* >& seq_table, Match* r_begin, Match* r_end, uint seqI, int64& gap_lend, int64& gap_rend ){
	// skip this sequence if it's undefined
	if( (r_end != NULL && r_end->Start( seqI ) == NO_MATCH) ||
		(r_begin != NULL && r_begin->Start( seqI ) == NO_MATCH) ){
		gap_lend = 0;
		gap_rend = 0;
		return true;
	}
			
	// determine the size of the gap
	gap_rend = r_end != NULL ? r_end->Start( seqI ) : seq_table[ seqI ]->length() + 1;
	gap_lend = r_begin != NULL ? r_begin->End( seqI ) + 1 : 1;
	if( gap_rend < 0 || gap_lend < 0 ){
		gap_rend = r_begin != NULL ? -r_begin->Start( seqI ) : seq_table[ seqI ]->length() + 1;
		gap_lend = r_end != NULL ? -r_end->Start( seqI ) + r_end->Length() : 1;
	}
	if( gap_rend <= 0 || gap_lend <= 0 ){
		// if either is still < 0 then there's a problem...
		genome::ErrorMsg( "Error constructing intervening coordinates" );
	}
	return true;
}

}

#endif // _GappedAligner_h_
