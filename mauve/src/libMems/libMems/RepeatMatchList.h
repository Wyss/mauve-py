/*******************************************************************************
 * $Id: MatchList.h,v 1.10 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _RepeatMatchList_h_
#define _RepeatMatchList_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <list>
#include "libMems/SortedMerList.h"
#include "libGenome/gnSequence.h"
#include "libMems/Match.h"
#include "libMems/MatchList.h"
//#include <valarray>

namespace mems {

// punt: Need to subclass AbstractMatchList, which can be a MatchList or RepeatMatchList

class RepeatMatchList : public MatchList {
public:
	RepeatMatchList();
	RepeatMatchList( const RepeatMatchList& ml );

	void LoadSequences( std::ostream* log_stream );
	void LoadSMLs( uint mer_size, std::ostream* log_stream );

	/**
	 * Reads a MatchList from an input stream
	 * Sequence and SML file names are read into the seq_filename
	 * and sml_filename vectors, but the actual files are not
	 * opened.  The calling function should load them after
	 * using this method.
	 * @param match_stream The input stream to read from
	 */
	void ReadList( std::istream& match_stream );

	/**
	 *  Writes a MatchList to the designated output stream
	 * @param match_stream The output stream to write to
	 */
	void WriteList( std::ostream& match_stream ) const;
		
	//vector<string> sml_filename;		/**< The file names of the sorted mer list for each sequence, may be empty or null */
	//vector<string> seq_filename;		/**< The file names of the sequence data, may be empty or null */
	//vector<SortedMerList*> sml_table;	/**< The sorted mer list associated with each sequence, may be empty or null */
	//vector<genome::gnSequence*> seq_table;		/**< The actual sequences associated with the matches stored in this list.  Should not be empty or null. */


protected:
	
};

}	// namespace mems

#endif 


