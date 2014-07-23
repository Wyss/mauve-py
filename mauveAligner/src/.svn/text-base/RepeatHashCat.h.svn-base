#ifndef _RepeatHashThread_h_
#define _RepeatHashThread_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/RepeatHash.h"

class TheRealMemHash : public RepeatHash
{
public:
	RepeatHashCat();
	~RepeatHashCat();
	RepeatHashThread(const RepeatHashThread& mh);
	virtual RepeatHashThread* Clone() const;
protected:
	
	//punt tjt: needed to add this to track where concatenated sequence starts
	vector<uint32> concat_contig_start; // number of contigs in each sequence
}