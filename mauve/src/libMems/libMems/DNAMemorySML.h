/*******************************************************************************
 * $Id: DNAMemorySML.h,v 1.3 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _DNAMemorySML_h_
#define _DNAMemorySML_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnSequence.h"
#include "libMems/MemorySML.h"

namespace mems {

/** The DNAMemorySML is an implementation of sorted mer lists which creates and
 *  stores the sorted mer list entirely in memory.  A DNAMemorySML consumes
 *  roughly 32 + alpha_bits bits of memory per character in the sequences.
 *  For unambiguous DNA sequences 4.25 bytes per base are required.
 *  The seed pattern for DNA SMLs must be palindromic
 */
class DNAMemorySML : public MemorySML
{
public:
	/** 
	 *  Create an empty DNAMemorySML
	 *  Creates an empty DNAMemorySML with the supplied translation
	 *  table and alphabet bit size.  Defaults to DNA settings
	 *  @param table The array used to translate characters into binary code
	 *  @param alpha_bits The number of bits each character consumes in binary
	 */
	DNAMemorySML(const uint8* table = SortedMerList::BasicDNATable(), const uint32 alpha_bits = DNA_ALPHA_BITS);
	DNAMemorySML(const DNAMemorySML& msa);
	DNAMemorySML(const SortedMerList& sa);
	DNAMemorySML& operator=(const DNAMemorySML& msa );
	DNAMemorySML* Clone() const;
	
	
	virtual uint64 GetMer(gnSeqI offset) const;
	virtual uint64 GetSeedMer( gnSeqI offset ) const;
	
protected:

	virtual void FillSML(const genome::gnSequence& seq, std::vector<bmer>& sml_array);

};

}

#endif   //_DNAMemorySML_h_
