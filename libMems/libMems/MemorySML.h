/*******************************************************************************
 * $Id: MemorySML.h,v 1.7 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _MemorySML_h_
#define _MemorySML_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnSequence.h"
#include "libMems/SortedMerList.h"

namespace mems {

/** The MemorySML is an implementation of sorted mer lists which creates and
 *  stores the sorted mer list entirely in memory.  A MemorySML consumes
 *  roughly 32 + alpha_bits bits of memory per character in the sequences.
 *  For unambiguous DNA sequences 4.25 bytes per base are required.
 */
class MemorySML : public SortedMerList
{
public:
	/** 
	 *  Create an empty MemorySML
	 *  Creates an empty MemorySML with the supplied translation
	 *  table and alphabet bit size.  Defaults to DNA settings
	 *  @param table The array used to translate characters into binary code
	 *  @param alpha_bits The number of bits each character consumes in binary
	 */
	MemorySML(const uint8* table = SortedMerList::BasicDNATable(), const uint32 alpha_bits = DNA_ALPHA_BITS);
	MemorySML(const MemorySML& msa);
	MemorySML& operator=(const MemorySML& msa );
	MemorySML* Clone() const;
	
	virtual void Clear();

	virtual void Create(const genome::gnSequence& seq, const uint64 seed);
	virtual boolean Read(std::vector<bmer>& readVector, gnSeqI size, gnSeqI offset = 0);
	virtual void Merge(SortedMerList& sa, SortedMerList& sa2);
	
	virtual bmer operator[](gnSeqI index);
	
protected:

//	virtual void FillSML(const gnSeqI seq_len, vector<gnSeqI>& sml_array);
	std::vector<smlSeqI_t> positions;

};

}

#endif   //_MemorySML_h_
