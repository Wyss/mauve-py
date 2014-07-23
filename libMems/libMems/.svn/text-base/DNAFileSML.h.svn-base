/*******************************************************************************
 * $Id: DNAFileSML.h,v 1.6 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _DNAFileSML_h_
#define _DNAFileSML_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/FileSML.h"

namespace mems {

/**
 *  The seed pattern for DNA SMLs must be palindromic
 */
class DNAFileSML : public FileSML
{
public:
	DNAFileSML();
	
	/** 
	 *  Load or create a DNAFileSML ()
	 *  Attempts to load a DNA sorted mer list from the named file if it exists.
	 *  If the given file does not exist it creates an empty DNAFileSML with 
	 *  the supplied translation table and alphabet bit size.
	 *  @param fname The name of the file to create.
	 *  @param table The array used to translate characters into binary code
	 *  @param alpha_bits The number of bits each character consumes in binary
	 */
	DNAFileSML(const std::string& fname, const uint8* table = SortedMerList::BasicDNATable(), const uint32 alpha_bits = DNA_ALPHA_BITS);
	DNAFileSML(const SortedMerList& sa);
	DNAFileSML& operator=(const DNAFileSML& msa );
	
	DNAFileSML* Clone() const;
	
	virtual uint64 GetMer(gnSeqI position) const;
	
	virtual uint32 FormatVersion();

	virtual uint64 GetSeedMer( gnSeqI offset ) const;

protected:
	virtual void FillSML(const genome::gnSequence& seq, std::vector<bmer>& sml_array);
	virtual uint32 CalculateMaxMerSize() const;
	virtual uint64 GetNeededMemory(gnSeqI len);
};

// version 3 was original DNAFileSML format
// version 4 was introduction of inexact seeds
// version 5 was fix in header struct for 64-bit seed size
inline
uint32 DNAFileSML::FormatVersion(){
	static uint32 f_version = 5;
	return f_version;
}

}

#endif   //_DNAFileSML_h_
