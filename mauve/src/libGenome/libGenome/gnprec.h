#pragma once


#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers

#ifndef _gnprec_h_
#define _gnprec_h_

//general headers
#include "libGenome/gnSetup.h"
#include "libGenome/gnDefs.h"
#include "libGenome/gnClone.h"
#include "libGenome/gnDebug.h"
#include "libGenome/gnVersion.h"
#include "libGenome/gnStringTools.h"
#include "libGenome/gnException.h"
#include "libGenome/gnExceptionCode.h"
//sequence headers
#include "libGenome/gnSequence.h"
#include "libGenome/gnDNASequence.h"
#include "libGenome/gnRNASequence.h"
#include "libGenome/gnProteinSequence.h"
//source headers
#include "libGenome/gnSourceFactory.h"
#include "libGenome/gnBaseSource.h"
#include "libGenome/gnFileSource.h"
#include "libGenome/gnFileContig.h"
#include "libGenome/gnABISource.h"
#include "libGenome/gnFASSource.h"
#include "libGenome/gnGBKSource.h"
#include "libGenome/gnRAWSource.h"
#include "libGenome/gnSEQSource.h"
#include "libGenome/gnDNXSource.h"
//filter and translate headers
#include "libGenome/gnBaseFilter.h"
#include "libGenome/gnFilter.h"
#include "libGenome/gnCompare.h"
#include "libGenome/gnTranslator.h"
#include "libGenome/gnFastTranslator.h"
#include "libGenome/gnPosSpecificTranslator.h"
//feature headers
#include "libGenome/gnBaseQualifier.h"
#include "libGenome/gnSourceQualifier.h"
#include "libGenome/gnStringQualifier.h"
#include "libGenome/gnLocation.h"
#include "libGenome/gnBaseHeader.h"
#include "libGenome/gnSourceHeader.h"
#include "libGenome/gnStringHeader.h"
#include "libGenome/gnBaseFeature.h"
#include "libGenome/gnFeature.h"
//spec headers
#include "libGenome/gnBaseSpec.h"
#include "libGenome/gnMultiSpec.h"
#include "libGenome/gnGenomeSpec.h"
#include "libGenome/gnFragmentSpec.h"
#include "libGenome/gnContigSpec.h"
#include "libGenome/gnSourceSpec.h"
#include "libGenome/gnStringSpec.h"

// Interval sequence tree
//#include "IntervalSequenceTree.h"

#endif  // _gnprec_h_
