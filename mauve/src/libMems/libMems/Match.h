/*******************************************************************************
 * $Id: GenericMatch.h,v 1.10 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _Match_h_
#define _Match_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnClone.h"
#include <iostream>
#include <set>
#include "libMems/UngappedLocalAlignment.h"
#include "libMems/SparseAbstractMatch.h"
#include "libMems/DenseAbstractMatch.h"
#include "libMems/HybridAbstractMatch.h"

namespace mems {

typedef UngappedLocalAlignment< HybridAbstractMatch<> > Match;

static uint seq_compare_start;


}

#endif // _Match_h_
