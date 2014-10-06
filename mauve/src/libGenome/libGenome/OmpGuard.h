/*******************************************************************************
 * $Id: OmpGuard.h,v 1.6 2004/02/27 23:08:55 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _OmpGuard_h_
#define _OmpGuard_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif


namespace genome {

#ifdef _OPENMP
/** This is a class for guard objects using OpenMP
*  It is adapted from the book
*  "Pattern-Oriented Software Architecture". */
class omp_guard {
public:
    /** Acquire the lock and store a pointer to it */
	omp_guard (omp_lock_t &lock){
	  lock_ = &lock;
      acquire ();
	}

	/** Set the lock explicitly */
	void acquire (){
		omp_set_lock (lock_);
		owner_ = true;
	}

	/** Release the lock explicitly (owner thread only!) */
	void release (){
		if (owner_) {
			owner_ = false;
			omp_unset_lock (lock_);
		};
	}

	/** Destruct guard object */
	~omp_guard (){ release (); }
 
private:
    omp_lock_t *lock_;  // pointer to our lock
    bool owner_;   // is this object the owner of the lock?
   
    // Disallow copies or assignment
    omp_guard (const omp_guard &);
    void operator= (const omp_guard &);
};

#else

// provide non-openmp stubs
typedef int omp_lock_t;
class omp_guard {
public:
	omp_guard (omp_lock_t &lock){}
	void acquire (){}
	void release (){}
 
private:
   
    // Disallow copies or assignment
    omp_guard (const omp_guard &);
    void operator= (const omp_guard &);
};

inline void omp_init_lock( omp_lock_t* dummy_lock ){}
inline void omp_destroy_lock( omp_lock_t* dummy_lock ){}

#endif	// ifdef _OPENMP


} // namespace genome


#endif	// _OmpGuard_h_

