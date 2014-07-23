/*
 *    This file is part of HMMoC 1.0, a hidden Markov model compiler.
 *    Copyright (C) 2006 by Gerton Lunter, Oxford University.
 *
 *    HMMoC is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    HMMOC is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with HMMoC; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
\*/
/*
 *
 * DPTable.h -- ordinary and sparse dynamic programming tables
 *
 * Gerton Lunter, 1 Oct 2006
 *
 * Modified for GCC 4.0.2 and Microsoft Visual Studio 2005 by Aaron Darling, 2007
 *
 */

#ifndef __dptable_h_
#define __dptable_h_


#include <map>
#include <cassert>


#ifdef __GNUC__
 #define HAVE_HASH_MAP
 #if __GNUC__ < 3
  #include <hash_map.h>
  namespace Sgi { using ::hash_map; }; // inherit globals
 #else
  #include <ext/hash_map>
  #if __GNUC_MINOR__ + __GNUC__ == 3
   namespace Sgi = std;               // GCC 3.0
  #else
   namespace Sgi = ::__gnu_cxx;       // GCC 3.1 and later
  #endif
 #endif
#else      // ...there are other compilers, right?
#ifdef _MSC_VER
// visual studio 2005 has no hash map.  older versions did.
#else
// default for all other compilers
#define HAVE_HASH_MAP
namespace Sgi = std;
#endif
#endif


using std::map;
#ifdef HAVE_HASH_MAP
using Sgi::hash_map;
#endif

// Define aliases for two maps: red-black trees, and hashes
// (GNU C++ does not define a hash function for long long int, so we have to define our own)

template<class key>
struct _hash {
  size_t operator()(long long x) const { return x; }
};

// typedefs can't be templated

template<class key, class value>
class treemap : public map<key,value> {};

#ifdef HAVE_HASH_MAP
template<class key, class value>
class hashmap : public hash_map<key,value,_hash<key> > {};
#endif

// select one of the maps (the hash map is faster and appears to use less memory)

#ifdef HAVE_HASH_MAP
#define _mymap hashmap
#else
#define _mymap treemap
#endif



// States are stored in a self-initializing array

template<class Real,int size> class States;

template<int size>
class States<double,size> {
private:
  double data[size];
public:
  enum { length = size };                                      // to know the size, just in case
  States() { for (int i=0; i<size; i++) data[i]=0; }           // initialization
  operator double* () { return data; }                         // cast to actual array
  operator const double* () const { return data; }
};

template<class Real, int size>
class States {
private:
  Real data[size];
public:
  enum { length = size };                                      // to know the size, just in case
  States() { for (int i=0; i<size; i++) data[i].clear(); }      // initialization
  operator Real* () { return data; }                           // cast to actual array
  operator const Real* () const { return data; }
};



// Define index types to serve as keys to the DP table position

template<int dim> class _index {};
template<> class _index<1> { public: typedef unsigned int t; };
template<> class _index<2> { public: typedef unsigned long t; };
template<> class _index<3> { public: typedef unsigned long long t; };
template<> class _index<4> { public: typedef unsigned long long t; };


//
// Base classes for a dynamic programming table
//
// DP tables provide the following methods:
//
// const States& read(...)   :   read access to state array
// States& write(...)        :   write access to state array
// void written()            :   signal that write access is finished
// void allocate(...)        :   inform table about its dimensions; must be called before read/write access
// void clear()              :   empties the table; keeps its dimensions
// void clear(int)           :   empties one column of a (folded) DP table
// void absolve()            :   ensures that table does not delete its data; another table with reference to the data,
//                               which is created by the default copy constructor, is now responsible.  Not allowed for
//                               folded tables
//

template<class States>
class _DPT {                            // base class, keeps track of the responsibility for the data
 private:
  void clear(int i) {}                  // placeholder, to allow dummy definition for non-folded tables
 protected:
  bool isInCharge;                      // true if this class' destructor destroys the data
 public:
  typedef States states_type;
  _DPT() : isInCharge(true) {}
  void absolve() { isInCharge=false; }  // take away the responsibility of destroying the data
  void written() {}                     // signal that we're done writing -- used by extensions
};

template<template<typename,int> class DPTable, class States, int dim>     // Wrapper for memory-efficient Fw/Bw/Baum-Welch
class _FoldedTable : public _DPT<States> {
 protected:
  DPTable<States,dim-1>* aTables[2];
 public:
  _FoldedTable() { aTables[0] = new DPTable<States,dim-1>(); aTables[1] = new DPTable<States,dim-1>(); }
  ~_FoldedTable() { assert(_DPT<States>::isInCharge); delete aTables[0]; delete aTables[1]; }        // do not allow data to be retained
  void clear(int i) { aTables[i%2]->clear(); }
};
  
template<class States, int dim>
class DPTable {};

template<class States, int dim>
class SparseDPTable {};

template<template<typename,int> class DPTable, class States, int dim>
class FoldedTable {};



// Explicit partial specializations for up to 4 spatial dimensions


template<template<typename,int> class DPTable, class States>
class FoldedTable<DPTable, States, 0> : public _FoldedTable<DPTable, States, 1> {
 public:
  void allocate() { this->aTables[0]->allocate(); };
  const States& read() const { return this->aTables[0]->read(); }
  States& write() { return this->aTables[0]->write(); }
  void written() { this->aTables[0]->written(); }
};


template<template<typename,int> class DPTable, class States>
class FoldedTable<DPTable, States, 1> : public _FoldedTable<DPTable, States, 1> {
  int z;
 public:
  void allocate(int a) { this->aTables[0]->allocate(); this->aTables[1]->allocate(); };
  const States& read(int a) const { return this->aTables[a%2]->read(); }
  States& write(int a) { return this->aTables[z=a%2]->write(); }
  void written() { this->aTables[z]->written(); }
};


template<template<typename,int> class DPTable, class States>
class FoldedTable<DPTable, States, 2> : public _FoldedTable<DPTable, States, 2> {
  int z;
 public:
  void allocate(int a, int b) { this->aTables[0]->allocate(a); this->aTables[1]->allocate(a); };
  const States& read(int a, int b) const { return this->aTables[b%2]->read(a); }
  States& write(int a, int b) { return this->aTables[z=b%2]->write(a); }
  void written() { this->aTables[z]->written(); }
};


template<template<typename,int> class DPTable, class States>
class FoldedTable<DPTable, States, 3> : public _FoldedTable<DPTable, States, 3> {
  int z;
 public:
  void allocate(int a, int b, int c) { this->aTables[0]->allocate(a,b); this->aTables[1]->allocate(a,b); };
  const States& read(int a, int b, int c) const { return this->aTables[c%2]->read(a,b); }
  States& write(int a, int b, int c) { return this->aTables[z=c%2]->write(a,b); }
  void written() { this->aTables[z]->written(); }
};


template<template<typename,int> class DPTable, class States>
class FoldedTable<DPTable, States, 4> : public _FoldedTable<DPTable, States, 4> {
  int z;
 public:
  void allocate(int a, int b, int c, int d) { this->aTables[0]->allocate(a,b,c); this->aTables[1]->allocate(a,b,c); };
  const States& read(int a, int b, int c, int d) const { return this->aTables[d%2]->read(a,b,c); }
  States& write(int a, int b, int c, int d) { return this->aTables[z=d%2]->write(a,b,c); }
  void written() { this->aTables[z]->written(); }
};


template<class States>
class DPTable<States,0> : public _DPT<States> {
private:
  States *pTable;
public:
  DPTable() { pTable = 0; };
  ~DPTable() { if (pTable && _DPT<States>::isInCharge) delete pTable; };
  void allocate() { pTable = new States(); };
  void clear() { delete pTable; allocate(); };
  const States& read() const { return *pTable; }
  States& write() { return *pTable; }
};


template<class States>
class DPTable<States,1> : public _DPT<States> {
private:
  States *pTable;
  int maxa;
public:
  DPTable() { pTable = 0; }
  ~DPTable() { if (pTable && _DPT<States>::isInCharge ) { delete[] pTable; } }
  void allocate(int a) { maxa = a; pTable = new States[a]; }
  void clear() { delete[] pTable; allocate(maxa); };
  const States& read(int a) const { return pTable[a]; }
  States& write(int a) { return pTable[a]; }
};


template<class States>
class DPTable<States,2> : public _DPT<States> {
private:
  States *pTable;
  int maxa, maxb;
public:
  DPTable() { pTable = 0; }
  ~DPTable() { if (pTable && _DPT<States>::isInCharge ) { delete[] pTable; } }
  void allocate(int a, int b) { maxa = a; maxb = b; pTable = new States[a*b]; }
  void clear() { delete[] pTable; allocate(maxa,maxb); };
  const States& read(int a, int b) const { return pTable[a+maxa*b]; }
  States& write(int a, int b) { return pTable[a+maxa*b]; }
};


template<class States>
class DPTable<States,3> : public _DPT<States> {
private:
  States *pTable;
  int maxa, maxb, maxc;
public:
  DPTable() { pTable = 0; }
  ~DPTable() { if (pTable && _DPT<States>::isInCharge ) { delete[] pTable; } }
  void allocate(int a, int b, int c) { maxa = a; maxb = b; maxc = c; pTable = new States[a*b*c]; }
  void clear() { delete[] pTable; allocate(maxa,maxb,maxc); };
  const States& read(int a, int b, int c) const { return pTable[a+maxa*(b+maxb*c)]; }
  States& write(int a, int b, int c) { return pTable[a+maxa*(b+maxb*c)]; }
};


template<class States>
class SparseDPTable<States,0> : public DPTable<States,0> {};


template<class States>
class SparseDPTable<States,1> : public _DPT<States> {
private:
  typedef _index<1>::t idx;
  _mymap<idx,States> &table;
  const States zero;
public:
  SparseDPTable() : table(*new _mymap<_index<1>::t,States>), zero() {};
  ~SparseDPTable() { if (_DPT<States>::isInCharge) delete &table; };
  void allocate(int a) {};
  void clear() { table.clear(); }
  States& write(int a) { return table[idx(a)]; }
  const States& read(int a) const {
    _mymap<_index<1>::t,char>::iterator iter2;
    typename _mymap<idx,States>::const_iterator iter = table.find(idx(a));
    if (iter == table.end()) return zero;
    return iter->second;
  }
};


template<class States>
class SparseDPTable<States,2> : public _DPT<States> {
private:
  typedef _index<2>::t idx;
  _mymap<idx,States> &table;
  idx maxa;
  const States zero;
public:
  SparseDPTable() : table(*new _mymap<idx,States>), zero() {};
  ~SparseDPTable() { if (_DPT<States>::isInCharge) delete &table; };
  void allocate(int a, int b) { maxa = a; };
  void clear() { table.clear(); }
  States& write(int a, int b) { return table[idx(a)+maxa*idx(b)]; }
  const States& read(int a, int b) const {
    typename _mymap<idx,States>::const_iterator iter = table.find(unsigned(a)+maxa*unsigned(b));
    if (iter == table.end()) return zero;
    return iter->second;
  }
};


template<class States>
class SparseDPTable<States,3> : public _DPT<States> {

private:
  typedef _index<3>::t idx;
  _mymap<idx,States> &table;
  idx maxa, maxb;
  const States zero;
public:
  SparseDPTable() : table(*new _mymap<idx,States>), zero() {};
  ~SparseDPTable() { if (_DPT<States>::isInCharge) delete &table; };
  void allocate(int a, int b, int c) { maxa = a; maxb = b; };
  void clear() { table.clear(); }
  States& write(int a, int b, int c) { return table[idx(a)+maxa*(idx(b)+maxb*idx(c))]; }
  const States& read(int a, int b, int c) const {
    typename _mymap<idx,States>::const_iterator iter = table.find(unsigned(a)+maxa*(unsigned(b)+maxb*unsigned(c)));
    if (iter == table.end()) return zero;
    return iter->second;
  }
};


template<class States>
class SparseDPTable<States,4> : public _DPT<States> {
private:
  typedef _index<4>::t idx;
  _mymap<idx,States> &table;
  idx maxa, maxb, maxc;
  const States zero;
public:
  SparseDPTable() : table(*new _mymap<idx,States>), zero() {};
  ~SparseDPTable() { if (_DPT<States>::isInCharge) delete &table; };
  void allocate(int a, int b, int c, int d) { maxa=a; maxb=b; maxc=c; }
  void clear() { table.clear(); }
  States& write(int a, int b, int c, int d) { return table[unsigned(a)+maxa*(unsigned(b)+maxb*(unsigned(c)+maxc*unsigned(d)))]; }
  const States& read(int a, int b, int c, int d) const {
    typename _mymap<idx,States>::const_iterator iter = table.find(unsigned(a)+maxa*(unsigned(b)+maxb*(unsigned(c)+maxc*unsigned(d))));
    if (iter == table.end()) return zero;
    return iter->second;
  }
};


#endif

