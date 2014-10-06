/*******************************************************************************
 * $Id: NumericMatrix.h,v 1.4 2004/03/01 02:40:08 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef _NumericMatrix_h_
#define _NumericMatrix_h_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libMems/Matrix.h"

template<class T>  // See section on templates for more
class NumericMatrix : public Matrix<T>
{
public:
	NumericMatrix(){};
	NumericMatrix(unsigned nrows, unsigned ncols);
   
      // Based on the Law Of The Big Three:
	~NumericMatrix();
	NumericMatrix(const NumericMatrix<T>& m);
	NumericMatrix<T>& operator= (const NumericMatrix<T>& m);
	
	// define some arithmetic operators 
	NumericMatrix<T>& operator+= (const NumericMatrix<T>& m);
	NumericMatrix<T>& operator-= (const NumericMatrix<T>& m);
	// not implemented
	NumericMatrix<T>& operator*= (const NumericMatrix<T>& m);
	NumericMatrix<T>& operator*= (const T& m);
	NumericMatrix<T>& operator/= (const NumericMatrix<T>& m);
	NumericMatrix<T>& operator/= (const T& m);

	// the following 5 are not implemented
	NumericMatrix<T>& operator+ (const NumericMatrix<T>& m ) const;
	const NumericMatrix<T>& operator- (const NumericMatrix<T>& m ) const;
	const NumericMatrix<T>& operator* (const NumericMatrix<T>& m ) const;
	const NumericMatrix<T>& operator* (const T& n) const;
	const NumericMatrix<T>& operator/ (const T& n) const;

};
   
template<class T>
inline NumericMatrix<T>::NumericMatrix(unsigned nrows, unsigned ncols)
	: Matrix<T>( nrows, ncols )
{
}
   
template<class T>
inline NumericMatrix<T>::NumericMatrix(const NumericMatrix<T>& m){
	*this = m;
}

template<class T>
inline NumericMatrix<T>& NumericMatrix<T>::operator= (const NumericMatrix<T>& m)
{
	Matrix<T>::operator=( m );
	return *this;
}

template<class T>
inline NumericMatrix<T>::~NumericMatrix()
{
}

template<class T>
inline
NumericMatrix<T>& NumericMatrix<T>::operator+= (const NumericMatrix<T>& m){
	// make sure matrix dimensions agree
	if (this->nrows_ != m.nrows_ || this->ncols_ != m.ncols_)
		throw typename Matrix<T>::BadSize();

	// do the arithmetic on each matrix entry
	for(unsigned i = 0; i < Matrix<T>::nrows_ * Matrix<T>::ncols_; i++ )
		this->data_[ i ] += m.data_[ i ];
	return *this;
}

template<class T>
inline
NumericMatrix<T>& NumericMatrix<T>::operator-= (const NumericMatrix<T>& m){
	// make sure matrix dimensions agree
	if (this->nrows_ != m.nrows_ || this->ncols_ != m.ncols_)
		throw typename Matrix<T>::BadSize();

	// do the arithmetic on each matrix entry
	for(unsigned i = 0; i < Matrix<T>::nrows_ * Matrix<T>::ncols_; i++ )
		this->data_[ i ] -= m.data_[ i ];
	return *this;
}

template<class T>
inline
NumericMatrix<T>& NumericMatrix<T>::operator*= (const NumericMatrix<T>& m){
	// make sure matrix dimensions agree
	if (this->ncols_ != m.nrows_)
		throw typename Matrix<T>::BadSize();
	// do a matrix multiply
	return *this;
}

template<class T>
inline
NumericMatrix<T>& NumericMatrix<T>::operator*= (const T& m){
	// do the arithmetic on each matrix entry
	for(unsigned i = 0; i < Matrix<T>::nrows_ * Matrix<T>::ncols_; i++ )
		this->data_[ i ] *= m;
	return *this;
}

template<class T>
inline
NumericMatrix<T>& NumericMatrix<T>::operator/= (const T& m){
	// do the arithmetic on each matrix entry
	for(unsigned i = 0; i < Matrix<T>::nrows_ * Matrix<T>::ncols_; i++ )
		this->data_[ i ] /= m;
	return *this;
}

template<class T>
inline
NumericMatrix<T>& NumericMatrix<T>::operator/= ( const NumericMatrix<T>& m ){
	// make sure matrix dimensions agree
	if (this->nrows_ != m.nrows_ || this->ncols_ != m.ncols_)
		throw typename Matrix<T>::BadSize();
	// do the arithmetic on each matrix entry
	for(unsigned i = 0; i < Matrix<T>::nrows_ * Matrix<T>::ncols_; i++ )
		this->data_[ i ] /= m.data_[ i ];
	return *this;
}

template<class T>
inline
NumericMatrix<T>& NumericMatrix<T>::operator+ (const NumericMatrix<T>& m) const {

}
template<class T>
inline
const NumericMatrix<T>& NumericMatrix<T>::operator- (const NumericMatrix<T>& m) const {

}
template<class T>
inline
const NumericMatrix<T>& NumericMatrix<T>::operator* (const NumericMatrix<T>& m) const {

}
template<class T>
inline
const NumericMatrix<T>& NumericMatrix<T>::operator* (const T& n) const {

}
template<class T>
inline
const NumericMatrix<T>& NumericMatrix<T>::operator/ (const T& n) const {

}


#endif // _NumericMatrix_h_
