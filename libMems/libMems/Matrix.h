/*******************************************************************************
 * $Id: Matrix.h,v 1.6 2004/02/27 23:08:55 darling Exp $
 * This file is copyright 2002-2007 Aaron Darling and authors listed in the AUTHORS file.
 * This file is licensed under the GPL.
 * Please see the file called COPYING for licensing details.
 * **************
 ******************************************************************************/

#ifndef __Matrix_h__
#define __Matrix_h__

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "libGenome/gnSetup.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <stdexcept>

template<class T>
class Matrix
{
public:
	Matrix();
	Matrix(unsigned nrows, unsigned ncols);
	// Throws a BadSize object if either size is zero
	class BadSize : public std::range_error{
	public:
		BadSize() : std::range_error( "Bad matrix size" ){}
	};

	// Based on the Law Of The Big Three:
	~Matrix();	
	Matrix(const Matrix<T>& m);
	Matrix<T>& operator= (const Matrix<T>& m);
   	// Access methods to get the (i,j) element:	
	T& operator() (unsigned i, unsigned j);
	const T& operator() (unsigned i, unsigned j) const;
	// These throw a BoundsViolation object if i or j is too big
	class BoundsViolation : public std::range_error{
	public:
		BoundsViolation() : std::range_error( "Index out of bounds" ){}
	 };
	// Support for initializing each matrix element to a value
	void init( const T& init_val );
	
	void print( std::ostream& os ) const;
	void read( std::istream& is );

	unsigned rows() const;
	unsigned cols() const;
protected:
	T* data_;
	unsigned nrows_, ncols_;
};
   
template<class T>
inline Matrix<T>::Matrix()
{
	data_ = NULL;
	nrows_ = 0;
	ncols_ = 0;
}

template<class T>
inline unsigned Matrix<T>::rows() const
{
	return nrows_;
}

template<class T>
inline unsigned Matrix<T>::cols() const
{
	return ncols_;
}

template<class T>
inline T& Matrix<T>::operator() (unsigned row, unsigned col)
{
	if (row >= nrows_ || col >= ncols_) 
		throw BoundsViolation();
	return data_[row*ncols_ + col];
}
   
template<class T>
inline const T& Matrix<T>::operator() (unsigned row, unsigned col) const
{
	if (row >= nrows_ || col >= ncols_) {
		std::cout << "debug me ";
		throw BoundsViolation();
	}
	return data_[row*ncols_ + col];
}
   
template<class T>
inline Matrix<T>::Matrix(unsigned nrows, unsigned ncols)
	: data_  (new T[nrows * ncols]),
	  nrows_ (nrows),
	  ncols_ (ncols)
{
}
template<class T>
inline Matrix<T>::Matrix(const Matrix<T>& m){
	*this = m;
}

template<class T>
inline Matrix<T>& Matrix<T>::operator=( const Matrix<T>& m )
{
	if( data_ != NULL )
		delete[] data_;
	data_ = new T[m.nrows_ * m.ncols_];
	nrows_ = m.nrows_;
	ncols_ = m.ncols_;
	memcpy( data_, m.data_, nrows_ * ncols_ * sizeof( T ) );
	return *this;
}

template<class T>
inline Matrix<T>::~Matrix()
{
	if( data_ != NULL )
		delete[] data_;
}

template<class T>
inline void Matrix<T>::init( const T& init_val )
{
	for( unsigned rowI = 0; rowI < nrows_; rowI++ )
		for( unsigned colI = 0; colI < ncols_; colI++ )
			data_[ rowI * ncols_ + colI ] = init_val;
}

template<class T>
inline void Matrix<T>::print( std::ostream& os ) const{
	for( unsigned rowI = 0; rowI < nrows_; rowI++ ){
		for( unsigned colI = 0; colI < ncols_; colI++ ){
			if( colI > 0 )
				os << '\t';
			os << data_[ rowI * ncols_ + colI ];
		}
		os << std::endl;
	}
}

template<class T>
inline void Matrix<T>::read( std::istream& is ){
	std::vector< std::string > lines;
	std::string cur_line;
	while( std::getline( is, cur_line ) )
		lines.push_back( cur_line );
		
	nrows_ = lines.size();
	// count ncols
	std::stringstream ss( lines[0] );
	ncols_ = 0;
	while( std::getline( ss, cur_line, '\t' ) )
		ncols_++;

	data_ = new T[nrows_ * ncols_];
	
	int valueI = 0;
	for( int lineI = 0; lineI < lines.size(); lineI++ ){
		ss = std::stringstream( lines[ lineI ] );
		std::getline( ss, cur_line, '\t' );
		std::stringstream type_stream( cur_line );
		type_stream >> data_[ valueI ];
		valueI++;
	}
}

#endif // __Matrix_h__
