/* $Id: matrix_wrapper.h,v 1.1 2008/01/25 11:47:49 ruthe Exp $ */

#ifndef _MATRIX_WRAPPER_H
#define _MATRIX_WRAPPER_H

// Following memberfunctions are required:
//
//   matrix()
//   matrix(int,int)
//   matrix(const matrix&)
//   operator=(const Matrix&M)
//   double& operator()(int,int) const
//	 double at(int r,int c) const
//   int rows() const
//   int columns() const
//   double det() const
//   matrix trans() const
//   matrix inv() const
//   matrix operator*(const matrix&) const
//   matrix operator*(double) const
//   valarray<double> to_vector() const

#include <valarray>
#include "matrix.h"

using namespace std;

class matrix;

// These function are implemented in matrix_wrapper.c++
matrix cof(const matrix& M);
double cof(const matrix& M,int r,int c);
valarray<double> to_vector(const matrix& M);
matrix to_matrix(const valarray<double>& v);


class matrix : public Matrix::matrix<double>
{
  public:

	matrix() : Matrix::matrix<double>()
		{}	
	matrix(int r,int c) : Matrix::matrix<double>(r,c)
		{}
	matrix(const matrix&M) : Matrix::matrix<double>(M)
		{}
	matrix(const Matrix::matrix<double>& M): Matrix::matrix<double>(M)
		{}
	
/* 	matrix& operator=(const matrix&M) */
/* 	{return Matrix::matrix<double>::operator=(M);} */
		
	inline int rows() const
		{return this->Matrix::matrix<double>::RowNo();}
	inline int columns() const
		{return this->Matrix::matrix<double>::ColNo();}
	inline double at(int r, int c) const
		{Matrix::matrix<double>* p=(Matrix::matrix<double>*)(this);
		return p->Matrix::matrix<double>::operator()(r,c);}
	
	inline matrix operator*(const matrix&M)
		{Matrix::matrix<double> ret;
		ret=*this;
		ret*=M;
		return ret;}

	inline matrix operator*(double t) const
		{Matrix::matrix<double> ret;
		ret=*this;
		ret*=t;
		return ret;}		

	inline double det() const
		{return this->Det();}
	inline matrix trans() const
		{return static_cast<matrix>(Matrix::operator~(*this));}
	inline matrix inv() const
		{return (1.0/det()) * this->Adj();}
	inline valarray<double> solve(const valarray<double>& v) const
		{return to_vector(this->Solve(to_matrix(v)));}
};

#endif // #ifndef _MATRIX_WRAPPER_H
