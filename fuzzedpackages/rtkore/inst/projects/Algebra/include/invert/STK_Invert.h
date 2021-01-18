/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff, Université Lille 1, Inria

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  stkpp::Algebra
 * created on: 9 août 2016
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Invert.h
 *  @brief In this file we implement inversion method for general matrices.
 **/

#ifndef STK_INVERT_H
#define STK_INVERT_H

#ifdef STKUSELAPACK
#include <Algebra/include/STK_lapack_Svd.h>
#else
#include <Algebra/include/STK_Svd.h>
#endif

namespace STK
{

namespace hidden
{
/** @ingroup hidden
 *  @brief Implementation of the inversion matrix method for general/square matrices
 **/
template<class Matrix, int Size_>
struct InvertMatrixImpl
{
  typedef typename Matrix::Type Type;
  /** @ingroup hidden
   *  @brief compute the inverse of the 1x1 matrix and store the result in inv.
   *  @param m, inv the matrices to invert and its inverse
   *  @return @c true if m is invertible, @c false otherwise
   **/
  static Type invertMatrix11( Matrix const& m, CArraySquare<Type, Size_>& inv)
  {
    const int iShift = m.beginRows(), jShift = m.beginCols();
    inv.resize(TRange<Size_>(0, 1));
    // cofactor (0, jShift+0) [0]
    inv(0, 0) = Type(1);
    // compute determinant
    Type det = m(iShift+0, jShift+0);
    if (det == Type(0)) return Type(0);
    // compute inverse matrix
    inv /= det;
    return det;
  }
  /** @ingroup hidden
   *  @brief compute the inverse of the 2x2 matrix m and store the result in inv.
   *  @param m, inv the matrices to invert and its inverse
   *  @return @c true if m is invertible, @c false otherwise
   **/
  static Type invertMatrix22( Matrix const& m, CArraySquare<Type, Size_>& inv)
  {
    const int iShift = m.beginRows(), jShift = m.beginCols();
    inv.resize(TRange<Size_>(0, 2));
    // cofactor (0, jShift+0) [0]
    inv(0, 0) =   m(iShift+1, jShift+1);
    // cofactor (1, jShift+0) [1]
    inv(1, 0) = - m(iShift+1, jShift+0);
    // cofactor (1, jShift+1) [2]
    inv(0, 1) = - m(iShift+0, jShift+1);
    // cofactor (1, jShift+1) [3]
    inv(1, 1) =   m(iShift+0, jShift+0);

    // determinant
    Type det = m(iShift+0, jShift+0) * inv(0, 0)
             + m(iShift+0, jShift+1) * inv(1, 0);
    if (det == Type(0)) return Type(0);
    // compute inverse matrix
    inv /= det;
    return det;
  }
  /** @ingroup hidden
   *  @brief compute the inverse of the 3x3 matrix m and store the result in inv.
   *  @param m, inv the matrices to invert and its inverse
   *  @return @c true if m is invertible, @c false otherwise
   **/
  static Type invertMatrix33( Matrix const& m, CArraySquare<Type, Size_>& inv)
  {
    const int iShift = m.beginRows(), jShift = m.beginCols();
    inv.resize(TRange<Size_>(0, 3));
    // cofactor
    inv(0, 0) = m(iShift+1, jShift+1) * m(iShift+2, jShift+2)
              - m(iShift+2, jShift+1) * m(iShift+1, jShift+2);
    inv(0, 1) = m(iShift+0, jShift+2) * m(iShift+2, jShift+1)
              - m(iShift+0, jShift+1) * m(iShift+2, jShift+2);
    inv(0, 2) = m(iShift+0, jShift+1) * m(iShift+1, jShift+2)
              - m(iShift+0, jShift+2) * m(iShift+1, jShift+1);
    inv(1, 0) = m(iShift+1, jShift+2) * m(iShift+2, jShift+0)
              - m(iShift+1, jShift+0) * m(iShift+2, jShift+2);
    inv(1, 1) = m(iShift+0, jShift+0) * m(iShift+2, jShift+2)
              - m(iShift+0, jShift+2) * m(iShift+2, jShift+0);
    inv(1, 2) = m(iShift+1, jShift+0) * m(iShift+0, jShift+2)
              - m(iShift+0, jShift+0) * m(iShift+1, jShift+2);
    inv(2, 0) = m(iShift+1, jShift+0) * m(iShift+2, jShift+1)
              - m(iShift+2, jShift+0) * m(iShift+1, jShift+1);
    inv(2, 1) = m(iShift+2, jShift+0) * m(iShift+0, jShift+1)
              - m(iShift+0, jShift+0) * m(iShift+2, jShift+1);
    inv(2, 2) = m(iShift+0, jShift+0) * m(iShift+1, jShift+1)
              - m(iShift+1, jShift+0) * m(iShift+0, jShift+1);

    // computes the inverse of a matrix m
    Type det = m(iShift+0, jShift+0) * inv(0, 0)
             + m(iShift+0, jShift+1) * inv(1, 0)
             + m(iShift+0, jShift+2) * inv(2, 0);
    if (det == Type(0)) return Type(0);
    // compute inverse matrix
    inv /= det;
    return det;
  }
  /** @ingroup hidden
   *  @brief compute the inverse of the 4x4 matrix m and store the result in inv.
   *  @param m, inv the matrices to invert and its inverse
   *  @return @c true if m is invertible, @c false otherwise
   **/
  static Type invertMatrix44( Matrix const& m, CArraySquare<Type, Size_>& inv)
  {
    const int iShift = m.beginRows(), jShift = m.beginCols();
    inv.resize(TRange<Size_>(0, 4));
    // cofactor (0, 0) [0]
    inv(0, 0) = m(iShift+1, jShift+1) * m(iShift+2, jShift+2) * m(iShift+3, jShift+3)
              - m(iShift+1, jShift+1) * m(iShift+3, jShift+2) * m(iShift+2, jShift+3)
              - m(iShift+1, jShift+2) * m(iShift+2, jShift+1) * m(iShift+3, jShift+3)
              + m(iShift+1, jShift+2) * m(iShift+3, jShift+1) * m(iShift+2, jShift+3)
              + m(iShift+1, jShift+3) * m(iShift+2, jShift+1) * m(iShift+3, jShift+2)
              - m(iShift+1, jShift+3) * m(iShift+3, jShift+1) * m(iShift+2, jShift+2);
    // cofactor (0, 1) [4]
    inv(0, 1) = -m(iShift+0, jShift+1) * m(iShift+2, jShift+2) * m(iShift+3, jShift+3)
              +  m(iShift+0, jShift+1) * m(iShift+3, jShift+2) * m(iShift+2, jShift+3)
              +  m(iShift+0, jShift+2) * m(iShift+2, jShift+1) * m(iShift+3, jShift+3)
              -  m(iShift+0, jShift+2) * m(iShift+3, jShift+1) * m(iShift+2, jShift+3)
              -  m(iShift+0, jShift+3) * m(iShift+2, jShift+1) * m(iShift+3, jShift+2)
              +  m(iShift+0, jShift+3) * m(iShift+3, jShift+1) * m(iShift+2, jShift+2);
    // cofactor (0, 2) [8]
    inv(0, 2) = m(iShift+0, jShift+1) * m(iShift+1, jShift+2) * m(iShift+3, jShift+3)
              - m(iShift+0, jShift+1) * m(iShift+3, jShift+2) * m(iShift+1, jShift+3)
              - m(iShift+0, jShift+2) * m(iShift+1, jShift+1) * m(iShift+3, jShift+3)
              + m(iShift+0, jShift+2) * m(iShift+3, jShift+1) * m(iShift+1, jShift+3)
              + m(iShift+0, jShift+3) * m(iShift+1, jShift+1) * m(iShift+3, jShift+2)
              - m(iShift+0, jShift+3) * m(iShift+3, jShift+1) * m(iShift+1, jShift+2);
    // cofactor (0, 3) [12]
    inv(0, 3) = -m(iShift+0, jShift+1) * m(iShift+1, jShift+2) * m(iShift+2, jShift+3)
              +  m(iShift+0, jShift+1) * m(iShift+2, jShift+2) * m(iShift+1, jShift+3)
              +  m(iShift+0, jShift+2) * m(iShift+1, jShift+1) * m(iShift+2, jShift+3)
              -  m(iShift+0, jShift+2) * m(iShift+2, jShift+1) * m(iShift+1, jShift+3)
              -  m(iShift+0, jShift+3) * m(iShift+1, jShift+1) * m(iShift+2, jShift+2)
              +  m(iShift+0, jShift+3) * m(iShift+2, jShift+1) * m(iShift+1, jShift+2);
    // cofactor (1, 0) [1]
    inv(1, 0) = -m(iShift+1, jShift+0) * m(iShift+2, jShift+2) * m(iShift+3, jShift+3)
              +  m(iShift+1, jShift+0) * m(iShift+3, jShift+2) * m(iShift+2, jShift+3)
              +  m(iShift+1, jShift+2) * m(iShift+2, jShift+0) * m(iShift+3, jShift+3)
              -  m(iShift+1, jShift+2) * m(iShift+3, jShift+0) * m(iShift+2, jShift+3)
              -  m(iShift+1, jShift+3) * m(iShift+2, jShift+0) * m(iShift+3, jShift+2)
              +  m(iShift+1, jShift+3) * m(iShift+3, jShift+0) * m(iShift+2, jShift+2);
    // cofactor (1, 1) [5]
    inv(1, 1) = m(iShift+0, jShift+0) * m(iShift+2, jShift+2) * m(iShift+3, jShift+3)
              - m(iShift+0, jShift+0) * m(iShift+3, jShift+2) * m(iShift+2, jShift+3)
              - m(iShift+0, jShift+2) * m(iShift+2, jShift+0) * m(iShift+3, jShift+3)
              + m(iShift+0, jShift+2) * m(iShift+3, jShift+0) * m(iShift+2, jShift+3)
              + m(iShift+0, jShift+3) * m(iShift+2, jShift+0) * m(iShift+3, jShift+2)
              - m(iShift+0, jShift+3) * m(iShift+3, jShift+0) * m(iShift+2, jShift+2);
    // cofactor (1, 2) [9]
    inv(1, 2) = -m(iShift+0, jShift+0) * m(iShift+1, jShift+2) * m(iShift+3, jShift+3)
              +  m(iShift+0, jShift+0) * m(iShift+3, jShift+2) * m(iShift+1, jShift+3)
              +  m(iShift+0, jShift+2) * m(iShift+1, jShift+0) * m(iShift+3, jShift+3)
              -  m(iShift+0, jShift+2) * m(iShift+3, jShift+0) * m(iShift+1, jShift+3)
              -  m(iShift+0, jShift+3) * m(iShift+1, jShift+0) * m(iShift+3, jShift+2)
              +  m(iShift+0, jShift+3) * m(iShift+3, jShift+0) * m(iShift+1, jShift+2);
    // cofactor (1, 3) [13]
    inv(1, 3) = m(iShift+0, jShift+0) * m(iShift+1, jShift+2) * m(iShift+2, jShift+3)
              - m(iShift+0, jShift+0) * m(iShift+2, jShift+2) * m(iShift+1, jShift+3)
              - m(iShift+0, jShift+2) * m(iShift+1, jShift+0) * m(iShift+2, jShift+3)
              + m(iShift+0, jShift+2) * m(iShift+2, jShift+0) * m(iShift+1, jShift+3)
              + m(iShift+0, jShift+3) * m(iShift+1, jShift+0) * m(iShift+2, jShift+2)
              - m(iShift+0, jShift+3) * m(iShift+2, jShift+0) * m(iShift+1, jShift+2);
    // cofactor (2, 0) [2]
    inv(2, 0) = m(iShift+1, jShift+0) * m(iShift+2, jShift+1) * m(iShift+3, jShift+3)
              - m(iShift+1, jShift+0) * m(iShift+3, jShift+1) * m(iShift+2, jShift+3)
              - m(iShift+1, jShift+1) * m(iShift+2, jShift+0) * m(iShift+3, jShift+3)
              + m(iShift+1, jShift+1) * m(iShift+3, jShift+0) * m(iShift+2, jShift+3)
              + m(iShift+1, jShift+3) * m(iShift+2, jShift+0) * m(iShift+3, jShift+1)
              - m(iShift+1, jShift+3) * m(iShift+3, jShift+0) * m(iShift+2, jShift+1);
    // cofactor (2, 1)
    inv(2, 1) = -m(iShift+0, jShift+0) * m(iShift+2, jShift+1) * m(iShift+3, jShift+3)
              +  m(iShift+0, jShift+0) * m(iShift+3, jShift+1) * m(iShift+2, jShift+3)
              +  m(iShift+0, jShift+1) * m(iShift+2, jShift+0) * m(iShift+3, jShift+3)
              -  m(iShift+0, jShift+1) * m(iShift+3, jShift+0) * m(iShift+2, jShift+3)
              -  m(iShift+0, jShift+3) * m(iShift+2, jShift+0) * m(iShift+3, jShift+1)
              +  m(iShift+0, jShift+3) * m(iShift+3, jShift+0) * m(iShift+2, jShift+1);
    // cofactor (2, 2)
    inv(2, 2) = m(iShift+0, jShift+0) * m(iShift+1, jShift+1) * m(iShift+3, jShift+3)
              - m(iShift+0, jShift+0) * m(iShift+3, jShift+1) * m(iShift+1, jShift+3)
              - m(iShift+0, jShift+1) * m(iShift+1, jShift+0) * m(iShift+3, jShift+3)
              + m(iShift+0, jShift+1) * m(iShift+3, jShift+0) * m(iShift+1, jShift+3)
              + m(iShift+0, jShift+3) * m(iShift+1, jShift+0) * m(iShift+3, jShift+1)
              - m(iShift+0, jShift+3) * m(iShift+3, jShift+0) * m(iShift+1, jShift+1);
    // cofactor (2, 3)
    inv(2, 3) = -m(iShift+0, jShift+0) * m(iShift+1, jShift+1) * m(iShift+2, jShift+3)
              +  m(iShift+0, jShift+0) * m(iShift+2, jShift+1) * m(iShift+1, jShift+3)
              +  m(iShift+0, jShift+1) * m(iShift+1, jShift+0) * m(iShift+2, jShift+3)
              -  m(iShift+0, jShift+1) * m(iShift+2, jShift+0) * m(iShift+1, jShift+3)
              -  m(iShift+0, jShift+3) * m(iShift+1, jShift+0) * m(iShift+2, jShift+1)
              +  m(iShift+0, jShift+3) * m(iShift+2, jShift+0) * m(iShift+1, jShift+1);
    // cofactor (3, 0)
    inv(3, 0) = -m(iShift+1, jShift+0) * m(iShift+2, jShift+1) * m(iShift+3, jShift+2)
              +  m(iShift+1, jShift+0) * m(iShift+3, jShift+1) * m(iShift+2, jShift+2)
              +  m(iShift+1, jShift+1) * m(iShift+2, jShift+0) * m(iShift+3, jShift+2)
              -  m(iShift+1, jShift+1) * m(iShift+3, jShift+0) * m(iShift+2, jShift+2)
              -  m(iShift+1, jShift+2) * m(iShift+2, jShift+0) * m(iShift+3, jShift+1)
              +  m(iShift+1, jShift+2) * m(iShift+3, jShift+0) * m(iShift+2, jShift+1);
    // cofactor (3, 1) [7]
    inv(3, 1) = m(iShift+0, jShift+0) * m(iShift+2, jShift+1) * m(iShift+3, jShift+2)
              - m(iShift+0, jShift+0) * m(iShift+3, jShift+1) * m(iShift+2, jShift+2)
              - m(iShift+0, jShift+1) * m(iShift+2, jShift+0) * m(iShift+3, jShift+2)
              + m(iShift+0, jShift+1) * m(iShift+3, jShift+0) * m(iShift+2, jShift+2)
              + m(iShift+0, jShift+2) * m(iShift+2, jShift+0) * m(iShift+3, jShift+1)
              - m(iShift+0, jShift+2) * m(iShift+3, jShift+0) * m(iShift+2, jShift+1);
    // cofactor (3, 2) [11]
    inv(3, 2) = -m(iShift+0, jShift+0) * m(iShift+1, jShift+1) * m(iShift+3, jShift+2)
              +  m(iShift+0, jShift+0) * m(iShift+3, jShift+1) * m(iShift+1, jShift+2)
              +  m(iShift+0, jShift+1) * m(iShift+1, jShift+0) * m(iShift+3, jShift+2)
              -  m(iShift+0, jShift+1) * m(iShift+3, jShift+0) * m(iShift+1, jShift+2)
              -  m(iShift+0, jShift+2) * m(iShift+1, jShift+0) * m(iShift+3, jShift+1)
              +  m(iShift+0, jShift+2) * m(iShift+3, jShift+0) * m(iShift+1, jShift+1);
    // cofactor (3, 3) [15]
    inv(3, 3) = m(iShift+0, jShift+0) * m(iShift+1, jShift+1) * m(iShift+2, jShift+2)
              - m(iShift+0, jShift+0) * m(iShift+2, jShift+1) * m(iShift+1, jShift+2)
              - m(iShift+0, jShift+1) * m(iShift+1, jShift+0) * m(iShift+2, jShift+2)
              + m(iShift+0, jShift+1) * m(iShift+2, jShift+0) * m(iShift+1, jShift+2)
              + m(iShift+0, jShift+2) * m(iShift+1, jShift+0) * m(iShift+2, jShift+1)
              - m(iShift+0, jShift+2) * m(iShift+2, jShift+0) * m(iShift+1, jShift+1);
    // compute determinant and inverse matrix
    Type det = m(iShift+0, jShift+0) * inv(0, 0)
             + m(iShift+1, jShift+0) * inv(0, 1)
             + m(iShift+2, jShift+0) * inv(0, 2)
             + m(iShift+3, jShift+0) * inv(0, 3);
    if (det == Type(0)) return Type(0);
    // compute inverse matrix
    inv /= det;
    return det;
  }

  /** @ingroup hidden
   *  @brief compute the inverse of the matrix m and store the result in inv.
   *  @note if the matrix is not invertible, the result will be a generalized inverse.
   *  @param m, inv the matrices to invert and its inverse
   *  @return The determinant value of m
   **/
  static Type invertMatrixXX( Matrix const& m, CArraySquare<Type, Size_>& inv)
  {
    inv.resize(TRange<Size_>(0, m.sizeRows()));
#ifdef STKUSELAPACK
      lapack::Svd decomp(m);
#else
      Svd<Matrix> decomp(m);
#endif
    if (!decomp.run()) return Type(0);
    // compute (generalized) inverse matrix
    decomp.ginv(inv);
    return decomp.det();
  }
};


/** @ingroup hidden utility class allowing to call the correct static function
 *  computing the inverse of a matrix. */
template<class Matrix, int Size_>
struct InvertMatrixDispatcher
{
  typedef typename Matrix::Type Type;
  inline static Type run( Matrix const& m, CArraySquare<Type, Size_>& inv)
  {
    switch (m.sizeRows())
    {
      case 1: return InvertMatrixImpl<Matrix, Size_>::invertMatrix11(m, inv);
      case 2: return InvertMatrixImpl<Matrix, Size_>::invertMatrix22(m, inv);
      case 3: return InvertMatrixImpl<Matrix, Size_>::invertMatrix33(m, inv);
      case 4: return InvertMatrixImpl<Matrix, Size_>::invertMatrix44(m, inv);
      default:return InvertMatrixImpl<Matrix, Size_>::invertMatrixXX(m, inv);
    }
  }
};


template<class Matrix>
struct InvertMatrixDispatcher<Matrix, 1>
{
  typedef typename Matrix::Type Type;
  inline static Type run( Matrix const& m, CArraySquare<Type, 1>& inv)
  { return InvertMatrixImpl<Matrix, 1>::invertMatrix11(m, inv);}
};

template<class Matrix>
struct InvertMatrixDispatcher<Matrix, 2>
{
  typedef typename Matrix::Type Type;
  inline static Type run( Matrix const& m, CArraySquare<Type, 2>& inv)
  { return InvertMatrixImpl<Matrix, 2>::invertMatrix22(m, inv);}
};

template<class Matrix>
struct InvertMatrixDispatcher<Matrix, 3>
{
  typedef typename Matrix::Type Type;
  inline static Type run( Matrix const& m, CArraySquare<Type, 3>& inv)
  { return InvertMatrixImpl<Matrix, 3>::invertMatrix33(m, inv);}
};

template<class Matrix>
struct InvertMatrixDispatcher<Matrix, 4>
{
  typedef typename Matrix::Type Type;
  inline static Type run( Matrix const& m, CArraySquare<Type, 4>& inv)
  { return InvertMatrixImpl<Matrix, 4>::invertMatrix44(m, inv);}
};

} // namespace hidden

} // namespace STK

#endif /* STK_INVERT_H */
