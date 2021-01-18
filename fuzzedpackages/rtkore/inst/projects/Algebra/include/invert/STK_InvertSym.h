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

/** @file STK_InvertSym.h
 *  @brief In this file we implement inversion method for symmetric matrices.
 **/

#ifndef STK_SYMINVERT_H
#define STK_SYMINVERT_H

#ifdef STKUSELAPACK
#include <Algebra/include/STK_lapack_SymEigen.h>
#else
#include <Algebra/include/STK_SymEigen.h>
#endif

namespace STK
{

namespace hidden
{
/** @ingroup hidden
 *  @brief Implementation of the computation of the inverse of a lower triangular matrix
 **/
template<class Matrix, int Size_>
struct InvertSymImpl
{
  typedef typename Matrix::Type Type;
  /** @ingroup hidden
   *  @brief compute the inverse of the symmetric matrix m of size 1x1 and store
   *  the result in inv.
   *  @param m, inv the matrices to invert and its inverse
   *  @return The determinant value of m
   **/
  static Type invertSymMatrix11( Matrix const& m, CArraySquare<Type, Size_>& inv)
  {
    const int iShift = m.beginRows(), jShift = m.beginCols();
    inv.resize(TRange<Size_>(0, 1));
    // cofactor (0,0) [0]
    inv(0, 0) = Type(1);
    // compute determinant
    Type det = m(iShift+0, jShift+0);
    if (det == Type(0)) return Type(0);
    // compute inverse matrix
    inv /= det;
    return det;
  }
  /** @ingroup hidden
   *  @brief compute the inverse of the symmetric 2x2 matrix m using only its lower part
   *  and store the result in inv.
   *
   *  @param m, inv the matrices to invert and its inverse
   *  @return The determinant value of m
   **/
  static Type invertSymMatrix22( Matrix const& m, CArraySquare<Type, Size_>& inv)
  {
    const int iShift = m.beginRows(), jShift = m.beginCols();
    inv.resize(TRange<Size_>(0, 2));
    // cofactor (0,0) [0]
    inv(0, 0) =   m(iShift+1, jShift+1);
    // cofactor (1,0) [1]
    inv(1, 0) = - m(iShift+1, jShift+0);
    // cofactor (1,1) [3]
    inv(1, 1) =   m(iShift+0, jShift+0);
    // symmetry
    inv(0, 1) = inv(1,0);
    // compute determinant
    Type det = m(iShift+0, jShift+0) * inv(0,0)
             + m(iShift+1, jShift+0) * inv(1,0);
    if (det == Type(0)) return Type(0);
    // inverse matrix
    inv /= det;
    return det;
  }
  /** @ingroup hidden
   *  @brief compute the inverse of the symmetric 3x3 matrix m using only its lower part
   *  and store the result in inv.
   *  @param m, inv the matrices to invert and its inverse
   *  @return The determinant value of m
   **/
  static Type invertSymMatrix33( Matrix const& m, CArraySquare<Type, Size_>& inv)
  {
    const int iShift = m.beginRows(), jShift = m.beginCols();
    inv.resize(TRange<Size_>(0, 3));
    // cofactor (0,0) [0]
    inv(0, 0) = m(iShift+1, jShift+1) * m(iShift+2, jShift+2)
              - m(iShift+2, jShift+1) * m(iShift+2, jShift+1);
    // cofactor (1,0) [1]
    inv(1, 0) = m(iShift+1, jShift+2) * m(iShift+2, jShift+0)
              - m(iShift+1, jShift+0) * m(iShift+2, jShift+2);
    inv(2, 0) = m(iShift+1, jShift+0) * m(iShift+2, jShift+1)
              - m(iShift+2, jShift+0) * m(iShift+1, jShift+1);
    inv(1, 1) = m(iShift+0, jShift+0) * m(iShift+2, jShift+2)
              - m(iShift+0, jShift+2) * m(iShift+2, jShift+0);
    inv(2, 1) = m(iShift+2, jShift+0) * m(iShift+1, jShift+0)
              - m(iShift+0, jShift+0) * m(iShift+2, jShift+1);
    inv(2, 2) = m(iShift+0, jShift+0) * m(iShift+1, jShift+1)
              - m(iShift+1, jShift+0) * m(iShift+0, jShift+1);
    // symmetric
    inv(0,1) = inv(1,0); inv(0,2) = inv(2,0);
    inv(1,2) = inv(2,1);
    // compute determinant
    Type det = m(iShift+0, jShift+0) * inv(0,0)
             + m(iShift+1, jShift+0) * inv(1,0)
             + m(iShift+2, jShift+0) * inv(2,0);
    if (det == Type(0)) return Type(0);
    // inverse matrix
    inv /= det;
    return det;
  }
  /** @ingroup hidden
   *  @brief compute the inverse of the symmetric 4x4 matrix m using only its lower part
   *  and store the result in inv.
   *  @param m, inv the matrices to invert and its inverse
   *  @return The determinant value of m
   **/
  static Type invertSymMatrix44( Matrix const& m, CArraySquare<Type, Size_>& inv)
  {
    const int iShift = m.beginRows(), jShift = m.beginCols();
    inv.resize(TRange<Size_>(0, 4));
    // cofactor (0,0) [0]
    inv(0,0) = m(iShift+1, jShift+1) * m(iShift+2, jShift+2) * m(iShift+3, jShift+3)
             - m(iShift+1, jShift+1) * m(iShift+3, jShift+2) * m(iShift+3, jShift+2)
             - m(iShift+2, jShift+1) * m(iShift+2, jShift+1) * m(iShift+3, jShift+3)
             + m(iShift+2, jShift+1) * m(iShift+3, jShift+1) * m(iShift+3, jShift+2)
             + m(iShift+3, jShift+1) * m(iShift+2, jShift+1) * m(iShift+3, jShift+2)
             - m(iShift+3, jShift+1) * m(iShift+3, jShift+1) * m(iShift+2, jShift+2);
    // cofactor (1,0) [1]
    inv(1,0) = -m(iShift+1, jShift+0) * m(iShift+2, jShift+2) * m(iShift+3, jShift+3)
             +  m(iShift+1, jShift+0) * m(iShift+3, jShift+2) * m(iShift+3, jShift+2)
             +  m(iShift+2, jShift+1) * m(iShift+2, jShift+0) * m(iShift+3, jShift+3)
             -  m(iShift+2, jShift+1) * m(iShift+3, jShift+0) * m(iShift+3, jShift+2)
             -  m(iShift+3, jShift+1) * m(iShift+2, jShift+0) * m(iShift+3, jShift+2)
             +  m(iShift+3, jShift+1) * m(iShift+3, jShift+0) * m(iShift+2, jShift+2);
    // cofactor (2,0) [2]
    inv(2,0) = m(iShift+1, jShift+0) * m(iShift+2, jShift+1) * m(iShift+3, jShift+3)
             - m(iShift+1, jShift+0) * m(iShift+3, jShift+1) * m(iShift+3, jShift+2)
             - m(iShift+1, jShift+1) * m(iShift+2, jShift+0) * m(iShift+3, jShift+3)
             + m(iShift+1, jShift+1) * m(iShift+3, jShift+0) * m(iShift+3, jShift+2)
             + m(iShift+3, jShift+1) * m(iShift+2, jShift+0) * m(iShift+3, jShift+1)
             - m(iShift+3, jShift+1) * m(iShift+3, jShift+0) * m(iShift+2, jShift+1);
    // cofactor (3,0) [3]
    inv(3,0) = -m(iShift+1, jShift+0) * m(iShift+2, jShift+1) * m(iShift+3, jShift+2)
             +  m(iShift+1, jShift+0) * m(iShift+3, jShift+1) * m(iShift+2, jShift+2)
             +  m(iShift+1, jShift+1) * m(iShift+2, jShift+0) * m(iShift+3, jShift+2)
             -  m(iShift+1, jShift+1) * m(iShift+3, jShift+0) * m(iShift+2, jShift+2)
             -  m(iShift+2, jShift+1) * m(iShift+2, jShift+0) * m(iShift+3, jShift+1)
             +  m(iShift+2, jShift+1) * m(iShift+3, jShift+0) * m(iShift+2, jShift+1);
    // cofactor (1,1) [5]
    inv(1,1) = m(iShift+0, jShift+0) * m(iShift+2, jShift+2) * m(iShift+3, jShift+3)
             - m(iShift+0, jShift+0) * m(iShift+3, jShift+2) * m(iShift+3, jShift+2)
             - m(iShift+2, jShift+0) * m(iShift+2, jShift+0) * m(iShift+3, jShift+3)
             + m(iShift+2, jShift+0) * m(iShift+3, jShift+0) * m(iShift+3, jShift+2)
             + m(iShift+3, jShift+0) * m(iShift+2, jShift+0) * m(iShift+3, jShift+2)
             - m(iShift+3, jShift+0) * m(iShift+3, jShift+0) * m(iShift+2, jShift+2);
    // cofactor (2,1) [6]
    inv(2,1) = -m(iShift+0, jShift+0) * m(iShift+2, jShift+1) * m(iShift+3, jShift+3)
             +  m(iShift+0, jShift+0) * m(iShift+3, jShift+1) * m(iShift+3, jShift+2)
             +  m(iShift+1, jShift+0) * m(iShift+2, jShift+0) * m(iShift+3, jShift+3)
             -  m(iShift+1, jShift+0) * m(iShift+3, jShift+0) * m(iShift+3, jShift+2)
             -  m(iShift+3, jShift+0) * m(iShift+2, jShift+0) * m(iShift+3, jShift+1)
             +  m(iShift+3, jShift+0) * m(iShift+3, jShift+0) * m(iShift+2, jShift+1);
    // cofactor (3,1) [7]
    inv(3,1) = m(iShift+0, jShift+0) * m(iShift+2, jShift+1) * m(iShift+3, jShift+2)
             - m(iShift+0, jShift+0) * m(iShift+3, jShift+1) * m(iShift+2, jShift+2)
             - m(iShift+1, jShift+0) * m(iShift+2, jShift+0) * m(iShift+3, jShift+2)
             + m(iShift+1, jShift+0) * m(iShift+3, jShift+0) * m(iShift+2, jShift+2)
             + m(iShift+2, jShift+0) * m(iShift+2, jShift+0) * m(iShift+3, jShift+1)
             - m(iShift+2, jShift+0) * m(iShift+3, jShift+0) * m(iShift+2, jShift+1);
    // cofactor (2,2) [10]
    inv(2,2) = m(iShift+0, jShift+0) * m(iShift+1, jShift+1) * m(iShift+3, jShift+3)
             - m(iShift+0, jShift+0) * m(iShift+3, jShift+1) * m(iShift+3, jShift+1)
             - m(iShift+1, jShift+0) * m(iShift+1, jShift+0) * m(iShift+3, jShift+3)
             + m(iShift+1, jShift+0) * m(iShift+3, jShift+0) * m(iShift+3, jShift+1)
             + m(iShift+3, jShift+0) * m(iShift+1, jShift+0) * m(iShift+3, jShift+1)
             - m(iShift+3, jShift+0) * m(iShift+3, jShift+0) * m(iShift+1, jShift+1);
    // cofactor (3,2) [11]
    inv(3,2) = -m(iShift+0, jShift+0) * m(iShift+1, jShift+1) * m(iShift+3, jShift+2)
             +  m(iShift+0, jShift+0) * m(iShift+3, jShift+1) * m(iShift+2, jShift+1)
             +  m(iShift+1, jShift+0) * m(iShift+1, jShift+0) * m(iShift+3, jShift+2)
             -  m(iShift+1, jShift+0) * m(iShift+3, jShift+0) * m(iShift+2, jShift+1)
             -  m(iShift+2, jShift+0) * m(iShift+1, jShift+0) * m(iShift+3, jShift+1)
             +  m(iShift+2, jShift+0) * m(iShift+3, jShift+0) * m(iShift+1, jShift+1);
    // cofactor (3,3) [15]
    inv(3,3) = m(iShift+0, jShift+0) * m(iShift+1, jShift+1) * m(iShift+2, jShift+2)
             - m(iShift+0, jShift+0) * m(iShift+2, jShift+1) * m(iShift+2, jShift+1)
             - m(iShift+1, jShift+0) * m(iShift+1, jShift+0) * m(iShift+2, jShift+2)
             + m(iShift+1, jShift+0) * m(iShift+2, jShift+0) * m(iShift+2, jShift+1)
             + m(iShift+2, jShift+0) * m(iShift+1, jShift+0) * m(iShift+2, jShift+1)
             - m(iShift+2, jShift+0) * m(iShift+2, jShift+0) * m(iShift+1, jShift+1);
    // symmetric
    inv(0,1) = inv(1,0); inv(0,2) = inv(2,0); inv(0,3) = inv(3,0);
    inv(1,2) = inv(2,1); inv(1,3) = inv(3,1);
    inv(2,3) = inv(3,2);
    // compute determinant
    Type det = m(iShift+0, jShift+0) * inv(0,0)
             + m(iShift+1, jShift+0) * inv(1,0)
             + m(iShift+2, jShift+0) * inv(2,0)
             + m(iShift+3, jShift+0) * inv(3,0);
    if (det == Type(0)) return Type(0);
    // inverse matrix
    inv /= det;
    return det;
  }
  /** @ingroup hidden
   *  @brief compute the inverse of the symmetric 4x4 matrix m using only its lower part
   *  and store the result in inv.
   *  @note if the matrix is not invertible, the result will be a generalized inverse.
   *  @param m, inv the matrices to invert and its inverse
   *  @return The determinant value of m
   **/
  static Type invertSymMatrixXX( Matrix const& m, CArraySquare<Type, Size_>& inv)
  {
#ifdef STKUSELAPACK
      lapack::SymEigen<Matrix> decomp(m);
#else
      SymEigen<Matrix> decomp(m);
#endif
    if (!decomp.run()) return Type(0);
    // compute (generalized) inverse matrix
    decomp.ginv(inv);
    return decomp.det();
  }

};

/** @ingroup hidden utility class allowing to call the correct static function
 *  computing the inverse of a symmetric matrix.
 **/
template<class Matrix, int Size_>
struct InvertSymMatrixDispatcher
{
  typedef typename Matrix::Type Type;
  inline static Type run( Matrix const& m, CArraySquare<Type, Size_>& inv)
  {
    switch (m.sizeRows())
    {
      case 1: return InvertSymImpl<Matrix, Size_>::invertSymMatrix11(m, inv);
      case 2: return InvertSymImpl<Matrix, Size_>::invertSymMatrix22(m, inv);
      case 3: return InvertSymImpl<Matrix, Size_>::invertSymMatrix33(m, inv);
      case 4: return InvertSymImpl<Matrix, Size_>::invertSymMatrix44(m, inv);
      default:return InvertSymImpl<Matrix, Size_>::invertSymMatrixXX(m, inv);
    }
  }
};

template<class Matrix>
struct InvertSymMatrixDispatcher<Matrix, 1>
{
  typedef typename Matrix::Type Type;
  inline static Type run( Matrix const& m, CArraySquare<Type, 1>& inv)
  { return InvertSymImpl<Matrix, 1>::invertSymMatrix11(m, inv);}
};

template<class Matrix>
struct InvertSymMatrixDispatcher<Matrix, 2>
{
  typedef typename Matrix::Type Type;
  inline static Type run( Matrix const& m, CArraySquare<Type, 2>& inv)
  { return InvertSymImpl<Matrix, 2>::invertSymMatrix22(m, inv);}
};

template<class Matrix>
struct InvertSymMatrixDispatcher<Matrix, 3>
{
  typedef typename Matrix::Type Type;
  inline static Type run( Matrix const& m, CArraySquare<Type, 3>& inv)
  { return InvertSymImpl<Matrix, 3>::invertSymMatrix33(m, inv);}
};

template<class Matrix>
struct InvertSymMatrixDispatcher<Matrix, 4>
{
  typedef typename Matrix::Type Type;
  inline static Type run( Matrix const& m, CArraySquare<Type, 4>& inv)
  { return InvertSymImpl<Matrix, 4>::invertSymMatrix44(m, inv);}
};

} // namespace hidden

} // namespace STK

#endif /* STK_SYMINVERT_H */
