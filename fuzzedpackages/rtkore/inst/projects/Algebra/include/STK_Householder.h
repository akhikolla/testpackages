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
 * Project:  Algebra
 * Purpose:  Define Householder methods for Algebra classes.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_Householder.h
 *  @brief In this file we define the Housholder methods used by the Algebra classes.
 **/

#ifndef STK_HOUSEHOLDER_H
#define STK_HOUSEHOLDER_H

#include <STKernel/include/STK_Misc.h>
#include <Arrays/include/STK_IArrayBase.h>

namespace STK
{
/** @ingroup Algebra
 *  @brief Compute the Householder vector v of a vector x.
 *
 *  Given a vector x, compute the vector v of the matrix of Householder
 *  \f$ P=I-2vv'/(v'v)  \f$ such that \f$ Px = v1 e_1 \f$.
 *  The vector v is of the form : \f$ (1,x_2/s,...,x_n/s)' \f$
 *  and is stored in x. The value 1 is skipped and
 *  \f$ \beta = -2/(v'v) \f$ is stored in front of v.
 *  The method return the value v1.
 *
 *  @param x the vector to rotate, it is overwritten by v
 **/
template <class Vector>
Real house(IArrayBase< Vector>& x)
{
  // compute L^{\infty} norm of X, declare result and norm2 of X
  Real scale  = x.normInf(), v1, norm2 = 0.0;
  if (scale)  // if not 0.0
  {
    norm2 = (x.asDerived().sub(_R(x.begin()+1,x.lastIdx()))/=scale).norm2();
  }
  // check if the lower part is significative
  if (norm2 < Arithmetic<Real>::epsilon())
  { v1 = x.front(); x.front() = 0.0; }
  else
  {
    Real s, aux1 = x.front()/scale;
    // compute v1 = P_v X and beta of the Householder vector
    v1 =  (norm2 = sign(aux1, sqrt(aux1*aux1+norm2))) * scale;
    // compute and save beta
    x.front() = (s = aux1-norm2)/norm2;
    // comp v and save it
    x.asDerived().sub(_R(x.begin()+1,x.lastIdx())) /= s;
  }
  return v1;
}

/** @ingroup Algebra
 *  @brief dot product with a Householder vector.
 *
 *  Scalar product of an 1D container with a Householder vector
 *  d = < x,v>. The first componant of a Householder vector is 1.0
 *
 *  @param x first vector
 *  @param v the Householder vector
 **/
template< class Lhs, class Rhs>
Real dotHouse( ExprBase< Lhs> const& x, ExprBase< Rhs> const& v)
{
  // compute the product
  Real sum = x[ v.begin()] /* *1.0 */;
  for (int i= v.begin()+1; i<v.end(); i++) sum += x[i] * v[i];
  // return <x,v>
  return(sum);
}

/** @ingroup Algebra
 *  @brief left multiplication by a Householder vector.
 *
 *  Perform a left multiplication of the matrix M with a Householder
 *  matrix \f$ H=I+beta vv' \f$. Overwrite M with HM.
 *
 *  @param M the matrix to multiply (input/output)
 *  @param v the Householder vector (input)
 **/
template < class Lhs, class Rhs>
void applyLeftHouseholderVector( ArrayBase<Lhs> const& M, ExprBase<Rhs> const& v)
{
  typedef typename hidden::Traits<Lhs>::Col ColVector;
  // get beta
  Real beta = v.front();
  if (beta)
  {
    // Multiplication of the cols by P=I+beta vv'
    for (int j=M.beginCols(); j<M.endCols(); j++)
    {
      // a ref on the jth column of M
      ColVector Mj(M.asDerived(), v.range(), j);
      // Computation of aux=beta* <v,M^j>
      Real aux =  dotHouse( Mj, v) * beta;
      // updating columns
      Mj.front() += aux;
      for(int i = v.begin()+1; i < v.end(); ++i) Mj[i] += v[i] * aux;
    }
  }
}

/** @ingroup Algebra
 *  @brief left multiplication by a Householder array.
 *
 *  Perform a left multiplication of the Array M by a Householder Matrix H.
 *  @c M <- HM with @c H = I + WZ'.
 *  The matrix @c H is represented as a product of elementary reflectors
 *  H = H(1) H(2) . . . H(k)
 *
 * @param M the matrix to multiply
 * @param H the matrix with the Householder vectors
 **/
template < class Lhs, class Rhs>
void applyLeftHouseholderArray( ArrayBase<Lhs> const& M, ArrayBase<Rhs> const& H)
{
  typedef typename hidden::Traits<Rhs>::Col ColVector;
  // compute the number of iterations
  int first = H.beginCols(), last = std::min( H.lastIdxCols(), H.lastIdxRows());
  // get range of the first Householder vector
  Range range_ve(last, H.lastIdxRows(), 0);
  // iterations
  for (int j=last; j>=first; j--)
  {
    // apply left Householder vector to M
    ColVector v(H.asDerived(), range_ve, j);
    applyLeftHouseholderVector(M, v);
    // decrease range of the Householder vector
    range_ve.decFirst(1);
  }
}

/** @ingroup Algebra
 *  @brief right multiplication by a Householder vector.
 *
 *  Perform a right multiplication of the matrix M with a Householder
 *  matrix \f$ H=I+beta vv' \f$. Overwrite M with MH.
 *
 *  @param M the matrix to multiply (input/output)
 *  @param v the Householder vector (input)
 **/
template < class Lhs, class Rhs>
void applyRightHouseholderVector( ArrayBase<Lhs> const& M, ExprBase<Rhs> const& v)
{
  // get beta
  Real beta = v.front();
  if (beta)
  {
    // Multiplication of the cols by P=I+beta vv'
    for (int i=M.beginRows(); i<M.endRows(); i++)
    {
      // a ref on the ith row of M
      typename hidden::Traits<Lhs>::Row Mi(M.asDerived(), v.range(), i);
      // Computation of aux=beta* <v,M_i>
      Real aux =  dotHouse( Mi, v) * beta;
      // updating column
      Mi.front() += aux;
      // Computation of M_i + beta <v,M_i>  v = M_i + aux v'
      for (int i=v.begin()+1; i<v.end(); i++) Mi[i] +=  v[i] * aux;
    }
  }
}

/** @ingroup Algebra
 *  @brief left multiplication by a Householder ArrayXX.
 *
 * Perform a right multiplication of the matrix M with a Householder
 * Marix H. M <- MP with H = I + WZ'. The Householder vectors are
 * stored in the rows of H.
 *
 * @param M the Array to multiply
 * @param H the Householder ArrayXX
 **/
template < class TContainer2D, class Rhs>
void applyRightHouseholderArray( ArrayBase<TContainer2D> const& M, ArrayBase<Rhs> const& H)
{
  // compute the number of iterations
  int first = H.beginCols(), last = std::min( H.lastIdxCols(), H.lastIdxRows());
  // get range of the first Householder vector
  Range range_ve(last, H.lastIdxRows());
  // iterations
  for (int j=last; j>=first; j--)
  {
    // apply left Householder vector to M
    typename Rhs::Col v(H.asDerived(), range_ve, j);
    applyRightHouseholderVector(M, v);
    // decrease range of the Householder vector
    range_ve.decFirst(1);
  }
}

} // namespace STK

#endif /*STK_HOUSEHOLDER_H*/
