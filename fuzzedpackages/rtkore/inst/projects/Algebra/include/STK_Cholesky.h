/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2017  Serge Iovleff, Université Lille 1, Inria

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
 * created on: 29 mars 2017
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Cholesky.h
 *  @brief In this file we implement the modified Cholesky decomposition.
 **/


#ifndef STK_CHOLESKY_H
#define STK_CHOLESKY_H

#include <Sdk/include/STK_Macros.h>
#include <Arrays/include/STK_Array2DDiagonal.h>
#include <Arrays/include/STK_Array2DLowerTriangular.h>
#include <Arrays/include/STK_CArraySquare.h>

namespace STK
{

/** @ingroup Algebra
 *  @brief Compute the Cholesky decomposition of a symmetric matrix.
 *  Compute a closely related variant of the classical Cholesky decomposition,
 *  the LDL decomposition \f$ \mathbf{A = L D L}^{*}\f$
 *  where L is a lower unit triangular matrix and D is a diagonal matrix.
 *
 *  \f[ D_{jj} = A_{jj} - \sum_{k=1}^{j-1} L_{jk}^2 D_{kk} \f]
 *  \f[ L_{ij} = \frac{1}{D_{jj}} \left( A_{ij} - \sum_{k=1}^{j-1} L_{ik} L_{jk} D_{kk} \right),
 *    \qquad\text{pour } i>j.\f]
 * @param A the square and symmetric matrix to decompose, only the lower part is used
 * @param D diagonal matrix of the decomposition
 * @param L lower triangular matrix of the decomposition
 * @return @c true if no error, @c false otherwise (results are misleading)
 */
template < class Lhs>
bool cholesky( ExprBase<Lhs> const& A
             , Array2DDiagonal<typename Lhs::Type>& D
             , Array2DLowerTriangular<typename Lhs::Type>& L)
{
  typedef typename Lhs::Type Type;
  if(A.rows() != A.cols())
  { STKRUNTIME_ERROR_NO_ARG(cholesky,A must be square symmetric matrix);}
  // start
  D.resize(A.rows(), A.cols());
  L.resize(A.rows(), A.cols());
  bool nozero = true;
  for (int j=A.beginCols(); j<A.endCols(); ++j )
  {
    Type sum1 = A(j,j);
    for (int k=A.beginCols(); k<j; ++k) { sum1 -= L(j,k) * L(j,k) * D[k];}
    D[j] = sum1;
    nozero &= (sum1!=0);
    L(j,j) = Type(1);
    for (int i=j+1; i<A.endRows(); ++i)
    {
      Type sum2 = A(i,j);
      for (int k=A.beginCols(); k<j; ++k) { sum2 -= L(i,k)*L(j,k)* D[k];}
      L(i,j) = sum1 ? sum2/sum1 : sum2;
    }
  }
  return nozero;
}


} // namespace STK

#endif /* STK_CHOLESKY_H */
