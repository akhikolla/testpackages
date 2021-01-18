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
 * created on: 1 avr. 2015
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_MultiLeastSquare.h
 *  @brief In this file we define and implement the class MultiLeastSQquare.
 **/


#ifndef STK_MULTILEASTSQUARE_H
#define STK_MULTILEASTSQUARE_H

#include "STK_ILeastSquare.h"
#include "STK_SymEigen.h"
#include "Arrays/include/STK_CArraySquare.h"

namespace STK
{
// forward declaration
template<class ArrayB, class ArrayA> class MultiLeastSquare;
namespace hidden
{
/** @ingroup hidden
 *  Specialization for the MultiLeastSquare class.
 **/
template<class ArrayB_, class ArrayA_>
struct AlgebraTraits< MultiLeastSquare<ArrayB_, ArrayA_> >
{
  typedef ArrayA_ ArrayA;
  typedef ArrayB_ ArrayB;
};

} // namespace hidden

/** @ingroup Algebra
 *  @brief The class MultiLeastSQquare solve the least square problem when
 *  the response @e b is multidimensional.
 *
 *  The class MultiLeastSquare allows to solve the least-square problem
 *  \f[
 *  \min_{x} \|b - A*x\|^2.
 *  \f]
 *  - A is an m-by-n matrix which may be rank-deficient.
 *  - B can be a vector or a Matrix.
 *
 *  It computes the minimum-norm solution to a real linear least squares
 *  problem: minimize 2-norm(| b - A*x |) using the singular value
 *  decomposition (SVD) of A.
 *  A is an M-by-N matrix which may be rank-deficient.
 *
 *  @sa ILeastSquare, lapack::MultiLeastSquare
 **/
template<class ArrayB, class ArrayA>
class MultiLeastSquare: public ILeastSquare<MultiLeastSquare<ArrayB, ArrayA> >
{
  public:
    typedef ILeastSquare<MultiLeastSquare<ArrayB, ArrayA> > Base;
    using Base::b_;
    using Base::a_;
    using Base::x_;
    using Base::rank_;
    /** @brief constructor
     *  @param b,a the left hand side and the right hand side of the least square problem.
     *  @param isBref,isAref are the left hand side and the right hand side references ?
     */
    inline MultiLeastSquare( ArrayB const& b, ArrayA const& a, bool isBref=false, bool isAref=false)
                          : Base(b, a, isBref, isAref) {}
    /** Destructor */
    inline virtual ~MultiLeastSquare() {};
    /** compute the multidimensional regression */
    bool runImpl();
    /** compute the weighted multidimensional regression */
    template<class Weights>
    bool runImpl(Weights const& weights);
};

template<class ArrayB, class ArrayA>
bool MultiLeastSquare<ArrayB, ArrayA>::runImpl()
{
  // compute a'a
  ArraySquareX prod = a_.transpose() * a_;
  // compute (a'a)^{-1}
  SymEigen<ArraySquareX> decomp(prod);
  decomp.run();
  rank_ = decomp.rank();
  // compute (a'a)^{-1}b'a
  if (a_.sizeRows() < b_.sizeCols())
  { x_ = (decomp.ginv(prod) * a_.transpose()) * b_;}
  else
  { x_ = decomp.ginv(prod) * (a_.transpose() * b_);}
  return true;
}

/* compute the weighted multidimensional regression */
template<class ArrayB, class ArrayA>
template<class Weights>
bool MultiLeastSquare<ArrayB, ArrayA>::runImpl(Weights const& weights)
{
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Weights);
  // compute a'a
  CSquareX prod = a_.transpose() * weights.diagonalize() * a_;
  // compute (a'a)^{-1}
  SymEigen<CSquareX> decomp(prod);
  decomp.run();
  rank_ = decomp.rank();
  // compute (a'a)^{-1}b'a
  if (a_.sizeRows() < b_.sizeCols())
  { x_ = (decomp.ginv(prod) * a_.transpose())  * weights.diagonalize() * b_;}
  else
  { x_ = decomp.ginv(prod) * (a_.transpose() * weights.diagonalize() * b_);}
  return true;
}

} // namespace STK
#endif /* STK_MULTILEASTSQUARE_H */
