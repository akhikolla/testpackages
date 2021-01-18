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

/** @file STK_ILeastSquare.h
 *  @brief In this file we define the interface class ILeastSquare.
 **/


#ifndef STK_ILEASTSQUARE_H
#define STK_ILEASTSQUARE_H

#include <Sdk/include/STK_IRunner.h>
#include "STK_Algebra_Util.h"

namespace STK
{
/** @ingroup Algebra
 *  @brief The class ILeastSquare is an interface class for the methods
 *  solving the least-square problem
 *  \f[
 *  \min_{x} \|b - A*x\|^2.
 *  \f]
 *  - A is an m-by-n matrix which may be rank-deficient.
 *  - B can be a vector or a Matrix.
 *
 *  The least-square problem can be solved using the Svd or the QR
 *  decomposition of A.
 **/
template<class Derived>
class ILeastSquare: public IRunnerBase, public IRecursiveTemplate<Derived>
{
  protected:
    typedef typename hidden::AlgebraTraits<Derived>::ArrayB ArrayB;
    typedef typename hidden::AlgebraTraits<Derived>::ArrayA ArrayA;
    /** @brief Default constructor
     *  @param b,a the left hand side and the right hand side of the least square problem.
     *  @param isBref,isAref are the left hand side and the right hand side references ?
     **/
    ILeastSquare( ArrayB const& b, ArrayA const& a, bool isBref=false, bool isAref=false)
                 : b_(b, isBref), a_(a, isAref), rank_(0)
    {
      if (b.beginRows() != b.beginCols())
        STKRUNTIME_ERROR_NO_ARG(ILeastSquare::ILeastSquare,Wrong data set: b.beginRows() must be equal to b.beginCols());
      if (a.beginRows() != a.beginCols())
        STKRUNTIME_ERROR_NO_ARG(ILeastSquare::ILeastSquare,Wrong data set: a.beginRows() must be equal to a.beginCols());
      if (a.beginRows() != b.beginRows())
        STKRUNTIME_ERROR_NO_ARG(ILeastSquare::ILeastSquare,Wrong data set: a.beginRows() must be equal to b.beginRows());
    }
    /** @brief template constructor
     *  @param b,a the left hand side and the right hand side of the least square
     *  problem.
     */
    template<class OtherArrayB, class OtherArrayA>
    ILeastSquare( ExprBase<OtherArrayB> const& b, ExprBase<OtherArrayA> const& a)
                 : b_(b), a_(a), rank_(0)
    {
      if (b.beginRows() != b.beginCols())
        STKRUNTIME_ERROR_NO_ARG(ILeastSquare::ILeastSquare,Wrong data set: b.beginRows() must be equal to b.beginCols());
      if (a.beginRows() != a.beginCols())
        STKRUNTIME_ERROR_NO_ARG(ILeastSquare::ILeastSquare,Wrong data set: a.beginRows() must be equal to a.beginCols());
      if (a.beginRows() != b.beginRows())
        STKRUNTIME_ERROR_NO_ARG(ILeastSquare::ILeastSquare,Wrong data set: a.beginRows() must be equal to b.beginRows());
    }

  public:
    /** Destructor */
    virtual ~ILeastSquare() {};
    // getters
    /** @return the rank of the matrix A */
    Real const& rank() const { return rank_;}
    /** @return the matrix A */
    inline ArrayA const& a() const { return a_;}
    /** @return the matrix B */
    inline ArrayB const& b() const { return b_;}
    /** @return the solution of the least-square problem */
    inline ArrayB const& x() const { return x_;}
    /** Compute the Least-Square solution.
     *  Delegate to Derived classes the concrete computation of the
     *  decomposition using @c runImpl method.
     *  @return @c true if the computation succeed, @c false otherwise
     **/
    virtual bool run();
    /** Compute the weighted Least-Square solution.
     *  Delegate to Derived classes the concrete computation of the
     *  decomposition using @c runImpl method.
     *  @return @c true if the computation succeed, @c false otherwise
     **/
    template<class VecWeights>
    bool run(VecWeights const& weights);

  protected:
    /** Array or vector of the left hand side. */
    ArrayB b_;
    /** Array of the right hand side */
    ArrayA a_;
    /** Array of the solution (a vector if b is a vector, a matrix otherwise) */
    ArrayB x_;
    /** rank of matrix A*/
    int rank_;
};

template<class Derived>
bool ILeastSquare<Derived>::run()
{
  if (a_.empty()||b_.empty()) { return true;}
  // compute Least Square
  return(this->asDerived().runImpl());
}

template<class Derived>
template<class Weights>
bool ILeastSquare<Derived>::run(Weights const& weights)
{
  if (a_.empty()||b_.empty()) { return true;}
  // compute Least Square
  return(this->asDerived().runImpl(weights));
}


} // namespace STK
#endif /* STK_ILEASTSQUARE_H */
