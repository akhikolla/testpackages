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
 * Project:  stkpp::STatisitK::Kernel
 * created on: 2 feb. 2016
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Kernel_Hamming.h
 *  @brief In this file we define the class and methods for computing a Hamming Kernel.
 **/


#ifndef STK_KERNEL_HAMMING_H
#define STK_KERNEL_HAMMING_H

#include "STK_Kernel_IKernelBase.h"
#include "STK_Stat_MultiFactor.h"

namespace STK
{

namespace Kernel
{
/** @ingroup Kernel
 * The Hamming Kernel is a kernel of the form
 * \f[
 * k(x,y) = \sum_{u\in D^p} \prod_{j=1}^p \lambda^{\delta(u_j,x_j)}\lambda^{\delta(u_j,y_j)}
 * \f]
 * where \f$ \lambda \in(0,1) \f$ represents the similarity index of the kernel.
 * It can be computed recursively.
 *
 * @note The creation of this Kernel trigger the computation of the factors
 * present in the data set.
 * @sa STK::Stat::MultiFactors
 */
template<class Array>
class Hamming: public IKernelBase<Array>
{
  public:
    typedef IKernelBase<Array> Base;
    typedef typename Array::Type Type;
    using Base::p_data_;
    using Base::gram_;
    using Base::hasRun_;

    /** constructor with a constant pointer on the data set
     *  @param p_data a pointer on a data set that will be "kernelized"
     *  @param lambda the size of the windows to use in the kernel
     **/
    Hamming( Array const* p_data, Real const& lambda= 1.)
          : Base(p_data), lambda_(lambda), diagElt_(1.), factors_(p_data)
    {
      if (!p_data)
      { STKRUNTIME_ERROR_NO_ARG(Hamming::Hamming(p_data,lambda),p_data is 0);}
      factors_.run();
      computeDiagonalElement();
    }
    /** constructor with a constant pointer on the data set
     *  @param data a reference on a data set that will be "kernelized"
     *  @param lambda the size of the windows to use in the kernel
     **/
    Hamming( Array const& data, Real const& lambda= 1.)
          : Base(data), lambda_(lambda), diagElt_(1.), factors_(data)
    {
      factors_.run();
      computeDiagonalElement();
    }
    /** constructor with an array of parameter.
      *  @param p_data a pointer on a data set that will be "kernelized"
      *  @param param array of parameter
      **/
     template<class Derived>
     Hamming( Array const* p_data, ExprBase<Derived> const& param)
           : Base(p_data), lambda_(param.empty() ? 1. : param.front())
            , diagElt_(1.), factors_(p_data)
     {
         if (!p_data)
         { STKRUNTIME_ERROR_NO_ARG(Hamming::Hamming(p_data,lambda),p_data is 0);}
         factors_.run();
         computeDiagonalElement();
     }
     /** constructor with a constant pointer on the data set
      *  @param data a reference on a data set that will be "kernelized"
      *  @param param array of parameter
      **/
     template<class Derived>
     Hamming( Array const& data, ExprBase<Derived> const& param)
           : Base(data), lambda_(param.empty() ? 1. : param.front())
            , diagElt_(1.), factors_(data)
     {
         factors_.run();
         computeDiagonalElement();
     }

    /** destructor */
    virtual ~Hamming() {}

    /** @return the lambda of the kernel */
    Real const& lambda() const {return lambda_;}
    /** @return the lambda of the kernel */
    Stat::MultiFactor<Array> const& factors() const {return factors_;}
    /** set the lambda of the kernel */
    void setLambda(Real const& lambda)
    {
      lambda_ = lambda;
      computeDiagonalElement();
    }
    /** Set parameter using an array
     *  @param param array of parameter
     **/
    template<class Derived>
    void setParam(  ExprBase<Derived> const& param)
    {
      lambda_ = (param.empty() ? 1. : param.front());
      computeDiagonalElement();
    }

    /** virtual method.
     *  @return diagonal value of the kernel for the ith individuals.
     *  @param i index of the individual
     **/
    virtual Real diag(int i) const;
    /** virtual method implementation.
     *  @return value of the kernel for the ith and jth individuals.
     *  @param i,j indexes of the individuals
     **/
    virtual Real comp(int i, int j) const;

  private:
    /** lambda of the kernel */
    Real lambda_;
    /** diagonal element of the kernel */
    Real diagElt_;
    /** factors of the values */
    Stat::MultiFactor<Array> factors_;
    /** Compute the diagonal element of the kernel*/
    void computeDiagonalElement();
};

/* Compute the diagonal element of the kernel*/
template<class Array>
void Hamming<Array>::computeDiagonalElement()
{
   diagElt_ = 1.;
   for(int j=factors_.nbLevels().begin(); j < factors_.nbLevels().end(); ++j)
   { diagElt_ *= lambda_*lambda_*(factors_.nbLevels()[j]-1) + 1;}
}

template<class Array>
inline Real Hamming<Array>::diag(int i) const { return diagElt_;}

template<class Array>
Real Hamming<Array>::comp(int i, int j) const
{
  typedef typename hidden::Traits<Array>::Row RowVector;
  if (hasRun_) return gram_(i,j);
  // create references on row i and j
  RowVector ind1(p_data_->row(i), true), ind2(p_data_->row(j), true);
  Real value = 1.;
  for(int j=factors_.nbLevels().begin(); j < factors_.nbLevels().end(); ++j)
  {
    value *= (ind1[j]==ind2[j]) ? lambda_*lambda_*(factors_.nbLevels()[j]-1) + 1.
                                : lambda_*(lambda_*(factors_.nbLevels()[j]-2) + 2.);
  }
  return value;
}

} // namespace Kernel

} // namespace STK

#endif /* STK_KERNEL_HAMMING_H */
