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
 * Project:  stkpp::
 * created on: 5 avr. 2015
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Kernel_Polynomial.h
 *  @brief In this file we define the class and methods for computing a Polynomial Kernel.
 **/


#ifndef STK_KERNEL_POLYNOMIAL_H
#define STK_KERNEL_POLYNOMIAL_H

#include "STK_Kernel_IKernelBase.h"

namespace STK
{

namespace Kernel
{
/** @ingroup Kernel
 * The Polynomial Kernel is a kernel of the form
 * \f[
 * k(x,y) = \left(<x-y>+c\right)^d
 * \f]
 * where @e c  represents the shift of the kernel (default is 0)
 * and @c d represents the degree.
 */
template<class Array>
class Polynomial: public IKernelBase<Array>
{
  public:
    typedef IKernelBase<Array> Base;
    typedef typename Array::Type Type;
    using Base::p_data_;
    using Base::gram_;
    using Base::hasRun_;

    /** Default constructor with the degree and the shift
     *  @param d degree of the polynomial
     *  @param shift the shift to use in the kernel
     **/
    Polynomial( Real const& d=2., Real const& shift= 0)
             : Base(0), d_(d), shift_(shift)
    {}
    /** constructor with a constant pointer on the data set
     *  @param p_data a pointer on a data set that will be "kernelized"
     *  @param d degree of the polynomial
     *  @param shift the shift to use in the kernel
     **/
    Polynomial( Array const* p_data, Real const& d=2., Real const& shift= 0)
             : Base(p_data), d_(d), shift_(shift)
    { if (d_ <= 0.)
      STKDOMAIN_ERROR_2ARG(Polynomial::Polynomial,shift,d,d must be>0);
    }
    /** constructor with a constant reference on the data set
     *  @param data a reference on a data set that will be "kernelized"
     *  @param d degree of the polynomial
     *  @param shift the shift to use in the kernel
     **/
    Polynomial( Array const& data, Real const& d=2., Real const& shift= 0.)
             : Base(data), d_(d), shift_(shift)
    { if (d_ <= 0.)
      STKDOMAIN_ERROR_2ARG(Polynomial::Polynomial,shift,d,d must be>0);
    }
    /** constructor with an array of parameter.
     *  @param p_data a pointer on a data set that will be "kernelized"
     *  @param param array of parameter
     **/
    template<class Derived>
    Polynomial( Array const* p_data, ExprBase<Derived> const& param)
             : Base(p_data)
              , d_(param.empty() ? 2. : param.front())
              , shift_(param.empty() ? 2. : param.elt(param.begin()+1))
    {}
    /** constructor with a constant pointer on the data set
     *  @param data a reference on a data set that will be "kernelized"
     *  @param param array of parameter
     **/
    template<class Derived>
    Polynomial( Array const& data, ExprBase<Derived> const& param)
             : Base(data)
              , d_(param.empty() ? 2. : param.front())
              , shift_(param.empty() ? 2. : param.elt(param.begin()+1))
    {}

    /** destructor */
    virtual ~Polynomial() {}
    /** @return the degree of the kernel */
    Real const& degree() const {return d_;}
    /** set the degree of the kernel */
    void setDegree(Real const& d) {d_ = d;}
    /** @return the shift of the kernel */
    Real const& shift() const {return shift_;}
    /** set the shift of the kernel */
    void setShift(Real const& shift) { shift_ = shift;}
    /** Set parameter using an array
     *  @param param array of parameter
     **/
    template<class Derived>
    void setParam(  ExprBase<Derived> const& param)
    { d_ = (param.empty() ? 2. : param.front());
      shift_ = (param.empty() ? 0. : param.elt(param.begin()+1));
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
    /** degree of the kernel */
    Real d_;
    /** shift of the kernel */
    Real shift_;
};

/* virtual method.
 *  @return diagonal value of the kernel for the ith individuals.
 *  @param i index of the individual
 **/
template<class Array>
inline Real Polynomial<Array>::diag(int i) const
{ return hasRun_ ? gram_(i,i)
                 :  std::pow(p_data_->row(i).norm2() + shift_, d_);
}

template<class Array>
inline Real Polynomial<Array>::comp(int i, int j) const
{ return hasRun_ ? gram_(i,j)
                 :  std::pow(p_data_->row(i).dot(p_data_->row(j)) + shift_, d_);}


} // namespace Kernel

} // namespace STK

#endif /* STK_KERNEL_POLYNOMIAL_H */
