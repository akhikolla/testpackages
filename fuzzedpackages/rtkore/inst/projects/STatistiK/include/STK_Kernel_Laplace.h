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

/** @file STK_Kernel_Laplace.h
 *  @brief In this file we define the class and methods for computing a Laplace Kernel.
 **/


#ifndef STK_KERNEL_LAPLACE_H
#define STK_KERNEL_LAPLACE_H

#include "STK_Kernel_IKernelBase.h"

namespace STK
{

namespace Kernel
{
/** @ingroup Kernel
 * The Laplace Kernel is a kernel of the form
 * \f[
 * k(x,y) = \exp\left(- \frac{\|x-y\|}{h} \right)
 * \f]
 * where @e h represents the bandwidth of the kernel.
 */
template<class Array>
class Laplace: public IKernelBase<Array>
{
  public:
    typedef IKernelBase<Array> Base;
    typedef typename Array::Type Type;
    using Base::p_data_;
    using Base::gram_;
    using Base::hasRun_;

    /** Default constructor with the width
     *  @param width the size of the windows to use in the kernel
     **/
    Laplace( Real const& width= 1.)
          : Base(0), width_(width)
    {}
    /** constructor with a constant pointer on the data set
     *  @param p_data a pointer on a data set that will be "kernelized"
     *  @param width the size of the windows to use in the kernel
     **/
    Laplace( Array const* p_data, Real const& width= 1.)
              : Base(p_data), width_(width) {}
    /** constructor with a constant pointer on the data set
     *  @param data a reference on a data set that will be "kernelized"
     *  @param width the size of the windows to use in the kernel
     **/
    Laplace( Array const& data, Real const& width= 1.)
              : Base(data), width_(width) {}
    /** constructor with an array of parameter.
     *  @param p_data a pointer on a data set that will be "kernelized"
     *  @param param array of parameter
     **/
    template<class Derived>
    Laplace( Array const* p_data, ExprBase<Derived> const& param)
          : Base(p_data), width_(param.empty() ? 1. : param.front())
    {}
    /** constructor with a constant pointer on the data set
     *  @param data a reference on a data set that will be "kernelized"
     *  @param param array of parameter
     **/
    template<class Derived>
    Laplace( Array const& data, ExprBase<Derived> const& param)
          : Base(data), width_(param.empty() ? 1. : param.front())
    {}

    /** destructor */
    virtual ~Laplace() {}

    /** @return the bandwidth of the kernel */
    Real const& width() const {return width_;}
    /** set the bandwidth of the kernel */
    void setWidth(Real const& width) {width_ = width;}
    /** Set parameter using an array
     *  @param param array of parameter
     **/
    template<class Derived>
    void setParam(  ExprBase<Derived> const& param)
    { width_ = (param.empty() ? 1. : param.front());}

    /** virtual method.
     *  @return diagonal value of the kernel for the ith individuals.
     *  @param i index of the individual
     **/
    virtual inline Real diag(int i) const {return 1.;};
    /** virtual method implementation.
     *  @return value of the kernel for the ith and jth individuals.
     *  @param i,j indexes of the individuals
     **/
    virtual Real comp(int i, int j) const;
    /** compute the value of the kernel for the given value
     *  @param v value
     *  @return the value of the kernel at v
     **/
    virtual Real value(Type const& v) const;

  private:
    /** bandwidth of the kernel */
    Real width_;
};

template<class Array>
inline Real Laplace<Array>::comp(int i, int j) const
{
  return hasRun_ ? gram_(i,j)
                 : std::exp(-(p_data_->row(i) - p_data_->row(j)).norm()/width_);
}

template<class Array>
inline Real Laplace<Array>::value(Type const& v) const
{ return std::exp(-std::abs(v)/width_);}

/** @ingroup Kernel
 * [Deprecated] The Exponential Kernel is a kernel of the form
 * \f[
 * k(x,y) = \exp\left(- \frac{\|x-y\|}{h} \right)
 * \f]
 * where @e h represents the bandwidth of the kernel.
 */
template<class Array>
class Exponential: public IKernelBase<Array>
{
  public:
    typedef IKernelBase<Array> Base;
    typedef typename Array::Type Type;
    using Base::p_data_;
    using Base::gram_;
    using Base::hasRun_;

    /** Default constructor with the width
     *  @param width the size of the windows to use in the kernel
     **/
    Exponential( Real const& width= 1.)
          : Base(0), width_(width)
    {}
    /** constructor with a constant pointer on the data set
     *  @param p_data a pointer on a data set that will be "kernelized"
     *  @param width the size of the windows to use in the kernel
     **/
    Exponential( Array const* p_data, Real const& width= 1.)
              : Base(p_data), width_(width) {}
    /** constructor with a constant pointer on the data set
     *  @param data a reference on a data set that will be "kernelized"
     *  @param width the size of the windows to use in the kernel
     **/
    Exponential( Array const& data, Real const& width= 1.)
              : Base(data), width_(width) {}
    /** constructor with an array of parameter.
     *  @param p_data a pointer on a data set that will be "kernelized"
     *  @param param array of parameter
     **/
    template<class Derived>
    Exponential( Array const* p_data, ExprBase<Derived> const& param)
          : Base(p_data), width_(param.empty() ? 1. : param.front())
    {}
    /** constructor with a constant pointer on the data set
     *  @param data a reference on a data set that will be "kernelized"
     *  @param param array of parameter
     **/
    template<class Derived>
    Exponential( Array const& data, ExprBase<Derived> const& param)
          : Base(data), width_(param.empty() ? 1. : param.front())
    {}

    /** destructor */
    virtual ~Exponential() {}

    /** @return the bandwidth of the kernel */
    Real const& width() const {return width_;}
    /** set the bandwidth of the kernel */
    void setWidth(Real const& width) {width_ = width;}
    /** Set parameter using an array
     *  @param param array of parameter
     **/
    template<class Derived>
    void setParam(  ExprBase<Derived> const& param)
    { width_ = (param.empty() ? 1. : param.front());}

    /** virtual method.
     *  @return diagonal value of the kernel for the ith individuals.
     *  @param i index of the individual
     **/
    virtual inline Real diag(int i) const {return 1.;};
    /** virtual method implementation.
     *  @return value of the kernel for the ith and jth individuals.
     *  @param i,j indexes of the individuals
     **/
    virtual Real comp(int i, int j) const;
    /** compute the value of the kernel for the given value
     *  @param v value
     *  @return the value of the kernel at v
     **/
    virtual Real value(Type const& v) const;

  private:
    /** bandwidth of the kernel */
    Real width_;
};

template<class Array>
inline Real Exponential<Array>::comp(int i, int j) const
{
  return hasRun_ ? gram_(i,j)
                 : std::exp(-(p_data_->row(i) - p_data_->row(j)).norm()/width_);
}

template<class Array>
inline Real Exponential<Array>::value(Type const& v) const
{ return std::exp(-std::abs(v)/width_);}
} // namespace Kernel

} // namespace STK

#endif /* STK_KERNEL_LAPLACE_H */
