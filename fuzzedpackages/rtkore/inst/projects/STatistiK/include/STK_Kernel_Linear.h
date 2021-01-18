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

/** @file STK_Kernel_Linear.h
 *  @brief In this file we define the class and methods for computing a Linear Kernel.
 **/


#ifndef STK_KERNEL_LINEAR_H
#define STK_KERNEL_LINEAR_H

#include "STK_Kernel_IKernelBase.h"

namespace STK
{

namespace Kernel
{
/** @ingroup Kernel
 * The Linear Kernel is a kernel of the form
 * \f[
 * k(x,y) = <x,y>.
 * \f]
 */
template<class Array>
class Linear: public IKernelBase<Array>
{
  public:
    typedef IKernelBase<Array> Base;
    typedef typename Array::Type Type;
    using Base::p_data_;
    using Base::gram_;
    using Base::hasRun_;

    /** Default constructor */
    Linear(): Base(0) {}
    /** constructor with a constant pointer on the data set
     *  @param p_data a pointer on a data set that will be "kernelized"
     **/
    Linear(Array const* p_data): Base(p_data) {}
    /** constructor with a constant pointer on the data set
     *  @param data a reference on a data set that will be "kernelized"
     **/
    Linear(Array const& data): Base(data) {}
    /** constructor with an array of parameter.
     *  @param p_data a pointer on a data set that will be "kernelized"
     *  @param param array of parameter
     **/
    template<class Derived>
    Linear( Array const* p_data, ExprBase<Derived> const& param): Base(p_data)
    {}
    /** constructor with a constant pointer on the data set
     *  @param data a reference on a data set that will be "kernelized"
     *  @param param array of parameter
     **/
    template<class Derived>
    Linear( Array const& data, ExprBase<Derived> const& param): Base(data)
    {}

    /** destructor */
    virtual ~Linear() {}
    /** Set parameter using an array
     *  @param param array of parameter
     **/
    template<class Derived>
    void setParam(  ExprBase<Derived> const& param) {}

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
};

/* virtual method.
 *  @return diagonal value of the kernel for the ith individuals.
 *  @param i index of the individual
 **/
template<class Array>
inline Real Linear<Array>::diag(int i) const
{ return hasRun_ ? gram_(i,i)
                 :  p_data_->row(i).norm2();
}

template<class Array>
Real Linear<Array>::comp(int i, int j) const
{ return hasRun_ ? gram_(i,j)
                 :  p_data_->row(i).dot(p_data_->row(j));}

} // namespace Kernel

} // namespace STK

#endif /* STK_KERNEL_LINEAR_H */
