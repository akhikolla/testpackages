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
 * Project:  stkpp::Stat::Kernel
 * created on: 5 avr. 2015
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Kernel_IKernelBase.h
 *  @brief In this file we define the Interface base class for computing a Kernels.
 **/


#ifndef STK_KERNEL_IKERNELBASE_H
#define STK_KERNEL_IKERNELBASE_H

#include "STK_Kernel_IKernel.h"

namespace STK
{
namespace Kernel
{
/** @ingroup Kernel
 *  Interface Base class for the kernels classes.
 */
template<class Array>
class IKernelBase: public IKernel
{
  public:
    typedef typename Array::Type Type;
    using IKernel::gram_;

    /** constructor with a constant pointer on the data set
     *  @param p_data a pointer on a data set that will be "kernelized"
     **/
    IKernelBase(Array const* p_data): IKernel(), p_data_(p_data) {}
    /** constructor with a constant reference on the data set
     *  @param data a reference on a data set that will be "kernelized"
     **/
    IKernelBase(Array const& data): IKernel(), p_data_(&data) {}
    /** copy constructor
     *  @param kernel kernel to copy
     **/
    IKernelBase(IKernelBase const& kernel): IKernel(kernel), p_data_(kernel.p_data_) {}
    /** destructor */
    virtual ~IKernelBase() {}

    /** @return the pointer on the data set */
    Array const* p_data() const { return p_data_;}

    /** compute Gram matrix
     *  @return @c true if the computation is successful, @c false otherwise */
    virtual bool run();
    /** @return the number of samples (the number of rows in the data set) */
    virtual int nbSample() const { return (p_data_) ? p_data_->sizeRows() : 0;}
    /** @return the number of variables (the number of columns in the data set) */
    virtual int nbVariable() const  { return (p_data_) ? p_data_->sizeCols() : 0;}

    /** compute the value of the kernel for the given value
     *  @param v value
     *  @return the value of the kernel at v
     **/
    virtual Real value(Type const& v) const
    { return 0;}

  protected:
    /** pointer on the data set */
    Array const* p_data_;
};

template<class Array>
bool IKernelBase<Array>::run()
{
  if(!p_data_) return false;
  gram_.resize(p_data_->rows());
  // upper part
  for (int j= gram_.begin(); j < gram_.end(); ++j)
  {
    { gram_(j,j) = this->diag(j);}
    for (int i= j+1; i < gram_.end(); ++i)
    { gram_(i,j) = this->comp(i,j);}
  }
  // lower part
  for (int j= gram_.begin(); j < gram_.end(); ++j)
  {
    for (int i= j+1; i < gram_.end(); ++i)
    { gram_(j,i) = gram_(i,j);}
  }
  this->hasRun_ = true;
  return true;
}


} // namespace Kernel

} // namespace STK

#endif /* STK_KERNEL_IKERNELBASE_H */
