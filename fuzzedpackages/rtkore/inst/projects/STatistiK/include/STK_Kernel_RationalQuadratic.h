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

/** @file STK_Kernel_RationalQuadratic.h
 *  @brief In this file we define the class and methods for computing a RationalQuadratic Kernel.
 **/


#ifndef STK_KERNEL_RATIONALQUADRATIC_H
#define STK_KERNEL_RATIONALQUADRATIC_H

#include "STK_Kernel_IKernelBase.h"

namespace STK
{

namespace Kernel
{
/** @ingroup Kernel
 * The RationalQuadratic Kernel is a kernel of the form
 * \f[
 * k(x,y) = 1 - \left( \frac{\|x-y\|^2}{\|x-y\|^2+h} \right)
 * \f]
 * where @e h represents the bandwidth of the kernel.
 */
template<class Array>
class RationalQuadratic: public IKernelBase<Array>
{
  public:
    typedef IKernelBase<Array> Base;
    typedef typename Array::Type Type;
    using Base::p_data_;
    using Base::gram_;
    using Base::hasRun_;

    /** Default constructor with the width
     *  @param shift the shift to use in the kernel
     **/
    RationalQuadratic( Real const& shift= 1.)
          : Base(0), shift_(shift)
    {}
   /** constructor with a constant pointer on the data set
     *  @param p_data a pointer on a data set that will be "kernelized"
     *  @param shift the shift to use in the kernel
     **/
    RationalQuadratic( Array const* p_data, Real const& shift= 1.)
                    : Base(p_data), shift_(shift)
    { if (shift_ == 0.)
        STKDOMAIN_ERROR_1ARG(RationalQuadratic::RationalQuadratic,shift,shift must be!=0);
    }
    /** constructor with a constant pointer on the data set
     *  @param data a reference on a data set that will be "kernelized"
     *  @param shift the size of the windows to use in the kernel
     **/
    RationalQuadratic( Array const& data, Real const& shift= 1.)
                    : Base(data),shift_(shift)
    { if (shift_ == 0.)
        STKDOMAIN_ERROR_1ARG(RationalQuadratic::RationalQuadratic,shift,shift must be!=0);
    }
    /** constructor with an array of parameter.
     *  @param p_data a pointer on a data set that will be "kernelized"
     *  @param param array of parameter
     **/
    template<class Derived>
    RationalQuadratic( Array const* p_data, ExprBase<Derived> const& param)
                    : Base(p_data), shift_(param.empty() ? 1. : param.front())
    {}
    /** constructor with a constant pointer on the data set
     *  @param data a reference on a data set that will be "kernelized"
     *  @param param array of parameter
     **/
    template<class Derived>
    RationalQuadratic( Array const& data, ExprBase<Derived> const& param)
                    : Base(data), shift_(param.empty() ? 1. : param.front())
    {}

    /** destructor */
    virtual ~RationalQuadratic() {}
    /** @return the shift of the kernel */
    Real const& shift() const {return shift_;}
    /** set the shift of the kernel */
    void setShift(Real const& shift) { shift_ = shift;}
    /** Set parameter using an array
     *  @param param array of parameter
     **/
    template<class Derived>
    void setParam(  ExprBase<Derived> const& param)
    { shift_ = (param.empty() ? 1. : param.front());}

    /** virtual method.
     *  @return diagonal value of the kernel for the ith individuals.
     *  @param i index of the individual
     **/
    virtual inline Real diag(int i) const {return 1.;};
    /** compute the kernel value between two individuals
     *  @param i,j index of the two individuals to compare using the kernel metric */
    virtual Real comp(int i, int j) const;
    /** compute the value of the kernel for the given value
     *  @param v value
     *  @return the value of the kernel at v
     **/
    virtual Real value(Type const& v) const;

  private:
    /** shift of the kernel */
    Real shift_;
};

template<class Array>
Real RationalQuadratic<Array>::comp(int i, int j) const
{
  if (hasRun_) return gram_(i,j);
  Real aux = (p_data_->row(i) - p_data_->row(j)).norm2();
  return 1. - aux/(aux + shift_);
}

template<class Array>
Real RationalQuadratic<Array>::value(Type const& v) const
{ return 1- v*v/(v*v + shift_);}

} // namespace Kernel

} // namespace STK

#endif /* STK_KERNEL_RATIONALQUADRATIC_H */
