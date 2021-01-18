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
 * Project:  stkpp::Arrays
 * created on: 27 sept. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ArrayBaseApplier.h
 *  @brief In this file we define the ArrayBaseApplier classes.
 **/

#ifndef STK_ARRAYBASEAPPLIER_H
#define STK_ARRAYBASEAPPLIER_H

namespace STK
{

/* Apply the Visitor @c visitor to the whole coefficients of the array.
  * The template parameter @c Visitor is the type of the visitor and provides
  * the following interface:
  * @code
  * struct MyVisitor {
  *   // called for all  coefficients
  *   void operator() (Type& value);
  * };
  * @endcode
  *
  * @note visitors offer automatic unrolling for small fixed size matrix.
  *
  * @sa setValue, setOnes(), setZeros()
  */
template<typename Derived>
template<typename Visitor>
inline void ArrayBase<Derived>::apply(Visitor& visitor)
{
  typedef typename hidden::VisitorSelector<Visitor, Derived>::Impl Impl;
  Impl::apply(this->asDerived(), visitor);
}

/*
 * start public function of the ArrayBase class using visitor
 */
/* @brief Set the value to all the Allocator */
template<typename Derived>
inline Derived& ArrayBase<Derived>::setValue(Type const& value)
{
  hidden::ValueApplier<Type> visitor(value);
  apply(visitor);
  return this->asDerived();
}

/* @brief Set the value one to all the Array */
template<typename Derived>
inline Derived& ArrayBase<Derived>::setOnes()
{
  hidden::ValueApplier<Type> visitor(Type(1));
  apply(visitor);
  return this->asDerived();
}

/* @brief Set the value one to all the Array */
template<typename Derived>
inline Derived& ArrayBase<Derived>::setZeros()
{
  hidden::ValueApplier<Type> visitor(Type(0));
  apply(visitor);
  return this->asDerived();
}

/* @brief Set the value one to all the Array */
template<typename Derived>
inline Derived& ArrayBase<Derived>::ones()
{ return setOnes();}

/* @brief Set the value one to all the Array */
template<typename Derived>
inline Derived& ArrayBase<Derived>::zeros()
{ return setZeros();}

template<typename Derived>
inline Derived& ArrayBase<Derived>::randUnif()
{
  hidden::RandUnifApplier<Type> visitor;
  apply(visitor);
  return this->asDerived();
}

/* set Gaussian random values to this */
template<typename Derived>
inline Derived& ArrayBase<Derived>::randGauss()
{
  hidden::RandGaussApplier<Type> visitor;
  apply(visitor);
  return this->asDerived();
}

/* set random values to this using a law given by the user */
template<typename Derived>
inline Derived& ArrayBase<Derived>::rand( Law::IUnivLaw<Type> const& law)
{
  hidden::RandApplier<Type> visitor(law);
  apply(visitor);
  return this->asDerived();
}


} // namespace STK

#endif /* STK_ARRAYBASEAPPLIER_H */
