/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff

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

/** @file STK_ExprBaseVisitor.h
 *  @brief In this file we define the Visitors for ExprBase.
 **/

#ifndef STK_EXPRBASEVISITOR_H
#define STK_EXPRBASEVISITOR_H

#include "visitors/STK_Visitors.h"
#include "visitors/STK_SlicingVisitors.h"

namespace STK
{

/** Run the visitor @a visitor to the whole coefficients of the array.
  * The template parameter @a Visitor is the type of the visitor and provides
  * the following interface:
  * @code
  * struct MyVisitor {
  *   // called for all  coefficients
  *   void operator() (Type const& value, Index i, Index j);
  * };
  * @endcode
  *
  * @note visitors offer automatic unrolling for small fixed size matrix.
  *
  * @sa minElt, maxElt
  */
template<typename Derived>
template<typename Visitor>
inline typename Visitor::ConstReturnType ExprBase<Derived>::visit(Visitor& visitor) const
{
  typedef typename hidden::VisitorSelector<Visitor, Derived>::Impl Impl;
  Impl::run(this->asDerived(), visitor);
  return visitor.result();
}

/* count the number of not-zero values in the expression */
template<typename Derived>
inline int ExprBase<Derived>::count() const
{
  hidden::CountVisitor<Type> visitor;
  return visit(visitor);
}

/* return true if one element is not zero */
template<typename Derived>
inline bool const ExprBase<Derived>::any() const
{
  hidden::AnyVisitor<Type> visitor;
  return visit(visitor);
}

/* count the number of not-zero values in the expression */
template<typename Derived>
inline bool const ExprBase<Derived>::all() const
{
  hidden::AllVisitor<Type> visitor;
  return visit(visitor);
}

/* count the number values in the expression */
template<typename Derived>
inline int ExprBase<Derived>::nbAvailableValues() const
{ return isFinite().count();}

// general min elt
template<typename Derived>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::minElt( int& row, int& col) const
{
  hidden::MinEltVisitor<Type> visitor;
  visit(visitor);
  row = visitor.row_;
  col = visitor.col_;
  return visitor.result();
}

template<typename Derived>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::minEltSafe( int& row, int& col) const
{
  typedef hidden::MinEltSafeVisitor<Type> Visitor;
  Visitor visitor;
  visit(visitor);
  row = visitor.row_;
  col = visitor.col_;
  return visitor.result();
}
// general max elt
template<typename Derived>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::maxElt( int& row, int& col) const
{
  typedef hidden::MaxEltVisitor<Type> Visitor;
  Visitor visitor;
  visit(visitor);
  row = visitor.row_;
  col = visitor.col_;
  return visitor.result();
}
//
template<typename Derived>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::maxEltSafe( int& row, int& col) const
{
  typedef hidden::MaxEltSafeVisitor<Type> Visitor;
  Visitor visitor;
  visit(visitor);
  row = visitor.row_;
  col = visitor.col_;
  return visitor.result();
}
// min elt with one index
template<typename Derived>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::minElt( int& idx) const
{
  typedef hidden::MinEltVisitor<Type> Visitor;
  Visitor visitor;
  visit(visitor);
  idx = hidden::GetIdx<Visitor, hidden::Traits<Derived>::structure_ >::idx(visitor);
  return visitor.result();
}
template<typename Derived>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::minEltSafe( int& idx) const
{
  typedef hidden::MinEltSafeVisitor<Type> Visitor;
  Visitor visitor;
  visit(visitor);
  idx = hidden::GetIdx<Visitor, hidden::Traits<Derived>::structure_ >::idx(visitor);
  return visitor.result();
}
// max elt with one index
template<typename Derived>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::maxElt( int& idx) const
{
  typedef hidden::MaxEltVisitor<Type> Visitor;
  Visitor visitor;
  visit(visitor);
  idx = hidden::GetIdx<Visitor, hidden::Traits<Derived>::structure_ >::idx(visitor);
  return visitor.result();
}
template<typename Derived>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::maxEltSafe( int& idx) const
{
  typedef hidden::MaxEltSafeVisitor<Type> Visitor;
  Visitor visitor;
  visit(visitor);
  idx = hidden::GetIdx<Visitor, hidden::Traits<Derived>::structure_ >::idx(visitor);
  return visitor.result();
}
// min without index
template<typename Derived>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::minElt() const
{
  typedef hidden::MinVisitor<Type> Visitor;
  Visitor visitor;
  return visit(visitor);
}
template<typename Derived>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::minEltSafe() const
{
  typedef hidden::MinSafeVisitor<Type> Visitor;
  Visitor visitor;
  return visit(visitor);
}
// max elt without index
template<typename Derived>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::maxElt() const
{
  typedef hidden::MaxVisitor<Type> Visitor;
  Visitor visitor;
  return visit(visitor);
}
template<typename Derived>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::maxEltSafe() const
{
  typedef hidden::MaxSafeVisitor<Type> Visitor;
  Visitor visitor;
  return visit(visitor);
}

/* sum the values of all the array */
template<typename Derived>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::sum() const
{
  hidden::SumVisitor<Type> visitor;
  return visit(visitor);
}
/* sum safely the values of all the array */
template<typename Derived>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::sumSafe() const
{
  hidden::SumVisitor<Type> visitor;
  return safe().visit(visitor);
}

//--------- Start result with return type
/* @return the norm of this*/
template<typename Derived>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::norm() const
{ return Type(std::sqrt(norm2()));}
template<typename Derived>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::normSafe() const
{ return static_cast<Type>(std::sqrt(safe().norm2()));}
/* @return the square norm of this*/
template<typename Derived>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::norm2() const
{ return square().sum();}
template<typename Derived>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::norm2Safe() const
{ return safe().square().sum();}
/* @return the norm of this*/
template<typename Derived>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::normInf() const
{ return abs().maxElt();}

/* average the values of all the array */
template<typename Derived>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::mean() const
{
  hidden::MeanVisitor<Type> visitor;
  return visit(visitor);
}
/* sum safely the values of all the array */
template<typename Derived>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::meanSafe() const
{
  hidden::MeanSafeVisitor<Type> visitor;
  return visit(visitor);
}

/* compute the variance of all the array */
template<typename Derived>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::variance() const
{
  Type mu = mean();
  return (*this-mu).square().mean();
}
/* compute the variance of all the array */
template<typename Derived>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::varianceSafe() const
{
  Type mu = meanSafe();
  return (*this-mu).square().meanSafe();
}
/* compute the variance with given mean of all the elements of this*/
template<typename Derived>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::variance(Type const& mean) const
{ return (*this-mean).square().mean();}
template<typename Derived>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::varianceSafe(Type const& mean) const
{ return (*this-mean).square().meanSafe();}

/* @return the weighted sum of all the elements of this using a Visitor*/
template<typename Derived>
template<typename Rhs>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::wsum(ExprBase<Rhs> const& weights) const
{
  STK_STATIC_ASSERT_POINT_OR_VECTOR_ONLY(Derived);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Rhs);
  return dot(weights);
}
template<typename Derived>
template<typename Rhs>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::wsumSafe(ExprBase<Rhs> const& weights) const
{
  STK_STATIC_ASSERT_POINT_OR_VECTOR_ONLY(Derived);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Rhs);
  return dotSafe(weights);
}
/* @return the norm of this*/
template<typename Derived>
template<typename Rhs>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::wnorm(ExprBase<Rhs> const& weights) const
{
  STK_STATIC_ASSERT_POINT_OR_VECTOR_ONLY(Derived);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Rhs);
  return static_cast<Type>(std::sqrt(wnorm2(weights)));
}
template<typename Derived>
template<typename Rhs>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::wnormSafe(ExprBase<Rhs> const& weights) const
{
  STK_STATIC_ASSERT_POINT_OR_VECTOR_ONLY(Derived);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Rhs);
  return static_cast<Type>(std::sqrt(wnorm2Safe(weights)));
}
/* @return the square norm of this*/
template<typename Derived>
template<typename Rhs>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::wnorm2(ExprBase<Rhs> const& weights) const
{
  STK_STATIC_ASSERT_POINT_OR_VECTOR_ONLY(Derived);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Rhs);
  return square().dot(weights);
}
template<typename Derived>
template<typename Rhs>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::wnorm2Safe(ExprBase<Rhs> const& weights) const
{
  STK_STATIC_ASSERT_POINT_OR_VECTOR_ONLY(Derived);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Rhs);
  return square().dotSafe(weights);
}

/* @return the weighted mean */
template<typename Derived>
template<typename Rhs>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::wmean(ExprBase<Rhs> const& weights) const
{
  STK_STATIC_ASSERT_POINT_OR_VECTOR_ONLY(Derived);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Rhs);
  Type size = weights.sum();
  if (size <= 0 || !STK::isFinite(size)) return Arithmetic<Type>::NA();
  return wsum(weights)/size;
}

template<typename Derived>
template<typename Rhs>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::wmeanSafe(ExprBase<Rhs> const& weights) const
{
  STK_STATIC_ASSERT_POINT_OR_VECTOR_ONLY(Derived);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Rhs);
  Type size = weights.sumSafe();
  if (size <= 0) return Arithmetic<Type>::NA();
  return wsumSafe(weights)/size;
}

/* @return the variance of all the elements of this using a Visitor*/
template<typename Derived>
template<typename Rhs>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::wvariance(ExprBase<Rhs> const& weights) const
{
  STK_STATIC_ASSERT_POINT_OR_VECTOR_ONLY(Derived);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Rhs);
  Type mean = wmean(weights);
  if (!STK::isFinite(mean)) return Arithmetic<Type>::NA();
  return (*this-mean).square().wmean(weights);
}
template<typename Derived>
template<typename Rhs>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::wvarianceSafe(ExprBase<Rhs> const& weights) const
{
  STK_STATIC_ASSERT_POINT_OR_VECTOR_ONLY(Derived);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Rhs);
  Type mean = wmeanSafe(weights);
  if (!STK::isFinite(mean)) return Arithmetic<Type>::NA();
  return (*this-mean).square().wmeanSafe(weights);
}

/* @return the variance with given mean of all the elements of this*/
template<typename Derived>
template<typename Rhs>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::wvariance(Type const& mean, ExprBase<Rhs> const& weights) const
{
  STK_STATIC_ASSERT_POINT_OR_VECTOR_ONLY(Derived);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Rhs);
  return (*this-mean).square().wmean(weights);
}
template<typename Derived>
template<typename Rhs>
inline typename hidden::Traits<Derived>::Type const ExprBase<Derived>::wvarianceSafe(Type const& mean, ExprBase<Rhs> const& weights) const
{
  STK_STATIC_ASSERT_POINT_OR_VECTOR_ONLY(Derived);
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Rhs);
  return (*this-mean).square().wmeanSafe(weights);
}


} // namespace STK

#endif /* STK_EXPRBASEVISITOR_H */
