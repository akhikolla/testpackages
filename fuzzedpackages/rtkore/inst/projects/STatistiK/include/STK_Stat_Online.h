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
 * Project:  stkpp::STatistiK::StatDesc
 * Purpose:  Compute elementary 1D statistics for all variables.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Stat_Online.h
 *  @brief This file contain the definition and implementation of the Online classes
 **/

#ifndef STK_STAT_ONLINE_H
#define STK_STAT_ONLINE_H

#include <STKernel/include/STK_Real.h>
#include <STKernel/include/STK_Range.h>
#include <Sdk/include/STK_StaticAssert.h>

namespace STK
{
namespace Stat
{

template < class Array, class Type> struct Online;
/** @ingroup StatDesc
 *  @brief Computation online of the statistics of a Real set of variable.
 **/
template < class Array>
struct Online<Array, Real>
{
  typedef typename Array::Type Type_;
  enum
  {
    value_ = hidden::isSame<Type_, Real>::value_
  };
  /** default constructor*/
  inline Online(): mean_(), variance_(), iter_(0)
  { STK_STATIC_ASSERT(value_, YOU_CANNOT_USED_THIS_TYPE_OF_DATA_WITH_THIS_OBJECT);
    release();
  }
  /** constructor for one dimensional arrays */
  inline Online(Range const& range): mean_(range), variance_(range), iter_(0)
  { STK_STATIC_ASSERT(value_, YOU_CANNOT_USED_THIS_TYPE_OF_DATA_WITH_THIS_OBJECT);
    release();
  }
  /** constructor for two dimensional arrays */
  inline Online(Range const& rows, Range const& cols): mean_(rows, cols), variance_(rows, cols), iter_(0)
  { STK_STATIC_ASSERT(value_, YOU_CANNOT_USED_THIS_TYPE_OF_DATA_WITH_THIS_OBJECT);
    release();
  }
  /** copy constructor */
  inline Online( Online const& stat): mean_(stat.mean_), variance_(stat.variance_), iter_(stat.iter_)
  {}
  /** @return the computed online mean */
  Array const& mean() const { return mean_;}
  /** @return the computed online variance */
  Array variance() const
  { return iter_ == 0 ? STK::Arithmetic<Real>::infinity() : variance_/iter_;}
  /** initialize one dimensional arrays */
  inline void resize(Range const& range)
  { mean_.resize(range) = 0.; variance_.resize(range) = 0.; iter_ =0;}
  /** initialize two dimensional arrays */
  inline void resize(Range const& rows, Range const& cols)
  { mean_.resize(rows, cols) = 0.; variance_.resize(rows, cols) = 0.; iter_ =0;}
  /** release the computed parameters */
  inline void release() { mean_ = 0.; variance_ = 0.; iter_ = 0;}
  /** update the parameters using the current estimated parameters
   *  @see https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm
   *  @param x the current value
   **/
  inline void update(Array const& x)
  {
    iter_++;
    Array delta = x - mean_;
    mean_ += delta/iter_;
    variance_ = variance_ + delta.prod(x - mean_);
  }
  /** overwrite the statistics with other.
   *  @param other the statistics to copy
   **/
  inline Online& operator=( Online const& other)
  {
    mean_ = other.mean_;
    variance_ = other.variance_;
    iter_  = other.iter_;
    return *this;
  }
  /** mean of the parameters */
  Array mean_;
  /** on line variance */
  Array variance_;
  /** number of stored values */
  int iter_;
};

/** @ingroup StatDesc
 *  @brief Computation online of the statistics of a Real set of variable.
 **/
template <>
struct Online<Real, Real>
{
  /** default constructor*/
  Online(): mean_(0.), variance_(0.), iter_(0) {}
  /** copy constructor */
  Online( Online const& stat): mean_(stat.mean_), variance_(stat.variance_), iter_(stat.iter_)
  {}
  /** @return the computed online mean */
  Real const& mean() const { return mean_;}
  /** @return the computed online variance */
  Real variance() const
  { return iter_ == 0 ? STK::Arithmetic<Real>::infinity() : variance_/iter_;}
  /** release the computed parameters */
  inline void release() { mean_ = 0.; variance_ = 0.; iter_ = 0;}
  /** update the parameters using the current estimated parameters
   *  @see https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm
   *  @param value the current value of the variable
   **/
  inline void update(Real const& value)
  {
    iter_++;
    Real delta = value - mean_;
    mean_ += delta/iter_;
    variance_ = variance_ + delta*(value - mean_);
  }
  /** overwrite the statistics with other.
   *  @param other the online statistics to copy
   **/
  inline Online& operator=( Online const& other)
  {
    mean_     = other.mean_;
    variance_ = other.variance_;
    iter_     = other.iter_;
    return *this;
  }
  /** on line mean of the variable */
  Real mean_;
  /** on line variance times n of the variable*/
  Real variance_;
  /** number of stored values */
  int iter_;
};

}  // namespace Stat

}  // namespace STK

#endif /*STK_STAT_ONLINE_H */
