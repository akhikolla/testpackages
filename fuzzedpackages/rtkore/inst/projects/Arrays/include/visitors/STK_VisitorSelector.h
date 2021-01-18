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

/** @file STK_VisitorSelector.h
 *  @brief In this file we implement the selector of the visitor and apply pattern.
 **/

#ifndef STK_VISITORSELECTOR_H
#define STK_VISITORSELECTOR_H

namespace STK
{

namespace hidden
{
// forward declaration
template<typename Visitor, typename Derived, int Structure_>
struct VisitorSelectorHelper;

/** @ingroup hidden
 *  @brief visitor selector. If @c Derived is a full two-dimensional array and
 *  the visitation can be unrolled,then we use directly the VisitorArrayImpl
 *  class in order to compute the result of the visitation.
 *
 *  Otherwise we delegate the search of the correct implementation to the
 *  VisitorSelectorHelper class.
 **/
template<typename Visitor, typename Derived>
struct VisitorSelector
{
  enum
  {
    structure_ = hidden::Traits<Derived>::structure_,
    orient_    = hidden::Traits<Derived>::orient_,
    sizeRows_  = hidden::Traits<Derived>::sizeRows_,
    sizeCols_  = hidden::Traits<Derived>::sizeCols_,
    storage_   = hidden::Traits<Derived>::storage_,
    unrollRows_ = (sizeRows_ < MaxUnrollSlice),
    unrollCols_ = (sizeCols_ < MaxUnrollSlice),
    is2D_       = (structure_ == (int)Arrays::array2D_ || structure_ == (int)Arrays::square_)
  };

  typedef typename VisitorSelectorHelper<Visitor, Derived, structure_>::Impl HelperImpl;
  typedef VisitorArrayUnrollImpl<Visitor, Derived, orient_, sizeRows_, sizeCols_> ArrayImpl;
  // the other cases will be take into account by the helper class
  typedef typename If<is2D_ && unrollRows_ && unrollCols_, ArrayImpl, HelperImpl >::Result Impl;
};


/** @ingroup hidden
 *  @brief Helper for the Visitor selector: allow to select the correct
 *  implementation to instantiate for arrays and square arrays when only one
 *  dimension (rows or columns but not both).
 *
 *  This selector will select the visitor method to use for 2D arrays (arrays_
 *  and square_).
 **/
template<typename Visitor, typename Derived, int Structure_>
struct VisitorSelectorHelper
{
  enum
  {
    orient_    = hidden::Traits<Derived>::orient_,
    sizeRows_  = hidden::Traits<Derived>::sizeRows_,
    sizeCols_  = hidden::Traits<Derived>::sizeCols_,
    unrollRows_ = (sizeRows_ < MaxUnroll),
    unrollCols_ = (sizeCols_ < MaxUnroll)
  };
  typedef VisitorArrayImpl<Visitor, Derived, sizeRows_, UnknownSize> RowUnrollImpl;
  typedef VisitorArrayImpl<Visitor, Derived, UnknownSize, sizeCols_> ColUnrollImpl;
  typedef VisitorArrayNoUnrollImpl<Visitor, Derived, orient_, UnknownSize, UnknownSize> NoUnrollImpl;
  typedef typename If<unrollCols_, ColUnrollImpl
                                 , typename If<unrollRows_, RowUnrollImpl, NoUnrollImpl>::Result >::Result Impl;
};

/** @ingroup hidden
 *  @brief specialization for the diagonal arrays */
template<typename Visitor, typename Derived>
struct VisitorSelectorHelper<Visitor, Derived, Arrays::diagonal_>
{
  enum
  {
    sizeRows_  = hidden::Traits<Derived>::sizeRows_,
    unrollRows_ = (sizeRows_ < MaxUnroll)
  };
  typedef VisitorDiagonalImpl<Visitor, Derived, sizeRows_> Unroll;
  typedef VisitorDiagonalImpl<Visitor, Derived, UnknownSize> Loop;

  typedef typename If<unrollRows_, Unroll, Loop  >::Result Impl;
};

/** @ingroup hidden
 *  @brief specialization for the vectors */
template<typename Visitor, typename Derived>
struct VisitorSelectorHelper<Visitor, Derived, Arrays::vector_>
{
  enum
  {
    sizeRows_  = hidden::Traits<Derived>::sizeRows_,
    unrollRows_ = (sizeRows_ < MaxUnroll)
  };
  typedef VisitorVectorImpl<Visitor, Derived, sizeRows_> Unroll;
  typedef VisitorVectorImpl<Visitor, Derived, UnknownSize> Loop;
  typedef typename If<unrollRows_, Unroll, Loop  >::Result Impl;
};

/** @ingroup hidden
 *  @brief specialization for the row vectors */
template<typename Visitor, typename Derived>
struct VisitorSelectorHelper<Visitor, Derived, Arrays::point_>
{
  enum
  {
    sizeCols_  = hidden::Traits<Derived>::sizeCols_,
    unrollCols_ = (sizeCols_ < MaxUnroll)
  };
  typedef VisitorPointImpl<Visitor, Derived, sizeCols_> Unroll;
  typedef VisitorPointImpl<Visitor, Derived, UnknownSize> Loop;

  typedef typename If<unrollCols_, Unroll, Loop  >::Result Impl;
};

/** @ingroup hidden
 *  @brief specialization for the upper triangular arrays */
template<typename Visitor, typename Derived>
struct VisitorSelectorHelper<Visitor, Derived, Arrays::upper_triangular_>
{
  enum
  { orient_    = hidden::Traits<Derived>::orient_};
  typedef VisitorUpperImpl<Visitor, Derived, orient_> Impl;
};

/** @ingroup hidden
 *  @brief specialization for the lower triangular arrays */
template<typename Visitor, typename Derived>
struct VisitorSelectorHelper<Visitor, Derived, Arrays::lower_triangular_>
{
  enum
  { orient_    = hidden::Traits<Derived>::orient_};
  typedef VisitorLowerImpl<Visitor, Derived, orient_> Impl;
};

/** @ingroup hidden
 *  @brief specialization for the numbers */
template<typename Visitor, typename Derived>
struct VisitorSelectorHelper<Visitor, Derived, Arrays::number_>
{ typedef VisitorNumberImpl<Visitor, Derived> Impl;};


} // namespace hidden

} // namespace STK

#endif /* STK_VISITORSIMPL_H */
