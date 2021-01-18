/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2017  Serge Iovleff, Université Lille 1, Inria

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
 * Project:  stkpp::STKernel
 * created on: Nov 22, 2017
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_STKernel_Util.h
 *  @brief In this file we define the utilities constant and classes for the STKernel project
 **/

#ifndef STK_STKERNEL_UTIL_H
#define STK_STKERNEL_UTIL_H

namespace STK
{

namespace hidden
{
/** @ingroup hidden, STKernel
 *  @brief The traits struct IteratorTraits must be specialized for any iterator
 *  derived from the base class IteratorBase.
 *  We use the type names defined by the STL for the iterator_traits class.
 *
 *  For example:
 *  @code
 *  template<typename Type>
 *  struct IteratorTraits
 *  {
 *    /// One of the iterator_tags types
 *    typedef std::random_access_iterator_tag  iterator_category;
 *    /// The type "pointed to" by the iterator.
 *    typedef Type        value_type;
 *    /// Distance between iterators is represented as this type.
 *    typedef int  difference_type;
 *    /// This type represents a pointer-to-value_type.
 *    typedef Type*   pointer;
 *    /// This type represents a reference-to-value_type.
 *    typedef Type& reference;
 *  };
 *  @endcode
 */
template <typename Derived> struct IteratorTraits;

} // namespace hidden

} // namespace STK


#endif /* STK_KERNEL_UTIL_H_ */
