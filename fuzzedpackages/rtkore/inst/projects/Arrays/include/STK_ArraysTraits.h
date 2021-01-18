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
 * created on: 14 déc. 2011
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ArraysTraits.h
 *  @brief In this file we define the main traits class we use for the
 *  STK++ Containers.
 **/


#ifndef STK_ARRAYSTRAITS_H
#define STK_ARRAYSTRAITS_H

namespace STK
{

namespace hidden
{
/** @ingroup hidden, Arrays
 *  @brief The traits struct Traits must be specialized for any
 *  container derived from the interface classes STK::ITContainer, STK::ITContainer1D
 *  and STK::ITContainer2D.
 *
 *  The enumerated values and type names defined in this trait struct are as
 *  follow (example taken from the STK::Array1D class)
 *  @code
 *  template<class Type_, int Size_>
 *  struct Traits< Array1D<Type_, Size_> >
 *  {
 *     typedef Array1D<Type_, 1> Row;
 *     typedef Array1D<Type_, Size_> Col;
 *     typedef Array1D<Type_, UnknownSize> SubRow;
 *     typedef Array1D<Type_, UnknownSize> SubCol;
 *     typedef Array1D<Type_, UnknownSize> SubArray;
 *     typedef Array1D<Type_, UnknownSize> SubVector;
 *
 *     typedef Type_ Type;
 *     typedef typename RemoveConst<Type>::Type const& ConstReturnType;
 *
 *     enum
 *     {
 *       structure_ = Arrays::vector_,
 *       orient_    = Arrays::by_col_,
 *       size_      = Size_,
 *       sizeCols_  = 1,
 *       sizeRows_  = Size_,
 *       storage_   = Arrays::dense_ // always dense
 *     };
 *
 *     // optional
 *     typedef RandomIterator1D<Array1D<Type_, Size_> > Iterator;
 *     typedef ConstRandomIterator1D<Array1D<Type_, Size_> > ConstIterator;
 *
 *     typedef std::reverse_iterator<Iterator> ReverseIterator;
 *     typedef std::reverse_iterator<ConstIterator> ConstReverseIterator;
 *  };
 *  @endcode
 *  @sa STK::Array1D, STK::Array2D, STK::CArray, STK::List1D
 */
template <typename Derived> struct Traits;

/** @ingroup hidden, Array
 *  @brief Traits class to get the correct returned Structure, Type, allocator,...
 *  of operator*. This traits struct is used by the functors classes operating
 *  on the STK::Array2D, STK::CArray,... classes.
 *  @note Impossible cases are tracked in ArrayByArrayProduct class.
 **/
template<typename Lhs, typename Rhs, int LStructure_, int RStructure_>
struct ProductTraits;

/** @ingroup hidden
 *  Utility class that will select the type of operator to apply.
 *  The result can be either a number if the data are in a vector or a point,
 *  or a vector if the data are in a matrix
 **/
template<typename Derived, template<class> class Functor>
struct FunctorTraits;

} // hidden

} // namespace STK

#endif /* STK_ARRAYSTRAITS_H */
