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
 * Project:  stkpp::
 * created on: Sep 23, 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_MetaTemplate.h
 *  @brief This file contains generic meta-programming classes which are not
 *  specifically related to STK++.
 **/


#ifndef STK_METATEMPLATE_H
#define STK_METATEMPLATE_H

#include <STKernel/include/STK_Constants.h>

namespace STK
{
namespace hidden
{

struct TrueType {  enum { value = 1 }; };
struct FalseType { enum { value = 0 }; };

template<bool Condition, typename Then, typename Else>
struct If;

template<typename Then, typename Else>
struct If <true, Then, Else> { typedef Then Result; };

template<typename Then, typename Else>
struct If <false, Then, Else> { typedef Else Result; };

/* C/C++ fundamental Types */
template<typename T>
struct IsArithmetic { enum { yes = false, no = true}; };

/* specializations */
template<> struct IsArithmetic<float>         { enum { yes = true, no = false }; };
template<> struct IsArithmetic<double>        { enum { yes = true, no = false }; };
template<> struct IsArithmetic<long double>   { enum { yes = true, no = false }; };

template<> struct IsArithmetic<bool>          { enum { yes = true, no = false }; };
template<> struct IsArithmetic<char>          { enum { yes = true, no = false }; };
template<> struct IsArithmetic<signed char>   { enum { yes = true, no = false }; };
template<> struct IsArithmetic<unsigned char> { enum { yes = true, no = false }; };
template<> struct IsArithmetic<signed short>  { enum { yes = true, no = false }; };
template<> struct IsArithmetic<unsigned short>{ enum { yes = true, no = false }; };
template<> struct IsArithmetic<signed int>    { enum { yes = true, no = false }; };
template<> struct IsArithmetic<unsigned int>  { enum { yes = true, no = false }; };
template<> struct IsArithmetic<signed long>   { enum { yes = true, no = false }; };
template<> struct IsArithmetic<unsigned long> { enum { yes = true, no = false }; };
//template<> struct IsArithmetic<signed long long>   { enum { yes = true, no = false }; };
//template<> struct IsArithmetic<unsigned long long> { enum { yes = true, no = false }; };

/* C/C++ fundamental Types */
template<typename Type>
struct IsInt { enum { yes = false, no = true}; };

template<> struct IsInt<bool>          { enum { yes = true, no = false }; };
template<> struct IsInt<char>          { enum { yes = true, no = false }; };
template<> struct IsInt<signed char>   { enum { yes = true, no = false }; };
template<> struct IsInt<unsigned char> { enum { yes = true, no = false }; };
template<> struct IsInt<signed short>  { enum { yes = true, no = false }; };
template<> struct IsInt<unsigned short>{ enum { yes = true, no = false }; };
template<> struct IsInt<signed int>    { enum { yes = true, no = false }; };
template<> struct IsInt<unsigned int>  { enum { yes = true, no = false }; };
template<> struct IsInt<signed long>   { enum { yes = true, no = false }; };
template<> struct IsInt<unsigned long> { enum { yes = true, no = false }; };
//template<> struct IsInt<signed long long>   { enum { yes = true, no = false }; };
//template<> struct IsInt<unsigned long long> { enum { yes = true, no = false }; };


// remove const and const& to typename
template<typename Type_> struct RemoveConst { typedef Type_ Type; };
template<typename Type_> struct RemoveConst<Type_ const>   { typedef typename RemoveConst<Type_>::Type Type; };
template<typename Type_> struct RemoveConst<Type_ const&>  { typedef typename RemoveConst<Type_>::Type Type; };



/** check if T and U are of the same type. */
template<typename T, typename U> struct isSame { enum { value_ = 0 }; };
template<typename T> struct isSame<T,T> { enum { value_ = 1 }; };

/** check if T and U are equal */
template<int M, int N> struct IsEqual { enum { value_ = (M==N) }; };

/** @ingroup hidden
  * Convenient struct to Promote the result Type of some binary functors.
  */
template<typename Type1, typename Type2>
struct Promote
{ typedef typename If<(sizeof(Type1) > sizeof(Type2)), Type1, Type2>::Result result_type;};
/** @ingroup STKernel
  * Specialization when we have the same type.
  */
template<typename Type>
struct Promote<Type, Type>
{ typedef Type result_type;};

/** @ingroup hidden
  * Convenient structure for computing the product of two template integer parameters
  * without overflow.
  */
template<int Size1, int Size2> struct ProductSizeRowsBySizeCols;

// easy part
template<>
struct ProductSizeRowsBySizeCols<1, 1>                    { enum { prod_ = 1};};
template<>
struct ProductSizeRowsBySizeCols<UnknownSize, 1>          { enum { prod_ = UnknownSize};};
template<>
struct ProductSizeRowsBySizeCols<1, UnknownSize>          { enum { prod_ = UnknownSize};};
template<>
struct ProductSizeRowsBySizeCols<UnknownSize, UnknownSize>{ enum { prod_ = UnknownSize};};
template<int Size1>
struct ProductSizeRowsBySizeCols<Size1, 1>                { enum { prod_ = Size1};};
template<int Size2>
struct ProductSizeRowsBySizeCols<1, Size2>                { enum { prod_ = Size2};};
template<int Size1>
struct ProductSizeRowsBySizeCols<Size1, UnknownSize>      { enum { prod_ = UnknownSize};};
template<int Size2>
struct ProductSizeRowsBySizeCols<UnknownSize, Size2>      { enum { prod_ = UnknownSize};};

template<bool isSmallSize1, bool isSmallSize2, int Size1, int Size2>
struct ProductSizesHelper { enum { prod_ =  UnknownSize};};
template<int Size1, int Size2>
struct ProductSizesHelper<false, false, Size1, Size2> { enum { prod_ =  Size1 * Size2};};

/** @ingroup hidden
 *  @brief Compute the product Size1*Size2
 *  trick as some compilers (g++ 4 and 5 at least) complain about overflow otherwise
 **/
template<int Size1, int Size2>
struct ProductSizeRowsBySizeCols { enum { prod_ = ProductSizesHelper< Size1 >= SqrtUnknownSize
                                                                    , Size2 >= SqrtUnknownSize
                                                                    , Size1
                                                                    , Size2>::prod_};};


}// namespace hidden

} // namespace STK

#endif /* STK_METATEMPLATE_H_ */
