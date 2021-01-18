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

    Contact: S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  stkpp::Functors
 * created on: 1 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Functors.h
 *  @brief In this file we implement the functors.
 **/


#ifndef STK_FUNCTORS_H
#define STK_FUNCTORS_H

#include <cmath>

namespace STK
{

namespace hidden
{

/** @ingroup hidden
 * The functor to use.  */
template< typename Functor, int NbParam_>
struct UsedFunctor;

/** @ingroup hidden
 *  produce the Nth parameter type of Functor.  If the type is note available
 *  because the parameter does not exits, generate an error using a "Void" class
 **/
template< typename Functor, int NbParam_>
struct FunctorParamTraits
{
  private:
    class Void
    {   class Private {};
      public:
        typedef Private Type;
    };
  public:
    typedef typename If< NbParam_ <= Functor::NbParam_
                       , UsedFunctor<Functor, NbParam_>
                       , Void>::Result::Type Type;
};

// specialization, until 4 parameters
template< typename Functor> struct UsedFunctor<Functor, 1>
{ typedef typename Functor::param1_type Type;};
template< typename Functor> struct UsedFunctor<Functor, 2>
{ typedef typename Functor::param2_type Type;};
template< typename Functor> struct UsedFunctor<Functor, 3>
{ typedef typename Functor::param3_type Type;};
template< typename Functor> struct UsedFunctor<Functor, 4>
{ typedef typename Functor::param4_type Type;};

/** @ingroup hidden
 *  produce the return type of the Functor.  */
template< typename Functor>
struct FunctorReturnTypeTraits
{ typedef typename Functor::result_type result_type;};


} // namespace hidden

/** @ingroup Functors
 *  @brief Template functor testing if a Char is an end of line.
 **/
struct TestEndOfLineOp
{
  enum { NbParam_ = 1 };
  typedef bool result_type;
  typedef Char param1_type ;
  inline TestEndOfLineOp(char* c): last_(c) {}

  inline result_type operator()(Char c) { *last_ = c; return (c == _T('\n'));}
  private:
    Char* last_;
};

//------------------------------------------------------------------------------
// binary functors:
//------------------------------------------------------------------------------
/** @ingroup Functors
  * @brief Template functor testing if a number is equal to an other number
  */
template<class Type1, class Type2 = Type1>
struct EqualOp
{
  enum { NbParam_ = 2 };
  typedef bool result_type;
  typedef typename hidden::RemoveConst<Type1>::Type param1_type ;
  typedef typename hidden::RemoveConst<Type2>::Type const& param2_type ;

  inline result_type operator()(param1_type const& value1, param2_type value2) const
  { return value1 == value2;}
};
/** @ingroup Functors
  * @brief Template functor testing if a number is not equal to an other number
  */
template<class Type1, class Type2 = Type1>
struct NotEqualOp
{
  enum { NbParam_ = 2 };
  typedef bool result_type;
  typedef typename hidden::RemoveConst<Type1>::Type param1_type ;
  typedef typename hidden::RemoveConst<Type2>::Type const& param2_type ;

  inline result_type operator()(param1_type const& value1, param2_type value2) const
  { return value1 != value2;}
};
/** @ingroup Functors
  * @brief Template functor testing if a number is less than an other number
  */
template<class Type1, class Type2 = Type1>
struct LessOp
{
  enum { NbParam_ = 2 };
  typedef bool result_type;
  typedef typename hidden::RemoveConst<Type1>::Type param1_type ;
  typedef typename hidden::RemoveConst<Type2>::Type const& param2_type ;

  inline result_type operator()(param1_type const& value1, param2_type value2) const
  { return value1 < value2;}
};
/** @ingroup Functors
  * @brief Template functor testing if a number is greater than an other number
  */
template<class Type1, class Type2 = Type1>
struct GreaterOp
{
  enum { NbParam_ = 2 };
  typedef bool result_type;
  typedef typename hidden::RemoveConst<Type1>::Type param1_type ;
  typedef typename hidden::RemoveConst<Type2>::Type const& param2_type ;

  inline result_type operator()(param1_type const& value1, param2_type value2) const
  { return value1 > value2;}
};
/** @ingroup Functors
  * @brief Template functor testing if a number is less or equal than an other number
  */
template<class Type1, class Type2 = Type1>
struct LeqOp
{
  enum { NbParam_ = 2 };
  typedef bool result_type;
  typedef typename hidden::RemoveConst<Type1>::Type param1_type ;
  typedef typename hidden::RemoveConst<Type2>::Type const& param2_type ;

  inline result_type operator()(param1_type const& value1, param2_type value2) const
  { return value1 <= value2;}
};
/** @ingroup Functors
  * @brief Template functor testing if a number is greater or equal than an other number
  */
template<class Type1, class Type2 = Type1>
struct GeqOp
{
  enum { NbParam_ = 2 };
  typedef bool result_type;
  typedef typename hidden::RemoveConst<Type1>::Type param1_type ;
  typedef typename hidden::RemoveConst<Type2>::Type const& param2_type ;

  inline result_type operator()(param1_type const& value1, param2_type value2) const
  { return value1 >= value2;}
};

/** @ingroup Functors
  * @brief Template functor which compute the sum of two numbers
  */
template<class Type1, class Type2 = Type1>
struct SumOp
{
  enum { NbParam_ = 2 };
  typedef typename hidden::Promote<Type1, Type2>::result_type result_type;
  typedef typename hidden::RemoveConst<Type1>::Type const& param1_type ;
  typedef typename hidden::RemoveConst<Type2>::Type const& param2_type ;

  inline result_type operator()(param1_type a, param2_type b) const { return a + b;}
};
/** @ingroup Functors
  * @brief Template functor which compute the difference of two numbers
  */
template<class Type1, class Type2 = Type1>
struct DifferenceOp
{
  enum { NbParam_ = 2 };
  typedef typename hidden::Promote<Type1, Type2>::result_type result_type;
  typedef typename hidden::RemoveConst<Type1>::Type const& param1_type ;
  typedef typename hidden::RemoveConst<Type2>::Type const& param2_type ;

  inline result_type operator()(param1_type a, param2_type b) const { return a - b;}
};
/** @ingroup Functors
  * @brief Template functor which compute the product of two numbers
  */
template<class Type1, class Type2 = Type1>
struct ProductOp
{
  enum { NbParam_ = 2 };
  typedef typename hidden::Promote<Type1, Type2>::result_type result_type;
  typedef typename hidden::RemoveConst<Type1>::Type const& param1_type ;
  typedef typename hidden::RemoveConst<Type2>::Type const& param2_type ;

  inline result_type operator()(param1_type a, param2_type b) const { return a * b;}
};
/** @ingroup Functors
  * @brief Template functor which compute the division of two numbers
  */
template<class Type1, class Type2 = Type1>
struct DivisionOp
{
  enum { NbParam_ = 2 };
  typedef typename hidden::Promote<Type1, Type2>::result_type result_type;
  typedef typename hidden::RemoveConst<Type1>::Type const& param1_type ;
  typedef typename hidden::RemoveConst<Type2>::Type const& param2_type ;

  inline result_type operator()(param1_type a, param2_type b) const { return a / b;}
};
/** @ingroup Functors
  * @brief Template functor which compute the modulo of two numbers
  */
template<class Type1, class Type2 = Type1>
struct ModuloOp
{
  enum { NbParam_ = 2 };
  typedef typename hidden::Promote<Type1, Type2>::result_type result_type;
  typedef typename hidden::RemoveConst<Type1>::Type const& param1_type ;
  typedef typename hidden::RemoveConst<Type2>::Type const& param2_type ;

  inline result_type operator()(param1_type a, param2_type b) const { return a % b;}
};
/** @ingroup Functors
  * @brief Template functor which compute the minimum of two numbers
  */
template<class Type1, class Type2 = Type1>
struct MinOp
{
  enum { NbParam_ = 2 };
  typedef typename hidden::Promote<Type1, Type2>::result_type result_type;
  typedef typename hidden::RemoveConst<Type1>::Type param1_type ;
  typedef typename hidden::RemoveConst<Type2>::Type const& param2_type ;

  inline result_type operator()(param1_type const&  value1, param2_type value2) const
  { return  (value1 < value2) ? value1: value2 ;}
};
/** @ingroup Functors
  * @brief Template functor which compute the maximum of two numbers
  */
template<class Type1, class Type2 = Type1>
struct MaxOp
{
  enum { NbParam_ = 2 };
  typedef typename hidden::Promote<Type1, Type2>::result_type result_type;
  typedef typename hidden::RemoveConst<Type1>::Type param1_type ;
  typedef typename hidden::RemoveConst<Type2>::Type const& param2_type ;

  inline result_type operator()(param1_type const&  value1, param2_type value2) const
  { return  (value1 < value2) ? value2: value1 ;}
};
/** @ingroup Functors
  * @brief Template functor which compute the logical and of two numbers
  */
template<class Type1, class Type2 = Type1>
struct LogicalAndOp
{
  enum { NbParam_ = 2 };
  typedef bool result_type;
  typedef typename hidden::RemoveConst<Type1>::Type const& param1_type ;
  typedef typename hidden::RemoveConst<Type2>::Type const& param2_type ;

  inline result_type operator()(param1_type a, param2_type b) const { return a && b;}
};
/** @ingroup Functors
  * @brief Template functor which compute the logical or of two numbers
  */
template<class Type1, class Type2 = Type1>
struct LogicalOrOp
{
  enum { NbParam_ = 2 };
  typedef bool result_type;
  typedef typename hidden::RemoveConst<Type1>::Type const& param1_type ;
  typedef typename hidden::RemoveConst<Type2>::Type const& param2_type ;

  inline result_type operator()(param1_type a, param2_type b) const { return a || b;}
};
/** @ingroup Functors
  * @brief Template functor which compute the bitwise and of two numbers
  */
template<class Type1, class Type2 = Type1>
struct BitwiseAndOp
{
  enum { NbParam_ = 2 };
  typedef typename hidden::Promote<Type1, Type2>::result_type result_type;
  typedef typename hidden::RemoveConst<Type1>::Type const& param1_type ;
  typedef typename hidden::RemoveConst<Type2>::Type const& param2_type ;

  inline result_type operator()(param1_type a, param2_type b) const { return a & b;}
};
/** @ingroup Functors
  * @brief Template functor which compute the bitwise or of two numbers
  */
template<class Type1, class Type2 = Type1>
struct BitwiseOrOp
{
  enum { NbParam_ = 2 };
  typedef typename hidden::Promote<Type1, Type2>::result_type result_type;
  typedef typename hidden::RemoveConst<Type1>::Type const& param1_type ;
  typedef typename hidden::RemoveConst<Type2>::Type const& param2_type ;

  inline result_type operator()(param1_type a, param2_type b) const { return a | b;}
};
/** @ingroup Functors
  * @brief Template functor which compute the bitwise xor of two numbers
  */
template<class Type1, class Type2 = Type1>
struct BitwiseXorOp
{
  enum { NbParam_ = 2 };
  typedef typename hidden::Promote<Type1, Type2>::result_type result_type;
  typedef typename hidden::RemoveConst<Type1>::Type const& param1_type ;
  typedef typename hidden::RemoveConst<Type2>::Type const& param2_type ;

  inline result_type operator()(param1_type a, param2_type b) const { return a ^ b;}
};

//------------------------------------------------------------------------------
// unary functors:
//------------------------------------------------------------------------------
/** @ingroup Functors
  * @brief Template functor testing if a number is a NA value
  */
template<class Type>
struct IsNaOp
{
  enum { NbParam_ = 1 };
  typedef bool result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline result_type operator()(param1_type const&  value) const
  { return Arithmetic<param1_type>::isNA(value);}
};
/** @ingroup Functors
  * @brief Template functor giving the not value
  */
template<class Type>
struct NegOp
{
  enum { NbParam_ = 1 };
  typedef bool result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline result_type operator()(param1_type const& value) const  { return !value;}
};
/** @ingroup Functors
  * @brief Template functor testing if a number is a finite value
  */
template<class Type>
struct IsFiniteOp
{
  enum { NbParam_ = 1 };
  typedef bool result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline result_type operator()(param1_type const& value) const
  { return Arithmetic<param1_type>::isFinite(value);}
};
/** @ingroup Functors
  * @brief Template functor testing if a number is an infinite value
  */
template<class Type>
struct IsInfiniteOp
{
  enum { NbParam_ = 1 };
  typedef bool result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type;

  inline result_type operator()(param1_type const& value1) const
  { return Arithmetic<param1_type>::isInfinite(value1);}
};
/** @ingroup Functors
  * @brief Template functor to test if a number is equal to a fixed value
  */
template<typename Type>
struct EqualWithOp
{
  enum { NbParam_ = 1 };
  typedef bool result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline EqualWithOp(EqualWithOp const& value): value_(value.value_) {}
  inline EqualWithOp(Type value): value_(value) {}

  inline result_type operator() (param1_type const& a) const { return a == value_;}
  private:
    param1_type const value_;
};
/** @ingroup Functors
  * @brief Template functor to test if a number is different than a fixed value
  */
template<typename Type>
struct NotEqualWithOp
{
  enum { NbParam_ = 1 };
  typedef bool result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline NotEqualWithOp(NotEqualWithOp const& value): value_(value.value_) {}
  inline NotEqualWithOp(Type value): value_(value) {}

  inline result_type operator() (param1_type const& a) const { return a != value_;}
  private:
    param1_type const value_;
};
/** @ingroup Functors
  * @brief Template functor to test if a number is less than a fixed value
  */
template<typename Type>
struct LessThanOp
{
  enum { NbParam_ = 1 };
  typedef bool result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type;

  inline LessThanOp(LessThanOp const& value): value_(value.value_) {}
  inline LessThanOp(Type value): value_(value) {}

  inline result_type operator() (param1_type const& a) const { return a < value_;}
  private:
    param1_type const value_;
};
/** @ingroup Functors
  * @brief Template functor to test if a number is greater than a fixed value
  */
template<typename Type>
struct GreaterThanOp
{
  enum { NbParam_ = 1 };
  typedef bool result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline GreaterThanOp(GreaterThanOp const& value): value_(value.value_) {}
  inline GreaterThanOp(Type value): value_(value) {}

  inline result_type operator() (param1_type const& a) const { return a > value_;}
  private:
    param1_type const value_;
};
/** @ingroup Functors
  * @brief Template functor to test if a number is less or equal than a fixed value
  */
template<typename Type>
struct LeqThanOp
{
  enum { NbParam_ = 1 };
  typedef bool result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline LeqThanOp(LeqThanOp const& value): value_(value.value_) {}
  inline LeqThanOp(Type value): value_(value) {}

  inline result_type operator() (param1_type const& a) const { return (a <= value_);}
  private:
    param1_type const value_;
};
/** @ingroup Functors
  * @brief Template functor to test if a number is greater or equal than a fixed value
  */
template<typename Type>
struct GeqThanOp
{
  enum { NbParam_ = 1 };
  typedef bool result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline GeqThanOp(GeqThanOp const& value): value_(value.value_) {}
  inline GeqThanOp(Type value): value_(value) {}

  inline result_type operator() (param1_type const& a) const { return a >= value_;}
  private:
    param1_type const value_;
};
/** @ingroup Functors
  * @brief Template functor to add a number to a fixed value
  */
template<typename Type>
struct SumWithOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline SumWithOp(SumWithOp const& value): value_(value.value_) {}
  inline SumWithOp(Type value): value_(value) {}

  inline result_type operator() (param1_type const& a) const { return a + value_;}
  private:
    param1_type const value_;
};
/** @ingroup Functors
  * @brief Template functor to add a number to a fixed value
  */
template<typename Type>
struct DifferenceWithOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline DifferenceWithOp(DifferenceWithOp const& value): value_(value.value_) {}
  inline DifferenceWithOp(Type value): value_(value) {}

  inline result_type operator() (param1_type const& a) const { return a - value_;}
  private:
    param1_type const value_;
};
/** @ingroup Functors
  * @brief Template functor to subtract a number to a fixed value
  */
template<typename Type>
struct SubstractToOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline SubstractToOp(SubstractToOp const& value): value_(value.value_) {}
  inline SubstractToOp(Type value): value_(value) {}

  inline result_type operator() (param1_type const& a) const { return value_ - a;}
  private:
    param1_type const value_;
};
/** @ingroup Functors
  * @brief Template functor computing the product of a number by a fixed value
  */
template<typename Type>
struct ProductWithOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline ProductWithOp(Type value): value_(value) {}
  inline ProductWithOp(ProductWithOp const& value): value_(value.value_) {}

  inline result_type operator() (param1_type const& a) const {return a * value_;}
  private:
    param1_type const value_;
};
/** @ingroup Functors
  * @brief Specialization for boolean
  */
template<>
struct ProductWithOp<bool>
{
  enum { NbParam_ = 1 };
  typedef bool const result_type;
  typedef bool const param1_type;

  inline ProductWithOp(bool value): value_(value) {}
  inline ProductWithOp( ProductWithOp const& value): value_(value.value_) {}

  inline result_type operator() (param1_type const& a) const {return a && value_;}
  private:
    bool const value_;
};
/** @ingroup Functors
 *  base class for DivisionWithOp functors
 **/
template<typename Type, bool IsInt>
struct DivisionWithBaseOp;

template<typename Type>
struct DivisionWithBaseOp<Type, false>
{
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;
  inline DivisionWithBaseOp(const DivisionWithBaseOp& value): value_(value.value_) {}
  inline DivisionWithBaseOp(param1_type const& value): value_( static_cast<Type>(1) / value) {}

  inline Type const operator() (param1_type const& a) const { return a * value_;}
  private:
    param1_type const value_;
};

template<typename Type>
struct DivisionWithBaseOp<Type, true>
{
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;
  inline DivisionWithBaseOp( DivisionWithBaseOp const& value): value_(value.value_) {}
  inline DivisionWithBaseOp( param1_type value): value_(value) {}

  inline Type const operator()( param1_type a) const { return a / value_;}
  private:
    param1_type const value_;
};
/** @ingroup Functors
  * @brief Template functor to divide a number by a fixed value
  *
  * This functor is used to implement the quotient of a matrix by
  * a number where the number type is not necessarily a floating point type.
  */
template<typename Type>
struct DivisionWithOp: public DivisionWithBaseOp<Type, hidden::IsInt<Type>::yes >
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline DivisionWithOp(param1_type const& value): DivisionWithBaseOp<Type, hidden::IsInt<Type>::yes >(value) {}
};
/** @ingroup Functors
  * @brief Template functor to compute the minimum with a fixed value
  */
template<typename Type>
struct ModuloWithOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline ModuloWithOp(Type const& value): value_(value) {}
  inline ModuloWithOp( ModuloWithOp const& value): value_(value.value_) {}

  inline result_type operator() (param1_type const& a) const {return a % value_;}
  private:
    param1_type const value_;
};
/** @ingroup Functors
  * @brief Template functor to compute the minimum between a fixed value
  */
template<typename Type>
struct MinWithOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline MinWithOp(Type const& value): value_(value) {}
  inline MinWithOp( MinWithOp const& value): value_(value.value_) {}

  inline result_type operator() (param1_type const& a) const {return std::min(a,value_);}
  private:
    param1_type const value_;
};
/** @ingroup Functors
  * @brief Template functor to compute the minimum with a fixed value
  */
template<typename Type>
struct MaxWithOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline MaxWithOp(Type const& value): value_(value) {}
  inline MaxWithOp( MaxWithOp const& value): value_(value.value_) {}

  inline result_type operator() (param1_type const& a) const {return std::max(a,value_);}
  private:
    param1_type const value_;
};
/** @ingroup Functors
  * @brief Template functor computing logical and with a fixed value
  */
template<typename Type>
struct LogicalAndWithOp
{
  enum { NbParam_ = 1 };
  typedef bool result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline LogicalAndWithOp(LogicalAndWithOp const& value): value_(value.value_) {}
  inline LogicalAndWithOp(Type value): value_(value) {}

  inline result_type operator() (param1_type const& a) const { return a && value_;}
  private:
    param1_type const value_;
};
/** @ingroup Functors
  * @brief Template functor computing logical and with a fixed value
  */
template<typename Type>
struct LogicalOrWithOp
{
  enum { NbParam_ = 1 };
  typedef bool result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline LogicalOrWithOp(LogicalOrWithOp const& value): value_(value.value_) {}
  inline LogicalOrWithOp(Type value): value_(value) {}

  inline result_type operator() (param1_type const& a) const { return a || value_;}
  private:
    param1_type const value_;
};
/** @ingroup Functors
  * @brief Template functor to compute the bitwise and with a fixed value
  */
template<typename Type>
struct BitwiseAndWithOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline BitwiseAndWithOp(BitwiseAndWithOp const& value): value_(value.value_) {}
  inline BitwiseAndWithOp(Type const& value): value_(value) {}

  inline result_type operator() (param1_type const& a) const {return a & value_;}
  private:
    param1_type const value_;
};
/** @ingroup Functors
  * @brief Template functor to compute the bitwise or with a fixed value
  */
template<typename Type>
struct BitwiseOrWithOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline BitwiseOrWithOp(BitwiseOrWithOp const& value): value_(value.value_) {}
  inline BitwiseOrWithOp(Type const& value): value_(value) {}

  inline result_type operator() (param1_type const& a) const {return a | value_;}
  private:
    param1_type const value_;
};
/** @ingroup Functors
  * @brief Template functor to compute the bitwise xor with a fixed value
  */
template<typename Type>
struct BitwiseXorWithOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline BitwiseXorWithOp(BitwiseXorWithOp const& value): value_(value.value_) {}
  inline BitwiseXorWithOp(Type const& value): value_(value) {}

  inline result_type operator() (param1_type const& a) const {return a ^ value_;}
  private:
    param1_type const value_;
};


// Functors without value_
/** @ingroup Functors
  * @brief Template functor which compute the opposite of a number
  */
template<class Type>
struct OppositeOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline result_type operator()(param1_type const& a) const { return -a;}
};
/** @ingroup Functors
  * @brief Template functor which return a default value if the value is NA
  */
template<class Type>
struct SafeOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline SafeOp(Type const value = Type()): value_(value) {}
  inline SafeOp( SafeOp const& value): value_(value.value_) {}

  inline result_type operator()(param1_type const& a) const
  { return Arithmetic<Type>::isFinite(a) ? a: value_;}
  private:
    Type const value_;
};
/** @ingroup Functors
  * @brief Template functor which compute safely the inverse of a number.
  * If the number is 0, return 0.
  */
template<typename Type>
struct SafeInverseOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline SafeInverseOp(Type const& tol = Arithmetic<Type>::epsilon()): tol_(tol) {}
  inline SafeInverseOp( SafeInverseOp const& value): tol_(value.tol_) {}

  inline result_type operator()(param1_type const& a) const
  { return (std::abs(a)>tol_) ? Type(1)/a: Type(0);}
  private:
    Type const tol_;
};
/** @ingroup Functors
  * @brief Template functor which compute the absolute value of a number
  */
template<class Type>
struct AbsOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline result_type operator()(param1_type const& a) const { return std::abs(a);}
};
/** @ingroup Functors
  * @brief Template functor which compute the exponential of a number
  */
template<class Type>
struct ExpOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline result_type operator()(param1_type const& a) const {return std::exp(a);}
};
/** @ingroup Functors
  * @brief Template functor which compute the logarithm of a number
  */
template<class Type>
struct LogOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline result_type operator()(param1_type const& a) const {return std::log(a);}
};
/** @ingroup Functors
  * @brief Template functor which compute the square root of a number
  */
template<typename Type> struct SqrtOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline result_type operator()(param1_type const& a) const { return std::sqrt(a);}
};
/** @ingroup Functors
  * @brief Template functor which compute the cosine of a number
  */
template<typename Type> struct CosOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline result_type operator() (param1_type const& a) const { return std::cos(a);}
};
/** @ingroup Functors
  * @brief Template functor which compute the sine of a number
  */
template<typename Type> struct SinOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline result_type operator() (param1_type const& a) const { return std::sin(a);}
};
/** @ingroup Functors
  * @brief Template functor which compute the tan of a number
  */
template<typename Type> struct TanOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline result_type operator() (param1_type const& a) const { return std::tan(a);}
};
/** @ingroup Functors
  * @brief Template functor which compute the arc cosine of a number
  */
template<typename Type> struct AcosOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline result_type operator() (param1_type const& a) const { return std::acos(a);}
};
/** @ingroup Functors
  * @brief Template functor which compute the arc sine of a number
  */
template<typename Type> struct AsinOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline result_type operator() (param1_type const& a) const { return std::asin(a);}
};
/** @ingroup Functors
  * @brief Template functor to raise a number to a power
  */
template<typename Type>
struct PowOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline PowOp(PowOp const& value): exponent_(value.exponent_) {}
  inline PowOp(Type exponent): exponent_(exponent) {}

  inline result_type operator()(param1_type const& a) const  { return std::pow(a, exponent_);}
  private:
    const Type exponent_;
};
/** @ingroup Functors
  * @brief Template functor which compute the inverse of a number
  */
template<typename Type>
struct InverseOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline result_type operator() (param1_type const& a) const { return Type(1)/a;}
};
/** @ingroup Functors
  * @brief Template functor which compute the square of a number
  */
template<typename Type>
struct SquareOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline result_type operator() (param1_type const& a) const { return a*a;}
};

/** @ingroup Functors
  * @brief Template functor which compute the cube of a number
  */
template<typename Type>
struct CubeOp
{
  enum { NbParam_ = 1 };
  typedef Type result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type;

  inline result_type operator() (param1_type const& a) const { return a*a*a;}
};
/** @ingroup Functors
  * @brief Template functor which cast a type to another type
  */
template<typename Type, typename OtherType>
struct CastOp
{
  enum { NbParam_ = 1 };
  typedef OtherType result_type;
  typedef typename hidden::RemoveConst<Type>::Type param1_type ;

  inline result_type operator()(param1_type const& a) const { return static_cast<OtherType>(a);}
};

} // namespace STK


#endif /* STK_FUNCTORS_H */
