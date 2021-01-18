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
 * created on: 21 jan. 2018
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_BinarySelector.h
 *  @brief In this file we implement the class allowing to disambiguate
 *  the use of an UnaryOperator or a BinaryOperator
 **/


#ifndef STK_BINARYSELECTOR_H
#define STK_BINARYSELECTOR_H


namespace STK
{

namespace hidden
{
/** @ingroup hidden
 *  @brief helper defining the return type of the selector involving expressions
 *  of the form expr + number and expr + other_expr
 **/
template<typename Lhs, typename Rhs, int operatorType_>
struct OperatorHelper;

/** @ingroup hidden
 *  @brief allow to disambiguate the case expr + number from the case expr + other_expr */
template<typename Functor, typename Lhs, typename Rhs, bool isNumber>
struct OperatorImpl;

template<typename Functor, typename Lhs, typename Rhs>
struct OperatorImpl<Functor, Lhs, Rhs, true>
{
  typedef typename hidden::Traits<Lhs>::Type Ltype;
  typedef typename hidden::Traits<Rhs>::Type Rtype;

  typedef UnaryOperator<Functor, Lhs> Result;

  inline static Result const run( ExprBase<Lhs> const& lhs, ExprBase<Rhs> const& rhs)
  { return Result(lhs.asDerived(), Functor(rhs.elt()) );}
};

template<typename Functor, typename Lhs, typename Rhs>
struct OperatorImpl<Functor, Lhs, Rhs, false>
{
  typedef typename hidden::Traits<Lhs>::Type LType;
  typedef typename hidden::Traits<Rhs>::Type RType;

  typedef BinaryOperator<Functor, Lhs, Rhs> Result;

  inline static Result const run( ExprBase<Lhs> const& lhs, ExprBase<Rhs> const& rhs)
  { return Result(lhs.asDerived(), rhs.asDerived());}
};

/** @ingroup hidden
 *  @brief allow to disambiguate the case expr + number from the case expr + other_expr */
template<typename Lhs, typename Rhs, typename UnaryFunctor, typename BinaryFunctor>
struct OperatorSelector
{
  enum
  {
    isUnary_ = int(hidden::Traits<Rhs>::structure_) == int(Arrays::number_)
  };
  typedef typename hidden::Traits<Lhs>::Type LType;
  typedef typename hidden::Traits<Rhs>::Type RType;

  typedef UnaryOperator<UnaryFunctor, Lhs> UnaryResult;
  typedef BinaryOperator<BinaryFunctor, Lhs, Rhs> BinaryResult;

  typedef typename If< bool(isUnary_), UnaryFunctor, BinaryFunctor>::Result Functor;
  typedef typename If< bool(isUnary_), UnaryResult, BinaryResult>::Result Result;

  inline static Result const run( ExprBase<Lhs> const& lhs, ExprBase<Rhs> const& rhs)
  { return OperatorImpl<Functor, Lhs, Rhs, bool(isUnary_)>::run( lhs, rhs);}

  inline static BinaryResult const binaryRun( ExprBase<Lhs> const& lhs, ExprBase<Rhs> const& rhs)
  { return BinaryResult(lhs.asDerived(), rhs.asDerived());}

  inline UnaryResult const unaryRun( ExprBase<Lhs> const& lhs, ExprBase<Rhs> const& rhs)
  { return UnaryResult(lhs.asDerived(), UnaryFunctor(lhs.elt()));}
};

/** @ingroup hidden
 *  @brief specialization for operator==
 **/
template<typename Lhs, typename Rhs>
struct OperatorHelper<Lhs, Rhs, Arrays::equalOp_>
{
  typedef OperatorSelector< Lhs
                          , Rhs
                          , EqualWithOp<typename hidden::Traits<Lhs>::Type>
                          , EqualOp<typename hidden::Traits<Rhs>::Type, typename hidden::Traits<Rhs>::Type>
                          > Selector;
  typedef typename Selector::Result Result;
};
/** @ingroup hidden
 *  @brief specialization for operator!=
 **/
template<typename Lhs, typename Rhs>
struct OperatorHelper<Lhs, Rhs, Arrays::notEqualOp_>
{
  typedef OperatorSelector< Lhs
                          , Rhs
                          , NotEqualWithOp<typename hidden::Traits<Lhs>::Type>
                          , NotEqualOp<typename hidden::Traits<Rhs>::Type, typename hidden::Traits<Rhs>::Type>
                          > Selector;
  typedef typename Selector::Result Result;
};
/** @ingroup hidden
 *  @brief specialization for operator>
 **/
template<typename Lhs, typename Rhs>
struct OperatorHelper<Lhs, Rhs, Arrays::greaterThanOp_>
{
  typedef OperatorSelector< Lhs
                          , Rhs
                          , GreaterThanOp<typename hidden::Traits<Lhs>::Type>
                          , GreaterOp<typename hidden::Traits<Rhs>::Type, typename hidden::Traits<Rhs>::Type>
                          > Selector;
  typedef typename Selector::Result Result;
};
/** @ingroup hidden
 *  @brief specialization for operator<
 **/
template<typename Lhs, typename Rhs>
struct OperatorHelper<Lhs, Rhs, Arrays::lessThanOp_>
{
  typedef OperatorSelector< Lhs
                          , Rhs
                          , LessThanOp<typename hidden::Traits<Lhs>::Type>
                          , LessOp<typename hidden::Traits<Rhs>::Type, typename hidden::Traits<Rhs>::Type>
                          > Selector;
  typedef typename Selector::Result Result;
};
/** @ingroup hidden
 *  @brief specialization for operator>=
 **/
template<typename Lhs, typename Rhs>
struct OperatorHelper<Lhs, Rhs, Arrays::greaterThanOrEqualOp_>
{
  typedef OperatorSelector< Lhs
                          , Rhs
                          , GeqThanOp<typename hidden::Traits<Lhs>::Type>
                          , GeqOp<typename hidden::Traits<Rhs>::Type, typename hidden::Traits<Rhs>::Type>
                          > Selector;
  typedef typename Selector::Result Result;
};
/** @ingroup hidden
 *  @brief specialization for operator<=
 **/
template<typename Lhs, typename Rhs>
struct OperatorHelper<Lhs, Rhs, Arrays::lessThanOrEqualOp_>
{
  typedef OperatorSelector< Lhs
                          , Rhs
                          , LeqThanOp<typename hidden::Traits<Lhs>::Type>
                          , LeqOp<typename hidden::Traits<Rhs>::Type, typename hidden::Traits<Rhs>::Type>
                          > Selector;
  typedef typename Selector::Result Result;
};
/** @ingroup hidden
 *  @brief specialization for operator+
 **/
template<typename Lhs, typename Rhs>
struct OperatorHelper<Lhs, Rhs, Arrays::sumOp_>
{
  typedef OperatorSelector< Lhs
                          , Rhs
                          , SumWithOp<typename hidden::Traits<Lhs>::Type>
                          , SumOp<typename hidden::Traits<Rhs>::Type, typename hidden::Traits<Rhs>::Type>
                          > Selector;
  typedef typename Selector::Result Result;
};
/** @ingroup hidden
 *  @brief specialization for operator-
 **/
template<typename Lhs, typename Rhs>
struct OperatorHelper<Lhs, Rhs, Arrays::differenceOp_>
{
 typedef OperatorSelector< Lhs
                         , Rhs
                         , DifferenceWithOp<typename hidden::Traits<Lhs>::Type>
                         , DifferenceOp<typename hidden::Traits<Rhs>::Type, typename hidden::Traits<Rhs>::Type>
                         > Selector;
 typedef typename Selector::Result Result;
};
/** @ingroup hidden
 *  @brief specialization for operator*
 **/
template<typename Lhs, typename Rhs>
struct OperatorHelper<Lhs, Rhs, Arrays::productOp_>
{
 typedef OperatorSelector< Lhs
                         , Rhs
                         , ProductWithOp<typename hidden::Traits<Lhs>::Type>
                         , ProductOp<typename hidden::Traits<Rhs>::Type, typename hidden::Traits<Rhs>::Type>
                         > Selector;
 typedef typename Selector::Result Result;
};
/** @ingroup hidden
 *  @brief specialization for operator/
 **/
template<typename Lhs, typename Rhs>
struct OperatorHelper<Lhs, Rhs, Arrays::divisionOp_>
{
 typedef OperatorSelector< Lhs
                         , Rhs
                         , DivisionWithOp<typename hidden::Traits<Lhs>::Type>
                         , DivisionOp<typename hidden::Traits<Rhs>::Type, typename hidden::Traits<Rhs>::Type>
                         > Selector;
 typedef typename Selector::Result Result;
};
/** @ingroup hidden
 *  @brief specialization for operator%
 **/
template<typename Lhs, typename Rhs>
struct OperatorHelper<Lhs, Rhs, Arrays::moduloOp_>
{
 typedef OperatorSelector< Lhs
                         , Rhs
                         , ModuloWithOp<typename hidden::Traits<Lhs>::Type>
                         , ModuloOp<typename hidden::Traits<Rhs>::Type, typename hidden::Traits<Rhs>::Type>
                         > Selector;
 typedef typename Selector::Result Result;
};
/** @ingroup hidden
 *  @brief specialization for operator min
 **/
template<typename Lhs, typename Rhs>
struct OperatorHelper<Lhs, Rhs, Arrays::minOp_>
{
 typedef OperatorSelector< Lhs
                         , Rhs
                         , MinWithOp<typename hidden::Traits<Lhs>::Type>
                         , MinOp<typename hidden::Traits<Rhs>::Type, typename hidden::Traits<Rhs>::Type>
                         > Selector;
 typedef typename Selector::Result Result;
};
/** @ingroup hidden
 *  @brief specialization for operator max
 **/
template<typename Lhs, typename Rhs>
struct OperatorHelper<Lhs, Rhs, Arrays::maxOp_>
{
 typedef OperatorSelector< Lhs
                         , Rhs
                         , MaxWithOp<typename hidden::Traits<Lhs>::Type>
                         , MaxOp<typename hidden::Traits<Rhs>::Type, typename hidden::Traits<Rhs>::Type>
                         > Selector;
 typedef typename Selector::Result Result;
};
/** @ingroup hidden
 *  @brief specialization for operator&&
 **/
template<typename Lhs, typename Rhs>
struct OperatorHelper<Lhs, Rhs, Arrays::logicalAndOp_>
{
 typedef OperatorSelector< Lhs
                         , Rhs
                         , LogicalAndWithOp<typename hidden::Traits<Lhs>::Type>
                         , LogicalAndOp<typename hidden::Traits<Rhs>::Type, typename hidden::Traits<Rhs>::Type>
                         > Selector;
 typedef typename Selector::Result Result;
};
/** @ingroup hidden
 *  @brief specialization for operator||
 **/
template<typename Lhs, typename Rhs>
struct OperatorHelper<Lhs, Rhs, Arrays::logicalOrOp_>
{
 typedef OperatorSelector< Lhs
                         , Rhs
                         , LogicalOrWithOp<typename hidden::Traits<Lhs>::Type>
                         , LogicalOrOp<typename hidden::Traits<Rhs>::Type, typename hidden::Traits<Rhs>::Type>
                         > Selector;
 typedef typename Selector::Result Result;
};
/** @ingroup hidden
 *  @brief specialization for operator&
 **/
template<typename Lhs, typename Rhs>
struct OperatorHelper<Lhs, Rhs, Arrays::bitwiseAndOp_>
{
 typedef OperatorSelector< Lhs
                         , Rhs
                         , BitwiseAndWithOp<typename hidden::Traits<Lhs>::Type>
                         , BitwiseAndOp<typename hidden::Traits<Rhs>::Type, typename hidden::Traits<Rhs>::Type>
                         > Selector;
 typedef typename Selector::Result Result;
};
/** @ingroup hidden
 *  @brief specialization for operator|
 **/
template<typename Lhs, typename Rhs>
struct OperatorHelper<Lhs, Rhs, Arrays::bitwiseOrOp_>
{
 typedef OperatorSelector< Lhs
                         , Rhs
                         , BitwiseOrWithOp<typename hidden::Traits<Lhs>::Type>
                         , BitwiseOrOp<typename hidden::Traits<Rhs>::Type, typename hidden::Traits<Rhs>::Type>
                         > Selector;
 typedef typename Selector::Result Result;
};
/** @ingroup hidden
 *  @brief specialization for operator^
 **/
template<typename Lhs, typename Rhs>
struct OperatorHelper<Lhs, Rhs, Arrays::bitwiseXorOp_>
{
 typedef OperatorSelector< Lhs
                         , Rhs
                         , BitwiseXorWithOp<typename hidden::Traits<Lhs>::Type>
                         , BitwiseXorOp<typename hidden::Traits<Rhs>::Type, typename hidden::Traits<Rhs>::Type>
                         > Selector;
 typedef typename Selector::Result Result;
};

} // namespace hidden


} // namespace STK

#endif /* STK_BINARYSELECTOR_H */
