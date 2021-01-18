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
 * created on: 13 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ExprBase.h
 *  @brief In this file we define the base class for Arrays and Expressions
 **/

#ifndef STK_EXPRBASEOPERATORS_H
#define STK_EXPRBASEOPERATORS_H

/// utility macro allowing to implement binary operators
#define IMPLEMENT_BINARY_OPERATOR(OPERATORNAME, BINARYOPERATORNAME) \
template<class Derived> \
template<typename Rhs> \
inline typename hidden::OperatorHelper<Derived, Rhs, Arrays::BINARYOPERATORNAME>::Result const \
ExprBase<Derived>::OPERATORNAME( ExprBase<Rhs> const& other) const \
{ return hidden::OperatorHelper<Derived, Rhs, Arrays::BINARYOPERATORNAME>::Selector::run(this->asDerived(), other.asDerived());}


namespace STK
{
IMPLEMENT_BINARY_OPERATOR(operator==,equalOp_)
IMPLEMENT_BINARY_OPERATOR(operator!=,notEqualOp_)
IMPLEMENT_BINARY_OPERATOR(operator>,greaterThanOp_)
IMPLEMENT_BINARY_OPERATOR(operator<,lessThanOp_)
IMPLEMENT_BINARY_OPERATOR(operator>=,greaterThanOrEqualOp_)
IMPLEMENT_BINARY_OPERATOR(operator<=,lessThanOrEqualOp_)

IMPLEMENT_BINARY_OPERATOR(operator+,sumOp_)
IMPLEMENT_BINARY_OPERATOR(operator-,differenceOp_)
IMPLEMENT_BINARY_OPERATOR(prod,productOp_)
IMPLEMENT_BINARY_OPERATOR(operator/,divisionOp_)
IMPLEMENT_BINARY_OPERATOR(operator%,moduloOp_)

IMPLEMENT_BINARY_OPERATOR(min,minOp_)
IMPLEMENT_BINARY_OPERATOR(max,maxOp_)

IMPLEMENT_BINARY_OPERATOR(operator&&,logicalAndOp_)
IMPLEMENT_BINARY_OPERATOR(operator||,logicalOrOp_)

IMPLEMENT_BINARY_OPERATOR(operator&,bitwiseAndOp_)
IMPLEMENT_BINARY_OPERATOR(operator|,bitwiseOrOp_)
IMPLEMENT_BINARY_OPERATOR(operator^,bitwiseXorOp_)

} // namespace STK

#undef IMPLEMENT_BINARY_OPERATOR

#endif /* STK_EXPRBASEOPERATORS_H */
