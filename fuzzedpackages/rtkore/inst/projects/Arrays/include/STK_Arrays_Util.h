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
 * created on: 17 févr. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Arrays_Util.h
 *  @brief In this file we define utilities functions and enum for the Array classes.
 **/


#ifndef STK_ARRAY_UTIL_H
#define STK_ARRAY_UTIL_H

#include <STKernel/include/STK_Range.h>

namespace STK
{

namespace Arrays
{

/** @ingroup Arrays
 *  Intrinsic dimension of the container : 0D, 1D, 2D, 3D or 4D. 0D is for scalar
 **/
enum Dimension
{
  _0D_ = 0, ///< a single scalar have no dimension
  _1D_ = 1,
  _2D_ = 2,
  _3D_ = 3,
  _4D_ = 4
};
/** @ingroup Arrays
 *  Define the Storage Orientation of the container
 **/
enum Orientation
{
  by_row_ =0,  ///< storage by row
  by_col_ =1   ///< storage by column
};

/** @ingroup Arrays
 *  Define the different type of Array that can be handle by STK++
 **/
enum Storage
{
  dense_ =1,  ///< dense matrix/vector/array/expression
  sparse_=0   ///< sparse matrix/vector/array/expression
};

/**  @ingroup Arrays
 *   structures of Arrays that can be handled by STK++
 **/
enum Structure
{
  array2D_ =0 ,        ///< general matrix/array/expression
  square_,             ///< square matrix/array/expression
  diagonal_,           ///< diagonal matrix/array/expression
  lower_triangular_,   ///< lower triangular matrix/array/expression
  upper_triangular_,   ///< upper triangular matrix/array/expression
  symmetric_,          ///< symmetric matrix/array/expression
  lower_symmetric_,    ///< lower symmetric matrix/array/expression
  upper_symmetric_,    ///< upper symmetric matrix/array/expression
  vector_,             ///< column oriented vector/array/expression
  point_,              ///< row oriented vector/array/expression
  number_,             ///< (1,1) matrix/vector/array/expression (like a number)
  expression_          ///< An expression that will be evaluated further
};

/** @ingroup Arrays
 *  @brief Kind of operators leading to a BinaryOperator or UnaryOperator.
 *
 *  This enum list the operators leading to binary or unary operators.
 *  If both the left-hand side and the right-hand side are expressions
 *  then we get a BinaryOperator while if the left-hand side is an expression
 *  and right hand side is a number then we get an UnaryOperator.
 *  For example, the return type of matrix1+matrix2 is a
 *  BinaryOperator. The return type of number+number is an UnaryOperator.
 **/
enum BinaryOperatorType
{
  equalOp_,              ///< operator==
  notEqualOp_,           ///< operator!=
  greaterThanOp_,        ///< operator>
  lessThanOp_,           ///< operator<
  greaterThanOrEqualOp_, ///< operator>=
  lessThanOrEqualOp_,    ///< operator<=

  sumOp_,                ///< operator+
  differenceOp_,         ///< operator-
  productOp_,            ///< operator*
  divisionOp_,           ///< operator/
  moduloOp_,             ///< operator%

  minOp_,                ///< min operator
  maxOp_,                ///< max operator

  logicalAndOp_,         ///< operator&&
  logicalOrOp_,          ///< operator||

  bitwiseAndOp_,         ///< operator&
  bitwiseOrOp_,          ///< operator|
  bitwiseXorOp_          ///< operator^
};

/** @ingroup Arrays
 *  Kind of operands in a BinaryOperator.
 **/
enum BinaryOpKind
{
  binary_op_0D_  = 0,     ///< both operand are number_
  binary_op_1D_  = 1,     ///< both operand are vector or point or diagonal
  binary_op_2D_  = 2,     ///< both operand are array2d or square
  binary_op_Diag_2D_,     ///< left operand is diagonal, right operand is 2D
  binary_op_Diag_UpTri_,
  binary_op_Diag_LowTri_,
  binary_op_Diag_Sym_,
  binary_op_Diag_UpSym_,
  binary_op_Diag_LowSym_,
  binary_op_2D_Diag_,     ///< left operand is 2D, right operand is diagonal
  binary_op_2D_UpTri_,
  binary_op_2D_LowTri_,
  binary_op_2D_Sym_,
  binary_op_2D_UpSym_,
  binary_op_2D_LowSym_,
  binary_op_UpTri_2D_,       ///< left operand is upper triangular, right operand is 2D
  binary_op_UpTri_Diag_,
  binary_op_UpTri_UpTri_,
  binary_op_UpTri_LowTri_,
  binary_op_UpTri_Sym_,
  binary_op_UpTri_UpSym_,
  binary_op_UpTri_LowSym_,
  binary_op_LowTri_2D_,      ///< left operand is lower triangular, right operand is 2D
  binary_op_LowTri_Diag_,
  binary_op_LowTri_UpTri_,
  binary_op_LowTri_LowTri_,
  binary_op_LowTri_Sym_,
  binary_op_LowTri_UpSym_,
  binary_op_LowTri_LowSym_,
  binary_op_Sym_2D_,    ///< left operand is symmetric, right operand is 2D
  binary_op_Sym_Diag_,
  binary_op_Sym_UpTri_,
  binary_op_Sym_LowTri_,
  binary_op_Sym_Sym_,
  binary_op_Sym_UpSym_,
  binary_op_Sym_LowSym_,
  binary_op_UpSym_2D_,    ///< left operand is upper symmetric, right operand is 2D
  binary_op_UpSym_Diag_,
  binary_op_UpSym_UpTri_,
  binary_op_UpSym_LowTri_,
  binary_op_UpSym_Sym_,
  binary_op_UpSym_UpSym_,
  binary_op_UpSym_LowSym_,
  binary_op_LowSym_2D_,   ///< left operand is lower symmetric, right operand is 2D
  binary_op_LowSym_Diag_,
  binary_op_LowSym_UpTri_,
  binary_op_LowSym_LowTri_,
  binary_op_LowSym_Sym_,
  binary_op_LowSym_UpSym_,
  binary_op_LowSym_LowSym_

};

/** @ingroup Arrays
 *  @brief Allow to disambiguate which rows() or cols() methods must be used
 *  which array to use when calling rows(), cols(), range() in Unary, Binary, Reshape,... operators ?
 **/
enum RangeOpUse
{
  useLhsSize_,      ///< use lhs.rows() (resp. lhs.cols())  in order to get rows() (resp. cols())
  useRhsSize_,      ///< use rhs.rows() (resp. rhs.cols())  in order to get rows() (resp. cols())
  useLhsOtherSize_, ///< use lhs.cols() in order to get rows() or lhs.rows() in order to get cols()
  useRhsOtherSize_, ///< use rhs.cols() in order to get rows() or rhs.rows() in order to get cols()
  useLhsRange_,     ///< use lhs.range() in order to get range()
  useRhsRange_      ///< use rhs.range() in order to get range()
};

/** @ingroup Arrays
 *  @return n+m, where n is the first number such that m < 2^n.
 *  @param m the size of the container
 **/
inline int evalSizeCapacity(int m)
{
  int n = 0;
  for (int k=1; k <= m; k <<= 1) {n++;}
  return(m+n);
}


/** @ingroup Arrays
 *  @return range of size n+m, where n is the first number such that m < 2^n.
 *  @param I the range of the container
 **/
template<int Size_>
TRange<Size_> evalRangeCapacity(TRange<Size_> const& I)
{
  int n = 0;
  for (int k=1; k <= I.size(); k <<= 1){ n++;}
  return TRange<Size_>(I.begin(),I.size() + n);
}

/** @ingroup Arrays
 *  convert an Arrays::Structure to a String.
 *  @param type the type of Structure we want to convert
 *  @return the string associated to this type.
 **/
inline std::string structureToString( Structure const& type)
{
  if (type == array2D_)          return String(_T("array2D"));
  if (type == square_)           return String(_T("square"));
  if (type == diagonal_)         return String(_T("diagonal"));
  if (type == lower_triangular_) return String(_T("lower_triangular"));
  if (type == upper_triangular_) return String(_T("upper_triangular"));
  if (type == symmetric_)        return String(_T("symmetric"));
  if (type == lower_symmetric_)  return String(_T("lower_symmetric"));
  if (type == upper_symmetric_)  return String(_T("upper_symmetric"));
  if (type == vector_)           return String(_T("vector"));
  if (type == point_)            return String(_T("point"));
  if (type == number_)           return String(_T("number"));
  if (type == expression_)       return String(_T("expression"));
  return String(_T("unknown"));
}


} // namespace Arrays

} // namespace STK



#endif /* STK_ARRAY_UTIL_H */
