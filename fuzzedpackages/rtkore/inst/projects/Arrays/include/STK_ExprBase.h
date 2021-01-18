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

#ifndef STK_EXPRBASE_H
#define STK_EXPRBASE_H

#include "STK_ITContainer.h"

/// utility macro allowing to define binary operators
#define DEFINE_BINARY_OPERATOR(OPERATORNAME, BINARYOPERATORNAME) \
template<typename Rhs> \
typename hidden::OperatorHelper<Derived, Rhs, Arrays::BINARYOPERATORNAME>::Result const \
OPERATORNAME( ExprBase<Rhs> const& other) const;

/// utility macro allowing to construct unary operators
#define MAKE_UNARY_OPERATOR_NOARG(FUNCTION, FUNCTOR) \
  inline UnaryOperator<FUNCTOR<Type>, Derived> FUNCTION() const \
  { return UnaryOperator<FUNCTOR<Type>, Derived>(this->asDerived()); }

/// utility macro allowing to construct unary operators with one argument
#define MAKE_UNARY_OPERATOR_1ARG(FUNCTION, FUNCTOR) \
  inline UnaryOperator<FUNCTOR<Type>, Derived> FUNCTION(Type const& value) const \
  { return UnaryOperator<FUNCTOR<Type>, Derived>(this->asDerived(), FUNCTOR<Type>(value)); }

/// utility macro allowing to construct reshape operators
#define MAKE_RESHAPE_OPERATOR(OPERATOR, SHAPE) \
  inline OPERATOR##Operator< Derived> const SHAPE() const \
  { return OPERATOR##Operator< Derived>(this->asDerived()); }


// forward declarations
namespace STK
{
template<class Derived> class ExprBase;
template<class Derived> class ArrayBase;
//<class Derived, class Rhs> struct  ProductProductType;
template<class Derived> class  ArrayInitializer;
namespace hidden
{
template<typename Lhs, typename Rhs, int LStructure_, int RStructure_> struct ProductSelector;
}

} // namespace STK

#include <STKernel/include/STK_Functors.h>
#include <STatistiK/include/STK_Law_Functors.h>

#include "products/STK_ProductOperators.h"
#include "products/STK_DotProduct.h"
#include "operators/STK_TransposeOperator.h"
#include "operators/STK_ReshapeOperators.h"
#include "operators/STK_UnaryOperators.h"
#include "operators/STK_BinaryOperators.h"
#include "operators/STK_BinarySelector.h"

namespace STK
{

/** @ingroup Arrays
 *  @brief base class for template evaluation expressions and visitors.
 *
 * This class is the base that is inherited by all objects (matrix, vector,
 * point and expression). The common API for these objects is contained in
 * this class.
 *
 *  @tparam Derived is the derived type, e.g., a matrix, vector, type or
 *  an expression.
 **/

template<class Derived>
class ExprBase: public ITContainer<Derived, hidden::Traits<Derived>::structure_>
{
  public:
    enum
    {
      structure_ = hidden::Traits<Derived>::structure_,
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };
    typedef ITContainer<Derived, structure_> Base;
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef typename hidden::Traits<Derived>::ConstReturnType ConstReturnType;

  protected:
    /** Default constructor */
    inline ExprBase(): Base() {}
    /** destructor */
    inline ~ExprBase() {}

  public:
    //--------------
    // Start Visitors
    /** @brief Visit the container using a constant visitor
     *  @param visitor the visitor to run
     **/
    template<typename Visitor>
    typename Visitor::ConstReturnType visit(Visitor& visitor) const;
    /** @brief compute the value of non-zero element in an expression.
     *  For example
     *  @code
     *  (a > 0).count();
     *  @endcode
     *  compute the value of positive element of the array @c a.
     *  @return the value of non-zero element in the expression.*/
    int count() const;
    /** @brief check if there is any non-zero element in an expression.
     * For example
     *  @code
     *  (a > 0).any();
     *  @endcode
     *  will return @c true if there exists positive values in the expression @c a.
     *  @return @c true if at least one element is not zero in the expression, @c false otherwise.*/
    bool const any() const;
    /** @brief check if all the elements in an expression are not zero.
     *  For example
     *  @code
     *  (a > 0).all();
     *  @endcode
     *  will return @c true if all the elements in the expression @c a are positive.
     *   @return @c true if all the elements are not zero in the expression, @c false otherwise.*/
    bool const all() const;
    /** @return the value of available values in the array (not count NA values).*/
    int nbAvailableValues() const;

    /** @return the minimum of all elements of this using a Visitor
      * and puts in (row, col) its location.
      * @sa maxElt(int,int), visitor(), minElt()
      */
    Type const minElt( int& row, int& col) const;
    /** @return the minimum of all elements of @c *this which are not NA values
     *  using a Visitor and puts in (row, col) its location.
      * @sa maxEltSafe(int, int), visitor(), minElt()
      */
    Type const minEltSafe( int& row, int& col) const;
    /** @return the maximum of all elements of this using a Visitor
      * and puts in (row, col) its location.
      * @sa minElt(), visitor()
      */
    Type const maxElt( int& row, int& col) const;
    /** @return the maximum of all elements of this using a Visitor
      * and puts in (row, col) its location.
      * @sa minElt(), visitor()
      */
    Type const maxEltSafe( int& row, int& col) const;
    /** @return the minimum of all elements of this using a Visitor
      * and puts in  pos its location.
      * @note Have to be used for col-vector, row-vector or diagonal matrix only.
      * @sa maxElt(), visitor()
      */
    Type const minElt( int& pos) const;
    /** @return the minimum of all elements of this using a Visitor
      * and puts in  pos its location.
      * @note Have to be used for col-vector, row-vector or diagonal matrix only.
      * @sa maxElt(), visitor()
      */
    Type const minEltSafe( int& pos) const;
    /** @return the maximum of all elements of this using a Visitor
      * and puts in pos its location.
      * @note Have to be used for col-vector, row-vector or diagonal matrix only.
      * @sa minElt(), visitor()
      */
    Type const maxElt( int& pos) const;
    /** @return the maximum of all elements of this using a Visitor
      * and puts in pos its location.
      * @note Have to be used for col-vector, row-vector or diagonal matrix only.
      * @sa minElt(), visitor()
      */
    Type const maxEltSafe( int& pos) const;
    /** @return the minimum of all elements of this using a Visitor.
      * @sa maxElt(row, col), visitor()
      */
    Type const minElt() const;
    /** @return the minimum of all elements of this using a Visitor
      * @sa maxElt(row, pos), visitor()
      */
    Type const minEltSafe() const;
    /** @return the maximum of all elements of this using a Visitor.
      * @sa minElt(row, col), visitor()
      */
    Type const maxElt() const;
    /** @return the maximum of all elements of this using a Visitor
      * @sa minElt(row, pos), visitor()
      */
    Type const maxEltSafe() const;

    /** @return the sum of all the elements of this using a Visitor*/
    Type const sum() const;
    /** @return safely the sum of all the elements of this */
    Type const sumSafe() const;

    /** @return the norm of this*/
    Type const norm() const;
    /** @return the norm of this*/
    Type const normSafe() const;
    /** @return the square norm of this*/
    Type const norm2() const;
    /** @return the square norm of this*/
    Type const norm2Safe() const;
    /** @return the infinite norm of this*/
    Type const normInf() const;

    /** @return the mean of all the elements of this*/
    Type const mean() const;
    /** @return safely the mean of all the elements of this*/
    Type const meanSafe() const;
    /** @return the variance of all the elements of this*/
    Type const variance() const;
    /** @return the variance of all the elements of this*/
    Type const varianceSafe() const;
    /** @return the variance with given mean of all the elements of this*/
    Type const variance(Type const& mean) const;
    /** @return safely the variance with given mean of all the elements of this*/
    Type const varianceSafe(Type const& mean) const;

    /** @return the weighted sum of all the elements of this using a Visitor
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    Type const wsum(ExprBase<Rhs> const& weights) const;
    /** @return the safely weighted sum of all the elements of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    Type const wsumSafe(ExprBase<Rhs> const& weights) const;
    /** @return the weighted norm of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    Type const wnorm(ExprBase<Rhs> const& weights) const;
    /** @return the weighted norm of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    Type const wnormSafe(ExprBase<Rhs> const& weights) const;
    /** @return the weighted square norm of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    Type const wnorm2(ExprBase<Rhs> const& weights) const;
    /** @return the weighted square norm of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    Type const wnorm2Safe(ExprBase<Rhs> const& weights) const;
    /** @return the weighted mean of all the elements of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    Type const wmean(ExprBase<Rhs> const& weights) const;
    /** @return the safely weighted mean of all the elements of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    Type const wmeanSafe(ExprBase<Rhs> const& weights) const;
    /** @return the weighted variance of all the elements of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    Type const wvariance(ExprBase<Rhs> const& weights) const;
    /** @return the weighted variance of all the elements of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    Type const wvarianceSafe(ExprBase<Rhs> const& weights) const;
    /** @return the variance with given mean of all the elements of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    Type const wvariance(Type const& mean, ExprBase<Rhs> const& weights) const;
    /** @return safely the variance with given mean of all the elements of this
     *  @note will only work with row-vectors or col-vectors
     **/
    template<typename Rhs>
    Type const wvarianceSafe(Type const& mean, ExprBase<Rhs> const& weights) const;
    // Visitors terminated
    //--------------------
    // BinaryOperators
    /** @return an expression with == operator of this and other.*/
    DEFINE_BINARY_OPERATOR(operator==,equalOp_)
    /** @return an expression with != operator of this and other.*/
    DEFINE_BINARY_OPERATOR(operator!=,notEqualOp_)
    /** @return an expression with > operator of this and other.*/
    DEFINE_BINARY_OPERATOR(operator>,greaterThanOp_)
    /** @return an expression with < operator of this and other*/
    DEFINE_BINARY_OPERATOR(operator<,lessThanOp_)
    /** @return an expression with >= operator of this and other.*/
    DEFINE_BINARY_OPERATOR(operator>=,greaterThanOrEqualOp_)
    /** @return an expression with <= operator of this and other.*/
    DEFINE_BINARY_OPERATOR(operator<=,lessThanOrEqualOp_)

    /** @return an expression with the addition of this and other. */
    DEFINE_BINARY_OPERATOR(operator+,sumOp_)
    /** @return an expression from the difference of this and  other. */
    DEFINE_BINARY_OPERATOR(operator-,differenceOp_)
    /** @return an expression with the product (element by element) of this and other.*/
    DEFINE_BINARY_OPERATOR(prod,productOp_)
    /** @return an expression with the quotient of this and other. */
    DEFINE_BINARY_OPERATOR(operator/,divisionOp_)
    /** @return an expression with the modulo of this and other. */
    DEFINE_BINARY_OPERATOR(operator%,moduloOp_)

    /** @return an expression with the min of this and other.*/
    DEFINE_BINARY_OPERATOR(min,minOp_)
    /** @return an expression with the max of this and other.*/
    DEFINE_BINARY_OPERATOR(max,maxOp_)

    /** @return an expression with the logical AND of this and other. */
    DEFINE_BINARY_OPERATOR(operator&&,logicalAndOp_)
    /** @return an expression with the logical OR of this and other. */
    DEFINE_BINARY_OPERATOR(operator||,logicalOrOp_)

    /** @return an expression with the bitwise AND of this and other. */
    DEFINE_BINARY_OPERATOR(operator&,bitwiseAndOp_)
    /** @return an expression with the bitwise OR of this and other. */
    DEFINE_BINARY_OPERATOR(operator|,bitwiseOrOp_)
    /** @return an expression with the bitwise XOR of this and other. */
    DEFINE_BINARY_OPERATOR(operator^,bitwiseXorOp_)
    // BinaryOperators terminated
    //--------------------
    // UnaryOperators
    /** @return an expression of the opposite of this */
    MAKE_UNARY_OPERATOR_NOARG(operator-, OppositeOp)
    /** @return which values of this is a NA value */
    MAKE_UNARY_OPERATOR_NOARG(isNA, IsNaOp)
    /** @return not of the values */
    MAKE_UNARY_OPERATOR_NOARG(neg, NegOp)
    /** @return which values of this are finite value */
    MAKE_UNARY_OPERATOR_NOARG(isFinite, IsFiniteOp)
    /** @return which values of this are infinite value */
    MAKE_UNARY_OPERATOR_NOARG(isInfinite, IsInfiniteOp)
    /** @return an expression of the absolute value of this */
    MAKE_UNARY_OPERATOR_NOARG(abs, AbsOp)
    /** @return an expression of the exponential of this. */
    MAKE_UNARY_OPERATOR_NOARG(exp, ExpOp)
    /** @return an expression of the logarithm of this. */
    MAKE_UNARY_OPERATOR_NOARG(log, LogOp)
    /** @return an expression of the square root of this. */
    MAKE_UNARY_OPERATOR_NOARG(sqrt, SqrtOp)
    /** @return an expression of the cosine of this. */
    MAKE_UNARY_OPERATOR_NOARG(cos, CosOp)
    /** @return an expression of the sine of this.  */
    MAKE_UNARY_OPERATOR_NOARG(sin, SinOp)
    /** @return an expression of the arc cosine of this. */
    MAKE_UNARY_OPERATOR_NOARG(acos, AcosOp)
    /** @return an expression of the arc sine of this. */
    MAKE_UNARY_OPERATOR_NOARG(asin, AsinOp)
    /** @return an expression of the tan of this. */
    MAKE_UNARY_OPERATOR_NOARG(tan, TanOp)
    /** @return an expression of the inverse of this.  */
    MAKE_UNARY_OPERATOR_NOARG(inverse, InverseOp)
    /** @return an expression of the square of this. */
    MAKE_UNARY_OPERATOR_NOARG(square, SquareOp)
    /** @return an expression of the cube of this. */
    MAKE_UNARY_OPERATOR_NOARG(cube, CubeOp)

    // boolean operations
    /** @return an expression of *this == value. */
    MAKE_UNARY_OPERATOR_1ARG(operator==, EqualWithOp)
    /** @return an expression of *this != value. */
    MAKE_UNARY_OPERATOR_1ARG(operator!=, NotEqualWithOp)
    /** @return an expression of *this > value. */
    MAKE_UNARY_OPERATOR_1ARG(operator>, GreaterThanOp)
    /** @return an expression of *this < value. */
    MAKE_UNARY_OPERATOR_1ARG(operator<, LessThanOp)
    /** @return an expression of *this <= value. */
    MAKE_UNARY_OPERATOR_1ARG(operator<=, LeqThanOp)
    /** @return the expression of *this >= value. */
    MAKE_UNARY_OPERATOR_1ARG(operator>=, GeqThanOp)

    /** @return an expression of the minimum of this and a value. */
    MAKE_UNARY_OPERATOR_1ARG(min, MinWithOp)
    /** @return an expression of the maximum of this and a value. */
    MAKE_UNARY_OPERATOR_1ARG(max, MaxWithOp)

    /** @return an expression of this + value */
    MAKE_UNARY_OPERATOR_1ARG(operator+, SumWithOp)
    /** @return an expression of this - value */
    MAKE_UNARY_OPERATOR_1ARG(operator-, DifferenceWithOp)
    /** @return an expression of this * value */
    MAKE_UNARY_OPERATOR_1ARG(operator*, ProductWithOp)
    /** @return an expression of this divided by a value */
    MAKE_UNARY_OPERATOR_1ARG(operator/, DivisionWithOp)
    /** @return an expression of this modulo a value */
    MAKE_UNARY_OPERATOR_1ARG(operator%, ModuloWithOp)

    /** @return a logical expression of this AND a value */
    MAKE_UNARY_OPERATOR_1ARG(operator&&, LogicalAndWithOp)
    /** @return a logical expression of this OR a value */
    MAKE_UNARY_OPERATOR_1ARG(operator||, LogicalOrWithOp)

    /** @return a bitwise expression of this AND a value */
    MAKE_UNARY_OPERATOR_1ARG(operator&, BitwiseAndWithOp)
    /** @return a bitwise expression of this OR a value */
    MAKE_UNARY_OPERATOR_1ARG(operator|, BitwiseOrWithOp)
    /** @return a bitwise expression of this OR a value */
    MAKE_UNARY_OPERATOR_1ARG(operator^, BitwiseXorWithOp)

    /** @return an expression of the power of this. */
    MAKE_UNARY_OPERATOR_1ARG(pow, PowOp)
    /** @return an expression of this divided by the value */
    MAKE_UNARY_OPERATOR_1ARG(safeInverse, SafeInverseOp)

    /** @return a safe value of this */
    inline UnaryOperator<SafeOp<Type>, Derived> const safe(Type const value = Type()) const
    { return UnaryOperator<SafeOp<Type>, Derived>(this->asDerived(), SafeOp<Type>(value)); }

    // friends
    /** @return an expression of value + this */
    inline friend UnaryOperator<SumWithOp<Type>, Derived> const
    operator+(Type const& value, ExprBase<Derived> const& other)
    { return other.asDerived() + value;}
    /** @return an expression of value - this */
    inline friend UnaryOperator<SubstractToOp<Type>, Derived > const
    operator-(Type const value, ExprBase<Derived> const& other)
    { return UnaryOperator<SubstractToOp<Type>, Derived>(other.asDerived(), SubstractToOp<Type>(value));}
    /** @return an expression of value * this */
    inline friend UnaryOperator< ProductWithOp<Type>, Derived> const
    operator*(Type const value, ExprBase<Derived> const& other)
    { return other.asDerived()*value; }

    /** @return an expression of *this with the  Type type casted to  OtherType. */
    template<typename OtherType>
    inline UnaryOperator<CastOp<Type, OtherType>, Derived> const cast() const
    { return UnaryOperator<CastOp<Type, OtherType>, Derived>(this->asDerived());}

    // pdf, lpdf, cdf, icdf
    /** compute pdf values to this using distribution law given by user */
    inline UnaryOperator<Law::PdfOp<Type>, Derived> pdf( Law::IUnivLaw<Type> const& law) const
    { return UnaryOperator<Law::PdfOp<Type>, Derived>(this->asDerived(), Law::PdfOp<Type>(law));}
    /** compute log-pdf values to this using distribution law given by user */
    inline UnaryOperator<Law::LogPdfOp<Type>, Derived> lpdf( Law::IUnivLaw<Type> const& law) const
    { return UnaryOperator<Law::LogPdfOp<Type>, Derived>(this->asDerived(), Law::LogPdfOp<Type>(law));}
    /** compute cumulative distribution function of this using distribution law given by user */
    inline UnaryOperator<Law::CdfOp<Type>, Derived> cdf( Law::IUnivLaw<Type> const& law) const
    { return UnaryOperator<Law::CdfOp<Type>, Derived>(this->asDerived(), Law::CdfOp<Type>(law));}
    /** compute log-cumulative distribution function of this using distribution law given by user */
    inline UnaryOperator<Law::LogCdfOp<Type>, Derived> lcdf( Law::IUnivLaw<Type> const& law) const
    { return UnaryOperator<Law::LogCdfOp<Type>, Derived>(this->asDerived(), Law::LogCdfOp<Type>(law));}
    /** compute complementary cumulative distribution function of this using distribution law given by user */
    inline UnaryOperator<Law::CdfcOp<Type>, Derived> cdfc( Law::IUnivLaw<Type> const& law) const
    { return UnaryOperator<Law::CdfcOp<Type>, Derived>(this->asDerived(), Law::CdfcOp<Type>(law));}
    /** compute complementary cumulative distribution function of this using distribution law given by user */
    inline UnaryOperator<Law::LogCdfcOp<Type>, Derived> lcdfc( Law::IUnivLaw<Type> const& law) const
    { return UnaryOperator<Law::LogCdfcOp<Type>, Derived>(this->asDerived(), Law::LogCdfcOp<Type>(law));}
    /** compute inverse cumulative distribution function using distribution law given by user */
    inline UnaryOperator<Law::IcdfOp<Type>, Derived> icdf( Law::IUnivLaw<Type> const& law) const
    { return UnaryOperator<Law::IcdfOp<Type>, Derived>(this->asDerived(), Law::IcdfOp<Type>(law));}

    // extension
    /** @return an expression of funct0 to this. */
    template< template<typename> class OtherOperator>
    inline UnaryOperator<OtherOperator<Type>, Derived> const funct0() const
    { return UnaryOperator<OtherOperator<Type>, Derived>(this->asDerived());}
    /** @return an expression of funct1 to this. */
    template< template<typename> class OtherOperator>
    inline UnaryOperator<OtherOperator<Type>, Derived> const funct1(Type const value) const
    { return UnaryOperator<OtherOperator<Type>, Derived>(this->asDerived(), OtherOperator<Type>(value));}

    // reshape operations
    /** @return the transposed expression of this. */
    MAKE_RESHAPE_OPERATOR(Transpose,transpose)
    /** @return this as a diagonal 1D expression (work only with vector/point/diagonal expressions). */
    MAKE_RESHAPE_OPERATOR(Diagonalize,diagonalize)
    /** @return the diagonal of this square expression (work only with square expressions). */
    MAKE_RESHAPE_OPERATOR(DiagonalGetter,getDiagonal)
    /** @return the upper triangular part of this expression. */
    MAKE_RESHAPE_OPERATOR(UpperTriangularize,upperTriangularize)
    /** @return the lower triangular part of this expression. */
    MAKE_RESHAPE_OPERATOR(LowerTriangularize,lowerTriangularize)
    /** @return this as a symmetric expression (work only with square expressions). */
    MAKE_RESHAPE_OPERATOR(Symmetrize,symmetrize)
    /** @return the upper part of this symmetric expression (work only with square expressions). */
    MAKE_RESHAPE_OPERATOR(UpperSymmetrize,upperSymmetrize)
    /** @return the lower part of this symmetric expression (work only with square expressions). */
    MAKE_RESHAPE_OPERATOR(LowerSymmetrize,lowerSymmetrize)

    // slice operators
    /** @return the sub-vector(I) of this. */
    template<int Size_>
    inline SubVectorOperator<Derived, Size_> const sub(TRange<Size_> const& I) const
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);
      return SubVectorOperator<Derived, Size_>(this->asDerived(), I);
    }

    /** @return the j-th column of this. */
    inline ColOperator<Derived> const col(int j) const
    { return ColOperator<Derived> (this->asDerived(), j);}
    /** @return the i-th row of this. */
    inline RowOperator<Derived> const row(int i) const
    { return RowOperator<Derived> (this->asDerived(), i);}
    /** @return the range J of columns of this. */
    template<int Size_>
    inline SubOperator<Derived, sizeRows_, Size_> const col(TRange<Size_> const& J) const
    { return SubOperator<Derived, sizeRows_, Size_> (this->asDerived(), this->rows(), J);}
    /** @return the range I of rows of this. */
    template<int Size_>
    inline SubOperator<Derived, Size_, sizeCols_> const row(TRange<Size_> const& I) const
    { return SubOperator<Derived, Size_, sizeCols_> (this->asDerived(), I, this->cols());}

    /** @return the sub-array(I,J) of this. */
    template<int SizeRows_, int SizeCols_>
    inline SubOperator<Derived, SizeRows_, SizeCols_> const sub(TRange<SizeRows_> const& I,TRange<SizeCols_> const& J) const
    {
      STK_STATIC_ASSERT_TWO_DIMENSIONS_ONLY(Derived);
      return SubOperator<Derived, SizeRows_, SizeCols_>(this->asDerived(), I, J);
    }

    // dot operations
    /** @returns the dot product of this with other.
     *  @sa norm2(), norm(), DotProduct
     */
    template<class Rhs>
    typename hidden::Promote<Type, typename Rhs::Type>::result_type const
    dot(ExprBase<Rhs> const& other) const;
    /** @returns the safe dot product of this with other.
      * @sa squaredNorm(), norm(), DotProduct
      */
    template<class Rhs>
    typename hidden::Promote<Type, typename Rhs::Type>::result_type const
    dotSafe(ExprBase<Rhs> const& other) const;

    /** @return the matrix multiplication of this with other.*/
    template<typename Rhs>
    typename hidden::ProductSelector<Derived, Rhs, hidden::Traits<Derived>::structure_, hidden::Traits<Rhs>::structure_>::ProductType const
    operator*( ExprBase<Rhs> const& other) const;
};

} // namespace STK

#undef DEFINE_BINARY_OPERATOR
#undef MAKE_UNARY_OPERATOR_NOARG
#undef MAKE_UNARY_OPERATOR_1ARG
#undef MAKE_RESHAPE_OPERATOR

#include "STK_ExprBaseVisitor.h"
#include "STK_ExprBaseDot.h"
#include "STK_ExprBaseProduct.h"
#include "STK_ExprBaseOperators.h"


#endif /* STK_EXPRBASE_H_ */
