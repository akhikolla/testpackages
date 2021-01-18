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

/** @file STK_ArrayBase.h
 *  @brief In this file we define the base class for Arrays. Derived ArrayBase can be
 *  a lhs
 **/

#ifndef STK_ARRAYBASE_H
#define STK_ARRAYBASE_H

#include "STK_ExprBase.h"
#include "STatistiK/include/STK_Law_IUnivLaw.h"

/// utility macro allowing to construct unary operators
#define MAKE_RESHAPE_OPERATOR(OPERATOR, SHAPE) \
  inline OPERATOR##Operator< Derived> const SHAPE() const \
  { return OPERATOR##Operator< Derived>(this->asDerived()); } \
  inline OPERATOR##Accessor< Derived> SHAPE() \
  { return OPERATOR##Accessor< Derived>(this->asDerived()); }

namespace STK
{

/** @ingroup Arrays
 *  @brief base class for template arrays.
 *
 * This class is the base that is inherited by all containers storing
 * values (matrix, vector, point). Expressions are not arrays. Any derived
 * class can be a lhs in an expression.
 *
 * The common API for these objects is contained in this class.
 *
 *  @tparam Derived is the derived type, e.g., a matrix, vector, point type or
 *  an expression.
 **/
template< class Derived>
class ArrayBase: public ExprBase<Derived>
{
  public:
    typedef ExprBase<Derived> Base;
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef typename hidden::Traits<Derived>::ConstReturnType ConstReturnType;

    enum
    {
      structure_ = hidden::Traits<Derived>::structure_,
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };

  protected:
    /** Default constructor. Default values are cols=(1:0) and rows=(1:0). */
    ArrayBase(): Base() {}
    /** destructor */
    ~ArrayBase() {}

  public:
    // start appliers
    /** Apply the Visitor @c visitor to the whole coefficients of the array.
      * The template parameter @c Visitor is the type of the visitor and provides
      * the following interface:
      * @code
      * struct MyVisitor {
      *   // called for all  coefficients
      *   void operator() (Type& value);
      * };
      * @endcode
      *
      * @note visitors offer automatic unrolling for small fixed size matrix.
      *
      * @sa setValue, setOnes(), setZeros()
      */
    template<typename Visitor>
    void apply(Visitor& visitor);
    /** set random values to this using a uniform law. @sa apply(), randGauss() */
    Derived& randUnif();
    /** set random values to this using a standard gaussian law.
     *  @sa randUnif(), rand(Law::IUnivLaw<Type> const& law), apply()
     **/
    Derived& randGauss();
    /** set random values to this using a distribution law given by the user.
     *  @sa randGauss(), randUnif(), apply() */
    Derived& rand( Law::IUnivLaw<Type> const& law);
    /** set one to this using a Visitor. @sa apply(), setValue(), setZeros() */
    Derived& setOnes();
    /** set zero to this using a Visitor. @sa apply(), setOnes(), setValue()*/
    Derived& setZeros();
    /** set one to this using a Visitor. @sa apply(), setValue(), zeros() */
    Derived& ones();
    /** set zero to this using a Visitor. @sa apply(), ones(), setValue()*/
    Derived& zeros();

    /** set a value to this container. @sa apply(), setOnes(), setZeros()
     *  @param value the value to set
     **/
    Derived& setValue(Type const& value);
    /** @return a copy of @c rhs inside @c this object.
     *  If the ranges of @c this and @c rhs are not exactly the same, the assign
     *  method will call the resize method on this.
     *
     *  @note If @c this is a reference, it cannot be resized and thus an
     *  exception will be thrown if range is not the same.
     **/
    template<class Rhs> Derived& assign(ExprBase<Rhs> const& rhs);

    /** @return the matrix or vector obtained by setting this constant*/
    Derived& operator=( Type const& value);
    /** @return the matrix or vector obtained by evaluating this expression */
    Derived& operator=( Derived const& rhs);
    /** @return the matrix or vector obtained by evaluating this expression */
    template<typename Rhs>
    Derived& operator=( ExprBase<Rhs> const& rhs);

    /** Add Rhs to this. */
    template<typename Rhs>
    Derived& operator+=( ExprBase<Rhs> const& other);
    /** subtract a Rhs to this. */
    template<typename Rhs>
    Derived& operator-=( ExprBase<Rhs> const& other);
    /** divide this by Rhs. */
    template<typename Rhs>
    Derived& operator/=( ExprBase<Rhs> const& other);
    /** Take modulo of this by Rhs. */
    template<typename Rhs>
    Derived& operator%=( ExprBase<Rhs> const& other);
    /** multiply this by Rhs. */
    template<typename Rhs>
    Derived& operator*=( ExprBase<Rhs> const& other);

    /** Adding a constant to this. */
    Derived& operator+=( Type const& other);
    /** Subtract a constant to this. */
    Derived& operator-=( Type const& other);
    /** product of this by a constant. */
    Derived& operator*=( Type const& other);
    /** dividing this by a constant. */
    Derived& operator/=( Type const& other);
    /** take modulo of this by a constant. */
    Derived& operator%=( Type const& other);

    /** overwrite @c this with @c src.
     *  @note this method does not take care of the possibility of overlapping
     *  @param rhs the right hand side to copy
     **/
    template<class Rhs>
    Derived& copy( ExprBase<Rhs> const& rhs);

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

    // slice operators and accessors
    /** @return the j-th column of this. */
    inline ColOperator<Derived> const col(int j) const
    { return ColOperator<Derived> (this->asDerived(), j);}
    /** @return the i-th row of this. */
    inline RowOperator<Derived> const row(int i) const
    { return RowOperator<Derived> (this->asDerived(), i);}
    /** @return the sub-vector(I) of this. */
    template<int Size_>
    inline SubVectorOperator<Derived, Size_> const sub(TRange<Size_> const& I) const
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);
      return SubVectorOperator<Derived, Size_>(this->asDerived(), I);
    }
    /** @return the sub-array(I,J) of this. */
    template<int SizeRows_, int SizeCols_>
    inline SubOperator<Derived, SizeRows_, SizeCols_> const sub(TRange<SizeRows_> const& I,TRange<SizeCols_> const& J) const
    {
      STK_STATIC_ASSERT_TWO_DIMENSIONS_ONLY(Derived);
      return SubOperator<Derived, SizeRows_, SizeCols_>(this->asDerived(), I, J);
    }
    /** @return the j-th column of this. */
    inline ColAccessor<Derived> col(int j)
    { return ColAccessor<Derived> (this->asDerived(), j);}
    /** @return the i-th row of this. */
    inline RowAccessor<Derived> row(int i)
    { return RowAccessor<Derived> (this->asDerived(), i);}
    /** @return the sub-vector(I) of this. */
    template<int Size_>
    inline SubVectorAccessor<Derived, Size_> sub(TRange<Size_> const& I)
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);
      return SubVectorAccessor<Derived, Size_>(this->asDerived(), I);
    }
    /** @return the sub-array(I,J) of this. */
    template<int SizeRows_, int SizeCols_>
    inline SubAccessor<Derived, SizeRows_, SizeCols_> sub(TRange<SizeRows_> const& I,TRange<SizeCols_> const& J)
    {
      STK_STATIC_ASSERT_TWO_DIMENSIONS_ONLY(Derived);
      return SubAccessor<Derived, SizeRows_, SizeCols_>(this->asDerived(), I, J);
    }

    /** Convenient operator to set the coefficients of a matrix.
     *
     * The coefficients must be provided in the row/column order and exactly
     * match the size of the matrix. Otherwise an exception is throwed.
     **/
    ArrayInitializer<Derived> operator<<(Type const& s);

    /** \sa operator<<(Type const&) */
    template<typename Rhs>
    ArrayInitializer<Derived> operator<<(ArrayBase<Rhs> const& other);
};

} // namespace STK

#undef MAKE_RESHAPE_OPERATOR

#endif /* STK_ARRAYBASE_H_ */
