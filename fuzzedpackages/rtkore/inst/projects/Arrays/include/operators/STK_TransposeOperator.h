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
 * created on: 17 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_TransposeOperator.h
 *  @brief In this file we implement the TransposeOperator class.
 **/

#ifndef STK_TRANSPOSEOPERATOR_H
#define STK_TRANSPOSEOPERATOR_H

#include "STK_SlicingOperators.h"

namespace STK
{

namespace hidden
{
/** @ingroup hidden
  * Helper Traits class for transposed operator. Default is structure_ unmodified
  * by transpose operator.
  **/
template<int Structure_> struct TransposeTraits
{ enum { structure_ = Structure_}; };
/** specialization for lower_triangular_ */
template<> struct TransposeTraits<Arrays::lower_triangular_>
{ enum { structure_ = Arrays::upper_triangular_}; };
/** specialization for upper_triangular_ */
template<> struct TransposeTraits<Arrays::upper_triangular_>
{ enum { structure_ = Arrays::lower_triangular_}; };
/** specialization for lower_triangular_ */
template<> struct TransposeTraits<Arrays::lower_symmetric_>
{ enum { structure_ = Arrays::upper_symmetric_}; };
/** specialization for upper_triangular_ */
template<> struct TransposeTraits<Arrays::upper_symmetric_>
{ enum { structure_ = Arrays::lower_symmetric_}; };
/** specialization for vector_ */
template<> struct TransposeTraits<Arrays::vector_>
{ enum { structure_ = Arrays::point_}; };
/** specialization for point_ */
template<> struct TransposeTraits<Arrays::point_>
{ enum { structure_ = Arrays::vector_}; };

} // namespace hidden

// forward declaration
template< typename Array> class TransposeOperator;
template< typename Array> class TransposeAccessor;

namespace hidden
{

/** @ingroup hidden
 *  @brief Traits class for the transposed operator
 */
template<typename Lhs>
struct Traits< TransposeOperator <Lhs> >
{
  enum
  {
    structure_ = TransposeTraits<Lhs::structure_>::structure_,
    orient_    = !Lhs::orient_,
    sizeRows_  = Lhs::sizeCols_,
    sizeCols_  = Lhs::sizeRows_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator<TransposeOperator < Lhs> > Row;
  typedef ColOperator<TransposeOperator < Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ConstReturnType ConstReturnType;
};

/** @ingroup hidden
 *  @brief Traits class for the transposed operator
 */
template<typename Lhs>
struct Traits< TransposeAccessor<Lhs> >
{
  enum
  {
    structure_ = TransposeTraits<Lhs::structure_>::structure_,
    orient_    = !Lhs::orient_,
    sizeRows_  = Lhs::sizeCols_,
    sizeCols_  = Lhs::sizeRows_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator<TransposeAccessor < Lhs> > Row;
  typedef ColOperator<TransposeAccessor < Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ConstReturnType ConstReturnType;
};

} // end namespace hidden


/** @ingroup Arrays
 *  @class TransposeOperator
  *
  * \brief Generic expression when an expression is transposed
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * transpose operator.
  *
  * This class represents an expression where a transpose operator is applied to
  * an expression. It is the return type of the transpose operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name TransposeOperator type explicitly.
  */
template< typename Lhs>
class TransposeOperator: public ExprBase< TransposeOperator<Lhs> >, public TRef<1>
{
  public:
    typedef ExprBase< TransposeOperator<Lhs> > Base;
    typedef typename hidden::Traits< TransposeOperator<Lhs> >::Type Type;
    typedef typename hidden::Traits< TransposeOperator<Lhs> >::ConstReturnType ConstReturnType;

    enum
    {
        structure_ = hidden::Traits< TransposeOperator<Lhs> >::structure_,
        orient_    = hidden::Traits< TransposeOperator<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< TransposeOperator<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< TransposeOperator<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< TransposeOperator<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;
    /** Constructor */
    inline TransposeOperator( Lhs const& lhs): Base(), lhs_(lhs) {}

    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs_.cols();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs_.rows();}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the element (i,j) of the transposed expression.
     *  @param i, j index of the row and of the column
     **/
    inline ConstReturnType elt2Impl(int i, int j) const { return lhs_.elt(j, i);}
    /** @return the element ith element of the transposed expression
     *  @param i index of the element to get
     **/
    inline ConstReturnType elt1Impl(int i) const { return lhs_.elt(i);}
    /** access to the element of the transposed expression */
    inline ConstReturnType elt0Impl() const { return lhs_.elt();}

  protected:
    Lhs const& lhs_;
};

/** @ingroup Arrays
 *  @class TransposeAccessor
  *
  * \brief Generic expression when an expression is transposed
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * transpose operator.
  *
  * This class represents an expression where a transpose operator is applied to
  * an expression. It is the return type of the transpose operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name TransposeAccessor type explicitly.
  */
template< typename Lhs>
class TransposeAccessor: public ArrayBase< TransposeAccessor<Lhs> >, public TRef<1>
{
  public:
    typedef ArrayBase< TransposeAccessor<Lhs> > Base;
    typedef typename hidden::Traits< TransposeOperator<Lhs> >::Type Type;
    typedef typename hidden::Traits< TransposeOperator<Lhs> >::ConstReturnType ConstReturnType;

    enum
    {
        structure_ = hidden::Traits< TransposeOperator<Lhs> >::structure_,
        orient_    = hidden::Traits< TransposeOperator<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< TransposeOperator<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< TransposeOperator<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< TransposeOperator<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;
    /** Constructor */
    inline TransposeAccessor( Lhs& lhs): Base(), lhs_(lhs) {}

    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs_.cols();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs_.rows();}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }

    /** @return the element (i,j) of the transposed expression.
     *  @param i, j index of the row and of the column
     **/
    inline ConstReturnType elt2Impl(int i, int j) const { return lhs_.elt(j, i);}
    /** @return the element ith element of the transposed expression
     *  @param i index of the element to get
     **/
    inline ConstReturnType elt1Impl(int i) const { return lhs_.elt(i);}
    /** access to the element of the transposed expression */
    inline ConstReturnType elt0Impl() const { return lhs_.elt();}
    /** @return the element (i,j) of the transposed expression.
     *  @param i, j index of the row and of the column
     **/
    inline Type& elt2Impl(int i, int j) { return lhs_.elt(j, i);}
    /** @return the element ith element of the transposed expression
     *  @param i index of the element to get
     **/
    inline Type& elt1Impl(int i) { return lhs_.elt(i);}
    /** access to the element of the transposed expression */
    inline Type& elt0Impl() { return lhs_.elt();}

  protected:
    Lhs& lhs_;
};

} // namespace STK

#endif /* STK_TRANSPOSEOPERATOR_H */
