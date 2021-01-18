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

/** @file STK_UnaryOperators.h
 *  @brief In this file we implement the UnaryOperator class.
 **/

#ifndef STK_UNARYOPERATORS_H
#define STK_UNARYOPERATORS_H

#include "STK_SlicingOperators.h"

namespace STK
{

// forward declaration
template<typename UnaryOp, typename Lhs> class UnaryOperator;

namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for the unary operators
 */
template<typename UnaryOp, typename Lhs>
struct Traits< UnaryOperator <UnaryOp, Lhs> >
{
  enum
  {
      structure_ = Traits<Lhs>::structure_,
      orient_    = Traits<Lhs>::orient_,
      sizeRows_  = Traits<Lhs>::sizeRows_,
      sizeCols_  = Traits<Lhs>::sizeCols_,
      storage_   = Traits<Lhs>::storage_
  };
  typedef RowOperator<UnaryOperator <UnaryOp, Lhs> > Row;
  typedef ColOperator<UnaryOperator <UnaryOp, Lhs> > Col;

  typedef typename UnaryOp::result_type Type;
  typedef typename UnaryOp::result_type ConstReturnType;
};

} // end namespace hidden

// forward declaration
template<typename UnaryOp, typename Lhs>
class UnaryOperatorBase;


/** @ingroup Arrays
  *
  * @brief Generic expression when unary operator is applied to an expression
  *
  * @tparam UnaryOp template functor implementing the operator
  * @tparam Lhs the type of the expression to which we are applying the unary operator
  *
  * This class represents an expression where a unary operator is applied to
  * an expression. It is the return type of all operations taking exactly 1
  * input expression, regardless of the presence of other inputs such as
  * numbers. For example, the operator* in the expression 3*matrix is
  * considered unary, because only the right-hand side is an expression, and its
  * return type is a specialization of UnaryOperator.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name UnaryOperator types explicitly.
  */
template<typename UnaryOp,  typename Lhs>
class UnaryOperator: public ExprBase< UnaryOperator<UnaryOp, Lhs> >, public TRef<1>
{
  public:
    typedef ExprBase<  UnaryOperator<UnaryOp, Lhs> > Base;
    typedef typename hidden::Traits< UnaryOperator<UnaryOp, Lhs> >::Type Type;
    typedef typename hidden::Traits< UnaryOperator<UnaryOp, Lhs> >::ConstReturnType ConstReturnType;

    enum
    {
        structure_ = hidden::Traits< UnaryOperator >::structure_,
        orient_    = hidden::Traits< UnaryOperator >::orient_,
        sizeRows_  = hidden::Traits< UnaryOperator >::sizeRows_,
        sizeCols_  = hidden::Traits< UnaryOperator >::sizeCols_,
        storage_   = hidden::Traits< UnaryOperator >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    inline UnaryOperator( Lhs const& lhs, UnaryOp const& functor = UnaryOp())
                       : Base(), lhs_(lhs), functor_(functor)
    {}
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs_.rows();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs_.cols();}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the functor representing the unary operation */
    inline UnaryOp const& functor() const { return functor_; }

    /** @return the element (i,j) of the operator.
     *  @param i index of the row
     *  @param j index of the column
     **/
    inline ConstReturnType elt2Impl(int i, int j) const { return functor_(lhs_.elt(i, j));}
    /** @return the ith element of the operator
     *  @param i index of the element to get
     **/
    inline ConstReturnType elt1Impl(int i) const { return functor_(lhs_.elt(i));}
    /** @return the element of the operator */
    inline ConstReturnType elt0Impl() const { return functor_(lhs_.elt());}

  protected:
    Lhs const& lhs_;
    UnaryOp const functor_;
};

} // namespace STK

#endif /* STK_UNARYOPERATORS_H */
