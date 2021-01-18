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

/** @file STK_BinaryOperators.h
 *  @brief In this file we implement the BinaryOperator class.
 **/


#ifndef STK_BINARYOPERATORS_H
#define STK_BINARYOPERATORS_H

#include "STK_SlicingOperators.h"
#include "STK_BinaryImpl.h"

#define EGAL(arg1, arg2) ((arg1::structure_ == int(Arrays::arg2)))

namespace STK
{

// forward declaration
template<typename FunctorOp, typename Lhs, typename Rhs>
class BinaryOperator;

namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for the BinaryOperator
 */
template<typename FunctorOp, typename Lhs, typename Rhs>
struct Traits< BinaryOperator<FunctorOp, Lhs, Rhs> >
{
  enum
  {
    // find the kind of binary operator and the Structure using helper class BinaryEltImpl
    binary_op_Kind_ = BinaryEltImpl<FunctorOp, Lhs, Rhs, Lhs::structure_, Rhs::structure_>::binary_op_Kind_,

    isLhs1D_ = EGAL(Lhs,vector_)||EGAL(Lhs,point_)||EGAL(Lhs,diagonal_),
    isRhs1D_ = EGAL(Rhs,vector_)||EGAL(Rhs,point_)||EGAL(Rhs,diagonal_),

    isRhs2D_ = EGAL(Rhs,array2D_)||EGAL(Rhs,square_)||EGAL(Rhs,diagonal_)
             ||EGAL(Rhs,lower_triangular_)||EGAL(Rhs,upper_triangular_)
             ||EGAL(Rhs,symmetric_)||EGAL(Rhs,lower_symmetric_)||EGAL(Rhs,upper_symmetric_),
    isLhs2D_ = EGAL(Lhs,array2D_)||EGAL(Lhs,square_)||EGAL(Lhs,diagonal_)
             ||EGAL(Lhs,lower_triangular_)||EGAL(Lhs,upper_triangular_)
             ||EGAL(Lhs,symmetric_)||EGAL(Lhs,lower_symmetric_)||EGAL(Lhs,upper_symmetric_),


    isRes0D_ = EGAL(Lhs,number_) && EGAL(Rhs,number_),
    isRes1D_ = (EGAL(Lhs,vector_)||EGAL(Lhs,point_)) && (EGAL(Rhs,vector_)||EGAL(Rhs,point_)),
    isRes2D_ = isLhs2D_ && isRhs2D_,

    is1D1D_  = isLhs1D_ && isRhs1D_,

    // get the structure from the helper class BinaryEltImpl
    structure_ = hidden::BinaryEltImpl<FunctorOp, Lhs, Rhs, Lhs::structure_, Rhs::structure_>::structure_,
    orient_    = Lhs::orient_,    // preserve the Lhs storage orientation. Could be optimized ?
    sizeRows_  = BinaryEltImpl<FunctorOp, Lhs, Rhs, Lhs::structure_, Rhs::structure_>::sizeRows_,
    sizeCols_  = BinaryEltImpl<FunctorOp, Lhs, Rhs, Lhs::structure_, Rhs::structure_>::sizeCols_,
    storage_   = (Lhs::storage_ == int(Arrays::dense_)) || (Rhs::storage_ == int(Arrays::dense_))
               ?  int(Arrays::dense_) : int(Arrays::sparse_),

    useForRows_    = BinaryEltImpl<FunctorOp, Lhs, Rhs, Lhs::structure_, Rhs::structure_>::useForRows_,
    useForCols_    = BinaryEltImpl<FunctorOp, Lhs, Rhs, Lhs::structure_, Rhs::structure_>::useForCols_
  };
  typedef RowOperator<BinaryOperator<FunctorOp, Lhs, Rhs> > Row;
  typedef ColOperator<BinaryOperator<FunctorOp, Lhs, Rhs> > Col;
  typedef typename FunctorOp::result_type Type;
  typedef typename FunctorOp::result_type ConstReturnType;
};

} // end namespace hidden


/** @ingroup Arrays
  * @brief Generic expression where a binary operator is applied to two expressions
  *
  * @tparam FunctorOp template functor implementing the binary operator
  * @tparam Lhs left-hand side type
  * @tparam Rhs right-hand side type
  *
  * This class represents an expression  where a binary operator is applied to
  * two expressions.
  * It is the return type of binary operators, by which we mean only those
  * binary operators where both the left-hand side and the right-hand side
  * are expressions. For example, the return type of matrix1+matrix2 is a
  * BinaryOperator. The return type of number+matrix is a unary operator.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name BinaryOperator types explicitly.
  **/

template<typename FunctorOp, typename Lhs, typename Rhs>
class BinaryOperator: public ExprBase< BinaryOperator<FunctorOp, Lhs, Rhs> >
                    , public TRef<1>
{
  public:
    typedef hidden::BinaryEltImpl<FunctorOp, Lhs, Rhs, Lhs::structure_, Rhs::structure_> EltImpl;
    typedef hidden::BinaryRowsImpl< Lhs, Rhs, EltImpl::sizeRows_, EltImpl::useForRows_ > RowsImpl;
    typedef hidden::BinaryColsImpl< Lhs, Rhs, EltImpl::sizeCols_, EltImpl::useForCols_ > ColsImpl;

    typedef typename hidden::Traits< BinaryOperator >::ConstReturnType ConstReturnType;
    typedef typename hidden::Traits<BinaryOperator >::Type Type;

    enum
    {
      isRes0D_   = hidden::Traits<BinaryOperator>::isRes0D_,
      isRes1D_   = hidden::Traits<BinaryOperator>::isRes1D_,
      isRes2D_   = hidden::Traits<BinaryOperator>::isRes2D_,
      is1D1D_    = hidden::Traits<BinaryOperator>::is1D1D_,

      structure_ = hidden::Traits<BinaryOperator>::structure_,
      orient_    = hidden::Traits<BinaryOperator>::orient_,
      sizeRows_  = hidden::Traits<BinaryOperator>::sizeRows_,
      sizeCols_  = hidden::Traits<BinaryOperator>::sizeCols_,
      storage_   = hidden::Traits<BinaryOperator>::storage_,

      // All the valid cases for binary operators
      isValid_ =( isRes0D_ || isRes1D_  || isRes2D_ || is1D1D_)
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** default constructor */

    inline BinaryOperator( Lhs const& lhs, Rhs const& rhs, FunctorOp const& func = FunctorOp())
                         : lhs_(lhs), rhs_(rhs), functor_(func)
    {
      // FIXME : not safe. Add more test in the 1D case at compile time (and runtime ?)
      STK_STATIC_ASSERT_BINARY_OPERATOR_MISMATCH( isValid_ );
      STK_STATIC_ASSERT_COLS_DIMENSIONS_MISMATCH(!( (int(Lhs::sizeCols_) != UnknownSize)
                                                &&  (int(Rhs::sizeCols_) != UnknownSize)
                                                &&  (int(Lhs::sizeCols_) != int(Rhs::sizeCols_))
                                                &&  (isRes2D_)
                                                 ));
      STK_STATIC_ASSERT_ROWS_DIMENSIONS_MISMATCH(!( (int(Lhs::sizeRows_) != UnknownSize)
                                                &&  (int(Rhs::sizeRows_) != UnknownSize)
                                                &&  (int(Lhs::sizeRows_) != int(Rhs::sizeRows_))
                                                &&  (isRes2D_)
                                                 ));
#ifdef STK_BOUNDS_CHECK
      if ((lhs.rows() != rhs.rows()) && (isRes2D_))
      { STKRUNTIME_ERROR_2ARG(BinaryOperator, lhs.rows(), rhs.rows(), Rows range mismatch in BinaryOperator);}
      if (( lhs.cols() != rhs.cols()) && (isRes2D_))
      { STKRUNTIME_ERROR_2ARG(BinaryOperator, lhs.cols(), rhs.cols(), Columns range mismatch in BinaryOperator);}
#endif
    }

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the right hand side nested expression */
    inline Rhs const& rhs() const { return rhs_; }
    /** @return the functor representing the binary operation */
    inline FunctorOp const& functor() const { return functor_; }

    /** @return element (i,j) */
    inline ConstReturnType elt2Impl(int i, int j) const { return EltImpl::elt2Impl(functor_, lhs_, rhs_, i, j);}
    /** @return element i */
    inline ConstReturnType elt1Impl(int i) const { return EltImpl::elt1Impl(functor_, lhs_, rhs_, i);}
    /** @return element */
    inline ConstReturnType elt0Impl() const { return EltImpl::elt0Impl(functor_, lhs_, rhs_);}
    /** @return range of the rows */
    inline RowRange const& rowsImpl() const { return RowsImpl::rowsImpl(lhs_, rhs_);}
    /** @return range of the columns */
    inline ColRange const& colsImpl() const { return ColsImpl::colsImpl(lhs_, rhs_);}

  protected:
    Lhs const& lhs_;
    Rhs const& rhs_;
    FunctorOp const functor_;
};


} // namespace STK

#undef EGAL

#endif /* STK_BINARYOPERATORS_H */
