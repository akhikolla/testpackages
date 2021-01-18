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
 * created on: 20 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_DotProduct.h
 *  @brief In this file we implement the DotProduct class.
 **/

#ifndef STK_DOTPRODUCT_H
#define STK_DOTPRODUCT_H

#include "../allocators/STK_CAllocator.h"

namespace STK
{


template<typename Lhs, typename Rhs> class DotProduct;

namespace hidden
{

/** @ingroup hidden
 *  @brief Traits class for the DotProduct class
 */
template< typename Lhs, typename Rhs>
struct Traits< DotProduct < Lhs, Rhs> >
{
  enum
  {
    structure_ = Arrays::number_,
    orient_    = Arrays::by_col_,
    sizeRows_  = 1,
    sizeCols_  = 1,
    storage_   = Arrays::dense_
  };
  typedef typename Promote< typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type const& ConstReturnType;
  typedef CAllocator<Type, sizeRows_, sizeCols_, (bool)orient_> Allocator;
};

} // end namespace hidden


/** @ingroup Arrays
  *
  * @brief Generic expression where a DotProduct is applied to two expressions.
  *
  * @tparam Lhs the type of the left-hand side
  * @tparam Rhs the type of the right-hand side
  *
  * This class represents an expression where a dot product is applied to
  * two expressions. The left hand side being a point_ (a row-oriented vector)
  * and the right hand side a vector_ (a column-oriented vector).
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name DotProduct types explicitly.
  **/
template<typename Lhs, typename Rhs>
class DotProduct: public ExprBase< DotProduct<Lhs, Rhs> >
                 , public TRef<1>
{
  public:
    typedef typename hidden::Traits<DotProduct>::Type Type;
    typedef typename hidden::Traits<DotProduct>::ConstReturnType ConstReturnType;
    typedef typename hidden::Traits<DotProduct>::Allocator Allocator;

    enum
    {
      structure_ = hidden::Traits<DotProduct>::structure_,
      orient_    = hidden::Traits<DotProduct>::orient_,
      sizeRows_  = hidden::Traits<DotProduct>::sizeRows_,
      sizeCols_  = hidden::Traits<DotProduct>::sizeCols_,
      storage_   = hidden::Traits<DotProduct>::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    inline DotProduct( const Lhs& lhs, const Rhs& rhs)
                     : lhs_(lhs), rhs_(rhs)
                     , result_()
    {
      STK_STATIC_ASSERT_POINT_ONLY(Lhs);
      STK_STATIC_ASSERT_VECTOR_ONLY(Rhs);
#ifdef STK_BOUNDS_CHECK
      if (lhs.range() != rhs.range())
      { STKRUNTIME_ERROR_NO_ARG(DotProduct:DotProduct, sizes mismatch between Lhs and Rhs);}
#endif
      result_.shift(lhs.beginRows(), rhs.beginCols());
      result_.elt() = lhs.dot(rhs);
    }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return result_.rows();}
    /** @return the columns range */
    inline ColRange const& colsImpl() const { return result_.cols();}

    /** access to the element */
    inline ConstReturnType elt2Impl(int i, int j) const
    {
#ifdef STK_BOUNDS_CHECK
      if ( (result_.beginRows() != i) || (result_.beginCols() !=j) )
      { STKRUNTIME_ERROR_2ARG(DotProduct:elt2Impl,i,j,wrong indexes);}
#endif
      return result_.elt(i, j);
    }
    /** access to the element */
    inline ConstReturnType elt0Impl() const { return result_.elt();}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the right hand side nested expression */
    inline Rhs const& rhs() const { return rhs_; }
    /** @return the result */
    inline Allocator const& result() const { return result_; }

  protected:
    Lhs const& lhs_;
    Rhs const& rhs_;

  private:
    Allocator result_;
};

}  // namespace STK

#endif /* STK_DOTPRODUCT_H */
