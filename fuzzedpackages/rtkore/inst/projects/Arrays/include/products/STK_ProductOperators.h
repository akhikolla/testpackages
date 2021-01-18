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

/** @file STK_ProductOperators.h
 *  @brief In this file we implement the ProductOperator class.
 **/

#ifndef STK_PRODUCTOPERATORS_H
#define STK_PRODUCTOPERATORS_H


#include <Sdk/include/STK_StaticAssert.h>
#include "../allocators/STK_CAllocator.h"
#include "../operators/STK_BinaryImpl.h"
#include "STK_ProductDispatcher.h"

#define EGAL(arg1, arg2) ((arg1::structure_ == int(Arrays::arg2)))

namespace STK
{

/** @ingroup Arrays
  * @class ArrayByDiagonalProduct
  * @brief Generic expression where a product operator is
  * applied to two expressions.
  *
  * @tparam Lhs the type of the left-hand side array expression
  * @tparam Rhs the type of the right-hand side diagonal expression
  *
  * This class represents an expression  where a product operator is applied to
  * two expressions. The left hand Side can be of any kind, the Left Hand Side
  * is a diagonal matrix or expression.
  **/
template<typename Lhs, typename Rhs> class ArrayByDiagonalProduct;

/** @ingroup Arrays
  * @class DiagonalByArrayProduct
  * @brief Generic expression where a product operator is
  * applied to two expressions.
  *
  * @tparam Lhs the type of the left-hand side, a diagonal expression
  * @tparam Rhs the type of the right-hand side, a vector expression
  *
  * This class represents an expression  where a product operator is applied to
  * two expressions. The left hand Side is a diagonal expression,
  * the Right Hand Side is a vector expression.
  **/
template<typename Lhs, typename Rhs> class DiagonalByArrayProduct;

/** @ingroup Arrays
  * @class VectorByPointProduct
  * @brief Generic expression where a product operator is
  * applied to two expressions.
  *
  * @tparam Lhs the type of the left-hand side, a vector expression
  * @tparam Rhs the type of the right-hand side, a row-vector expression
  *
  * This class represents an expression  where a product operator is applied to
  * two expressions. The left hand Side is a vector expression,
  * the Right Hand Side is a row-vector expression.
  **/
template<typename Lhs, typename Rhs> class VectorByPointProduct;

/** @ingroup Arrays
  * @class PointByArrayProduct
  * @brief Generic expression where a product operator is
  * applied to two expressions.
  *
  * @tparam Lhs the type of the left-hand side
  * @tparam Rhs the type of the right-hand side
  *
  * This class represents an expression  where a product operator is applied to
  * two expressions. The left hand Side is a point, the Right Hand Side
  * is any compatible expression or matrix.
  **/
template<typename Lhs, typename Rhs> class PointByArrayProduct;

/** @ingroup Arrays
  * @class ArrayByVectorProduct
  * @brief Generic expression where a product operator is
  * applied to two expressions.
  *
  * @tparam Lhs the type of the left-hand side
  * @tparam Rhs the type of the right-hand side
  *
  * This class represents an expression  where a product operator is applied to
  * two expressions. The left hand Side can be of any kind, the Right Hand Side
  * is a vector or vectorial expression.
  **/
template<typename Lhs, typename Rhs> class ArrayByVectorProduct;

/** @ingroup Arrays
  * @class ArrayByArrayProduct
  * @brief Generic expression where a product operator is
  * applied to two expressions.
  *
  * @tparam Lhs the type of the left-hand side
  * @tparam Rhs the type of the right-hand side
  *
  * This class represents an expression  where a product operator is applied to
  * two expressions.
  * It is the return type of product operators, by which we mean only those
  * product operators where both the left-hand side and the right-hand side
  * are expressions. For example, the return type of matrix1*matrix2 is a
  * ArrayByArrayProduct.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name ArrayByArrayProduct types explicitly.
  **/
template<typename Lhs, typename Rhs> class ArrayByArrayProduct;

namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for array2d by diagonal product
 */
template< typename Lhs, typename Rhs>
struct Traits< ArrayByDiagonalProduct < Lhs, Rhs> >
{
  enum
  {
    structure_ = Lhs::structure_,
    sizeRows_  = Lhs::sizeRows_ ,
    sizeCols_  = Rhs::sizeCols_,
    orient_    = Lhs::orient_,
    storage_   = Lhs::storage_
  };
  typedef typename Promote<typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type ConstReturnType;
};

} // namespace hidden

/** @ingroup Arrays
  * @brief Generic expression where a product operator is applied to two expressions
  *
  * @tparam Lhs left-hand side type
  * @tparam Rhs right-hand side type
  *
  * This class represents an expression where a product operator is applied to
  * two expressions.
  *
  * It is the return type of product operator when the left hand side expression
  * is an array and the right and side expression is a diagonal expression.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name ArrayByDiagonalProduct types explicitly.
  **/
template<typename Lhs, typename Rhs>
class ArrayByDiagonalProduct: public ExprBase< ArrayByDiagonalProduct<Lhs, Rhs> >
                            , public TRef<1>
{
  public:
    typedef ExprBase< ArrayByDiagonalProduct<Lhs, Rhs> > Base;
    typedef typename hidden::Traits<ArrayByDiagonalProduct>::Type Type;
    typedef typename hidden::Traits<ArrayByDiagonalProduct>::ConstReturnType ConstReturnType;

    enum
    {
      structure_ = hidden::Traits<ArrayByDiagonalProduct>::structure_,
      orient_    = hidden::Traits<ArrayByDiagonalProduct>::orient_,
      sizeRows_  = hidden::Traits<ArrayByDiagonalProduct>::sizeRows_,
      sizeCols_  = hidden::Traits<ArrayByDiagonalProduct>::sizeCols_,
      storage_   = hidden::Traits<ArrayByDiagonalProduct>::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    inline ArrayByDiagonalProduct( const Lhs& lhs, const Rhs& rhs)
                                 : Base(), lhs_(lhs), rhs_(rhs)
    {
      if (lhs.cols() != rhs.rows())
      { STKRUNTIME_ERROR_NO_ARG(ArrayByDiagonalProduct, sizes mismatch);}
    }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs_.rows();}
    /** @return the columns range */
    inline ColRange const& colsImpl() const { return rhs_.cols();}

    /** @return the element (i,j) */
    inline ConstReturnType elt2Impl(int i, int j) const { return lhs_.elt(i,j)*rhs_.elt(j);}
    /** @return ith element
     *  @note make sense if lhs_ is a diagonal matrix, in this case this operator is also a
     *  diagonal matrix.
     **/
    inline ConstReturnType elt1Impl(int i) const { return lhs_.elt(i)*rhs_.elt(i);}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the right hand side nested expression */
    inline Rhs const& rhs() const { return rhs_; }

  protected:
    Lhs const& lhs_;
    Rhs const& rhs_;
};

namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for the DiagonalByArrayProduct class
 */
template< typename Lhs, typename Rhs>
struct Traits< DiagonalByArrayProduct < Lhs, Rhs> >
{
  enum
  {
    structure_ = Rhs::structure_,
    sizeRows_  = Lhs::sizeRows_,
    sizeCols_  = Rhs::sizeCols_,
    orient_    = Rhs::orient_,
    storage_   = Rhs::storage_
  };
  typedef typename Promote<typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type ConstReturnType;  // not a reference as it is a temporary
};

} // namespace hidden

/** @ingroup Arrays
  * @brief Generic expression where a product operator is applied to two expressions
  *
  * @tparam Lhs left-hand side diagonal expression
  * @tparam Rhs right-hand side array expression
  *
  * This class represents an expression  where a product operator is applied to
  * two expressions.
  *
  * It is the return type of product operator when the left hand side expression
  * is a diagonal expression and the right and side expression is an array.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name ArrayByDiagonalProduct types explicitly.
  *
  * @note the product is evaluated when element is needed
  *
  **/
template<typename Lhs, typename Rhs>
class DiagonalByArrayProduct: public ExprBase< DiagonalByArrayProduct<Lhs, Rhs> >
                            , public TRef<1>
{
  public:
    typedef ExprBase< DiagonalByArrayProduct<Lhs, Rhs> > Base;
    typedef typename hidden::Traits<DiagonalByArrayProduct>::Type Type;
    typedef typename hidden::Traits<DiagonalByArrayProduct>::ConstReturnType ConstReturnType;

    enum
    {
      structure_ = hidden::Traits<DiagonalByArrayProduct>::structure_,
      orient_    = hidden::Traits<DiagonalByArrayProduct>::orient_,
      sizeRows_  = hidden::Traits<DiagonalByArrayProduct>::sizeRows_,
      sizeCols_  = hidden::Traits<DiagonalByArrayProduct>::sizeCols_,
      storage_   = hidden::Traits<DiagonalByArrayProduct>::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    inline DiagonalByArrayProduct( const Lhs& lhs, const Rhs& rhs)
                                : Base(), lhs_(lhs), rhs_(rhs)
    {
      if (lhs.cols() != rhs.rows())
      { STKRUNTIME_ERROR_NO_ARG(DiagonalByArrayProduct, sizes mismatch);}
    }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs_.rows();}
    /** @return the columns range */
    inline ColRange const& colsImpl() const { return rhs_.cols();}

    /** @return element (i,j) */
    inline ConstReturnType elt2Impl(int i, int j) const { return lhs_.elt(i)*rhs_.elt(i,j);}
    /** @return ith element
     *  @note make sense if rhs_ is a diagonal matrix, in this case this operator is also a
     *  diagonal matrix.
     **/
    inline ConstReturnType elt1Impl(int i) const { return lhs_.elt(i)*rhs_.elt(i);}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the right hand side nested expression */
    inline Rhs const& rhs() const { return rhs_; }

  protected:
    Lhs const& lhs_;
    Rhs const& rhs_;
};


namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for the VectorByPointProduct class
 *  @note the product is evaluated when element is needed
 */
template< typename Lhs, typename Rhs>
struct Traits< VectorByPointProduct < Lhs, Rhs> >
{
  enum
  {
    structure_ = Arrays::array2D_,
    sizeRows_  = Traits<Lhs>::sizeRows_,
    sizeCols_  = Traits<Rhs>::sizeCols_,
    orient_    = Arrays::by_col_,
    storage_   = ( Traits<Lhs>::storage_ == int(Arrays::dense_)) || (Traits<Rhs>::storage_ == int(Arrays::dense_))
                 ? int(Arrays::dense_) : int(Arrays::sparse_)
  };
  typedef typename Promote<typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type ConstReturnType; // not a reference
};

} // namespace hidden

/** @ingroup Arrays
  * @brief Generic expression where a product operator is applied to two expressions
  *
  * @tparam Lhs left-hand side type
  * @tparam Rhs right-hand side type
  *
  * This class represents an expression  where a product operator is applied to
  * two expressions.
  *
  * It is the return type of product operator when the left hand side expression
  * is a (column) vector and the right and side expression is a point (a row vector).
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name VectorByPointProduct types explicitly.
  **/
template<typename Lhs, typename Rhs>

class VectorByPointProduct: public ExprBase< VectorByPointProduct<Lhs, Rhs> >
                          , public TRef<1>
{
  public:
    typedef ExprBase< VectorByPointProduct<Lhs, Rhs> > Base;
    typedef typename hidden::Traits<VectorByPointProduct>::Type Type;
    typedef typename hidden::Traits<VectorByPointProduct>::ConstReturnType ConstReturnType;

    enum
    {
      structure_ = hidden::Traits<VectorByPointProduct>::structure_,
      orient_    = hidden::Traits<VectorByPointProduct>::orient_,
      sizeRows_  = hidden::Traits<VectorByPointProduct>::sizeRows_,
      sizeCols_  = hidden::Traits<VectorByPointProduct>::sizeCols_,
      storage_   = hidden::Traits<VectorByPointProduct>::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** constructor */
    VectorByPointProduct( const Lhs& lhs, const Rhs& rhs): Base(), lhs_(lhs), rhs_(rhs){}
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs_.rows();}
    /** @return the columns range */
    inline ColRange const& colsImpl() const { return rhs_.cols();}

    /** @return the element (i,j) */
    inline ConstReturnType elt2Impl(int i, int j) const { return lhs_.elt(i)*rhs_.elt(j);}
    /** @return the element */
    inline ConstReturnType elt0Impl() const { return lhs_.elt()*rhs_.elt();}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the right hand side expression */
    inline Rhs const& rhs() const { return rhs_; }

  protected:
    Lhs const& lhs_;
    Rhs const& rhs_;
};



namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for the PointByArrayProduct
 */
template< typename Lhs, typename Rhs>
struct Traits< PointByArrayProduct < Lhs, Rhs> >
{
  enum
   {
     structure_ = Arrays::point_,
     sizeRows_  = 1,
     sizeCols_  = Traits<Rhs>::sizeCols_,
     orient_    = Arrays::by_row_,
     storage_   = ( Traits<Lhs>::storage_ == int(Arrays::dense_)) || (Traits<Rhs>::storage_ == int(Arrays::dense_))
                  ? int(Arrays::dense_) : int(Arrays::sparse_)
   };

  typedef typename hidden::Promote< typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type const& ConstReturnType;

  typedef CAllocator<Type, sizeRows_, sizeCols_, orient_> Allocator;
};

} // namespace hidden

/** @ingroup Arrays
  * @brief Generic expression where a product operator is applied to two expressions
  *
  * @tparam Lhs left-hand side type
  * @tparam Rhs right-hand side type
  *
  * This class represents an expression  where a product operator is applied to
  * two expressions.
  *
  * It is the return type of product operator when the left hand side expression
  * is a point (a  row vector) and the right and side expression is an array.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name PointByArrayProduct types explicitly.
  **/
template<typename Lhs, typename Rhs>
class PointByArrayProduct: public ExprBase< PointByArrayProduct<Lhs, Rhs> >
                         , public TRef<1>
{
  public:
    typedef ExprBase< PointByArrayProduct<Lhs, Rhs> > Base;
    typedef typename hidden::Traits<PointByArrayProduct>::Type Type;
    typedef typename hidden::Traits<PointByArrayProduct>::ConstReturnType ConstReturnType;
    typedef typename hidden::Traits<PointByArrayProduct>::Allocator Allocator;

    enum
    {
      structure_ = hidden::Traits<PointByArrayProduct>::structure_,
      orient_    = hidden::Traits<PointByArrayProduct>::orient_,
      sizeRows_  = hidden::Traits<PointByArrayProduct>::sizeRows_,
      sizeCols_  = hidden::Traits<PointByArrayProduct>::sizeCols_,
      storage_   = hidden::Traits<PointByArrayProduct>::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    inline PointByArrayProduct( const Lhs& lhs, const Rhs& rhs)
                             : Base(), lhs_(lhs), rhs_(rhs)
                             , result_(1, rhs.sizeCols(), Type(0))
    {
      if (lhs.range() != rhs.rows())
      { STKRUNTIME_ERROR_2ARG(PointByArrayProduct, lhs.range(), rhs.rows(), sizes mismatch);}
      result_.shift(lhs_.beginRows(), rhs_.beginCols());
      hidden::ProductDispatcher<Lhs, Rhs, Allocator>::run(lhs, rhs, result_);
    }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return result_.rows();}
    /** @return the columns range */
    inline ColRange const& colsImpl() const { return result_.cols();}

    /** @return the element (i,j) */
    inline ConstReturnType elt2Impl(int i, int j) const { return result_.elt(i, j);}
    /** @return the ith element */
    inline ConstReturnType elt1Impl(int i) const { return result_.elt(i);}
    /** @return the element */
    inline ConstReturnType elt0Impl() const { return result_.elt();}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the right hand side nested expression */
    inline Rhs const& rhs() const { return rhs_; }
    /** @return the right hand side nested expression */
    Allocator const& result() const { return result_; }

  protected:
    Lhs const& lhs_;
    Rhs const& rhs_;

  private:
    Allocator result_;
};

namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for the ArrayByVectorProduct class
 */
template< typename Lhs, typename Rhs>
struct Traits< ArrayByVectorProduct < Lhs, Rhs> >
{
  enum
   {
     structure_ = Arrays::vector_,
     orient_    = Arrays::by_col_,
     sizeRows_  = Traits<Lhs>::sizeRows_,
     sizeCols_  = 1,
     storage_   = ( Traits<Lhs>::storage_ == int(Arrays::dense_)) || (Traits<Rhs>::storage_ == int(Arrays::dense_))
                  ? int(Arrays::dense_) : int(Arrays::sparse_)
   };

  typedef typename hidden::Promote< typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type const& ConstReturnType;

  typedef CAllocator<Type, sizeRows_, sizeCols_, orient_> Allocator;
};

} // nmaespace hidden

/** @ingroup Arrays
  * @brief Generic expression where a product operator is applied to two expressions
  *
  * @tparam Lhs left-hand side type
  * @tparam Rhs right-hand side type
  *
  * This class represents an expression  where a product operator is applied to
  * two expressions.
  *
  * It is the return type of product operator when the left hand side expression
  * is an array and the right and side expression is a (column) vector.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name ArrayByVectorProduct types explicitly.
  **/
template<typename Lhs, typename Rhs>
class ArrayByVectorProduct: public ExprBase< ArrayByVectorProduct<Lhs, Rhs> >
                           , public TRef<1>
{
  public:
    typedef ExprBase< ArrayByVectorProduct<Lhs, Rhs> > Base;
    typedef typename hidden::Traits<ArrayByVectorProduct>::Type Type;
    typedef typename hidden::Traits<ArrayByVectorProduct>::ConstReturnType ConstReturnType;
    typedef typename hidden::Traits<ArrayByVectorProduct>::Allocator Allocator;
    enum
    {
      structure_ = hidden::Traits<ArrayByVectorProduct>::structure_,
      orient_    = hidden::Traits<ArrayByVectorProduct>::orient_,
      sizeRows_  = hidden::Traits<ArrayByVectorProduct>::sizeRows_,
      sizeCols_  = hidden::Traits<ArrayByVectorProduct>::sizeCols_,
      storage_   = hidden::Traits<ArrayByVectorProduct>::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    ArrayByVectorProduct( const Lhs& lhs, const Rhs& rhs)
                        : Base(), lhs_(lhs), rhs_(rhs)
                        , result_(lhs.sizeRows(), 1, Type(0)) // create result
    {
      if (lhs.cols() != rhs.range())
      { STKRUNTIME_ERROR_NO_ARG(ArrayByVectorProduct, sizes mismatch);}
      result_.shift(lhs_.beginRows(), rhs_.beginCols());
      hidden::ProductDispatcher<Lhs, Rhs, Allocator>::run(lhs, rhs, result_);
    }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return result_.rows();}
    /** @return the columns range */
    inline ColRange const& colsImpl() const { return result_.cols();}

    /** @return the element (i,j) */
    inline ConstReturnType elt2Impl(int i, int j) const { return result_.elt(i, j);}
    /** @return the ith element */
    inline ConstReturnType elt1Impl(int i) const { return result_.elt(i);}
    /** @return the element */
    inline ConstReturnType elt0Impl() const { return result_.elt();}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the right hand side nested expression */
    inline Rhs const& rhs() const { return rhs_; }
    /** @return the right hand side nested expression */
    inline Allocator const& result() const { return result_; }

  protected:
    Lhs const& lhs_;
    Rhs const& rhs_;

  private:
    Allocator result_;
};


namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for the ArrayByArrayProduct class
 */
template< typename Lhs, typename Rhs>
struct Traits< ArrayByArrayProduct<Lhs, Rhs> >
{
  // delegate to ProductTraits all the stuff
  typedef ProductTraits<Lhs, Rhs, Traits<Lhs>::structure_, Traits<Rhs>::structure_> Base;

  enum
  {
    lhs_structure_  = Traits<Lhs>::structure_,
    rhs_structure_  = Traits<Rhs>::structure_,
    structure_      = Base::structure_,
    sizeRows_       = Base::sizeRows_,
    sizeCols_       = Base::sizeCols_,
    orient_         = Base::orient_,
    storage_        = Base::storage_,
    useForRows_     = Base::useForRows_,
    useForCols_     = Base::useForCols_
  };

  typedef typename Base::Type Type;
  typedef typename Base::ConstReturnType ConstReturnType;
  typedef typename Base::Allocator Allocator;
};

} // namespace hidden

/** @ingroup Arrays
  * @brief Generic expression where a product operator is applied to two expressions
  *
  * @tparam Lhs left-hand side type
  * @tparam Rhs right-hand side type
  *
  * This class represents an expression  where a product operator is applied to
  * two expressions.
  *
  * It is the return type of product operator when the left hand side expression
  * is an array and the right and side expression is an array.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name ArrayByArrayProduct types explicitly.
  **/
template<typename Lhs, typename Rhs>
class ArrayByArrayProduct: public ExprBase< ArrayByArrayProduct<Lhs, Rhs> >, public TRef<1>
{
  public:
    enum
    {
      // All the valid cases for ArrayByArray operator
      isValid_ = (  EGAL(Lhs,array2D_)||EGAL(Lhs,square_)
                  ||EGAL(Lhs,lower_triangular_)||EGAL(Lhs,upper_triangular_)
                  ||EGAL(Lhs,symmetric_)||EGAL(Lhs,lower_symmetric_)||EGAL(Lhs,upper_symmetric_)
                 )
                 &&
                 (  EGAL(Rhs,array2D_)||EGAL(Rhs,square_)
                  ||EGAL(Rhs,lower_triangular_)||EGAL(Rhs,upper_triangular_)
                  ||EGAL(Rhs,symmetric_)||EGAL(Rhs,lower_symmetric_)||EGAL(Rhs,upper_symmetric_)
                 ),
      lhs_structure_ = hidden::Traits< ArrayByArrayProduct >::lhs_structure_,
      rhs_structure_ = hidden::Traits< ArrayByArrayProduct >::rhs_structure_,
      structure_     = hidden::Traits< ArrayByArrayProduct >::structure_,
      orient_        = hidden::Traits< ArrayByArrayProduct >::orient_,
      sizeRows_      = hidden::Traits< ArrayByArrayProduct >::sizeRows_,
      sizeCols_      = hidden::Traits< ArrayByArrayProduct >::sizeCols_,
      storage_       = hidden::Traits< ArrayByArrayProduct >::storage_,
      useForRows_    = hidden::Traits< ArrayByArrayProduct >::useForRows_,
      useForCols_    = hidden::Traits< ArrayByArrayProduct >::useForCols_,
    };

    typedef ExprBase< ArrayByArrayProduct<Lhs, Rhs> > Base;
    typedef typename hidden::Traits< ArrayByArrayProduct>::Type Type;
    typedef typename hidden::Traits< ArrayByArrayProduct>::ConstReturnType ConstReturnType;

    typedef typename hidden::Traits< ArrayByArrayProduct < Lhs, Rhs> >::Allocator Allocator;
    typedef hidden::ProductDispatcher<Lhs, Rhs, Allocator, lhs_structure_, rhs_structure_> Dispatcher;

    typedef hidden::BinaryRowsImpl< Lhs, Rhs, sizeRows_, useForRows_ > RowsImpl;
    typedef hidden::BinaryColsImpl< Lhs, Rhs, sizeCols_, useForCols_ > ColsImpl;

    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;
    /** Constructor. Compute the result of the product and store it for further use in result_ */
    ArrayByArrayProduct( Lhs const& lhs, Rhs const& rhs);

    /** @return range of the rows */
    inline RowRange const& rowsImpl() const { return RowsImpl::rowsImpl(lhs_, rhs_);}
    /** @return range of the columns */
    inline ColRange const& colsImpl() const { return ColsImpl::colsImpl(lhs_, rhs_);}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the right hand side nested expression */
    inline Rhs const& rhs() const { return rhs_; }
    /** @return the result */
    inline Allocator const& result() const { return result_; }

    /** @return element (i,j) */
    inline ConstReturnType elt2Impl(int i, int j) const { return result_.elt(i,j);}
    /** @return ith element */
    inline ConstReturnType elt1Impl(int i) const  { return result_.elt(i);}
    /** @return number */
    inline ConstReturnType elt0Impl() const { return result_.elt();}

  protected:
    Lhs const& lhs_;
    Rhs const& rhs_;

  private:
    Allocator result_;
};

template<typename Lhs, typename Rhs>
ArrayByArrayProduct<Lhs,Rhs>::ArrayByArrayProduct( Lhs const& lhs, Rhs const& rhs)
                                                 : Base(), lhs_(lhs), rhs_(rhs)
                                                 , result_(lhs.sizeRows(), rhs.sizeCols(), Type(0))
{
  STK_STATIC_ASSERT_PRODUCT_OPERATOR_MISMATCH( isValid_ );
  if (lhs.cols() != rhs.rows())
  { STKRUNTIME_ERROR_NO_ARG(ArrayByArrayProduct,sizes mismatch for 2D array);}
  result_.shift(lhs_.beginRows(), rhs_.beginCols());
  // general sizes
  Dispatcher::run(lhs, rhs, result_);
}

}  // namespace STK

#undef EGAL

#endif /* STK_PRODUCTOPERATORS_H */
