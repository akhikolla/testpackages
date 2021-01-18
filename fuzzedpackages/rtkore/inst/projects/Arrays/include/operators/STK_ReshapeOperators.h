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

/** @file STK_ReshapeOperators.h
 *  @brief In this file we implement the DiagonalGetterOperator and DiagonalizeOperator classes.
 **/

#ifndef STK_RESHAPEOPERATORS_H
#define STK_RESHAPEOPERATORS_H


#include "Sdk/include/STK_StaticAssert.h"
#include "STK_SlicingOperators.h"

#define EGAL(arg1, arg2) ((arg1::structure_ == int(Arrays::arg2)))

namespace STK
{

// forward declaration
template< typename Lhs> class DiagonalizeOperator;
template< typename Lhs> class DiagonalGetterOperator;
template< typename Lhs> class UpperTriangularizeOperator;
template< typename Lhs> class LowerTriangularizeOperator;
template< typename Lhs> class SymmetrizeOperator;
template< typename Lhs> class UpperSymmetrizeOperator;
template< typename Lhs> class LowerSymmetrizeOperator;

template< typename Lhs> class DiagonalizeAccessor;
template< typename Lhs> class DiagonalGetterAccessor;
template< typename Lhs> class UpperTriangularizeAccessor;
template< typename Lhs> class LowerTriangularizeAccessor;
template< typename Lhs> class SymmetrizeAccessor;
template< typename Lhs> class UpperSymmetrizeAccessor;
template< typename Lhs> class LowerSymmetrizeAccessor;

namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for DiagonalizeOperator operator
 */
template<typename Lhs>
struct Traits< DiagonalizeOperator <Lhs> >
{
  enum
  {
    structure_ = Arrays::diagonal_,
    orient_    = Lhs::orient_,
    sizeRows_  = (Lhs::structure_== int(Arrays::diagonal_) )
                 ? Lhs::sizeRows_ != UnknownSize ? Lhs::sizeRows_ : Lhs::sizeCols_
                 : Lhs::structure_== int(Arrays::point_) ? int(Lhs::sizeCols_) : int(Lhs::sizeRows_),
    sizeCols_  = (Lhs::structure_== int(Arrays::diagonal_) )
                 ? Lhs::sizeRows_ != UnknownSize ? Lhs::sizeRows_ : Lhs::sizeCols_
                 : Lhs::structure_== int(Arrays::point_) ?  Lhs::sizeCols_ : Lhs::sizeRows_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator<DiagonalizeOperator < Lhs> > Row;
  typedef ColOperator<DiagonalizeOperator < Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ConstReturnType ConstReturnType;
};
/** @ingroup hidden
 *  @brief Traits class for DiagonalizeAccessor operator
 */
template<typename Lhs>
struct Traits< DiagonalizeAccessor <Lhs> >
{
  enum
  {
    structure_ = Arrays::diagonal_,
    orient_    = Lhs::orient_,
    sizeRows_  = (Lhs::structure_== int(Arrays::diagonal_) )
                 ? Lhs::sizeRows_ != UnknownSize ? Lhs::sizeRows_ : Lhs::sizeCols_
                 : Lhs::structure_== int(Arrays::point_) ? int(Lhs::sizeCols_) : int(Lhs::sizeRows_),
    sizeCols_  = (Lhs::structure_== int(Arrays::diagonal_) )
                 ? Lhs::sizeRows_ != UnknownSize ? Lhs::sizeRows_ : Lhs::sizeCols_
                 : Lhs::structure_== int(Arrays::point_) ?  Lhs::sizeCols_ : Lhs::sizeRows_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator<DiagonalizeAccessor < Lhs> > Row;
  typedef ColOperator<DiagonalizeAccessor < Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ConstReturnType ConstReturnType;
};

} // namespace hidden

/** @ingroup Arrays
 *  @class DiagonalizeOperator
  *
  * @brief Generic expression when a one dimensional vector/point/idagonal
  * expression is "diagonalized".
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * DiagonalizeOperator operator.
  *
  * This class represents an expression where a DiagonalizeOperator operator is
  * applied to a vector/point/diagonal expression. It is the return type of the
  * diagonalize() const operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name DiagonalizeOperator type explicitly.
  */
template< typename Lhs>
class DiagonalizeOperator: public ExprBase< DiagonalizeOperator< Lhs> >, public TRef<1>
{
  public:
    typedef ExprBase< DiagonalizeOperator< Lhs> > Base;
    typedef typename hidden::Traits< DiagonalizeOperator<Lhs> >::Type Type;
    typedef typename hidden::Traits< DiagonalizeOperator<Lhs> >::ConstReturnType ConstReturnType;
    enum
    {
        structure_ = hidden::Traits< DiagonalizeOperator<Lhs> >::structure_,
        orient_    = hidden::Traits< DiagonalizeOperator<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< DiagonalizeOperator<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< DiagonalizeOperator<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< DiagonalizeOperator<Lhs> >::storage_,
        // this is safe as we can use DiagonalizeOperator only on 1D container
        size_      = (sizeRows_ != UnknownSize) ? sizeRows_ : sizeCols_
    };
    /** Type of the Range for the rows */
    typedef TRange<size_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<size_> ColRange;
    /** Constructor */
    inline DiagonalizeOperator( Lhs const& lhs): Base(), lhs_(lhs)
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Lhs);
    }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs().range();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs().range();}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j indexes of the element to get
     **/
    inline ConstReturnType elt2Impl(int i, int j) const
    {
#ifdef STK_BOUNDS_CHECK
      if (i != j) { STKOUT_OF_RANGE_2ARG(DiagonalizeOperator::elt, i, j, i != j);}
#endif
      return (lhs_.elt(i));
    }
    /** @return a constant reference on the ith element of the expression
     *  @param i index of the element to get
     **/
    inline ConstReturnType elt1Impl(int i) const { return (lhs_.elt(i));}
    /** @return a constant reference on the element of the expression */
    inline ConstReturnType elt0Impl() const { return (lhs_.elt());}

  protected:
    Lhs const& lhs_;
};

/** @ingroup Arrays
 *  @class DiagonalizeAccessor
  *
  * @brief Generic expression when a one dimensional vector/point/idagonal
  * expression is "diagonalized".
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * DiagonalizeAccessor operator.
  *
  * This class represents an expression where a DiagonalizeAccessor operator is
  * applied to a vector/point/diagonal array. It is the return type of the
  * diagonalizeize() operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name DiagonalizeAccessor type explicitly.
  */
template< typename Lhs>
class DiagonalizeAccessor: public ArrayBase< DiagonalizeAccessor< Lhs> >, public TRef<1>
{
  public:
    typedef ArrayBase< DiagonalizeAccessor< Lhs> > Base;
    typedef typename hidden::Traits< DiagonalizeAccessor<Lhs> >::Type Type;
    typedef typename hidden::Traits< DiagonalizeAccessor<Lhs> >::ConstReturnType ConstReturnType;
    enum
    {
        structure_ = hidden::Traits< DiagonalizeAccessor<Lhs> >::structure_,
        orient_    = hidden::Traits< DiagonalizeAccessor<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< DiagonalizeAccessor<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< DiagonalizeAccessor<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< DiagonalizeAccessor<Lhs> >::storage_,
        // this is safe as we can use DiagonalizeOperator only on 1D container
        size_      = (sizeRows_ != UnknownSize) ? sizeRows_ : sizeCols_
    };
    /** Type of the Range for the rows */
    typedef TRange<size_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<size_> ColRange;
    /** Constructor */
    inline DiagonalizeAccessor(Lhs& lhs): Base(), lhs_(lhs)
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Lhs);
    }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs().range();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs().range();}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j indexes of the element to get
     **/
    inline ConstReturnType elt2Impl(int i, int j) const { return (lhs_.elt(i, j));}
    /** @return a constant reference on the ith element of the expression
     *  @param i index of the element to get
     **/
    inline ConstReturnType elt1Impl(int i) const { return (lhs_.elt(i));}
    /** @return a constant reference on the element of the expression */
    inline ConstReturnType elt0Impl() const { return (lhs_.elt());}

    /** @return a reference on the (i,j) element of the expression.
     *  @param i, j indexes of the element to get
     **/
    inline Type& elt2Impl(int i, int j)
    {
#ifdef STK_BOUNDS_CHECK
      if (i != j) { STKOUT_OF_RANGE_2ARG(DiagonalizeAccessor::elt, i, j, i != j);}
#endif
      return (lhs_.elt(i));
    }
    /** @return a reference on the ith element of the expression
     *  @param i index of the element to get
     **/
    inline Type& elt1Impl(int i) { return (lhs_.elt(i));}
    /** @return a reference on the element of the expression */
    inline Type& elt0Impl() { return (lhs_.elt());}

  protected:
    Lhs& lhs_;
};


namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for DiagonalGetterOperator operator
 */
template<typename Lhs>
struct Traits< DiagonalGetterOperator <Lhs> >
{
  enum
  {
    structure_ = Arrays::diagonal_,
    orient_    = Lhs::orient_,
    sizeRows_  = ((Lhs::sizeRows_ < Lhs::sizeCols_)) ?  Lhs::sizeRows_ : Lhs::sizeCols_,
    sizeCols_  = sizeRows_,
    storage_   = Lhs::storage_,
    isValid_   = ( hidden::Traits<Lhs>::structure_==(int)Arrays::array2D_)
                ||(hidden::Traits<Lhs>::structure_==(int)Arrays::square_)
                ||(hidden::Traits<Lhs>::structure_==(int)Arrays::diagonal_)
                ||(hidden::Traits<Lhs>::structure_==(int)Arrays::lower_triangular_)
                ||(hidden::Traits<Lhs>::structure_==(int)Arrays::upper_triangular_)
                ||(hidden::Traits<Lhs>::structure_==(int)Arrays::symmetric_)
                ||(hidden::Traits<Lhs>::structure_==(int)Arrays::lower_symmetric_)
                ||(hidden::Traits<Lhs>::structure_==(int)Arrays::upper_symmetric_)
                ||(hidden::Traits<Lhs>::structure_==(int)Arrays::number_),
    use_       = (sizeRows_ != UnknownSize) ? Arrays::useLhsSize_ : Arrays::useLhsOtherSize_,
    size_      = (sizeRows_ != UnknownSize) ? sizeRows_ : sizeCols_
  };
  typedef RowOperator<DiagonalGetterOperator < Lhs> > Row;
  typedef ColOperator<DiagonalGetterOperator < Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ConstReturnType ConstReturnType;
};

/** @ingroup hidden
 *  @brief Traits class for DiagonalGetterAccessor operator
 */
template<typename Lhs>
struct Traits< DiagonalGetterAccessor <Lhs> >
{
  enum
  {
    structure_ = Arrays::diagonal_,
    orient_    = Lhs::orient_,
    sizeRows_  = ((Lhs::sizeRows_ < Lhs::sizeCols_)) ?  Lhs::sizeRows_ : Lhs::sizeCols_,
    sizeCols_  = sizeRows_,
    storage_   = Lhs::storage_,
    isValid_   = ( hidden::Traits<Lhs>::structure_==(int)Arrays::array2D_)
                ||(hidden::Traits<Lhs>::structure_==(int)Arrays::square_)
                ||(hidden::Traits<Lhs>::structure_==(int)Arrays::diagonal_)
                ||(hidden::Traits<Lhs>::structure_==(int)Arrays::lower_triangular_)
                ||(hidden::Traits<Lhs>::structure_==(int)Arrays::upper_triangular_)
                ||(hidden::Traits<Lhs>::structure_==(int)Arrays::symmetric_)
                ||(hidden::Traits<Lhs>::structure_==(int)Arrays::lower_symmetric_)
                ||(hidden::Traits<Lhs>::structure_==(int)Arrays::upper_symmetric_)
                ||(hidden::Traits<Lhs>::structure_==(int)Arrays::number_),
    use_       = (sizeRows_ != UnknownSize) ? Arrays::useLhsSize_ : Arrays::useLhsOtherSize_,
    size_      = (sizeRows_ != UnknownSize) ? sizeRows_ : sizeCols_
  };
  typedef RowOperator<DiagonalGetterAccessor < Lhs> > Row;
  typedef ColOperator<DiagonalGetterAccessor < Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ConstReturnType ConstReturnType;
};

/** @ingroup hidden
 *  Allow to disambiguate the call to range to use */
template<typename Lhs, int Size_, int use_>
struct DiagonalRangeImpl
{
  /** Type of the Range for the rows */
  typedef TRange<Size_> RangeType;
  /**  @return the range of the rows */
  inline static RangeType const& rangeImpl(Lhs const& lhs) { return lhs.range();}
};

template<typename Lhs, int Size_>
struct DiagonalRangeImpl<Lhs, Size_, Arrays::useLhsSize_>
{
  /** Type of the Range for the rows */
  typedef TRange<Size_> RangeType;
  /**  @return the range of the rows */
  inline static RangeType const& rangeImpl(Lhs const& lhs) { return lhs.rows();}
};

template<typename Lhs, int Size_>
struct DiagonalRangeImpl<Lhs, Size_, Arrays::useLhsOtherSize_>
{
  /** Type of the Range for the rows */
  typedef TRange<Size_> RangeType;
  /**  @return the range of the rows */
  inline static RangeType const& rangeImpl(Lhs const& lhs) { return lhs.cols();}
};


} // namespace hidden


/** @ingroup Arrays
 *  @class DiagonalGetterOperator
  *
  * @brief Generic expression when we want to get the diagonal of a
  * two-dimensional square expression.
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * DiagonalGetterOperator operator.
  *
  * This class represents an expression where a diagonal operator is applied to
  * an expression. It is the return type of the diagonal operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name DiagonalGetterOperator type explicitly.
  */
template< typename Lhs>
class DiagonalGetterOperator: public ExprBase< DiagonalGetterOperator< Lhs> >, public TRef<1>
{
  public:
    typedef ExprBase< DiagonalGetterOperator< Lhs> > Base;
    typedef typename hidden::Traits< DiagonalGetterOperator<Lhs> >::Type Type;
    typedef typename hidden::Traits< DiagonalGetterOperator<Lhs> >::ConstReturnType ConstReturnType;

    enum
    {
        structure_ = hidden::Traits< DiagonalGetterOperator<Lhs> >::structure_,
        orient_    = hidden::Traits< DiagonalGetterOperator<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< DiagonalGetterOperator<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< DiagonalGetterOperator<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< DiagonalGetterOperator<Lhs> >::storage_,
        isValid_   = hidden::Traits< DiagonalGetterOperator<Lhs> >::isValid_,
        use_       = hidden::Traits< DiagonalGetterOperator<Lhs> >::use_,
        size_      = hidden::Traits< DiagonalGetterOperator<Lhs> >::size_
    };
    typedef hidden::DiagonalRangeImpl<Lhs, size_, use_> RangeImpl;
    /** Type of the Range for the rows */
    typedef TRange<size_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<size_> ColRange;

    /** Constructor */
    inline DiagonalGetterOperator( Lhs const& lhs): Base(), lhs_(lhs)
    {
      STK_STATIC_ASSERT(isValid_,YOU_TRIED_CALLING_A_MATRIX_METHOD_ON_A_VECTOR);
#ifdef STK_BOUNDS_CHECK
      if (lhs.rows()!=lhs.cols())
        STKRUNTIME_ERROR_NO_ARG(DiagonalGetterOperatorBase,lhs.rows()!=lhs.cols());
#endif
    }
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return RangeImpl::rangeImpl(lhs_);}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return RangeImpl::rangeImpl(lhs_);}

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j index of the row and of the column
     **/
    inline ConstReturnType elt2Impl(int i, int j) const { return (lhs_.elt(i, j));}
    /** @return a constant reference on the ith element of the expression
     *  @param i index of the element to get
     **/
    inline ConstReturnType elt1Impl(int i) const { return (lhs_.elt(i,i));}
    /** @return a constant reference on the element of the expression */
    inline ConstReturnType elt0Impl() const { return (lhs_.elt());}

  protected:
    Lhs const& lhs_;
};

/** @ingroup Arrays
 *  @class DiagonalGetterAccessor
  *
  * @brief Generic expression when we want to get the diagonal of a
  * two-dimensional square expression.
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * DiagonalGetterOperator operator.
  *
  * This class represents an expression where a diagonal operator is applied to
  * an expression. It is the return type of the diagonal operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name DiagonalGetterAccessor type explicitly.
  */
template< typename Lhs>
class DiagonalGetterAccessor: public ArrayBase< DiagonalGetterAccessor< Lhs> >, public TRef<1>
{
  public:
    typedef ArrayBase< DiagonalGetterAccessor< Lhs> > Base;
    typedef typename hidden::Traits< DiagonalGetterAccessor<Lhs> >::Type Type;
    typedef typename hidden::Traits< DiagonalGetterAccessor<Lhs> >::ConstReturnType ConstReturnType;

    enum
    {
        structure_ = hidden::Traits< DiagonalGetterAccessor<Lhs> >::structure_,
        orient_    = hidden::Traits< DiagonalGetterAccessor<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< DiagonalGetterAccessor<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< DiagonalGetterAccessor<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< DiagonalGetterAccessor<Lhs> >::storage_,
        isValid_   = hidden::Traits< DiagonalGetterAccessor<Lhs> >::isValid_,
        use_       = hidden::Traits< DiagonalGetterAccessor<Lhs> >::use_,
        size_      = hidden::Traits< DiagonalGetterAccessor<Lhs> >::size_
    };
    typedef hidden::DiagonalRangeImpl<Lhs, size_, use_> RangeImpl;
    /** Type of the Range for the rows */
    typedef TRange<size_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<size_> ColRange;

    /** Constructor */
    inline DiagonalGetterAccessor( Lhs& lhs): Base(), lhs_(lhs)
    {
      STK_STATIC_ASSERT(isValid_,YOU_TRIED_CALLING_A_MATRIX_METHOD_ON_A_VECTOR);
#ifdef STK_BOUNDS_CHECK
      if (lhs.rows()!=lhs.cols())
        STKRUNTIME_ERROR_NO_ARG(DiagonalGetterOperatorBase,lhs.rows()!=lhs.cols());
#endif
    }
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return RangeImpl::rangeImpl(lhs_);}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return RangeImpl::rangeImpl(lhs_);}

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j index of the row and of the column
     **/
    inline ConstReturnType elt2Impl(int i, int j) const { return (lhs_.elt(i, j));}
    /** @return a constant reference on the ith element of the expression
     *  @param i index of the element to get
     **/
    inline ConstReturnType elt1Impl(int i) const { return (lhs_.elt(i,i));}
    /** @return a constant reference on the element of the expression */
    inline ConstReturnType elt0Impl() const { return (lhs_.elt());}

    /** @return a reference on the (i,j) element of the expression.
     *  @param i, j index of the row and of the column
     **/
    inline Type& elt2Impl(int i, int j) { return (lhs_.elt(i, j));}
    /** @return a reference on the ith element of the expression
     *  @param i index of the element to get
     **/
    inline Type& elt1Impl(int i) { return (lhs_.elt(i,i));}
    /** @return a reference on the element of the expression */
    inline Type& elt0Impl() { return (lhs_.elt());}

  protected:
    Lhs& lhs_;
};

namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for UpperTriangularizeOperator operator
 */
template<typename Lhs>
struct Traits< UpperTriangularizeOperator<Lhs> >
{
  enum
  {
    structure_ = Arrays::upper_triangular_,
    orient_    = Lhs::orient_,
    sizeRows_  = Lhs::sizeRows_,
    sizeCols_  = Lhs::sizeCols_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator< UpperTriangularizeOperator< Lhs> > Row;
  typedef ColOperator< UpperTriangularizeOperator< Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ConstReturnType ConstReturnType;
};

/** @ingroup hidden
 *  @brief Traits class for UpperTriangularizeAccessor operator
 */
template<typename Lhs>
struct Traits< UpperTriangularizeAccessor<Lhs> >
{
  enum
  {
    structure_ = Arrays::upper_triangular_,
    orient_    = Lhs::orient_,
    sizeRows_  = Lhs::sizeRows_,
    sizeCols_  = Lhs::sizeCols_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator< UpperTriangularizeAccessor< Lhs> > Row;
  typedef ColOperator< UpperTriangularizeAccessor< Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ConstReturnType ConstReturnType;
};

} // namespace hidden

/** @ingroup Arrays
 *  @class UpperTriangularizeOperator
  *
  * @brief Generic expression when we want to get the upper-part of a
  * two-dimensional expression.
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * UpperTriangularizeOperator operator.
  *
  * This class represents an expression where an UpperTriangularizeOperator
  * operator is applied to an expression. It is the return type of the
  * upperTriangularize operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name UpperTriangularizeOperator type explicitly.
  */
template< typename Lhs>
class UpperTriangularizeOperator: public ExprBase< UpperTriangularizeOperator< Lhs> >, public TRef<1>
{
  public:
    typedef ExprBase< UpperTriangularizeOperator< Lhs> > Base;
    typedef typename hidden::Traits< UpperTriangularizeOperator<Lhs> >::Type Type;
    typedef typename hidden::Traits< UpperTriangularizeOperator<Lhs> >::ConstReturnType ConstReturnType;

    enum
    {
        structure_ = hidden::Traits< UpperTriangularizeOperator<Lhs> >::structure_,
        orient_    = hidden::Traits< UpperTriangularizeOperator<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< UpperTriangularizeOperator<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< UpperTriangularizeOperator<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< UpperTriangularizeOperator<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** Constructor */
    inline UpperTriangularizeOperator( Lhs const& lhs): Base(), lhs_(lhs) {}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs().rows();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs().cols();}

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j index of the row and of the column
     **/
    inline ConstReturnType elt2Impl(int i, int j) const
    {
#ifdef STK_BOUNDS_CHECK
      if (j<i)
        STKRUNTIME_ERROR_2ARG(UpperTriangularizeOperator::elt2Impl,i,j,use of the lower part);
#endif
      return (lhs_.elt(i, j));
    }

  protected:
    Lhs const& lhs_;
};

/** @ingroup Arrays
 *  @class UpperTriangularizeOperator
  *
  * @brief Generic expression when we want to get the upper-part of a
  * two-dimensional expression.
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * UpperTriangularizeAccessor operator.
  *
  * This class represents an expression where an UpperTriangularizeOperator
  * operator is applied to an expression. It is the return type of the
  * upperTriangularize operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name UpperTriangularizeAccessor type explicitly.
  */
template< typename Lhs>
class UpperTriangularizeAccessor: public ArrayBase< UpperTriangularizeAccessor< Lhs> >, public TRef<1>
{
  public:
    typedef ArrayBase< UpperTriangularizeAccessor< Lhs> > Base;
    typedef typename hidden::Traits< UpperTriangularizeAccessor<Lhs> >::Type Type;
    typedef typename hidden::Traits< UpperTriangularizeAccessor<Lhs> >::ConstReturnType ConstReturnType;

    enum
    {
        structure_ = hidden::Traits< UpperTriangularizeAccessor<Lhs> >::structure_,
        orient_    = hidden::Traits< UpperTriangularizeAccessor<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< UpperTriangularizeAccessor<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< UpperTriangularizeAccessor<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< UpperTriangularizeAccessor<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** Constructor */
    inline UpperTriangularizeAccessor( Lhs& lhs): Base(), lhs_(lhs) {}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs().rows();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs().cols();}

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j index of the row and of the column
     **/
    inline ConstReturnType elt2Impl(int i, int j) const
    {
#ifdef STK_BOUNDS_CHECK
      if (j<i)
        STKRUNTIME_ERROR_2ARG(UpperTriangularizeOperator::elt2Impl,i,j,use of the lower part);
#endif
      return (lhs_.elt(i, j));
    }

    /** @return a reference on the (i,j) element of the expression.
     *  @param i, j indexes of the element to get
     **/
    inline Type& elt2Impl(int i, int j) {
#ifdef STK_BOUNDS_CHECK
      if (i>j)
        STKRUNTIME_ERROR_2ARG(UpperTriangularizeAccessor::elt2Impl,i,j,use of the lower part);
#endif
      return (lhs_.elt(i, j));
    }

  protected:
    Lhs& lhs_;
};

namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for UpperTriangularizeOperator operator
 */
template<typename Lhs>
struct Traits< LowerTriangularizeOperator<Lhs> >
{
  enum
  {
    structure_ = Arrays::lower_triangular_,
    orient_    = Lhs::orient_,
    sizeRows_  = Lhs::sizeRows_,
    sizeCols_  = Lhs::sizeCols_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator< LowerTriangularizeOperator< Lhs> > Row;
  typedef ColOperator< LowerTriangularizeOperator< Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ConstReturnType ConstReturnType;
};
/** @ingroup hidden
 *  @brief Traits class for LowerTriangularizeAccessor operator
 */
template<typename Lhs>
struct Traits< LowerTriangularizeAccessor<Lhs> >
{
  enum
  {
    structure_ = Arrays::lower_triangular_,
    orient_    = Lhs::orient_,
    sizeRows_  = Lhs::sizeRows_,
    sizeCols_  = Lhs::sizeCols_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator< LowerTriangularizeAccessor< Lhs> > Row;
  typedef ColOperator< LowerTriangularizeAccessor< Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ConstReturnType ConstReturnType;
};

} // namespace hidden

/** @ingroup Arrays
 *  @class LowerTriangularizeOperator
  *
  * @brief Generic expression when we want to get the lower-part of a
  * two-dimensional expression.
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * LowerTriangularizeOperator operator.
  *
  * This class represents an expression where an LowerTriangularizeOperator
  * operator is applied to an expression. It is the return type of the
  * lowerTriangularize operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name LowerTriangularizeOperator type explicitly.
  */
template< typename Lhs>
class LowerTriangularizeOperator: public ExprBase< LowerTriangularizeOperator< Lhs> >, public TRef<1>
{
  public:
    typedef ExprBase< LowerTriangularizeOperator< Lhs> > Base;
    typedef typename hidden::Traits< LowerTriangularizeOperator<Lhs> >::Type Type;
    typedef typename hidden::Traits< LowerTriangularizeOperator<Lhs> >::ConstReturnType ConstReturnType;

    enum
    {
        structure_ = hidden::Traits< LowerTriangularizeOperator<Lhs> >::structure_,
        orient_    = hidden::Traits< LowerTriangularizeOperator<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< LowerTriangularizeOperator<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< LowerTriangularizeOperator<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< LowerTriangularizeOperator<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** Constructor */
    inline LowerTriangularizeOperator( Lhs const& lhs): Base(), lhs_(lhs) {}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs().rows();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs().cols();}

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j index of the row and of the column
     **/
    inline ConstReturnType elt2Impl(int i, int j) const
    {
#ifdef STK_BOUNDS_CHECK
      if (j>i)
        STKRUNTIME_ERROR_2ARG(LowerTriangularizeOperator::elt2Impl,i,j,use of the upper part);
#endif
      return (lhs_.elt(i, j));
    }

  protected:
    Lhs const& lhs_;
};


/** @ingroup Arrays
 *  @class LowerTriangularizeAccessor
  *
  * @brief Generic expression when we want to get the lower-part of a
  * two-dimensional expression.
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * LowerTriangularizeAccessor operator.
  *
  * This class represents an expression where an LowerTriangularizeOperator
  * operator is applied to an expression. It is the return type of the
  * lowerTriangularize operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name LowerTriangularizeAccessor type explicitly.
  */
template< typename Lhs>
class LowerTriangularizeAccessor: public ArrayBase< LowerTriangularizeAccessor< Lhs> >, public TRef<1>
{
  public:
    typedef ArrayBase< LowerTriangularizeAccessor< Lhs> > Base;
    typedef typename hidden::Traits< LowerTriangularizeAccessor<Lhs> >::Type Type;
    typedef typename hidden::Traits< LowerTriangularizeAccessor<Lhs> >::ConstReturnType ConstReturnType;

    enum
    {
        structure_ = hidden::Traits< LowerTriangularizeAccessor<Lhs> >::structure_,
        orient_    = hidden::Traits< LowerTriangularizeAccessor<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< LowerTriangularizeAccessor<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< LowerTriangularizeAccessor<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< LowerTriangularizeAccessor<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** Constructor */
    inline LowerTriangularizeAccessor( Lhs& lhs): Base(), lhs_(lhs) {}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs().rows();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs().cols();}

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j index of the row and of the column
     **/
    inline ConstReturnType elt2Impl(int i, int j) const
    {
#ifdef STK_BOUNDS_CHECK
      if (j>i)
        STKRUNTIME_ERROR_2ARG(LowerTriangularizeOperator::elt2Impl,i,j,use of the upper part);
#endif
      return (lhs_.elt(i, j));
    }

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j indexes of the element to get
     **/
    inline Type& elt2Impl(int i, int j)
    {
#ifdef STK_BOUNDS_CHECK
      if (j>i)
        STKRUNTIME_ERROR_2ARG(LowerTriangularizeAccessor::elt2Impl,i,j,use of the upper part);
#endif
      return (lhs_.elt(i, j));
    }

  protected:
    Lhs& lhs_;
};

namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for SymmetrizeOperator operator
 */
template<typename Lhs>
struct Traits< SymmetrizeOperator<Lhs> >
{
  enum
  {
    structure_ = Arrays::symmetric_,
    orient_    = Lhs::orient_,
    sizeRows_  = Lhs::sizeRows_,
    sizeCols_  = Lhs::sizeCols_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator< SymmetrizeOperator< Lhs> > Row;
  typedef ColOperator< SymmetrizeOperator< Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ConstReturnType ConstReturnType;
};
/** @ingroup hidden
 *  @brief Traits class for SymmetrizeAccessor operator
 */
template<typename Lhs>
struct Traits< SymmetrizeAccessor<Lhs> >
{
  enum
  {
    structure_ = Arrays::symmetric_,
    orient_    = Lhs::orient_,
    sizeRows_  = Lhs::sizeRows_,
    sizeCols_  = Lhs::sizeCols_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator< SymmetrizeAccessor< Lhs> > Row;
  typedef ColOperator< SymmetrizeAccessor< Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ConstReturnType ConstReturnType;
};

} // namespace hidden

/** @ingroup Arrays
 *  @class SymmetrizeOperator
  *
  * @brief Generic expression when we want to get the upper-part of a
  * two-dimensional symmetric expression.
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * SymmetrizeOperator operator.
  *
  * This class represents an expression where an SymmetrizeOperator
  * operator is applied to an expression. It is the return type of the
  * symmetrize operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name SymmetrizeOperator type explicitly.
  */
template< typename Lhs>
class SymmetrizeOperator: public ExprBase< SymmetrizeOperator< Lhs> >, public TRef<1>
{
  public:
    typedef ExprBase< SymmetrizeOperator< Lhs> > Base;
    typedef typename hidden::Traits< SymmetrizeOperator<Lhs> >::Type Type;
    typedef typename hidden::Traits< SymmetrizeOperator<Lhs> >::ConstReturnType ConstReturnType;

    enum
    {
        structure_ = hidden::Traits< SymmetrizeOperator<Lhs> >::structure_,
        orient_    = hidden::Traits< SymmetrizeOperator<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< SymmetrizeOperator<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< SymmetrizeOperator<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< SymmetrizeOperator<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** Constructor */
    inline SymmetrizeOperator( Lhs const& lhs): Base(), lhs_(lhs)
    { STK_STATIC_ASSERT_TWO_DIMENSIONS_ONLY(Lhs);}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs().rows();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs().cols();}

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j index of the row and of the column
     **/
    inline ConstReturnType elt2Impl(int i, int j) const
    { return (lhs_.elt(i, j));}

  protected:
    Lhs const& lhs_;
};

/** @ingroup Arrays
 *  @class SymmetrizeAccessor
  *
  * @brief Generic expression when we want to get the upper-part of a
  * two-dimensional symmetric expression.
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * SymmetrizeAccessor operator.
  *
  * This class represents an expression where an SymmetrizeOperator
  * operator is applied to an expression. It is the return type of the
  * symmetrize operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name SymmetrizeAccessor type explicitly.
  */
template< typename Lhs>
class SymmetrizeAccessor: public ArrayBase< SymmetrizeAccessor< Lhs> >, public TRef<1>
{
  public:
    typedef ArrayBase< SymmetrizeAccessor< Lhs> > Base;
    typedef typename hidden::Traits< SymmetrizeAccessor<Lhs> >::Type Type;
    typedef typename hidden::Traits< SymmetrizeAccessor<Lhs> >::ConstReturnType ConstReturnType;

    enum
    {
        structure_ = hidden::Traits< SymmetrizeAccessor<Lhs> >::structure_,
        orient_    = hidden::Traits< SymmetrizeAccessor<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< SymmetrizeAccessor<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< SymmetrizeAccessor<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< SymmetrizeAccessor<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** Constructor */
    inline SymmetrizeAccessor( Lhs& lhs): Base(), lhs_(lhs)
    { STK_STATIC_ASSERT_TWO_DIMENSIONS_ONLY(Lhs);}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs().rows();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs().cols();}

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j index of the row and of the column
     **/
    inline ConstReturnType elt2Impl(int i, int j) const
    { return (lhs_.elt(i, j));}

    /** @return a reference on the (i,j) element of the expression.
     *  @param i, j indexes of the element to get
     **/
    inline Type& elt2Impl(int i, int j) { return (lhs_.elt(i, j));}

  protected:
    Lhs& lhs_;
};

namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for UpperSymmetrizeOperator operator
 */
template<typename Lhs>
struct Traits< UpperSymmetrizeOperator<Lhs> >
{
  enum
  {
    structure_ = Arrays::upper_symmetric_,
    orient_    = Lhs::orient_,
    sizeRows_  = Lhs::sizeRows_,
    sizeCols_  = Lhs::sizeCols_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator< UpperSymmetrizeOperator< Lhs> > Row;
  typedef ColOperator< UpperSymmetrizeOperator< Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ConstReturnType ConstReturnType;
};
/** @ingroup hidden
 *  @brief Traits class for UpperSymmetrizeAccessor operator
 */
template<typename Lhs>
struct Traits< UpperSymmetrizeAccessor<Lhs> >
{
  enum
  {
    structure_ = Arrays::upper_symmetric_,
    orient_    = Lhs::orient_,
    sizeRows_  = Lhs::sizeRows_,
    sizeCols_  = Lhs::sizeCols_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator< UpperSymmetrizeAccessor< Lhs> > Row;
  typedef ColOperator< UpperSymmetrizeAccessor< Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ConstReturnType ConstReturnType;
};

} // namespace hidden

/** @ingroup Arrays
 *  @class UpperSymmetrizeOperator
  *
  * @brief Generic expression when we want to get the upper-part of a
  * two-dimensional symmetric expression.
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * UpperSymmetrizeOperator operator.
  *
  * This class represents an expression where an UpperSymmetrizeOperator
  * operator is applied to an expression. It is the return type of the
  * upperSymmetrize operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name UpperSymmetrizeOperator type explicitly.
  */
template< typename Lhs>
class UpperSymmetrizeOperator: public ExprBase< UpperSymmetrizeOperator< Lhs> >, public TRef<1>
{
  public:
    typedef ExprBase< UpperSymmetrizeOperator< Lhs> > Base;
    typedef typename hidden::Traits< UpperSymmetrizeOperator<Lhs> >::Type Type;
    typedef typename hidden::Traits< UpperSymmetrizeOperator<Lhs> >::ConstReturnType ConstReturnType;

    enum
    {
        structure_ = hidden::Traits< UpperSymmetrizeOperator<Lhs> >::structure_,
        orient_    = hidden::Traits< UpperSymmetrizeOperator<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< UpperSymmetrizeOperator<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< UpperSymmetrizeOperator<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< UpperSymmetrizeOperator<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** Constructor */
    inline UpperSymmetrizeOperator( Lhs const& lhs): Base(), lhs_(lhs)
    { STK_STATIC_ASSERT_TWO_DIMENSIONS_ONLY(Lhs);}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs().rows();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs().cols();}

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j index of the row and of the column
     **/
    inline ConstReturnType elt2Impl(int i, int j) const
    { return ((j<i) ? lhs_.elt(j, i) : lhs_.elt(i, j));}

  protected:
    Lhs const& lhs_;
};

/** @ingroup Arrays
 *  @class UpperSymmetrizeAccessor
  *
  * @brief Generic expression when we want to get the upper-part of a
  * two-dimensional symmetric expression.
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * UpperSymmetrizeAccessor operator.
  *
  * This class represents an expression where an UpperSymmetrizeOperator
  * operator is applied to an expression. It is the return type of the
  * upperSymmetrize operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name UpperSymmetrizeAccessor type explicitly.
  */
template< typename Lhs>
class UpperSymmetrizeAccessor: public ArrayBase< UpperSymmetrizeAccessor< Lhs> >, public TRef<1>
{
  public:
    typedef ArrayBase< UpperSymmetrizeAccessor< Lhs> > Base;
    typedef typename hidden::Traits< UpperSymmetrizeAccessor<Lhs> >::Type Type;
    typedef typename hidden::Traits< UpperSymmetrizeAccessor<Lhs> >::ConstReturnType ConstReturnType;

    enum
    {
        structure_ = hidden::Traits< UpperSymmetrizeAccessor<Lhs> >::structure_,
        orient_    = hidden::Traits< UpperSymmetrizeAccessor<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< UpperSymmetrizeAccessor<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< UpperSymmetrizeAccessor<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< UpperSymmetrizeAccessor<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** Constructor */
    inline UpperSymmetrizeAccessor( Lhs& lhs): Base(), lhs_(lhs)
    { STK_STATIC_ASSERT_TWO_DIMENSIONS_ONLY(Lhs);}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs().rows();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs().cols();}

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j index of the row and of the column
     **/
    inline ConstReturnType elt2Impl(int i, int j) const
    { return ((j<i) ? lhs_.elt(j, i) : lhs_.elt(i, j));}

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j indexes of the element to get
     **/
    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j indexes of the element to get
     **/
    inline Type& elt2Impl(int i, int j)
    { return ((j<i) ? lhs_.elt(j, i) : lhs_.elt(i, j));}

  protected:
    Lhs& lhs_;
};

namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for LowerSymmetrizeOperator operator
 */
template<typename Lhs>
struct Traits< LowerSymmetrizeOperator<Lhs> >
{
  enum
  {
    structure_ = Arrays::lower_symmetric_,
    orient_    = Lhs::orient_,
    sizeRows_  = Lhs::sizeRows_,
    sizeCols_  = Lhs::sizeCols_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator< LowerSymmetrizeOperator< Lhs> > Row;
  typedef ColOperator< LowerSymmetrizeOperator< Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ConstReturnType ConstReturnType;
};
/** @ingroup hidden
 *  @brief Traits class for LowerSymmetrizeAccessor operator
 */
template<typename Lhs>
struct Traits< LowerSymmetrizeAccessor<Lhs> >
{
  enum
  {
    structure_ = Arrays::lower_symmetric_,
    orient_    = Lhs::orient_,
    sizeRows_  = Lhs::sizeRows_,
    sizeCols_  = Lhs::sizeCols_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator< LowerSymmetrizeAccessor< Lhs> > Row;
  typedef ColOperator< LowerSymmetrizeAccessor< Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ConstReturnType ConstReturnType;
};

} // end namespace hidden

/** @ingroup Arrays
 *  @class LowerSymmetrizeOperator
  *
  * @brief Generic expression when we want to get the lower-part of a
  * two-dimensional symmetric expression.
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * LowerSymmetrizeOperator operator.
  *
  * This class represents an expression where an LowerTriangularizeOperator
  * operator is applied to an expression. It is the return type of the
  * lowerSymmetrize operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name LowerSymmetrizeOperator type explicitly.
  */
template< typename Lhs>
class LowerSymmetrizeOperator: public ExprBase< LowerSymmetrizeOperator< Lhs> >, public TRef<1>
{
  public:
    typedef ExprBase< LowerSymmetrizeOperator< Lhs> > Base;
    typedef typename hidden::Traits< LowerSymmetrizeOperator<Lhs> >::Type Type;
    typedef typename hidden::Traits< LowerSymmetrizeOperator<Lhs> >::ConstReturnType ConstReturnType;

    enum
    {
        structure_ = hidden::Traits< LowerSymmetrizeOperator<Lhs> >::structure_,
        orient_    = hidden::Traits< LowerSymmetrizeOperator<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< LowerSymmetrizeOperator<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< LowerSymmetrizeOperator<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< LowerSymmetrizeOperator<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** Constructor */
    inline LowerSymmetrizeOperator( Lhs const& lhs): Base(), lhs_(lhs)
    { STK_STATIC_ASSERT_TWO_DIMENSIONS_ONLY(Lhs);}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs().rows();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs().cols();}

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j index of the row and of the column
     **/
    inline ConstReturnType elt2Impl(int i, int j) const
    { return ((j>i) ? lhs_.elt(j, i) : lhs_.elt(i, j));}

  protected:
    Lhs const& lhs_;
};


/** @ingroup Arrays
 *  @class LowerSymmetrizeAccessor
  *
  * @brief Generic expression when we want to get the lower-part of a
  * two-dimensional symmetric expression.
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * LowerSymmetrizeAccessor operator.
  *
  * This class represents an expression where an LowerTriangularizeOperator
  * operator is applied to an expression. It is the return type of the
  * lowerSymmetrize operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name LowerSymmetrizeAccessor type explicitly.
  */
template< typename Lhs>
class LowerSymmetrizeAccessor: public ArrayBase< LowerSymmetrizeAccessor< Lhs> >, public TRef<1>
{
  public:
    typedef ArrayBase< LowerSymmetrizeAccessor< Lhs> > Base;
    typedef typename hidden::Traits< LowerSymmetrizeAccessor<Lhs> >::Type Type;
    typedef typename hidden::Traits< LowerSymmetrizeAccessor<Lhs> >::ConstReturnType ConstReturnType;

    enum
    {
        structure_ = hidden::Traits< LowerSymmetrizeAccessor<Lhs> >::structure_,
        orient_    = hidden::Traits< LowerSymmetrizeAccessor<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< LowerSymmetrizeAccessor<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< LowerSymmetrizeAccessor<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< LowerSymmetrizeAccessor<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** Constructor */
    inline LowerSymmetrizeAccessor( Lhs& lhs): Base(), lhs_(lhs)
    { STK_STATIC_ASSERT_TWO_DIMENSIONS_ONLY(Lhs);}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs().rows();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs().cols();}

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j index of the row and of the column
     **/
    inline ConstReturnType elt2Impl(int i, int j) const
    { return ((j>i) ? lhs_.elt(j, i) : lhs_.elt(i, j));}

    /** @return a reference on the (i,j) element of the expression.
     *  @param i, j indexes of the element to get
     **/
    inline Type& elt2Impl(int i, int j)
    { return ((j>i) ? lhs_.elt(j, i) : lhs_.elt(i, j));}

  protected:
    Lhs& lhs_;
};

} // namespace STK


#undef EGAL

#endif /* STK_RESHAPEOPERATORS_H */
