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
 * created on: 21 nov. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_SlicingOperators.h
 *  @brief In this file we implement the RowOperator, ColOperator, SubVectorOperator and SubOperator classes.
 **/

#ifndef STK_SLICINGOPERATORS_H
#define STK_SLICINGOPERATORS_H

namespace STK
{

// forward declaration
template< typename Array> class RowOperator;
template< typename Array> class ColOperator;
// only for vectors/points/diagonal expressions
template< typename Array, int Size_> class SubVectorOperator;
// only for array expressions
template< typename Array, int SizeRows_, int SizeCols_> class SubOperator;

template< typename Array> class RowAccessor;
template< typename Array> class ColAccessor;
// only for vectors/points/diagonal expressions
template< typename Array, int Size_> class SubVectorAccessor;
// only for array expressions
template< typename Array, int SizeRows_, int SizeCols_> class SubAccessor;

namespace hidden
{

/** @ingroup hidden
 *  @brief Traits class for the row operator
 */
template<typename Lhs>
struct Traits< RowOperator <Lhs> >
{
  enum
  {
    structure_ = ( Lhs::structure_ != int(Arrays::vector_) && Lhs::structure_ != int(Arrays::number_) )
                 ? Arrays::point_ : Arrays::number_,
    orient_    = Lhs::orient_,
    sizeRows_  = 1,
    sizeCols_  = Lhs::sizeCols_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator<RowOperator < Lhs> > Row;
  typedef ColOperator<RowOperator < Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ConstReturnType ConstReturnType;
};

/** @ingroup hidden
 *  @brief Traits class for the row accessor
 */
template<typename Lhs>
struct Traits< RowAccessor <Lhs> >
{
  enum
  {
    structure_ = ( Lhs::structure_ != int(Arrays::vector_) && Lhs::structure_ != int(Arrays::number_) )
                 ? Arrays::point_ : Arrays::number_,
    orient_    = Lhs::orient_,
    sizeRows_  = 1,
    sizeCols_  = Lhs::sizeCols_,
    storage_   = Lhs::storage_
  };
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ConstReturnType ConstReturnType;
};

} // namespace hidden
/** @ingroup Arrays
  * @class RowOperator
  *
  * \brief Generic expression when the row of an expression is accessed
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * row operator.
  *
  * This class represents an expression where a row operator is applied to
  * an expression. It is the return type of the row operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name RowOperator type explicitly.
  */
template< typename Lhs>
class RowOperator: public ExprBase< RowOperator< Lhs> >, public TRef<1>
{
  public:
    typedef ExprBase< RowOperator< Lhs> > Base;
    typedef typename hidden::Traits< RowOperator< Lhs> >::Type Type;
    typedef typename hidden::Traits< RowOperator< Lhs> >::ConstReturnType ConstReturnType;
    enum
    {
        structure_ = hidden::Traits< RowOperator<Lhs> >::structure_,
        orient_    = hidden::Traits< RowOperator<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< RowOperator<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< RowOperator<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< RowOperator<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange; // will not be used

    /** constructor */
    inline RowOperator( Lhs const& lhs, int i): Base(), lhs_(lhs), i_(i), rows_(i_, 1), cols_(lhs_.rangeColsInRow(i_)) {}
//    /** Copy constructor */
//    inline RowOperator( RowOperator const& row, bool ref = true)
//                      : Base(), lhs_(row.lhs_), i_(row.i_), rows_(row.rows_), cols_(row.cols_) {}
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return rows_;}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return cols_;}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_;}
    /** @return the element (i,j)
     *  @param i, j index of the row and column
     **/
    inline ConstReturnType elt2Impl(int i, int j) const
    {
#ifdef STK_BOUNDS_CHECK
      if (i != i_) { STKRUNTIME_ERROR_2ARG(RowOperatorBase::elt2Impl,i,j,row index is not valid);}
#endif
      return (lhs().elt(i_, j));
    }
    /** @return the element jth element
     *  @param j index of the jth element
     **/
    inline ConstReturnType elt1Impl(int j) const { return (lhs().elt(i_, j));}
    /** accesses to the element */
    inline ConstReturnType elt0Impl() const { return (lhs().elt());}

  protected:
    Lhs const& lhs_;
    int i_;
    const RowRange  rows_;
    const ColRange  cols_;
};

/** @ingroup Arrays
  * @class RowAccessor
  *
  * \brief Generic expression when the row of an expression is accessed
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * row accessor.
  *
  * This class represents an expression where a row operator is applied to
  * an expression. It is the return type of the row operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name RowOperator type explicitly.
  */
template< typename Lhs>
class RowAccessor: public ExprBase< RowAccessor< Lhs> >, public TRef<1>
{
  public:
    typedef ExprBase< RowAccessor< Lhs> > Base;
    typedef typename hidden::Traits< RowAccessor< Lhs> >::Type Type;
    typedef typename hidden::Traits< RowAccessor< Lhs> >::ConstReturnType ConstReturnType;
    enum
    {
        structure_ = hidden::Traits< RowAccessor<Lhs> >::structure_,
        orient_    = hidden::Traits< RowAccessor<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< RowAccessor<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< RowAccessor<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< RowAccessor<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange; // will not be used

    /** constructor */
    inline RowAccessor( Lhs& lhs, int i): Base(), lhs_(lhs), i_(i), rows_(i_, 1), cols_(lhs_.rangeColsInRow(i_)) {}
//    /** Copy constructor */
//    inline RowAccessor( RowAccessor& row, bool ref = true)
//                      : Base(), lhs_(row.lhs_), i_(row.i_), rows_(row.rows_), cols_(row.cols_) {}
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return rows_;}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return cols_;}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_;}
    /** @return the element (i,j)
     *  @param i, j index of the row and column
     **/
    inline ConstReturnType elt2Impl(int i, int j) const
    {
#ifdef STK_BOUNDS_CHECK
      if (i != i_) { STKRUNTIME_ERROR_2ARG(RowOperatorBase::elt2Impl,i,j,row index is not valid);}
#endif
      return (lhs().elt(i_, j));
    }
    /** @return the element jth element
     *  @param j index of the jth element
     **/
    inline ConstReturnType elt1Impl(int j) const { return (lhs().elt(i_, j));}
    /** accesses to the element */
    inline ConstReturnType elt0Impl() const { return (lhs().elt());}

    /** @return the element (i,j)
     *  @param i, j index of the row and column
     **/
    inline Type& elt2Impl(int i, int j)
    {
#ifdef STK_BOUNDS_CHECK
      if (i != i_) { STKRUNTIME_ERROR_2ARG(RowOperatorBase::elt2Impl,i,j,row index is not valid);}
#endif
      return (lhs().elt(i_, j));
    }
    /** @return the element jth element
     *  @param j index of the jth element
     **/
    inline Type& elt1Impl(int j)
    {  return (lhs().elt(i_, j));}
    /** accesses to the element */
    inline Type& elt0Impl() { return (lhs().elt());}

  protected:
    Lhs& lhs_;
    int i_;
    const RowRange  rows_;
    const ColRange  cols_;
};


namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for the column operator
 */
template<typename Lhs>
struct Traits< ColOperator <Lhs> >
{
  enum
  {
    structure_ = ( Lhs::structure_ != int(Arrays::point_) && Lhs::structure_ != int(Arrays::number_) )
               ? Arrays::vector_ : Arrays::number_,
    orient_    = Lhs::orient_,
    sizeRows_  = Lhs::sizeRows_,
    sizeCols_  = 1,
    storage_   = Lhs::storage_
  };
  typedef RowOperator<ColOperator < Lhs> > Row;
  typedef ColOperator<ColOperator < Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ConstReturnType ConstReturnType;
};

/** @ingroup hidden
 *  @brief Traits class for the column operator
 */
template<typename Lhs>
struct Traits< ColAccessor <Lhs> >
{
  enum
  {
    structure_ = ( Lhs::structure_ != int(Arrays::point_) && Lhs::structure_ != int(Arrays::number_) )
               ? Arrays::vector_ : Arrays::number_,
    orient_    = Lhs::orient_,
    sizeRows_  = Lhs::sizeRows_,
    sizeCols_  = 1,
    storage_   = Lhs::storage_
  };
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ConstReturnType ConstReturnType;
};

} // namespace hidden

/** @ingroup Arrays
  * @class ColOperator
  *
  * @brief Generic expression when the column of an expression is accessed
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * column operator.
  *
  * This class represents an expression where a column accessor is applied to
  * an expression. It is the return type of the col operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name ColOperator type explicitly.
  */
template< typename Lhs>
class ColOperator: public ExprBase< ColOperator< Lhs> >, public TRef<1>
{
  public:
    typedef ExprBase< ColOperator< Lhs> > Base;
    typedef typename hidden::Traits< ColOperator<Lhs> >::Type Type;
    typedef typename hidden::Traits< ColOperator<Lhs> >::ConstReturnType ConstReturnType;
    enum
    {
        structure_ = hidden::Traits< ColOperator<Lhs> >::structure_,
        orient_    = hidden::Traits< ColOperator<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< ColOperator<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< ColOperator<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< ColOperator<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** constructor */
    inline ColOperator( Lhs const& lhs, int j)
                      : Base(), lhs_(lhs), j_(j), rows_(lhs_.rangeRowsInCol(j_)), cols_(j_,1) {}
//    /** Copy constructor */
//    inline ColOperator( ColOperator const& col, bool ref = true)
//                      : Base(), lhs_(col.lhs_), j_(col.j_), rows_(col.rows_), cols_(col.cols_) {}
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return rows_;}
    /** @return the columns range */
    inline ColRange const& colsImpl() const { return cols_;}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the ith element of the column
     *  @param i,j indexes of the element
     **/
    inline ConstReturnType elt2Impl(int i, int j) const
    {
#ifdef STK_BOUNDS_CHECK
      if (j != j_) { STKRUNTIME_ERROR_2ARG(ColOperatorBase::elt2Impl,i,j,column index is not valid);}
#endif
      return (lhs_.elt(i, j_));
    }
    /** @return the element ith element
     *  @param i index of the element to get
     **/
    inline ConstReturnType elt1Impl(int i) const { return (lhs_.elt(i, j_));}
    /** access to the element */
    inline ConstReturnType elt0Impl() const { return (lhs_.elt(j_));}

  protected:
    Lhs const& lhs_;
    const int j_;
    const RowRange rows_;
    const ColRange cols_;
};

/** @ingroup Arrays
  * @class ColAccessor
  *
  * @brief Generic expression when the column of an expression is accessed
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * column accessor.
  *
  * This class represents an expression where a column accessor is applied to
  * an expression. It is the return type of the column operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name ColOperator type explicitly.
  */
template< typename Lhs>
class ColAccessor: public ExprBase< ColAccessor< Lhs> >, public TRef<1>
{
  public:
    typedef ExprBase< ColAccessor< Lhs> > Base;
    typedef typename hidden::Traits< ColAccessor<Lhs> >::Type Type;
    typedef typename hidden::Traits< ColAccessor<Lhs> >::ConstReturnType ConstReturnType;
    enum
    {
        structure_ = hidden::Traits< ColAccessor<Lhs> >::structure_,
        orient_    = hidden::Traits< ColAccessor<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< ColAccessor<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< ColAccessor<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< ColAccessor<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** constructor */
    inline ColAccessor( Lhs& lhs, int j)
                      : Base(), lhs_(lhs), j_(j), rows_(lhs_.rangeRowsInCol(j_)), cols_(j_,1) {}
//    /** Copy constructor */
//    inline ColAccessor( ColAccessor& col, bool ref = true)
//                      : Base(), lhs_(col.lhs_), j_(col.j_), rows_(col.rows_), cols_(col.cols_) {}
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return rows_;}
    /** @return the columns range */
    inline ColRange const& colsImpl() const { return cols_;}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the ith element of the column
     *  @param i,j indexes of the element
     **/
    inline ConstReturnType elt2Impl(int i, int j) const
    {
#ifdef STK_BOUNDS_CHECK
      if (j != j_) { STKRUNTIME_ERROR_2ARG(ColOperatorBase::elt2Impl,i,j,column index is not valid);}
#endif
      return (lhs_.elt(i, j_));
    }
    /** @return the element ith element
     *  @param i index of the element to get
     **/
    inline ConstReturnType elt1Impl(int i) const { return (lhs_.elt(i, j_));}
    /** access to the element */
    inline ConstReturnType elt0Impl() const { return (lhs_.elt(j_));}

    /** @return the ith element of the column
     *  @param i,j indexes of the element
     **/
    inline Type& elt2Impl(int i, int j)
    {
#ifdef STK_BOUNDS_CHECK
      if (j != j_) { STKRUNTIME_ERROR_2ARG(ColOperatorBase::elt2Impl,i,j,column index is not valid);}
#endif
      return (lhs_.elt(i, j_));
    }
    /** @return the element ith element
     *  @param i index of the element to get
     **/
    inline Type& elt1Impl(int i) { return (lhs_.elt(i, j_));}
    /** access to the element */
    inline Type& elt0Impl() { return (lhs_.elt(j_));}

  protected:
    Lhs& lhs_;
    const int j_;
    const RowRange rows_;
    const ColRange cols_;
};


namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for the sub vector operator
 */
template<typename Lhs, int Size_>
struct Traits< SubVectorOperator <Lhs, Size_> >
{
  enum
  {
    structure_ = Lhs::structure_,
    orient_    = Lhs::orient_,
    sizeRows_  = (structure_ == int(Arrays::point_)) ? 1 : Size_,
    sizeCols_  = (structure_ == int(Arrays::vector_)) ? 1 : Size_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator<SubVectorOperator <Lhs, Size_> > Row;
  typedef ColOperator<SubVectorOperator <Lhs, Size_> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ConstReturnType ConstReturnType;
};

/** @ingroup hidden
 *  @brief Traits class for the sub vector accessor
 */
template<typename Lhs, int Size_>
struct Traits< SubVectorAccessor <Lhs, Size_> >
{
  enum
  {
    structure_ = Lhs::structure_,
    orient_    = Lhs::orient_,
    sizeRows_  = (structure_ == int(Arrays::point_)) ? 1 : Size_,
    sizeCols_  = (structure_ == int(Arrays::vector_)) ? 1 : Size_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator<SubVectorOperator <Lhs, Size_> > Row;
  typedef ColOperator<SubVectorOperator <Lhs, Size_> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ConstReturnType ConstReturnType;
};

} // namespace hidden

// forward declaration
template< typename Lhs, int Size_, int Structure_>
class SubVectorOperatorBase;
template< typename Lhs, int Size_, int Structure_>
class SubVectorAccessorBase;

/** @ingroup Arrays
  * @class SubVectorOperator
  *
  * @brief Generic expression when the sub-part of an expression is accessed
  * (specialization for vectors)
  *
  * @tparam Lhs the type of the array or expression to which we are
  * applying the sub operator.
  *
  * This class represents an expression where a subVectorOperator is applied to
  * an array expression. It is the return type of the sub(Range I) operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name SubVectorOperator type explicitly.
  */
template< typename Lhs, int Size_>
class SubVectorOperator: public SubVectorOperatorBase< Lhs, Size_
                                                     , hidden::Traits< SubVectorOperator<Lhs, Size_> >::structure_
                                                     >
                       , public TRef<1>
{
  public:
    typedef typename hidden::Traits< SubVectorOperator< Lhs, Size_> >::Type Type;
    typedef typename hidden::Traits< SubVectorOperator< Lhs, Size_> >::ConstReturnType ConstReturnType;
    enum
    {
      structure_ = hidden::Traits< SubVectorOperator<Lhs, Size_> >::structure_,
      orient_    = hidden::Traits< SubVectorOperator<Lhs, Size_> >::orient_,
      sizeRows_  = hidden::Traits< SubVectorOperator<Lhs, Size_> >::sizeRows_,
      sizeCols_  = hidden::Traits< SubVectorOperator<Lhs, Size_> >::sizeCols_,
      storage_   = hidden::Traits< SubVectorOperator<Lhs, Size_> >::storage_
    };
    typedef SubVectorOperatorBase<Lhs, Size_, structure_> Base;

    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** constructor */
    inline SubVectorOperator( Lhs const& lhs, TRange<Size_> const& I): Base(I), lhs_(lhs) {}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_;}
    /** @return element (i,j)
     *  @param i,j row and column indexes
     **/
    inline ConstReturnType elt2Impl(int i, int j) const { return (lhs_.elt2Impl(i, j));}
    /** @return i-th element
     *  @param i element index
     **/
    inline ConstReturnType elt1Impl(int i) const { return (lhs_.elt1Impl(i));}
    /** accesses to the element */
    inline ConstReturnType elt0Impl() const { return (lhs_.elt());}

  protected:
    /** @return element (i,j)
     *  @param i,j row and column indexes
     **/
    inline Type& elt2Impl(int i, int j) { return (lhs_.elt2Impl(i, j));}
    /** @return i-th element
     *  @param i element index
     **/
    inline Type& elt1Impl(int i) { return (lhs_.elt1Impl(i));}
    /** accesses to the element */
    inline Type& elt0Impl() { return (lhs_.elt());}

  protected:
    Lhs const& lhs_;
};


/** @ingroup Arrays
 *  @class SubVectorOperatorBase
 *  Specialization for point_
 **/
template< typename Lhs, int Size_>
class SubVectorOperatorBase<Lhs, Size_, Arrays::point_ >: public ExprBase< SubVectorOperator<Lhs, Size_> >
{
  public:
    typedef SubVectorOperator<Lhs, Size_> Derived;
    typedef ExprBase< Derived > Base;
    enum
    {
      sizeRows_  = hidden::Traits< Derived >::sizeRows_,
      sizeCols_  = hidden::Traits< Derived >::sizeCols_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;
    /** constructor */
    inline SubVectorOperatorBase( ColRange const& J): Base(), cols_(J) {}
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return this->asDerived().lhs().rows();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return cols_;}

  protected:
    ColRange cols_;
};

/** @ingroup Arrays
 *  @class SubVectorOperatorBase
 *  Specialization for vector_
 **/
template< typename Lhs, int Size_>
class SubVectorOperatorBase<Lhs, Size_, Arrays::vector_ >: public ExprBase< SubVectorOperator<Lhs, Size_> >
{
  public:
    typedef SubVectorOperator<Lhs, Size_> Derived;
    typedef ExprBase< Derived > Base;
    enum
    {
      sizeRows_  = hidden::Traits< Derived >::sizeRows_,
      sizeCols_  = hidden::Traits< Derived >::sizeCols_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;
    /** constructor */
    inline SubVectorOperatorBase( RowRange const& I): Base(), rows_(I) {}
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return rows_;}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return this->asDerived().lhs().cols();}

  protected:
    RowRange rows_;
};

/** @ingroup Arrays
 *  @class SubVectorOperatorBase
 *  Specialization for diagonal_
 **/
template< typename Lhs, int Size_>
class SubVectorOperatorBase<Lhs, Size_, Arrays::diagonal_ >: public ExprBase< SubVectorOperator<Lhs, Size_> >
{
    public:
      typedef SubVectorOperator<Lhs, Size_> Derived;
      typedef ExprBase< Derived > Base;
      enum
      {
        sizeRows_  = hidden::Traits< Derived >::sizeRows_,
        sizeCols_  = hidden::Traits< Derived >::sizeCols_
      };
      /** Type of the Range for the rows */
      typedef TRange<sizeRows_> RowRange;
      /** Type of the Range for the columns */
      typedef TRange<sizeCols_> ColRange;
      /** constructor */
      inline SubVectorOperatorBase( RowRange const& I): Base(), range_(I) {}
      /**  @return the range of the rows */
      inline RowRange const& rowsImpl() const { return range_;}
      /** @return the range of the Columns */
      inline ColRange const& colsImpl() const { return range_;}

  protected:
    RowRange range_;
};

/** @ingroup Arrays
  * @class SubVectorOperator
  *
  * @brief Generic expression when the sub-part of an expression is accessed
  * (specialization for vectors)
  *
  * @tparam Lhs the type of the array or expression to which we are
  * applying the sub operator.
  *
  * This class represents an expression where a subVectorAccessor is applied to
  * an array expression. It is the return type of the @c sub(Range I) operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name SubVectorOperator type explicitly.
  */
template< typename Lhs, int Size_>
class SubVectorAccessor: public SubVectorAccessorBase< Lhs, Size_
                                                     , hidden::Traits< SubVectorAccessor<Lhs, Size_> >::structure_
                                                     >
                       , public TRef<1>
{
  public:
    typedef typename hidden::Traits< SubVectorAccessor< Lhs, Size_> >::Type Type;
    typedef typename hidden::Traits< SubVectorAccessor< Lhs, Size_> >::ConstReturnType ConstReturnType;
    enum
    {
      structure_ = hidden::Traits< SubVectorAccessor<Lhs, Size_> >::structure_,
      orient_    = hidden::Traits< SubVectorAccessor<Lhs, Size_> >::orient_,
      sizeRows_  = hidden::Traits< SubVectorAccessor<Lhs, Size_> >::sizeRows_,
      sizeCols_  = hidden::Traits< SubVectorAccessor<Lhs, Size_> >::sizeCols_,
      storage_   = hidden::Traits< SubVectorAccessor<Lhs, Size_> >::storage_
    };
    typedef SubVectorOperatorBase<Lhs, Size_, structure_> Base;

    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** constructor */
    inline SubVectorAccessor( Lhs& lhs, TRange<Size_> const& I): Base(I), lhs_(lhs) {}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_;}
    /** @return element (i,j)
     *  @param i,j row and column indexes
     **/
    inline ConstReturnType elt2Impl(int i, int j) const { return (lhs_.elt2Impl(i, j));}
    /** @return i-th element
     *  @param i element index
     **/
    inline ConstReturnType elt1Impl(int i) const { return (lhs_.elt1Impl(i));}
    /** accesses to the element */
    inline ConstReturnType elt0Impl() const { return (lhs_.elt());}

    /** @return element (i,j)
     *  @param i,j row and column indexes
     **/
    inline Type& elt2Impl(int i, int j) { return (lhs_.elt2Impl(i, j));}
    /** @return i-th element
     *  @param i element index
     **/
    inline Type& elt1Impl(int i) { return (lhs_.elt1Impl(i));}
    /** accesses to the element */
    inline Type& elt0Impl() { return (lhs_.elt());}

  protected:
    Lhs& lhs_;
};


/** @ingroup Arrays
 *  @class SubVectorOperatorBase
 *  Specialization for point_
 **/
template< typename Lhs, int Size_>
class SubVectorAccessorBase<Lhs, Size_, Arrays::point_>: public ExprBase< SubVectorAccessor<Lhs, Size_> >
{
  public:
    typedef SubVectorAccessor<Lhs, Size_> Derived;
    typedef ExprBase< Derived > Base;
    enum
    {
      sizeRows_  = hidden::Traits< Derived >::sizeRows_,
      sizeCols_  = hidden::Traits< Derived >::sizeCols_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;
    /** constructor */
    inline SubVectorAccessorBase( ColRange const& J): Base(), cols_(J) {}
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return this->asDerived().lhs().rows();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return cols_;}

  protected:
    ColRange cols_;
};

/** @ingroup Arrays
 *  @class SubVectorOperatorBase
 *  Specialization for vector_
 **/
template< typename Lhs, int Size_>
class SubVectorAccessorBase<Lhs, Size_, Arrays::vector_ >: public ExprBase< SubVectorAccessor<Lhs, Size_> >
{
  public:
    typedef SubVectorAccessor<Lhs, Size_> Derived;
    typedef ExprBase< Derived > Base;
    enum
    {
      sizeRows_  = hidden::Traits< Derived >::sizeRows_,
      sizeCols_  = hidden::Traits< Derived >::sizeCols_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;
    /** constructor */
    inline SubVectorAccessorBase( RowRange const& I): Base(), rows_(I) {}
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return rows_;}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return this->asDerived().lhs().cols();}

  protected:
    RowRange rows_;
};

/** @ingroup Arrays
 *  @class SubVectorOperatorBase
 *  Specialization for diagonal_
 **/
template< typename Lhs, int Size_>
class SubVectorAccessorBase<Lhs, Size_, Arrays::diagonal_ >: public ExprBase< SubVectorAccessor<Lhs, Size_> >
{
    public:
      typedef SubVectorAccessor<Lhs, Size_> Derived;
      typedef ExprBase< Derived > Base;
      enum
      {
        sizeRows_  = hidden::Traits< Derived >::sizeRows_,
        sizeCols_  = hidden::Traits< Derived >::sizeCols_
      };
      /** Type of the Range for the rows */
      typedef TRange<sizeRows_> RowRange;
      /** Type of the Range for the columns */
      typedef TRange<sizeCols_> ColRange;
      /** constructor */
      inline SubVectorAccessorBase( RowRange const& I): Base(), range_(I) {}
      /**  @return the range of the rows */
      inline RowRange const& rowsImpl() const { return range_;}
      /** @return the range of the Columns */
      inline ColRange const& colsImpl() const { return range_;}

  protected:
    RowRange range_;
};



namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for the sub operator
 */
template<typename Lhs, int SizeRows_, int SizeCols_>
struct Traits< SubOperator <Lhs, SizeRows_, SizeCols_> >
{
  enum
  {
    structure_ = Lhs::structure_,
    orient_    = Lhs::orient_,
    sizeRows_  = SizeRows_,
    sizeCols_  = SizeCols_,
    storage_   = Lhs::storage_
  };
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ConstReturnType ConstReturnType;
};
/** @ingroup hidden
 *  @brief Traits class for the sub accessor
 */
template<typename Lhs, int SizeRows_, int SizeCols_>
struct Traits< SubAccessor <Lhs, SizeRows_, SizeCols_> >
{
  enum
  {
    structure_ = Lhs::structure_,
    orient_    = Lhs::orient_,
    sizeRows_  = SizeRows_,
    sizeCols_  = SizeCols_,
    storage_   = Lhs::storage_
  };
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ConstReturnType ConstReturnType;
};

} // namespace hidden

/** @ingroup Arrays
  * @class SubOperator
  *
  * @brief Generic expression when the sub-part of an expression is accessed
  *
  * @tparam Lhs the type of the array or expression to which we are
  * applying the sub operator.
  *
  * This class represents an expression where a subOperator is applied to
  * an array expression. It is the return type of the @c sub(I,J) operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name SubOperator type explicitly.
  */
template< typename Lhs, int SizeRows_, int SizeCols_>
class SubOperator: public ExprBase< SubOperator<Lhs, SizeRows_, SizeCols_> >
                 , public TRef<1>
{
  public:
    typedef typename hidden::Traits< SubOperator<Lhs, SizeRows_, SizeCols_> >::Type Type;
    typedef typename hidden::Traits< SubOperator<Lhs, SizeRows_, SizeCols_> >::ConstReturnType ConstReturnType;
    enum
    {
      structure_ = hidden::Traits< SubOperator<Lhs, SizeRows_, SizeCols_> >::structure_,
      orient_    = hidden::Traits< SubOperator<Lhs, SizeRows_, SizeCols_> >::orient_,
      sizeRows_  = hidden::Traits< SubOperator<Lhs, SizeRows_, SizeCols_> >::sizeRows_,
      sizeCols_  = hidden::Traits< SubOperator<Lhs, SizeRows_, SizeCols_> >::sizeCols_,
      storage_   = hidden::Traits< SubOperator<Lhs, SizeRows_, SizeCols_> >::storage_
    };
    typedef ExprBase< SubOperator<Lhs, SizeRows_, SizeCols_> > Base;

    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** constructor */
    inline SubOperator( Lhs const& lhs, RowRange const& I, ColRange const& J)
                      : Base(), lhs_(lhs), rows_(I), cols_(J) {}
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return rows_;}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return cols_;}
   /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_;}

    /** @return element (i,j)
     *  @param i,j row and column indexes
     **/
    inline ConstReturnType elt2Impl(int i, int j) const { return (lhs_.elt2Impl(i, j));}
    /** @return i-th element
     *  @param i element index
     **/
    inline ConstReturnType elt1Impl(int i) const { return (lhs_.elt1Impl(i));}
    /** accesses to the element */
    inline ConstReturnType elt0Impl() const { return (lhs_.elt());}

  protected:
    Lhs const& lhs_;
    RowRange rows_;
    ColRange cols_;
};

/** @ingroup Arrays
* @class SubAccessor
*
* @brief Generic expression when the sub-part of an expression is accessed
*
* @tparam Lhs the type of the array or expression to which we are
* applying the sub operator.
*
* This class represents an expression where a subVectorOperator is applied to
* an array expression. It is the return type of the @c sub(I, J) operation.
*
* Most of the time, this is the only way that it is used, so you typically
* don't have to name SubOperator type explicitly.
*/
template< typename Lhs, int SizeRows_, int SizeCols_>
class SubAccessor: public ExprBase< SubAccessor<Lhs, SizeRows_, SizeCols_> >
               , public TRef<1>
{
  public:
    typedef typename hidden::Traits< SubAccessor<Lhs, SizeRows_, SizeCols_> >::Type Type;
    typedef typename hidden::Traits< SubAccessor<Lhs, SizeRows_, SizeCols_> >::ConstReturnType ConstReturnType;
    enum
    {
      structure_ = hidden::Traits< SubAccessor<Lhs, SizeRows_, SizeCols_> >::structure_,
      orient_    = hidden::Traits< SubAccessor<Lhs, SizeRows_, SizeCols_> >::orient_,
      sizeRows_  = hidden::Traits< SubAccessor<Lhs, SizeRows_, SizeCols_> >::sizeRows_,
      sizeCols_  = hidden::Traits< SubAccessor<Lhs, SizeRows_, SizeCols_> >::sizeCols_,
      storage_   = hidden::Traits< SubAccessor<Lhs, SizeRows_, SizeCols_> >::storage_
    };
    typedef ExprBase< SubAccessor<Lhs, SizeRows_, SizeCols_> > Base;

    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** constructor */
    inline SubAccessor( Lhs& lhs, RowRange const& I, ColRange const& J)
                      : Base(), lhs_(lhs), rows_(I), cols_(J) {}
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return rows_;}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return cols_;}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_;}

    /** @return element (i,j)
     *  @param i,j row and column indexes
     **/
    inline ConstReturnType elt2Impl(int i, int j) const { return (lhs_.elt2Impl(i, j));}
    /** @return i-th element
     *  @param i element index
     **/
    inline ConstReturnType elt1Impl(int i) const { return (lhs_.elt1Impl(i));}
    /** accesses to the element */
    inline ConstReturnType elt0Impl() const { return (lhs_.elt());}

    /** @return element (i,j)
     *  @param i,j row and column indexes
     **/
    inline Type& elt2Impl(int i, int j) { return (lhs_.elt2Impl(i, j));}
    /** @return i-th element
     *  @param i element index
     **/
    inline Type& elt1Impl(int i) { return (lhs_.elt1Impl(i));}
    /** accesses to the element */
    inline Type& elt0Impl() { return (lhs_.elt());}

  protected:
    Lhs& lhs_;
    RowRange rows_;
    ColRange cols_;
};


} // namespace STK

#endif /* STK_SLICINGOPERATORS_H */
