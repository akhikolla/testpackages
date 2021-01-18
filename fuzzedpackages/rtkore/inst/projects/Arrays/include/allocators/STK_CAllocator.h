/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016 Serge Iovleff

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
 * created on: 10 août 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_CAllocator.h
 *  @brief In this file we define the CAllocator class.
 **/

#ifndef STK_CALLOCATOR_H
#define STK_CALLOCATOR_H

#include "../STK_ITContainer2D.h"
#include "STK_MemAllocator.h"

namespace STK
{

/** @ingroup Arrays
 *  @brief Allocator for dense Array classes.
 *  The data are stored in two dimensions.
 *  It can be the columns or the rows allocator of any dense container.
 */
template<typename Type, int SizeRows_, int SizeCols_, bool Orient_>
class CAllocator;

namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of Traits class for CAllocator.
 */
template< typename Type_, int SizeRows_, int SizeCols_, bool Orient_>
struct Traits< CAllocator<Type_, SizeRows_, SizeCols_, Orient_> >
{
    typedef Type_  Type;
    typedef typename RemoveConst<Type_>::Type const& ConstReturnType;

    enum
    {
      structure_ = (SizeRows_==1) ? (SizeCols_== 1) ? Arrays::number_ : Arrays::point_
                                  : (SizeCols_== 1) ? Arrays::vector_
                                                    : (SizeRows_ == SizeCols_ && SizeRows_!= UnknownSize) ?
                                                      Arrays::square_ : Arrays::array2D_,
      orient_    = Orient_,
      sizeRows_  = SizeRows_,
      sizeCols_  = SizeCols_,
      sizeProd_  = ProductSizeRowsBySizeCols<SizeRows_, SizeCols_>::prod_,
      storage_   = Arrays::dense_ // always dense
    };

    typedef CAllocator<Type_, 1, SizeCols_, Orient_> Row;
    typedef CAllocator<Type_, SizeRows_, 1, Orient_> Col;

    /** Type of the base allocator allocating data */
    typedef MemAllocator<Type_, sizeProd_> Allocator;
};

} // namespace hidden

// forward declaration
template<class Derived> class CAllocatorBase;
template<class Derived, bool Orient_> class OrientedCAllocator;
template<class Derived, int SizeRows_, int SizeCols_> class StructuredCAllocator;

/** Base class for CAllocators */
template<class Derived>
class CAllocatorBase: public ITContainer2D<Derived>
{
  public:
    typedef typename hidden::Traits< Derived >::Row Row;
    typedef typename hidden::Traits< Derived >::Col Col;
    typedef typename hidden::Traits< Derived >::Type Type;
    typedef typename hidden::Traits< Derived >::ConstReturnType ConstReturnType;

    enum
    {
      structure_ = hidden::Traits< Derived >::structure_,
      orient_    = hidden::Traits< Derived >::orient_,
      sizeRows_  = hidden::Traits< Derived >::sizeRows_,
      sizeCols_  = hidden::Traits< Derived >::sizeCols_,
      storage_   = hidden::Traits< Derived >::storage_
    };
    typedef ITContainer2D< Derived > Base;
    typedef typename hidden::Traits< Derived >::Allocator Allocator;

    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

  protected:
    /** default constructor */
    CAllocatorBase(): Base() {}
    /** constructor */
    CAllocatorBase( RowRange const& I, ColRange const& J): Base(I, J) {}
    /** copy constructor */
    CAllocatorBase( CAllocatorBase const& A): Base(A){}

  public:
    /** @return the range of the effectively stored elements in the column. */
    RowRange const& rangeRowsInCol(int) const { return this->rows();}
    /** @return the range of the effectively stored elements in the row. */
    ColRange const& rangeColsInRow(int) const { return this->cols();}

    /** Access to the ith row of the Allocator.
     *  @param i index of the row
     *  @return a reference on the ith row
     **/
    inline Row row(int i) const
    {
#ifdef STK_BOUNDS_CHECK
      if (this->beginRows() > i) { STKOUT_OF_RANGE_1ARG(CAllocatorBase::row, i, beginRows() > i);}
      if (this->endRows() <= i)  { STKOUT_OF_RANGE_1ARG(CAllocatorBase::row, i, endRows() <= i);}
#endif
      return Row(this->asDerived(), TRange<1>(i,1), this->cols());
    }
    /** Access to the row (i,J) of the Allocator.
     *  @param i,J index of the row and range of the columns
     *  @return a reference on the ith row
     **/
    template<int Size_>
    inline CAllocator<Type, 1, Size_, orient_> row(int i, TRange<Size_> const& J) const
    {
#ifdef STK_BOUNDS_CHECK
      if (this->beginRows() > i) { STKOUT_OF_RANGE_2ARG(CAllocatorBase::row, i, J, beginRows() > i);}
      if (this->endRows() <= i)  { STKOUT_OF_RANGE_2ARG(CAllocatorBase::row, i, J, endRows() <= i);}
      if (this->beginCols() > J.begin()) { STKOUT_OF_RANGE_2ARG(CAllocatorBase::row, i, J, beginCols() > J.begin());}
      if (this->endCols() < J.end())  { STKOUT_OF_RANGE_2ARG(CAllocatorBase::row, i, J, endCols() < J.end());}
#endif
      return CAllocator<Type, 1, Size_, orient_>(this->asDerived(), TRange<1>(i,1), J);
    }

    /** Access to the jth column of the Allocator.
     *  @param j index of the column
     *  @return a reference on the jth column
     **/
    inline Col col(int j) const
    {
#ifdef STK_BOUNDS_CHECK
      if (this->beginCols() > j) { STKOUT_OF_RANGE_1ARG(CAllocatorBase::col, j, beginCols() > j);}
      if (this->endCols() <= j)  { STKOUT_OF_RANGE_1ARG(CAllocatorBase::col, j, endCols() <= j);}
#endif
      return Col(this->asDerived(), this->rows(), TRange<1>(j,1));
    }
    /** Access to the column (I,j) of the Allocator.
     *  @param I,j range of the rows and index of the column
     *  @return a reference on the jth column
     **/
    template<int Size_>
    inline CAllocator<Type, Size_, 1, orient_> col(TRange<Size_> const& I, int j) const
    {
#ifdef STK_BOUNDS_CHECK
      if (this->beginCols() > j)         { STKOUT_OF_RANGE_2ARG(CAllocatorBase::col, I, j, beginCols() > j);}
      if (this->endCols() <= j)          { STKOUT_OF_RANGE_2ARG(CAllocatorBase::col, I, j, endCols() <= j);}
      if (this->beginRows() > I.begin()) { STKOUT_OF_RANGE_2ARG(CAllocatorBase::col, I, j, beginRows() > I.begin());}
      if (this->endRows() < I.end())     { STKOUT_OF_RANGE_2ARG(CAllocatorBase::col, I, j, endRows() < I.end());}
#endif
      return CAllocator<Type, Size_, 1, orient_>(this->asDerived(), I, TRange<1>(j,1));
    }

    /** Access to the sub-part (I,J) of the Allocator.
      *  @param I,J range of the rows and columns
      *  @return a reference on a sub-part of the Allocator
      **/
     template<int OtherRows_, int OtherCols_>
     inline CAllocator<Type, OtherRows_, OtherCols_, orient_>
     sub(TRange<OtherRows_> const& I, TRange<OtherCols_> const& J) const
     {
#ifdef STK_BOUNDS_CHECK
      if (this->beginRows() > I.begin()) { STKOUT_OF_RANGE_2ARG(CAllocatorBase::sub, I, J, beginRows() > I.begin());}
      if (this->endRows() < I.end())     { STKOUT_OF_RANGE_2ARG(CAllocatorBase::sub, I, J, endRows() < I.end());}
      if (this->beginCols() > J.begin()) { STKOUT_OF_RANGE_2ARG(CAllocatorBase::sub, I, J, beginCols() > j);}
      if (this->endCols() < J.end())     { STKOUT_OF_RANGE_2ARG(CAllocatorBase::sub, I, J, endCols() <= j);}
#endif
       return CAllocator<Type, OtherRows_, OtherCols_, orient_>(this->asDerived(), I, J);
     }


    /** @param pos1, pos2 position of the first and second columns to swap */
    void swapCols(int pos1, int pos2)
    {
#ifdef STK_BOUNDS_CHECK
      if (this->beginCols() > pos1)
      { STKOUT_OF_RANGE_2ARG(CAllocatorBase::swapCols,pos1, pos2,beginCols() > pos1);}
      if (this->endCols() <= pos1)
      { STKOUT_OF_RANGE_2ARG(CAllocatorBase::swapCols,pos1, pos2,endCols() <= pos1);}
      if (this->beginCols() > pos2)
      { STKOUT_OF_RANGE_2ARG(CAllocatorBase::swapCols,pos1, pos2,beginCols() > pos2);}
      if (this->endCols() <= pos2)
      { STKOUT_OF_RANGE_2ARG(CAllocatorBase::swapCols,pos1, pos2,endCols()<= pos2);}
#endif
      for (int i=this->beginRows(); i< this->endRows(); ++i)
      { std::swap(this->asDerived().elt2Impl(i, pos1),this->asDerived().elt2Impl(i, pos2));}
    }
    /** @param pos1, pos2 position of the first and second rows to swap */
    void swapRows(int pos1, int pos2)
    {
#ifdef STK_BOUNDS_CHECK
      if (this->beginRows() > pos1)
      { STKOUT_OF_RANGE_2ARG(CAllocatorBase::swapRows,pos1, pos2,beginRows() > pos1);}
      if (this->endRows() <= pos1)
      { STKOUT_OF_RANGE_2ARG(CAllocatorBase::swapRows,pos1, pos2,endRows() <= pos1);}
      if (this->beginRows() > pos2)
      { STKOUT_OF_RANGE_2ARG(CAllocatorBase::swapRows,pos1, pos2,beginRows() > pos2);}
      if (this->endRows() <= pos2)
      { STKOUT_OF_RANGE_2ARG(CAllocatorBase::swapRows,pos1, pos2,endRows() <= pos2);}
#endif
      for (int j=this->beginCols(); j< this->endCols(); ++j)
      { std::swap(this->asDerived().elt2Impl(pos1, j),this->asDerived().elt2Impl(pos2, j));}
    }
};


/**  @ingroup Arrays
 *   @brief Specialization for column-oriented Allocators.
 **/
template<class Derived>
class OrientedCAllocator<Derived, Arrays::by_col_>: public CAllocatorBase<Derived>
{
  public:
    typedef typename hidden::Traits< Derived >::Row Row;
    typedef typename hidden::Traits< Derived >::Col Col;
    typedef typename hidden::Traits< Derived >::Type Type;
    typedef typename hidden::Traits< Derived >::ConstReturnType ConstReturnType;

    enum
    {
      structure_ = hidden::Traits<Derived>::structure_,
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };
    typedef CAllocatorBase<Derived> Base;
    typedef typename hidden::Traits<Derived>::Allocator Allocator;

    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

  protected:
    /** default constructor */
    OrientedCAllocator( RowRange const& I, ColRange const& J)
                      : Base(I, J), ldx_(I.size()), allocator_(prod(I, J))
    {}
    /** copy constructor */
    OrientedCAllocator( OrientedCAllocator const& A, bool ref)
                     : Base(A), ldx_(A.ldx()), allocator_(A.allocator_, ref)
    { if (!ref) allocator_.assign(A.allocator_);}
    /** Reference constructor */
    template<class OtherDerived>
    inline OrientedCAllocator( OrientedCAllocator<OtherDerived, Arrays::by_col_> const& A
                             , RowRange const& I, ColRange const& J)
                             : Base(I, J), ldx_(A.ldx()), allocator_(A.allocator(), true)  {}
    /** wrapper constructor for 0 based C-Array*/
    OrientedCAllocator( Type* const& q, int nbRow, int nbCol)
                      : Base(Range(0,nbRow), Range(0,nbCol)), ldx_(nbRow)
                      , allocator_(q, nbRow * nbCol, true)
    {}
    /** destructor */
    ~OrientedCAllocator() {}

  public:
    /** @return is this own its data ? */
    inline bool isRef() const { return allocator_.isRef();}
    /** @return a reference on the main pointer*/
    inline Type*& p_data() { return allocator_.p_data();}
    /** @return a constant reference on the main pointer*/
    inline Type* const& p_data() const { return allocator_.p_data();}
    /**  @return a constant reference on the memory manager */
    inline Allocator const& allocator() const { return allocator_;}
    /** @return the leading dimension of the allocator */
    inline int ldx() const { return ldx_;}

    /** @return a constant reference on the element (i,j) of the Allocator.
     *  @param i, j indexes of the element
     **/
    inline ConstReturnType elt2Impl(int i, int j) const { return p_data()[j*ldx_ + i];}
    /** @return a reference on the element (i,j) of the Allocator.
     *  @param i, j indexes of the element
     **/
    inline Type& elt2Impl(int i, int j) { return p_data()[j*ldx_ + i];}
    /** set a value to this allocator.
     *  @param v the value to set
     **/
    void setValue(Type const& v)
    {
      for (int j= this->beginCols(); j < this->endCols(); ++j)
        for (int i = this->beginRows(); i < this->endRows(); ++i)
        { this->elt(i, j) = v;}
    }
    /** move T to this without copying data*/
    Derived& move(OrientedCAllocator const& T)
    {
      if (this->asPtrDerived() != T.asPtrDerived())
      {
        allocator_.move(T.allocator_);
        this->asDerived().moveImpl(T.asDerived()); // move any other data fields
        this->setRanges(T.rows(), T.cols());
        ldx_ = T.ldx_;
      }
      return this->asDerived();
    }

  protected:
    /** leading dimension of the data set */
    int ldx_;
    /** manager of the memory */
    Allocator allocator_;
    /** do nothing by default */
    void moveImpl(OrientedCAllocator const& T) {}
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    void exchange(OrientedCAllocator &T)
    {
      allocator_.exchange(T.allocator_);
      Base::exchange(T);
      std::swap(ldx_, T.ldx_);
    }
    /** @brief Compute the range of the 1D Allocator when we want to
     *  allocate a 2D array with  range I for the rows and range J for the
     *  columns.
     *  @param I,J the range of the rows and columns
     *  @return The range of the 1D allocator
     **/
    static Range prod(Range const& I, Range const& J)
    { return Range(I.size()*J.begin()+I.begin(), I.size()*J.size()); }
    /** return the increment to apply to a zero based pointer corresponding to
     *  the actual first row and first column indexes. */
    int shiftInc(int beginRows, int beginCols)
    { return ldx_*beginCols+beginRows; }
    /** set the index corresponding to the actual size of the allocator. */
    void setSizedIdx() {ldx_ = this->sizeRows();}
};

/**  @ingroup Arrays
 *   @brief Specialization for row-oriented Allocators.*/
template<class Derived>
class OrientedCAllocator<Derived, Arrays::by_row_>: public CAllocatorBase<Derived>
{
  public:
    typedef typename hidden::Traits< Derived >::Row Row;
    typedef typename hidden::Traits< Derived >::Col Col;
    typedef typename hidden::Traits< Derived >::Type Type;
    typedef typename hidden::Traits< Derived >::ConstReturnType ConstReturnType;

    enum
    {
      structure_ = hidden::Traits<Derived>::structure_,
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };
    typedef CAllocatorBase<Derived> Base;
    typedef typename hidden::Traits<Derived>::Allocator Allocator;

    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

  protected:
    /** constructor with specified ranges */
    OrientedCAllocator( RowRange const& I, ColRange const& J)
                      : Base(I, J), ldx_(J.size()), allocator_(prod(I, J))
    {}
    /** copy constructor */
    OrientedCAllocator( OrientedCAllocator const& A, bool ref)
                      : Base(A), ldx_(A.ldx_), allocator_(A.allocator_, ref)
    { if (!ref) allocator_.assign(A.allocator_);}
    /** Reference constructor */
    template<class OtherDerived>
    inline OrientedCAllocator( OrientedCAllocator<OtherDerived, Arrays::by_row_> const& A
                             , RowRange const& I, ColRange const& J)
                             : Base(I, J), ldx_(A.ldx()), allocator_(A.allocator(), true)
    {}
    /** wrapper constructor for 0 based C-Array*/
    OrientedCAllocator( Type* const& q, int nbRow, int nbCol)
                            : Base(Range(0,nbRow), Range(0,nbCol)), ldx_(nbCol)
                             , allocator_(q, nbRow * nbCol, true)
    {}
    /** destructor */
    ~OrientedCAllocator() {}

  public:
    /** @return is this own its data ? */
    inline bool isRef() const { return allocator_.isRef();}
    /** @return a reference on the main pointer*/
    inline Type*& p_data() { return allocator_.p_data();}
    /** @return a constant reference on the main pointer*/
    inline Type* const& p_data() const { return allocator_.p_data();}
    /**  @return a constant reference on the memory manager */
    inline Allocator const& allocator() const { return allocator_;}
    /** @return the leading dimension of the allocator */
    inline int ldx() const { return ldx_;}
    /** @return a constant reference on the element (i,j) of the Allocator.
     *  @param i,j indexes of the element
     **/
    inline ConstReturnType elt2Impl(int i, int j) const  { return p_data()[i*ldx_ + j];}
    /** @return a reference on the element (i,j) of the Allocator.
     *  @param i,j indexes of the element
     **/
    inline Type& elt2Impl(int i, int j) { return p_data()[i*ldx_ + j];}
    /** set a value to this container.
     *  @param v the value to set
     **/
    void setValue(Type const& v)
    {
      for (int i = this->beginRows(); i < this->endRows(); ++i)
        for (int j= this->beginCols(); j < this->endCols(); ++j)
        { this->elt(i, j) = v;}
    }
    /** move T to this without copying data*/
    Derived& move(OrientedCAllocator const& T)
    {
      if (this->asPtrDerived() != T.asPtrDerived())
      {
        allocator_.move(T.allocator_);
        this->asDerived().moveImpl(T.asDerived()); // move any other data fields
        this->setRanges(T.rows(), T.cols());
        ldx_ = T.ldx_;
      }
      return this->asDerived();
    }

  protected:
    /** leading dimension of the data set */
    int ldx_;
    /** manager of the memory */
    Allocator allocator_;
    /** exchange T with this.
     *  @param T the container to move
     **/
    void exchange(OrientedCAllocator &T)
    {
      allocator_.exchange(T.allocator_);
      Base::exchange(T);
      std::swap(ldx_, T.ldx_);
    }
    /** @brief Compute the range of the 1D Allocator when we want to
     *  allocate a 2D array with I indexes in the first dimension and J indexes
     *  in the second dimension.
     *  @param I,J range of the rows and columns
     *  @return The range of the 1D allocator
     **/
    static Range prod(Range const& I, Range const& J)
    { return Range(J.size()*I.begin()+J.begin(), I.size()*J.size());}
    /** return the increment corresponding to the actual first row an column. */
    int shiftInc(int beginRows, int beginCols)
    { return ldx_*beginRows+beginCols; }
    /** set the index corresponding to the actual size of the allocator. */
    void setSizedIdx() { ldx_ = this->sizeCols();}
};

/** @ingroup Arrays
 *  @brief  Base class for the general by_col_ structured case.
 **/
template<class Derived, int SizeRows_, int SizeCols_>
class StructuredCAllocator: public OrientedCAllocator<Derived, hidden::Traits<Derived>::orient_>
{
  public:
    typedef typename hidden::Traits< Derived >::Row Row;
    typedef typename hidden::Traits< Derived >::Col Col;
    typedef typename hidden::Traits< Derived >::Type Type;
    typedef typename hidden::Traits< Derived >::ConstReturnType ConstReturnType;

    enum
    {
      structure_ = hidden::Traits< Derived >::structure_,
      orient_    = hidden::Traits< Derived >::orient_,
      sizeRows_  = hidden::Traits< Derived >::sizeRows_,
      sizeCols_  = hidden::Traits< Derived >::sizeCols_,
      storage_   = hidden::Traits< Derived >::storage_
    };
    typedef OrientedCAllocator<Derived, orient_> Base;
    typedef IContainer2D<sizeRows_, sizeCols_> LowBase;
    typedef typename hidden::Traits<Derived>::Allocator Allocator;

    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

  protected:
    /** Default constructor */
    StructuredCAllocator( RowRange const& I, ColRange const& J): Base(I, J) {}
    /** copy constructor */
    StructuredCAllocator( StructuredCAllocator const& A, bool ref): Base(A, ref) {}
    /** Reference constructor */
    template<class OtherDerived, int OtherRows_, int OtherCols_>
    inline StructuredCAllocator( StructuredCAllocator<OtherDerived, OtherRows_, OtherCols_> const& A
                               , RowRange const& I, ColRange const& J)
                               : Base(A, I, J)
    {}
    /** wrapper constructor for 0 based C-Array*/
    StructuredCAllocator( Type* const& q, int nbRow, int nbCol): Base(q, nbRow, nbCol) {}

  public:
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    void exchange(StructuredCAllocator &T) { Base::exchange(T);}
    /** shift the first indexes of the allocator
     *  @param beginRows, beginCols indexes to set
     **/
    void shift2Impl(int beginRows, int beginCols)
    {
      if ((beginRows == this->beginRows())&&(beginCols == this->beginCols())) return;
      LowBase::shift(beginRows, beginCols);
      this->allocator_.shift(this->shiftInc(beginRows, beginCols));
    }
    /** shift the first indexes of the allocator (for square and diagonal matrices).
     *  @param begin the index of the first row and column
     **/
    void shift1Impl(int begin) { shift2Impl(begin, begin);}
    /** resize the allocator (for square and diagonal matrices).
     *  @param size the size of the point
     **/
    Derived& resize1Impl(int size)
    { return this->asDerived().resize2Impl(size, size);}
};


/** @ingroup Arrays
 *  @brief specialization for the point_ case (derived class is a RowVector).
 **/
template<class Derived, int SizeCols_>
class StructuredCAllocator<Derived, 1, SizeCols_>
   : public OrientedCAllocator<Derived, hidden::Traits<Derived>::orient_>
{
  public:
    typedef typename hidden::Traits< Derived >::Row Row;
    typedef typename hidden::Traits< Derived >::Col Col;
    typedef typename hidden::Traits< Derived >::Type Type;
    typedef typename hidden::Traits< Derived >::ConstReturnType ConstReturnType;

    enum
    {
      structure_ = hidden::Traits< Derived >::structure_,
      orient_    = hidden::Traits< Derived >::orient_,
      sizeRows_  = hidden::Traits< Derived >::sizeRows_,
      sizeCols_  = hidden::Traits< Derived >::sizeCols_,
      storage_   = hidden::Traits< Derived >::storage_
    };
    typedef OrientedCAllocator<Derived, orient_> Base;
    typedef IContainer2D<sizeRows_, sizeCols_> LowBase;
    typedef typename hidden::Traits< Derived >::Allocator Allocator;

    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

  protected:
    /** Default constructor */
    StructuredCAllocator( RowRange const& I, ColRange const& J): Base(I, J), row_(I.begin()) {}
    /** copy constructor */
    StructuredCAllocator( StructuredCAllocator const& A, bool ref): Base(A, ref), row_(A.row_) {}
    /** Reference constructor */
    template<class OtherDerived, int OtherRows_, int OtherCols_>
    inline StructuredCAllocator( StructuredCAllocator<OtherDerived, OtherRows_, OtherCols_> const& A
                               , RowRange const& I, ColRange const& J)
                               : Base(A, I, J), row_(I.begin())
    {}
    /** wrapper constructor for 0 based C-Array*/
    StructuredCAllocator( Type* const& q, int , int nbCol): Base(q, 1, nbCol), row_(0) {}
  public:
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    void exchange(StructuredCAllocator &T)
    { Base::exchange(T);
      std::swap(row_, T.row_);
    }

    /** @return the index of the first element */
    inline int begin() const { return this->beginCols();}
    /**  @return the ending index of the elements */
    inline int end() const { return this->endCols();}
    /**  @return the size of the allocator */
    inline int size() const { return this->sizeCols();}

    /** @return a constant reference on the element (i,j) of the Allocator.
     *  @param j index of the column
     **/
    inline ConstReturnType elt1Impl( int j) const { return this->elt2Impl(row_,j);}
    /** @return a reference on the element (i,j) of the Allocator.
     *  @param j index of the columns
     **/
    inline Type& elt1Impl( int j) { return this->elt2Impl(row_,j);}
    /** shift the first indexes of the allocator
     *  @param beginRows, beginCols indexes to set
     **/
    void shift2Impl(int beginRows, int beginCols)
    {
      if ((beginRows == this->beginRows())&&(beginCols == this->beginCols())) return;
      LowBase::shift(beginRows, beginCols);
      this->allocator_.shift(this->shiftInc(beginRows, beginCols));
      row_ = beginRows;
    }
    /** shift the first column index of the allocator.
     *  @param beginCols the index of the first column */
    void shift1Impl(int beginCols)
    { shift2Impl(beginCols, beginCols);
      row_ = beginCols;
    }
    /** resize the allocator.
     *  @param sizeCols the size of the point
     **/
    Derived& resize1Impl(int sizeCols)
    { this->asDerived().resize2Impl(1, sizeCols);
      row_ = this->beginRows();
      return this->asDerived();
    }
     /** @return a sub-vector in the specified range of the Allocator.
     *  @param J range of the sub-vector
     **/
    template< int Size_>
    CAllocator<Type, 1, Size_, orient_> subVector( TRange<Size_> const& J) const
    { return Base::template row<Size_>(row_, J);}
    /** move T to this.
     *  @param T the container to move
     **/
    void moveImpl(StructuredCAllocator const& T)
    { row_ = T.row_;}

  private:
    /** row of the point (needed when this is a reference) */
    int row_;
};

/** @ingroup Arrays
 *  @brief specialization for the vector_ case.
 **/
template<class Derived, int SizeRows_>
class StructuredCAllocator<Derived, SizeRows_, 1>
   : public OrientedCAllocator<Derived, hidden::Traits<Derived>::orient_>
{
  public:
    typedef typename hidden::Traits< Derived >::Row Row;
    typedef typename hidden::Traits< Derived >::Col Col;
    typedef typename hidden::Traits< Derived >::Type Type;
    typedef typename hidden::Traits< Derived >::ConstReturnType ConstReturnType;

    enum
    {
      structure_ = hidden::Traits< Derived >::structure_,
      orient_    = hidden::Traits< Derived >::orient_,
      sizeRows_  = hidden::Traits< Derived >::sizeRows_,
      sizeCols_  = hidden::Traits< Derived >::sizeCols_,
      storage_   = hidden::Traits< Derived >::storage_
    };
    typedef OrientedCAllocator<Derived, orient_> Base;
    typedef IContainer2D<sizeRows_, sizeCols_> LowBase;
    typedef typename hidden::Traits<Derived>::Allocator Allocator;

    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

  protected:
    /** Default constructor */

    StructuredCAllocator( RowRange const& I, ColRange const& J): Base(I, J), col_(J.begin()) {}
    /** copy constructor */
    StructuredCAllocator( StructuredCAllocator const& A, bool ref): Base(A, ref), col_(A.col_) {}
    /** Reference constructor */
    template<class OtherDerived, int OtherRows_, int OtherCols_>
    inline StructuredCAllocator( StructuredCAllocator<OtherDerived, OtherRows_, OtherCols_> const& A
                               , RowRange const& I, ColRange const& J)
                               : Base(A, I, J)
                               , col_(J.begin())
    {}
    /** wrapper constructor for 0 based C-Array*/
    StructuredCAllocator( Type* const& q, int nbRow, int): Base(q, nbRow, 1), col_(0) {}
  public:
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    void exchange(StructuredCAllocator &T)
    { Base::exchange(T);
      std::swap(col_, T.col_);
    }
    /** @return the index of the first element */
    inline int begin() const { return this->beginRows();}
    /**  @return the ending index of the elements */
    inline int end() const { return this->endRows();}
    /**  @return the size of the allocator */
    inline int size() const { return this->sizeRows();}

    /** @return a constant reference on the element (i,j) of the Allocator.
     *  @param i index of the row
     **/
    inline ConstReturnType elt1Impl( int i) const { return this->elt2Impl(i, col_);}
    /** @return a reference on the element (i,j) of the Allocator.
     *  @param i index of the row
     **/
    inline Type& elt1Impl( int i) { return this->elt2Impl(i, col_);}
    /** shift the first indexes of the allocator
     *  @param beginRows, beginCols indexes to set
     **/
    void shift2Impl(int beginRows, int beginCols)
    {
      if ((beginRows == this->beginRows())&&(beginCols == this->beginCols())) return;
      LowBase::shift(beginRows, beginCols);
      this->allocator_.shift(this->shiftInc(beginRows, beginCols));
      col_ = beginCols;
    }
    /** shift the first row index of the allocator.
     *  @param beginRows the index of the first row
     **/
    void shift1Impl(int beginRows)
    { shift2Impl(beginRows, beginRows);
      col_ = beginRows;
    }
    /** resize the allocator.
     *  @param sizeRow the size of the vector
     **/
    Derived& resize1Impl(int sizeRow)
    {
      this->asDerived().resize2Impl(sizeRow, 1);
      col_ = this->beginCols();
      return this->asDerived();
    }
    /** @return a sub-vector in the specified range of the Allocator.
     *  @param I range of the sub-vector
     **/
     template<int Size_>
     inline CAllocator<Type, Size_, 1, orient_> subVector(TRange<Size_> const& I) const
     { return Base::template col<Size_>(I, col_);}
    /** move T to this.
     *  @param T the container to move
     **/
    void moveImpl(StructuredCAllocator const& T)
    { col_ = T.col_;}

  private:
    int col_;
};


/** @ingroup Arrays
 *  @brief specialization for the number_ case.
 **/
template<class Derived>
class StructuredCAllocator<Derived, 1, 1>
   : public OrientedCAllocator<Derived, hidden::Traits<Derived>::orient_>
{
  public:
    typedef typename hidden::Traits< Derived >::Row Row;
    typedef typename hidden::Traits< Derived >::Col Col;
    typedef typename hidden::Traits< Derived >::Type Type;
    typedef typename hidden::Traits< Derived >::ConstReturnType ConstReturnType;

    enum
    {
      structure_ = hidden::Traits< Derived >::structure_,
      orient_    = hidden::Traits< Derived >::orient_,
      sizeRows_  = hidden::Traits< Derived >::sizeRows_,
      sizeCols_  = hidden::Traits< Derived >::sizeCols_,
      storage_   = hidden::Traits< Derived >::storage_
    };
    typedef OrientedCAllocator<Derived, orient_> Base;
    typedef IContainer2D<sizeRows_, sizeCols_> LowBase;
    typedef typename hidden::Traits<Derived>::Allocator Allocator;

    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

  protected:
    /** Default constructor */
    StructuredCAllocator( RowRange const& I, ColRange const& J)
                        : Base(I, J), row_(I.begin()), col_(J.begin())
    {}
    /** copy constructor */
    StructuredCAllocator( StructuredCAllocator const& A, bool ref)
                        : Base(A, ref), row_(A.row_), col_(A.col_)
    {}
    /** Reference constructor */
    template<class OtherDerived, int OtherRows_, int OtherCols_>
    inline StructuredCAllocator( StructuredCAllocator<OtherDerived, OtherRows_, OtherCols_> const& A
                               , RowRange const& I, ColRange const& J)
                               : Base(A, I, J)
                               , row_(I.begin()), col_(J.begin())
    {}
    /** wrapper constructor for 0 based C-Array*/
    StructuredCAllocator( Type* const& q, int , int): Base(q, 1, 1), row_(0), col_(0) {}
    /** destructor */
    ~StructuredCAllocator() {}
    /** move T to this.
     *  @param T the container to move
     **/
    StructuredCAllocator& move(StructuredCAllocator const& T)
    { row_ = T.row_; col_ = T.col_; return *this;}

  public:
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    void exchange(StructuredCAllocator &T)
    { Base::exchange(T);
      std::swap(row_, T.row_);
      std::swap(col_, T.col_);
    }
    /**  @return the size of the allocator */
    inline int size() const { return 1;}

    /** @return a constant reference on the element of the Allocator. */
    inline ConstReturnType elt0Impl() const { return this->elt2Impl(row_, col_);}
    /** @return a reference on the element of the Allocator. */
    inline Type& elt0Impl() { return this->elt2Impl(row_, col_);}
    /** @return a constant reference on the element of the Allocator. */
    inline ConstReturnType elt1Impl(int) const { return this->elt2Impl(row_, col_);}
    /** @return a reference on the element of the Allocator. */
    inline Type& elt1Impl(int) { return this->elt2Impl(row_, col_);}

    /** shift the first indexes of the allocator
     *  @param beginRows, beginCols indexes to set
     **/
    void shift2Impl(int beginRows, int beginCols)
    {
      if ((beginRows == this->beginRows())&&(beginCols == this->beginCols())) return;
      LowBase::shift(beginRows, beginCols);
      this->allocator_.shift(this->shiftInc(beginRows, beginCols));
      row_ = beginRows;
      col_ = beginCols;
   }
    /** shift the first indexes of the allocator.
     *  @param begin the index of the first row and column
     **/
    void shift1Impl(int begin)
    {
      this->shift2Impl(begin, begin);
      row_ = begin; col_ = begin;
    }
    /** resize the allocator (for square and diagonal matrices). */
    Derived& resize1Impl(int) { return this->asDerived().resize2Impl(1, 1);}
    /** move T to this.
     *  @param T the container to move
     **/
    void moveImpl(StructuredCAllocator const& T) { row_ = T.row_; col_ = T.col_;}

  private:
    int row_;
    int col_;
};

/** @ingroup Arrays
 *  @brief Allocator for the dense CArray classes.
 *  The size of the Allocator is known in both dimension
 */
template<typename Type_, int SizeRows_, int SizeCols_, bool Orient_>
class CAllocator: public StructuredCAllocator<CAllocator<Type_, SizeRows_, SizeCols_, Orient_>, SizeRows_, SizeCols_>
{
  public:
    typedef typename hidden::Traits< CAllocator >::Row Row;
    typedef typename hidden::Traits< CAllocator >::Col Col;
    typedef typename hidden::Traits< CAllocator >::Type Type;
    typedef typename hidden::Traits< CAllocator >::ConstReturnType ConstReturnType;

    enum
    {
      structure_ = hidden::Traits< CAllocator >::structure_,
      orient_    = hidden::Traits< CAllocator >::orient_,
      sizeRows_  = hidden::Traits< CAllocator >::sizeRows_,
      sizeCols_  = hidden::Traits< CAllocator >::sizeCols_,
      storage_   = hidden::Traits< CAllocator >::storage_
    };
    typedef MemAllocator<Type, SizeRows_* SizeCols_> Allocator;
    typedef StructuredCAllocator<CAllocator, SizeRows_, SizeCols_> Base;
    typedef IContainer2D<SizeRows_, SizeCols_> LowBase;

    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    using Base::allocator_;
    /** default constructor */
    CAllocator(): Base(SizeRows_, SizeCols_) {}
    CAllocator( int, int): Base(SizeRows_, SizeCols_) {}
    CAllocator( int, int, Type const& v): Base(SizeRows_, SizeCols_) {Base::Base::setValue(v);}
    CAllocator( CAllocator const& A, bool ref = true): Base(A, ref)
    { if (!ref) { allocator_.assign(A.allocator_);} }

    template< int OtherRows_, int OtherCols_>
    inline CAllocator( CAllocator<Type, OtherRows_, OtherCols_, Orient_> const& A
                     , RowRange const& I, ColRange const& J)
                     : Base(A, I, J)
    {}
    /** wrapper constructor for 0 based C-Array*/
    CAllocator( Type* const& q, int , int ): Base(q, SizeRows_, SizeCols_) {}
    /** destructor */
    ~CAllocator() {}
    /**  clear allocated memories */
    void clear() {}

    /** resize this */
    CAllocator& resize2Impl( int, int) { return *this;}
    void realloc(int, int) {}
};

/** @brief Specialized Allocator for the dense Arrays classes.
 *  The sizes of the columns and of the rows are unknown. The Orientation is
 *  either by rows or by column.
 */
template<typename Type_, bool Orient_>
class CAllocator<Type_, UnknownSize, UnknownSize, Orient_>
    : public StructuredCAllocator<CAllocator<Type_, UnknownSize, UnknownSize, Orient_>, UnknownSize, UnknownSize>
{
  public:
    typedef typename hidden::Traits< CAllocator >::Row Row;
    typedef typename hidden::Traits< CAllocator >::Col Col;
    typedef typename hidden::Traits< CAllocator >::Type Type;
    typedef typename hidden::Traits< CAllocator >::ConstReturnType ConstReturnType;

    enum
    {
      structure_ = hidden::Traits<CAllocator>::structure_,
      orient_    = hidden::Traits<CAllocator>::orient_,
      sizeRows_  = hidden::Traits<CAllocator>::sizeRows_,
      sizeCols_  = hidden::Traits<CAllocator>::sizeCols_,
      storage_   = hidden::Traits<CAllocator>::storage_
    };
    typedef MemAllocator<Type, UnknownSize> Allocator;
    typedef StructuredCAllocator<CAllocator, UnknownSize, UnknownSize> Base;

    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    using Base::allocator_;

    /** Default constructor */
    CAllocator(): Base(0, 0)  {}
    /** Constructor with specified size.
     *  @param sizeRows, sizeCols size of the rows and columns
     **/
    CAllocator( int sizeRows, int sizeCols): Base(sizeRows, sizeCols) {}
    /** Constructor with specified size and specified value.
     *  @param sizeRows, sizeCols size of the rows and columns
     *  @param v the initial value
     **/
    CAllocator( int sizeRows, int sizeCols, Type const& v): Base(sizeRows, sizeCols)
    { Base::Base::setValue(v);}
    /** Copy or wrapper constructor.
     *  @param A : the array to copy
     *  @param ref : is this a wrapper of A ?
     **/
    CAllocator( CAllocator const& A, bool ref = true): Base(A, ref)
    { if (!ref) { allocator_.assign(A.allocator_);}}
    /** Wrapper constructor. This become a reference on (some part of) the Allocator A.
     *  @param A original allocator
     *  @param I,J range of the rows and columns to wrap.
     **/
    template< int OtherRows_, int OtherCols_>
    inline CAllocator( CAllocator<Type, OtherRows_, OtherCols_, Orient_> const& A
                     , RowRange const& I, ColRange const& J)
                     : Base(A, I, J)
    {}
    /** wrapper constructor for 0 based C-Array*/
    CAllocator( Type* const& q, int nbRow, int nbCol): Base(q, nbRow, nbCol) {}
    /** Destructor */
    ~CAllocator() {}
    /**  clear allocated memories */
    void clear() { allocator_.free(); this->setRanges();}
    /** resize CAllocator */
    CAllocator& resize2Impl( int sizeRows, int sizeCols)
    {
      if (this->sizeRows() == sizeRows && this->sizeCols() == sizeCols) return *this;
     // check size
     if ((sizeRows <= 0)||(sizeCols<=0))
     {
       // free any allocated memory if this is not a reference
       allocator_.free();
       this->setRanges(sizeRows, sizeCols);
       this->setSizedIdx();
       return *this;
     }
     // allocate
     allocator_.malloc(this->prod(sizeRows, sizeCols));
     this->setRanges(sizeRows, sizeCols);
     this->setSizedIdx();
     return *this;
    }
    /** @brief function for memory reallocation.
     *  The function assume you want to increase or reduce the size without
     *  modification of the bases ranges.
     *  @param sizeRows, sizeCols size of the rows and columns
    **/
    void realloc(int sizeRows, int sizeCols)
    {
      if ((sizeRows == this->sizeRows())&&(sizeCols == this->sizeCols())) return;
      // create a copy the original data set
      CAllocator copy;
      this->exchange(copy);
      try
      {
        // create new container
        resize2Impl(sizeRows, sizeCols);
        this->shift2Impl(copy.beginRows(), copy.beginCols());
        // copy data
        const int endRow = std::min(copy.endRows(), this->endRows());
        const int endCol = std::min(copy.endCols(), this->endCols());
        for (int j= this->beginCols(); j < endCol; ++j)
          for (int i = this->beginRows(); i< endRow; ++i)
          { this->elt(i, j) = copy.elt(i, j);}
      }
      catch (std::bad_alloc & error)  // if an alloc error occur
      {
        this->exchange(copy); // restore the original container
        STKRUNTIME_ERROR_2ARG(Allocator::realloc, sizeRows, sizeCols, memory allocation failed);
      }
    }
};

/** @brief Specialized Allocator for the dense Arrays classes.
 *  The number of rows is known, the number of columns unknown
 **/
template<typename Type_, int SizeRows_, bool Orient_>
class CAllocator<Type_, SizeRows_, UnknownSize, Orient_>
     : public StructuredCAllocator<CAllocator<Type_, SizeRows_, UnknownSize, Orient_>, SizeRows_, UnknownSize>
{
  public:
    typedef typename hidden::Traits< CAllocator >::Row Row;
    typedef typename hidden::Traits< CAllocator >::Col Col;
    typedef typename hidden::Traits< CAllocator >::Type Type;
    typedef typename hidden::Traits< CAllocator >::ConstReturnType ConstReturnType;

    enum
    {
      structure_  = hidden::Traits< CAllocator >::structure_,
      orient_     = hidden::Traits< CAllocator >::orient_,
      sizeRows_   = hidden::Traits< CAllocator >::sizeRows_,
      sizeCols_   = hidden::Traits< CAllocator >::sizeCols_,
      storage_    = hidden::Traits< CAllocator >::storage_
    };
    typedef StructuredCAllocator<CAllocator, SizeRows_, UnknownSize> Base;
    typedef MemAllocator<Type, UnknownSize> Allocator;

    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    using Base::allocator_;


    CAllocator(): Base(sizeRows_, 0) {}
    CAllocator( int, int sizeCols): Base(sizeRows_, sizeCols) {}
    CAllocator( int, int sizeCols, Type const& v): Base(sizeRows_, sizeCols) { Base::Base::setValue(v);}
    CAllocator( CAllocator const& A, bool ref = true): Base(A, ref) { if (!ref) { allocator_.assign(A.allocator_);}}

    template< int OtherRows_, int OtherCols_>
    inline CAllocator( CAllocator<Type, OtherRows_, OtherCols_, Orient_> const& A
                     , RowRange const& I, ColRange const& J)
                     : Base(A, I, J)
    {}
    /** wrapper constructor for 0 based C-Array*/
    CAllocator( Type* const& q, int , int nbCol): Base(q, sizeRows_, nbCol) {}
    /** destructor */
    ~CAllocator() {}
    /** clear allocated memories. */
    void clear() { allocator_.free(); this->setCols();}
    /** resize allocator */
    CAllocator& resize2Impl( int, int sizeCols)
    {
      if (this->sizeCols() == sizeCols) return *this;
      // check size
      if (sizeCols<=0)
      {
        // free any allocated memory if it is not a reference
        allocator_.free();
        this->setRanges(SizeRows_, sizeCols);
        this->setSizedIdx();
        return *this;
      }
      // allocate
      allocator_.malloc(this->prod(SizeRows_, sizeCols));
      this->setRanges(SizeRows_, sizeCols);
      this->setSizedIdx();
      return *this;
    }
    void realloc(int, int sizeCols)
    {
      CAllocator copy;// create a copy the original data set
      this->exchange(copy);
      try
      {
        // create new container
        resize2Impl(SizeRows_, sizeCols);
        this->shift2Impl(copy.beginRows(), copy.beginCols());
        // copy data
        const int endCol = std::min(copy.endCols(), this->endCols());
        for (int j= this->beginCols(); j < endCol; ++j)
          for (int i = this->beginRows(); i < this->endRows(); ++i)
        { this->elt(i, j) = copy.elt(i, j);}

      }
      catch (std::bad_alloc const& error)  // if an alloc error occur
      {
        this->exchange(copy); // restore the original container
        STKRUNTIME_ERROR_2ARG(Allocator::realloc, sizeRows_, sizeCols, memory allocation failed);
      }
    }
};

/** @brief Specialized Allocator for the dense Arrays classes.
 *  The sizes of the columns is known, the number of rows is unknown
 */
template<typename Type_, bool Orient_, int SizeCols_>
class CAllocator<Type_, UnknownSize, SizeCols_, Orient_>
     : public StructuredCAllocator<CAllocator<Type_, UnknownSize, SizeCols_, Orient_>, UnknownSize, SizeCols_>
{
  public:
    typedef typename hidden::Traits< CAllocator >::Row Row;
    typedef typename hidden::Traits< CAllocator >::Col Col;
    typedef typename hidden::Traits< CAllocator >::Type Type;
    typedef typename hidden::Traits< CAllocator >::ConstReturnType ConstReturnType;

    enum
    {
      structure_  = hidden::Traits< CAllocator >::structure_,
      orient_     = hidden::Traits< CAllocator >::orient_,
      sizeRows_   = hidden::Traits< CAllocator >::sizeRows_,
      sizeCols_   = hidden::Traits< CAllocator >::sizeCols_,
      storage_    = hidden::Traits< CAllocator >::storage_
    };

    typedef MemAllocator<Type, UnknownSize> Allocator;
    typedef StructuredCAllocator<CAllocator, UnknownSize, SizeCols_> Base;
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    using Base::allocator_;

    CAllocator(): Base(0, SizeCols_) {}
    CAllocator( int sizeRows, int): Base(sizeRows, SizeCols_) {}

    CAllocator( int sizeRows, int, Type const& v): Base(sizeRows, SizeCols_){ this->asDerived().setValue(v);}
    CAllocator( CAllocator const& A, bool ref = true): Base(A, ref)
    { if (!ref) { allocator_.assign(A.allocator_);} }

    /** wrap other allocator */
    template< int OtherRows_, int OtherCols_>
    inline CAllocator( CAllocator<Type, OtherRows_, OtherCols_, Orient_> const& A
                     , RowRange const& I, ColRange const& J)
                     : Base(A, I, J)
    {}
    /** wrapper constructor for 0 based C-Array*/
    CAllocator( Type* const& q, int nbRow, int ): Base(q, nbRow, SizeCols_) {}
    /** destructor */
    ~CAllocator() {}
    /**  clear allocated memories. */
    void clear() { allocator_.free(); this->setRows();}
    /** resize the rows */
    CAllocator& resize2Impl( int sizeRows, int)
    {
      if (this->sizeRows() == sizeRows) return *this;
      if (sizeRows <= 0)// check if
      {
        // free any allocated memory if it is not a reference
        allocator_.free();
        this->setRanges(sizeRows, SizeCols_);
        this->setSizedIdx();
        return *this;
     }
     // allocate
     allocator_.malloc(this->prod(sizeRows,  SizeCols_));
     this->setRanges(sizeRows,  SizeCols_);
     this->setSizedIdx();
     return *this;
   }
   void realloc(int sizeRows, int)
   {
     // create a copy the original data set
     CAllocator copy;
     this->exchange(copy);
     try
     {
       // create new container
       resize2Impl(sizeRows, SizeCols_);
       this->shift2Impl(copy.beginRows(), copy.beginCols());
       // copy data
       const int endRow = std::min(copy.endRows(), this->endRows());
       for (int j= this->beginCols(); j < this->endCols(); ++j)
         for (int i = this->beginRows(); i< endRow; ++i)
         { this->elt(i, j) = copy.elt(i, j);}
     }
     catch (std::bad_alloc & error)  // if an alloc error occur
     {
       this->exchange(copy); // restore the original container
       STKRUNTIME_ERROR_2ARG(Allocator::realloc, sizeRows, SizeCols_, memory allocation failed);
     }
   }
};

} // namespace STK

#endif /* STK_CALLOCATOR_H */
