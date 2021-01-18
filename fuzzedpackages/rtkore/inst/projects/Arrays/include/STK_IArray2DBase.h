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
 * Purpose:  Define the Interface for the Array classes.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_IArray2DBase.h
 *  @brief Interface base class for the Array2D classes, this is an internal
 *  header file, included by other containers library headers.
 *
 *  You should not attempt to use it directly but rather used one of the
 *  derived class like Array2D, except if you want to create your own
 *  Container Class.
 **/

#ifndef STK_IARRAY2DBASE_H
#define STK_IARRAY2DBASE_H

#include "STK_IArrayBase.h"

#include "STK_ArrayBaseApplier.h"
#include "STK_ArrayBaseAssign.h"
#include "STK_ArrayBaseInitializer.h"

#include "STK_Array1D.h"

namespace STK
{
namespace hidden
{
// forward declaration
//template< typename Type> class Array2D;
//template< typename Type> class Array2DSquare;
//template< typename Type> class Array2DDiagonal;
//template< typename Type> class Array2DNumber;
//template< typename Type> class Array2DPoint;
//template< typename Type> class Array2DVector;
//template< typename Type> class Array2DLowerTriangular;
//template< typename Type> class Array2DUpperTriangular;

///** @ingroup hidden
// *  @brief The traits struct Slice2D must be specialized allow to disambiguate return
// *  type of the col/row/sub operators for Array2D family
// *
// *  @sa STK::IArray2D
// */
//template<typename Derived, int SizeRows, int SizeCols>
//struct Slice2D
//{
//    typedef typename Traits< Derived >::Type Type;
//    enum
//    {
//      structure_ = Traits<Derived>::structure_,
//      orient_    = Traits<Derived>::orient_,
//      sizeRows_  = Traits<Derived>::sizeRows_,
//      sizeCols_  = Traits<Derived>::sizeCols_,
//      storage_   = Traits<Derived>::storage_,
//
//      isNumber_ = (SizeRows==SizeCols)&&(SizeRows==1),
//      isVector_ = (SizeCols == 1)&&(!isNumber_),
//      isPoint_  = (SizeRows == 1)&&(!isNumber_),
//      isSquare_ = (SizeRows==SizeCols)&&(SizeRows!=UnknownSize)&&(!isNumber_),
//      isArray_  = (SizeRows!=SizeCols)&&(!isNumber_)
//    };
//   typedef typename If< (isNumber_), Array2DNumber<Type>
//                      , typename If< isVector_, Array2DVector<Type>
//                          , typename If< isPoint_, Array2DPoint<Type>
//                              , typename If< isSquare_, Array2DSquare<Type>
//                                           , Array2D<Type>
//                                           >::Result
//                                       >::Result
//                                   >::Result
//                      >::Result Result;
//};
//
///** @ingroup hidden
// *  helper allowing to disambiguate SubVector access
// **/
//template<typename Derived, int Size>
//struct Slice2DDispatcher
//{
//    enum
//    { structure_ = Traits<Derived>::structure_};
//    typedef typename If< structure_ == (int)Arrays::vector_
//                       , typename Slice2D<Derived, Size, 1>::Result
//                       , typename Slice2D<Derived, 1, Size>::Result
//                       >::Result Result;
//};

}

/** @ingroup Arrays
 *  @brief template interface base class for two-dimensional arrays.
 *
 * A IArray2DBase is an interface class for two-dimensional Arrays
 * stored in columns and having flexible dimensions. It is possible
 * to add, remove easily columns and rows in Derived class.
 *
 * Each column has a Range stored in the array @c rangeCols_ and a
 * capacity stored in the array @c availableRows_. It should be worth
 * noting that we should have
 * @code
 *   (rangeCols_[j].size() <= availableRows_[j]) == true;
 *   (rangeCols_[j].isIn(this->rows()) == true;
 * @endcode
 *
 * @tparam PTRCOL is the type of the ptr of column in a two-dimensional
 * array: for exemple @c TYPE*, @c Array1D<TYPE>*, @c DBACCESS*....
 * @tparam Derived is the name of the class that implements @c IArray2DBase.
 **/
template < class PTRCOL, class Derived>
class IArray2DBase: protected IContainer2D<hidden::Traits<Derived>::sizeRows_, hidden::Traits<Derived>::sizeCols_>
                  , public IArrayBase<Derived>
{
   // needed by merge
   template <class OTHERPTRCOL, class OtherDerived>
   friend    class IArray2DBase;

  public:
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef typename hidden::Traits<Derived>::ConstReturnType ConstReturnType;

    typedef typename hidden::Traits<Derived>::Row Row;
    typedef typename hidden::Traits<Derived>::Col Col;
    typedef typename hidden::Traits<Derived>::SubRow SubRow;
    typedef typename hidden::Traits<Derived>::SubCol SubCol;
    typedef typename hidden::Traits<Derived>::SubArray SubArray;
    // for 1D container
    typedef typename hidden::Traits<Derived>::SubVector SubVector;

    enum
    {
      structure_ = hidden::Traits<Derived>::structure_,
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;
    /** Type for the IContainer2D base Class. */
    typedef IContainer2D<sizeRows_, sizeCols_ > Base2D;
    /** Type for the Base Class. */
    typedef MemAllocator<PTRCOL, sizeCols_> Allocator;
    /** type of the Base Container Class. */
    typedef IArrayBase<Derived> Base;

  protected:
    /** Default constructor */
    IArray2DBase(): Base2D(), Base(), allocator_()
                  , availableRows_(), rangeCols_()
                  , availableCols_(0), capacityByCols_(0)
    { mallocCols(this->cols());}
    /** constructor with specified ranges
     *  @param I,J rows and columns range
     **/
    IArray2DBase( Range const& I, Range const& J)
                : Base2D(I, J), Base(), allocator_()
                , availableRows_(), rangeCols_()
                , availableCols_(0), capacityByCols_(0)
    { mallocCols(this->cols());}
    /** Copy constructor If we want to wrap T, the main ptr will be wrapped
     *  in MemAllocator class. If we want to copy  T, Allocator is
     *  initialized to default values.
     *  @note bug correction, we have to use a copy of T.rangeCols_. in case
     *  we are using the code
     *  @code
     *  Array2DVector<TYPE> Dref(D.sub(J), true)
     *  @endcode
     *  we would get an error.
     *  @param T the container to copy
     *  @param ref true if we wrap T
     **/
    IArray2DBase( IArray2DBase const& T, bool ref =false)
               : Base2D(T), Base()
                , allocator_(T.allocator_, ref)
                , availableRows_(T.availableRows_, ref)
                , rangeCols_(T.rangeCols_) //  have to be created again, in case T is a temporary
                , availableCols_(T.availableCols_), capacityByCols_(T.capacityByCols_)
    { if (!ref) mallocCols(this->cols());}
    /** constructor by reference, ref_=1.
     *  @param T the container to wrap
     *  @param I,J rows and columns to wrap
     **/
    template<class OtherDerived>
    IArray2DBase( IArray2DBase<PTRCOL, OtherDerived> const& T, Range const& I, Range const& J)
               : Base2D(I, J), Base()
                , allocator_(T.allocator(), J)         // a reference
                , availableRows_(T.availableRows(), J) // a reference
                , rangeCols_(T.rangeCols())            // have to be created again
                , availableCols_(J.size()), capacityByCols_(I.size())
    {
      for (int j=J.begin(); j<J.end(); j++)
      { rangeCols_[j] = inf(I, T.rangeCols()[j]);}
    }
    /** Wrapper constructor. We get a reference of the data.
     *  @param q pointer on data
     *  @param I,J range of the rows and columns to wrap
     **/
    IArray2DBase( PTRCOL* q, Range const& I, Range const& J)
                : Base2D(I, J), Base(), allocator_(q, J, true)
                , availableRows_(J, I.size()), rangeCols_(J, I)
                , availableCols_(I.size()), capacityByCols_(J.size())
    {}
    /** destructor. Allocated horizontal memory (the array with the pointers
     *  on the columns) is liberated by the Allocator.
     **/
    ~IArray2DBase() {}

  public:
    /** @return the Vertical range */
    inline RowRange const& rowsImpl() const { return Base2D::rows();}
    /** @return the Vertical range */
    inline RowRange const& rows() const { return Base2D::rows();}
    /** @return the index of the first row */
    inline int beginRows() const { return Base2D::beginRows();}
    /** @return the ending index of the rows */
    inline int endRows() const { return Base2D::endRows();}
    /** @return the number of rows */
    inline int sizeRows() const { return Base2D::sizeRows();}

    /**@return the Horizontal range */
    inline ColRange const& colsImpl() const { return Base2D::cols();}
    /**@return the Horizontal range */
    inline ColRange const& cols() const { return Base2D::cols();}
    /** @return the index of the first column */
    inline int beginCols() const { return Base2D::beginCols();}
    /**  @return the ending index of columns */
    inline int endCols() const { return Base2D::endCols();}
    /** @return the number of columns */
    inline int sizeCols() const { return Base2D::sizeCols();}

    /**  @return the index of the last column */
    inline int lastIdxCols() const { return Base2D::lastIdxCols();}
    /** @return the index of the last row */
    inline int lastIdxRows() const { return Base2D::lastIdxRows();}

    /**  @return @c true if the container is empty, @c false otherwise */
    bool empty() const { return Base2D::empty();}

    /** access to an element.
     *  @note Take care that @c PTRCOL can be accessed using @c operator[]
     *  @param i,j indexes of the element to get
     *  @return a reference on the (i,j) element
     **/
    inline Type& elt2Impl( int i, int j) { return this->data(j)[i];}
    /** constant access to an element.
     *  @note Take care that @c PTRCOL can be accessed using @c operator[]
     *  @param i,j indexes of the element to get
     *  @return a constant reference on the (i,j) element
     **/
    inline ConstReturnType elt2Impl( int i, int j) const { return this->data(j)[i];}
    /** access to a part of a column.
     *  @param j index of the column
     *  @return a reference in the range I of the column j of this
     **/
    inline Col col( int j) const
    { return Col( this->asDerived(), this->rangeRowsInCol(j), j);}
    /** access to a part of a column.
     *  @param I range of the rows
     *  @param j index of the col
     *  @return a reference in the range I of the column j of this
     **/
    inline SubCol col(Range const& I, int j) const
    { return SubCol( this->asDerived(), inf(I, this->rangeRowsInCol(j)), j);}
    /** access to many columns.
     *  @param J range of the index of the cols
     *  @return a 2D array containing the Container in the Horizontal range @c J
     **/
    inline SubArray col(Range const& J) const
    { return SubArray( this->asDerived(), this->rows(), J);}
    /** access to a part of a row.
     *  @param i index of the row
     *  @return a reference of the row i.
     **/
    inline Row row( int i) const
    { return Row( this->asDerived(), this->rangeColsInRow(i), i);}
    /** access to a part of a row.
     *  @param i index of the row
     *  @param J range of the columns
     *  @return a reference of the row i.
     **/
    inline SubRow row(int i, Range const& J) const
    { return SubRow( this->asDerived(), inf(J, this->rangeColsInRow(i)), i);}
    /** access to many rows.
     *  @param I range of the index of the rows
     *  @return a 2D array containing the Container in the vertical range @c I
     **/
    inline SubArray row(Range const& I) const
    { return SubArray(this->asDerived(), I, this->cols());}
    /** @return  many elements.
     *  @param J Range of the elements
     **/
    inline SubVector sub(Range const& J) const
    { return SubVector(this->asDerived(), J);}
    /** access to a sub-array.
     *  @param I,J range of the rows and of the columns
     **/
    inline SubArray sub(Range const& I, Range const& J) const
    { return SubArray(this->asDerived(), I, J);}

    // overloaded operators
    /** @return a constant reference on the element (i,j) of the 2D container.
     *  @param i,j row and column indexes
     **/
    inline ConstReturnType operator()(int i, int j) const
    {
#ifdef STK_BOUNDS_CHECK
       if (this->beginRows() > i) { STKOUT_OF_RANGE_2ARG(IArrayBase::elt, i, j, beginRows() > i);}
       if (this->endRows() <= i)  { STKOUT_OF_RANGE_2ARG(IArrayBase::elt, i, j, endRows() <= i);}
       if (this->beginCols() > j) { STKOUT_OF_RANGE_2ARG(IArrayBase::elt, i, j, beginCols() > j);}
       if (this->endCols() <= j)  { STKOUT_OF_RANGE_2ARG(IArrayBase::elt, i, j, endCols() <= j);}
#endif
      return this->elt(i,j);}
    /** @return a reference on the element (i,j) of the 2D container.
     *  @param i, j indexes of the element to get
     **/
    inline Type& operator()(int i, int j)
    {
#ifdef STK_BOUNDS_CHECK
       if (this->beginRows() > i) { STKOUT_OF_RANGE_2ARG(IArray2DBase::elt, i, j, beginRows() > i);}
       if (this->endRows() <= i)  { STKOUT_OF_RANGE_2ARG(IArray2DBase::elt, i, j, endRows() <= i);}
       if (this->beginCols() > j) { STKOUT_OF_RANGE_2ARG(IArray2DBase::elt, i, j, beginCols() > j);}
       if (this->endCols() <= j)  { STKOUT_OF_RANGE_2ARG(IArray2DBase::elt, i, j, endCols() <= j);}
#endif
      return this->elt(i,j);
    }
    /** @return the ith element
     *  @param i index of the element to get
     **/
    inline ConstReturnType operator[](int i) const
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);
#ifdef STK_BOUNDS_CHECK
      if (this->asDerived().begin() > i) { STKOUT_OF_RANGE_1ARG(IArra2DyBase::elt, i, begin() > i);}
      if (this->asDerived().end() <= i)  { STKOUT_OF_RANGE_1ARG(IArray2DBase::elt, i, end() <= i);}
#endif
      return this->elt(i);
    }
    /** @return a reference on the ith element
     *  @param i index of the element to get
     **/
    inline Type& operator[](int i)
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);
#ifdef STK_BOUNDS_CHECK
      if (this->asDerived().begin() > i) { STKOUT_OF_RANGE_1ARG(IArray2DBase::elt, i, begin() > i);}
      if (this->asDerived().end() <= i)  { STKOUT_OF_RANGE_1ARG(IArray2DBase::elt, i, end() <= i);}
#endif
      return this->elt(i);
    }
    /** @return a constant reference on the number */
    inline ConstReturnType operator()() const { return this->elt();}
    /** @return the number */
    inline Type& operator()() { return this->elt();}
    // overloaded operators for sub-arrays/vectors
    /** @return the sub-vector in given range
     *  @param I range to get
     **/
    inline SubVector operator[](Range const& I) const
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);
      return this->sub(I);
    }
    /** @param I range of the index of the rows
     *  @param j index of the column
     *  @return a Vertical container containing the column @c j of this
     *  in the range @c I
     **/
    inline SubCol operator()(Range const& I, int j) const { return this->col(I, j);}
    /** @param i index of the row
     *  @param J range of the columns
     *  @return an Horizontal container containing the row @c i of this
     *  in the range @c J
     **/
    inline SubRow operator()(int i, Range const& J) const { return this->row(i, J);}
    /** @param I,J range of the rows and of the columns
     *  @return a 2D container containing this in the range @c I, @c J
     **/
    inline SubArray operator()(Range const& I, Range const& J) const { return this->sub(I, J);}

    /** @return a constant pointer on the j-th column of the container
     *  @param j the index of the column
     **/
    inline PTRCOL const& data(int j) const { return allocator_.data(j);}
    /** @return a constant pointer on the main pointer of the container */
    inline PTRCOL* const& p_data() const { return allocator_.p_data();}
    /**  @return @c true if the container is empty, @c false otherwise */
    bool isRef() const { return allocator_.isRef();}
    /** @return the column j.
     *  @param j index of the column
     **/
    SubCol atCol(int j) const
    {
      if (this->beginCols() > j)
      { STKOUT_OF_RANGE_1ARG(IArray2DBase::atCol, j, beginCols() > j);}
      if (this->endCols() <= j)
      { STKOUT_OF_RANGE_1ARG(IArray2DBase::atCol, j, endCols() <= j);}
      return this->asDerived().col(j);
    }
    /** @return the row i.
     *  @param i the index of the row
     **/
    Row atRow(int i) const
    {
      if (this->beginRows() > i)
      { STKOUT_OF_RANGE_1ARG(IArray2DBase::atRow, i, beginRows() > i);}
      if (this->endRows() <= i)
      { STKOUT_OF_RANGE_1ARG(IArray2DBase::at, i, lastIdxRows() < i);}
      return this->asDerived().row(i);
    }

    /** @return the allocator. */
    inline Allocator const& allocator() const { return allocator_;}
    /** @return the maximum possible number of columns without reallocation. */
    int availableCols() const { return availableCols_;}
    /** @return the maximum possible number of rows without reallocation for all Cols. */
    const Array1D<int>& availableRows() const { return availableRows_;}
    /** @return the capacity to used in a column. */
    int capacityByCols() const { return capacityByCols_;}
    /** @return the capacity of the column @c col.
     *  @param col index of the column we want the capacity
     **/
    int capacityCol(int col) const { return availableRows_[col];}
    /** @return the range of each columns. */
    Array1D<Range> const& rangeCols() const { return rangeCols_;}
    /** @return the range of a column.
     *  @param col index of the column we want the range
     **/
    Range const rangeCol(int col) const { return rangeCols_[col];}
    /** Reserve a certain amount of columns
     *  @param sizeCols the size to reserve.
     **/
    void reserveCols(int sizeCols);
    /** New beginning index for the columns of the object.
     *  @param cbeg the index of the first column to set
     **/
    void shiftBeginCols(int cbeg)
    {
      // if there is something to do
      if (cbeg == this->beginCols()) return;
      // is this structure just a pointer?
      if (this->isRef())
      { STKRUNTIME_ERROR_1ARG(IArray2DBase::shiftBeginCols,cbeg,cannot operate on references);}
      // shift beginCols()
      allocator_.shift(cbeg);   // translate data
      availableRows_.shift(cbeg);   // tranlate availableRows_
      rangeCols_.shift(cbeg);       // translate rangeCols_
      Base2D::shiftBeginCols(cbeg); // adjust dimensions
    }
    /** Swapping two columns.
     *  @param pos1, pos2 positions of the columns to swap
     **/
    void swapCols(int pos1, int pos2)
    {
#ifdef STK_BOUNDS_CHECK
      if (this->beginCols() > pos1)
      { STKOUT_OF_RANGE_2ARG(IArray2D::swapCols,pos1, pos2,beginCols() >pos1);}
      if (this->endCols() <= pos1)
      { STKOUT_OF_RANGE_2ARG(IArray2D::swapCols,pos1, pos2,endCols() <= pos1);}
      if (this->beginCols() > pos2)
      { STKOUT_OF_RANGE_2ARG(IArray2D::swapCols,pos1, pos2,beginCols() >pos2);}
      if (this->endCols() <= pos2)
      { STKOUT_OF_RANGE_2ARG(IArray2D::swapCols,pos1, pos2,endCols() <=pos2);}
#endif
      allocator_.swap(pos1, pos2);      // swap allocator part
      std::swap(availableRows_[pos1], availableRows_[pos2]);      // swap availableRows_
      std::swap(rangeCols_[pos1], rangeCols_[pos2]);      // swap rangeCols_
    }
    /** exchange this container with T.
     *  @param T the container to exchange with this
     **/
    void exchange(IArray2DBase &T)
    {
      // swap MemAllocator part
      allocator_.exchange(T.allocator_);
      // swap IContainer2D part
      Base2D::exchange(T);
      // swap this part
      std::swap(availableCols_, T.availableCols_);
      std::swap(capacityByCols_, T.capacityByCols_);
      availableRows_.exchange(T.availableRows_);
      rangeCols_.exchange(T.rangeCols_);
    }
    /** swap two elements: only for vectors and points
     * @param i,j indexes of the elemet to swap
     **/
    void swap(int i, int  j) { std::swap(this->elt(i), this->elt(j)); }
    /** Append the container @c other to @c this without copying the data
     *  explicitly. The column of @c other are appended to this and
     *  @c other will become a reference container. Observe that the @c const
     *  keyword is not respected in this method: but it is useful to
     *  define this method even for constant objects. The data in itself are not
     *  altered, the Array2D become a reference on its own data.
     *  @param other the container to merge with this
     **/
    template<class Other>
    void merge(IArray2DBase<PTRCOL, Other> const& other)
    {
      //checks
      if (this->isRef())
      { STKRUNTIME_ERROR_NO_ARG(IArray2DBase::merge(other),*this is a reference.);}
      if (other.isRef())
      { STKRUNTIME_ERROR_NO_ARG(IArray2DBase::merge(other),other is a reference.);}
      // if there is no columns, we can safely modify the vertical range
      if (this->sizeCols() <= 0) this->setRows(other.rows());
      if (this->rows() != other.rows())
      { STKRUNTIME_ERROR_NO_ARG(IArray2DBase::merge(other),this->rows() != other.rows());}
      // break const reference
      IArray2DBase<PTRCOL, Other>& Tref = const_cast<IArray2DBase<PTRCOL, Other>&>(other);
      // compute horizontal range of the container after insertion
      Range cols(this->cols());
      // compute first index of the first column added
      const int first = cols.end();
      // reallocate memory for the columns
      cols.incLast(Tref.sizeCols());
      reallocCols(cols);
      Base2D::setCols(cols);
      // align other range
      Tref.shiftBeginCols(first); // easiest like that
      // copy data from other
      for (int j=first; j< cols.end(); j++) { transferColumn(Tref, j, j);}
      // delete and set view on the data
      Tref.allocator().free();
      Tref.allocator().setPtr(allocator_.p_data(), Tref.cols(), true);
    }
    /** Append the vector @c other to @c this without copying the data
     *  explicitly. @c other is appended to this and
     *  @c other will become a reference container. The data in itself are not
     *  altered, the Array1D become a reference on its own data.
     *  @param other the container to merge with this
     **/
    template<class Other>
    void merge(IArray1D<Other> const& other)
    {
      // checks
      if (this->isRef())
      { STKRUNTIME_ERROR_NO_ARG(IArray2DBase::merge(IArray1D),*this is a reference.);}
      if (other.isRef())
      { STKRUNTIME_ERROR_NO_ARG(IArray2DBase::merge(IArray1D),other is a reference.);}
      // if there is no columns, we can safely modify the vertical range
      if (this->sizeCols() <= 0) this->setRows(other.range());
      if (this->rows() != other.range())
      { STKRUNTIME_ERROR_NO_ARG(IArray2DBase::merge(IArray1D),this->rows() != other.range());}
      // compute horizontal range of the container after insertion
      Range cols(this->cols());
      // reallocate memory for the columns
      cols.incLast(1);
      this->reallocCols(cols);
      this->setCols(cols);
      // set column
      data(cols.lastIdx()) = other.p_data();
      availableRows_[cols.lastIdx()] = other.capacity();
      rangeCols_[cols.lastIdx()] = other.range();
      // set other as reference
      other.setRef(true);
    }

  protected:
    /** allocator of the column data set */
    Allocator allocator_;
    /** capacity of the columns of the container (for each column: number of
     *  available Rows without reallocation in this column)
    **/
    Array1D<int> availableRows_;
    /** range of the index of the columns of the container. **/
    Array1D<Range> rangeCols_;

    /** set the maximum possible number of columns without reallocation.
     *  @param capacity the maximum number of columns
     **/
    void setAvailableCols(int capacity = 0) { availableCols_ = capacity;}
    /** set the default capacity of a column.
     *  @param capacity the capacity to used
     **/
    void setCapacityByCols(int capacity) { capacityByCols_ = capacity;}
    /** @return the allocator. */
    inline Allocator& allocator() { return allocator_;}
    /** @return a pointer on the j-th column of the container
     *  @param j the index of the column
     **/
    inline PTRCOL& data(int j) { return allocator_.data(j);}
    /** move T to this
     *  @param T the container to move
     **/
    void move(Derived const& T)
    {
      allocator_.move(T.allocator_); // T become a reference
      // move this part
      availableRows_.move(T.availableRows_);
      rangeCols_.move(T.rangeCols_);
      availableCols_ = T.availableCols_;
      // Set IContainer2D part
      this->setCols(T.cols());
      this->setRows(T.rows());
    }
    /** Transfer the column pos2 of the container T to the column
     *  pos1 of this. Set the column pos2 in T to a default value.
     *  The column pos1 should not exists or should be deleted previously
     *  otherwise user will experiment a memory leak.
     *
     *  @param T the container with the column to transfer
     *  @param pos1 index of the column to initialize
     *  @param pos2 the column in the container T to transfer in this
     **/
    template<class Other>
    void transferColumn( IArray2DBase<PTRCOL, Other>& T, int pos1, int pos2)
    {
      // copy column pos2 of T in pos1 of this
      data(pos1) = T.data(pos2);
      // set availableRows_
      availableRows_[pos1] = T.availableRows_[pos2];
      // set rangeCols_
      rangeCols_[pos1] = T.rangeCols_[pos2];
      // set column of T to default
      T.setDefaultCol(pos2);
    }
    /** Method for memory allocation and initialization of the horizontal
     *  range of the container.
     *  The vertical range is not set in this method. If an
     *  error occur, we set the cols_ of the container to default.
     *  @param J horizontal range
     **/
    void mallocCols(Range const& J)
    {
      // compute the size necessary (can be 0)
      int size= Arrays::evalSizeCapacity(J.size());
      // try to allocate memory
      try
      {
        availableRows_.resize(J);        // resize availableRows_ and rangeCols
        rangeCols_.resize(J);
        allocator_.malloc(Range(J.begin(), size)); // allocate memory for the columns
        setAvailableCols(size);
      }
      catch (runtime_error & error)   // if an error occur
      {
        setAvailableCols();     // set default capacity (0)
        Base2D::setCols();      // set default range
        availableRows_.clear(); // clear this->availableRows_
        rangeCols_.clear();     // clear this->rangeCols_
        throw error;            // throw the error
      }
    }
    /** Method for memory reallocation and initialization of the horizontal
     *  range of the container.
     *  The vertical range is not set in this method. If an
     *  error occur, we set the cols_ of the container to default.
     *  @param J horizontal range
     **/
    void reallocCols(Range const& J)
    {
      // compute the necessary size (can be 0)
      int size= Arrays::evalSizeCapacity(J.size());
      // try to allocate memory
      try
      {
        allocator_.realloc(Range(J.begin(), size)); // reallocate memory for the columns
        availableRows_.resize(J);        // initialize this->availableRows_
        rangeCols_.resize(J);        // initialize this->rangeCols_
      }
      catch (runtime_error const& error)   // if an error occur
      {
        setAvailableCols();       // set default capacity (0)
        Base2D::setCols();        // set default range
        availableRows_.clear();   // clear this->availableRows_
        rangeCols_.clear();       // clear this->rangeCols_
        throw error;              // throw the error
      }
      // set new capacity if no error occur
      setAvailableCols(size);
    }
    /** Horizontal Memory deallocation.
     *  This method clear all allocated memory. The range of the columns
     *  is set to (firstCol_:firstCol_-1). The range of the Rows remain
     *  unmodified. If there is allocated memory for the columns, it
     *  should be liberated prior to this method.
     **/
    void freeRows()
    {
      // Nothing to do for reference
      if (this->isRef()) return;
      // free memory allocated in MemAllocator
      allocator_.free();
      setAvailableCols(0);
      Base2D::setCols(Range(this->beginCols(), 0));
      // clear arrays
      availableRows_.clear();
      rangeCols_.clear();
    }

  private:
    /** Horizontal capacity of the container (number of available
     *  columns without reallocation)
     **/
    int availableCols_;
    /** default capacity of a column */
    int capacityByCols_;
    /** set the default parameters and dimension to a column of the container.
     *  @param col the position of the column to initialize to a default value.
     *  @note if data is allocated, it will be lost
     **/
    void setDefaultCol(int col)
    {
      data(col) = 0;              // set column of T to default
      availableRows_[col] = 0;    // set availableRows_
      rangeCols_[col] = Range();  // set rangeCols_
    }
};

/** Reserve a certain amount of columns
 *  @param sizeCols the size to reserve.
 **/
template< class PTRCOL, class Derived>
void IArray2DBase< PTRCOL, Derived>::reserveCols(int sizeCols)
{
  if (availableCols_ >= sizeCols) return;
  Range J(this->beginCols(), sizeCols);
  // try to allocate memory
  try
  {
    // re-allocate memory for the columns
    allocator_.realloc(J);
    availableRows_.resize(J);
    rangeCols_.resize(J);
  }
  catch (runtime_error const& error)   // if an error occur
  {
    setAvailableCols();        // set default capacity (0)
    Base2D::setCols();        // set default range
    this->availableRows_.clear();  // clear this->availableRows_
    this->rangeCols_.clear();      // clear this->rangeCols_
    throw error;                   // throw the error
  }
  // set new capacity if no error occur
  setAvailableCols(sizeCols);
}
} // namespace STK

#endif
// STK_ITARRAY2DBASE_H
