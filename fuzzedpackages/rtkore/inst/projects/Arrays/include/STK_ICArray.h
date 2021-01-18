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
 * Project:  stkpp::Array
 * created on: 10 août 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ICArray.h
 *  @brief Interface base class for the CArray, this is an internal header file,
 *  included by other Containers library headers.
 *
 *  You should not attempt to use it directly but rather used one of the
 *  derived class like CArray, except if you want to create your own
 *  Container Class.
 **/

#ifndef STK_ICARRAY_H
#define STK_ICARRAY_H

#include "STK_IArrayBase.h"

#include "STK_ArrayBaseApplier.h"
#include "STK_ArrayBaseAssign.h"
#include "STK_ArrayBaseInitializer.h"

namespace STK
{
// forward declaration of all CArray classes
template< typename Type, int SizeRows_, int SizeCols_, bool Orient_> class CArray;
template< typename Type, int Size_, bool Orient_> class CArraySquare;
template< typename Type, int SizeCols_, bool Orient_> class CArrayPoint;
template< typename Type, int SizeRows_, bool Orient_> class CArrayVector;
template< typename Type, bool Orient_> class CArrayNumber;

namespace hidden
{
/** @ingroup hidden
 *  @brief The traits struct CSlice must be specialized allow to disambiguate return
 *  type of the col/row/sub operators for CArray family
 *
 *  @sa STK::ICArray
 */
template<typename Derived, int SizeRows, int SizeCols>
struct CSlice
{
    typedef typename Traits< Derived >::Type Type;
    enum
    {
      structure_ = Traits<Derived>::structure_,
      orient_    = Traits<Derived>::orient_,
      sizeRows_  = Traits<Derived>::sizeRows_,
      sizeCols_  = Traits<Derived>::sizeCols_,
      storage_   = Traits<Derived>::storage_,

      isNumber_ = (SizeRows==SizeCols)&&(SizeRows==1),
      isVector_ = (SizeCols == 1)&&(!isNumber_),
      isPoint_  = (SizeRows == 1)&&(!isNumber_),
      isSquare_ = (SizeRows==SizeCols)&&(SizeRows!=UnknownSize)&&(!isNumber_),
      isArray_  = (SizeRows!=SizeCols)&&(!isNumber_) // not used
    };
   typedef typename If< (isNumber_), CArrayNumber<Type, orient_>
                      , typename If< isVector_, CArrayVector<Type, SizeRows, orient_>
                          , typename If< isPoint_, CArrayPoint<Type, SizeCols, orient_>
                              , typename If< isSquare_, CArraySquare<Type, SizeRows, orient_>
                                           , CArray<Type, SizeRows, SizeCols, orient_>
                                           >::Result
                                       >::Result
                                   >::Result
                      >::Result Result;
};

/** @ingroup hidden
 *  helper allowing to disambiguate SubVector access
 **/
template<typename Derived, int Size>
struct CSliceDispatcher
{
    enum
    { structure_ = Traits<Derived>::structure_};
    typedef typename If< structure_ == (int)Arrays::vector_
                       , typename CSlice<Derived, Size, 1>::Result
                       , typename CSlice<Derived, 1, Size>::Result
                       >::Result Result;
};

}
/** @ingroup Arrays
  * @class ICArray
  *s
  * @brief Interface class for CArray, CArrayPoint, CArrayVector, CArraySquare, CArrayNumber.
  *
  * This class is the base that is inherited by all objects (matrix, vector,
  * point) which are not expression and stored as CArrays. The common API for
  * these objects is contained in this class.
  *
  * This is essentially a wrapper of a CAllocator
  *
  * @tparam Derived is the derived type, e.g., a matrix type.
  * @sa CAllocator, CArray, CArrayPoint, CArrayVector, CArraySquare, CArrayNumber
  */
template<class Derived>

class ICArray: public IArrayBase<Derived>
{
  public:
    typedef IArrayBase<Derived> Base;
    typedef typename hidden::Traits< Derived >::Allocator Allocator;

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
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

  protected:
    /** allocator of the memory  */
    Allocator allocator_;
    /** default constructor. */
    ICArray(): Base(), allocator_() {}
    /** constructor with specified sizes.
     *  @param sizeRows,sizeCols size of the rows and columns
     **/
    ICArray( int sizeRows, int sizeCols): Base(), allocator_(sizeRows, sizeCols) {}
    /** constructor with specified sizes and value.
     *  @param sizeRows,sizeCols size of the rows and columns
     *  @param value the value to set
     **/
    ICArray( int sizeRows, int sizeCols, Type const& value): Base(), allocator_(sizeRows, sizeCols, value) {}
    /** copy or wrapper constructor.
     *  @param T size of the rows
     *  @param ref is this owning its own data ?
     **/
    ICArray( Derived const& T, bool ref = false): Base(), allocator_(T.allocator_, ref) {}
    /** wrapper constructor for 0 based C-Array.
     *  @param q pointer on the array
     *  @param sizeRows,sizeCols size of the rows and columns
     **/
    ICArray( Type* const& q, int sizeRows, int sizeCols): Base(), allocator_(q, sizeRows, sizeCols){}
    /** constructor by reference, ref_=1.
     *  @param allocator the allocator to wrap
     *  @param I,J range of the rows and columns to wrap
     **/
    template<class OtherAllocator>
    inline ICArray( ITContainer2D<OtherAllocator> const& allocator, Range const& I, Range const& J)
                 : Base(), allocator_(allocator.asDerived(), I, J)
    {}
    /** constructor by reference, ref_=1.
     *  @param allocator with the data
     **/
    template< class OtherAllocator>
    inline ICArray( ITContainer2D<OtherAllocator> const& allocator)
                  : Base(), allocator_(allocator.asDerived(), true)
    {}
    /**  destructor */
    ~ICArray() {}
  public:
    /** @return the Horizontal range */
    inline ColRange const& colsImpl() const { return allocator_.cols();};
    /** @return the Vertical range */
    inline RowRange const& rowsImpl() const { return allocator_.rows();}

    /** clear all allocated memory . */
    void clear() { allocator_.clear();}

    /** @return @c true if the container is empty, @c false otherwise */
    bool empty() const { return allocator_.empty();}
    /** @return @c true if *this is reference container, @c false otherwise */
    bool isRef() const { return allocator_.isRef();}

    /** Get a constant reference on the main allocator. */
    inline Allocator const& allocator() const { return allocator_;}
    /** @return a constant reference on the main pointer. */
    inline Type* const& p_data() const { return allocator_.p_data();}

    /** implement the const element accessor */
    inline Type& elt2Impl( int i, int j) { return allocator_.elt(i, j);}
    /** implement the writable element accessor */
    inline ConstReturnType elt2Impl( int i, int j) const { return allocator_.elt(i, j);}

    /** implement the const element accessor for vector/point/diagonal arrays*/
    inline Type& elt1Impl( int j) { return allocator_.elt(j);}
    /** implement the writable element accessor for vector/point/diagonal arrays*/
    inline ConstReturnType elt1Impl( int j) const{ return allocator_.elt(j);}

    /** implement the const element accessor for number arrays*/
    inline Type& elt0Impl() { return allocator_.elt();}
    /** implement the writable element accessor for number arrays*/
    inline ConstReturnType elt0Impl() const { return allocator_.elt();}

    // overloaded operators
    /** @return a reference on the element (i,j) of the 2D container.
     *  @param i, j indexes of the element to get
     **/
    inline Type& operator()(int i, int j)
    {
#ifdef STK_BOUNDS_CHECK
       if (this->beginRows() > i) { STKOUT_OF_RANGE_2ARG(ICArray::operator(), i, j, beginRows() > i);}
       if (this->endRows() <= i)  { STKOUT_OF_RANGE_2ARG(ICArray::operator(), i, j, endRows() <= i);}
       if (this->beginCols() > j) { STKOUT_OF_RANGE_2ARG(ICArray::operator(), i, j, beginCols() > j);}
       if (this->endCols() <= j)  { STKOUT_OF_RANGE_2ARG(ICArray::operator(), i, j, endCols() <= j);}
#endif
      return this->elt(i,j);
    }
    /** @return a constant reference on the element (i,j) of the 2D container.
     *  @param i,j row and column indexes
     **/
    inline ConstReturnType operator()(int i, int j) const
    {
#ifdef STK_BOUNDS_CHECK
       if (this->beginRows() > i) { STKOUT_OF_RANGE_2ARG(ICArray::operator(), i, j, beginRows() > i);}
       if (this->endRows() <= i)  { STKOUT_OF_RANGE_2ARG(ICArray::operator(), i, j, endRows() <= i);}
       if (this->beginCols() > j) { STKOUT_OF_RANGE_2ARG(ICArray::operator(), i, j, beginCols() > j);}
       if (this->endCols() <= j)  { STKOUT_OF_RANGE_2ARG(ICArray::operator(), i, j, endCols() <= j);}
#endif
      return this->elt(i,j);
    }

    /** @return a reference on the ith element
     *  @param i index of the element to get
     **/
    inline Type& operator[](int i)
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);
#ifdef STK_BOUNDS_CHECK
      if (this->asDerived().begin() > i) { STKOUT_OF_RANGE_1ARG(ICArray::operator[], i, begin() > i);}
      if (this->asDerived().end() <= i)  { STKOUT_OF_RANGE_1ARG(ICArray::operator[], i, end() <= i);}
#endif
      return this->elt(i);
    }
    /** @return the ith element
     *  @param i index of the element to get
     **/
    inline ConstReturnType operator[](int i) const
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);
#ifdef STK_BOUNDS_CHECK
      if (this->asDerived().begin() > i) { STKOUT_OF_RANGE_1ARG(ICArray::operator[], i, begin() > i);}
      if (this->asDerived().end() <= i)  { STKOUT_OF_RANGE_1ARG(ICArray::operator[], i, end() <= i);}
#endif
      return this->elt(i);
    }

    /** @return the number */
    inline Type& operator()() { return this->elt();}
    /** @return a constant reference on the number */
    inline ConstReturnType operator()() const { return this->elt();}

    // row operators
    /** implement the row operator using a reference on the row of the allocator
     *  @param i index of the row to reference
     **/
    inline typename hidden::CSlice<Derived, 1, sizeCols_>::Result row(int i) const
    { return typename hidden::CSlice<Derived, 1, sizeCols_>::Result( allocator_.row(i));}
    /** implement the row operator using a reference on the row of the allocator
     *  @param i,J index of the row and range of the columns to reference
     **/
    template<int Size_>
    inline typename hidden::CSlice<Derived, 1, Size_>::Result row(int i, TRange<Size_> const& J) const
    { return typename hidden::CSlice<Derived, 1, Size_>::Result( allocator_.row( i, J));}
    /** @param i,J index of the row and range of the columns
     *  @return an Horizontal container referencing row @c i in range @c J
     **/
    template<int Size_>
    inline typename hidden::CSlice<Derived, 1, Size_>::Result operator()(int i, TRange<Size_> const& J) const
    { return row(i, J);}
    /** implement the row operator using a reference on a range of rows of the allocator
     *  @param I range of the rows to reference
     **/
    template<int Size_>
    inline typename hidden::CSlice<Derived, Size_, sizeCols_>::Result row(TRange<Size_> const& I) const
    { return typename hidden::CSlice<Derived, Size_, sizeCols_>::Result( allocator_.sub(I, this->cols()));}

    // col operators
    /** implement the col operator using a reference on the column of the allocator
     * @param j index of the column to reference
     **/
    inline typename hidden::CSlice<Derived, sizeRows_, 1>::Result col(int j) const
    { return typename hidden::CSlice<Derived, sizeRows_, 1>::Result( allocator_.col(j));}
    /** implement the col operator using a reference on the column of the allocator
     * @param I,j range of the rows and index of the column to reference
     **/
    template<int Size_ >
    inline typename hidden::CSlice<Derived, Size_, 1>::Result col(TRange<Size_> const& I, int j) const
    { return typename hidden::CSlice<Derived, Size_, 1>::Result( allocator_.col( I, j));}
    /** implement the col operator using a reference on a range of columns of the allocator
     *  @param J range of columns to reference
     **/
    template<int Size_>
    inline typename hidden::CSlice<Derived, sizeRows_, Size_>::Result col(TRange<Size_> const& J) const
    { return typename hidden::CSlice<Derived, sizeRows_, Size_>::Result( allocator_.sub(this->rows(), J));}
    /** @param I,j range of the rows and index of the column to reference
     *  @return a Vertical container containing the column @c j of this in the range @c I
     **/
    template<int Size_>
    inline typename hidden::CSlice<Derived, Size_, 1>::Result operator()(TRange<Size_> const& I, int j) const
    { return col(I, j);}

    // sub operators
    /** implement the sub operator for 1D arrays using a reference on the row/column of the allocator
     *  @param J range to get
     **/
    template<int Size>
    inline typename hidden::CSliceDispatcher<Derived, Size>::Result sub( TRange<Size> const& J) const
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);
      return typename hidden::CSliceDispatcher<Derived, Size>::Result( allocator_.subVector(J));
    }
    /** @return the sub-vector in given range
     *  @param I range to get
     **/
    template<int Size>
    inline typename hidden::CSliceDispatcher<Derived, Size>::Result operator[](TRange<Size> const& I) const
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);
      return this->sub(I);
    }
    /** implement the sub operator for 2D arrays using references on a range of rows and columns
     *  of the allocator
     *  @param I,J range of the rows and columns to reference
     **/
    template<int OtherRows_, int OtherCols_>
    inline typename hidden::CSlice<Derived, OtherRows_, OtherCols_>::Result sub(TRange<OtherRows_> const& I, TRange<OtherCols_> const& J) const
    { return typename hidden::CSlice<Derived, OtherRows_, OtherCols_>::Result(allocator_.sub(I, J));}
    /** @param I,J range of the rows and columns
     *  @return a 2D container containing this in the range @c I, @c J
     **/
    template<int OtherRows_, int OtherCols_>
    inline typename hidden::CSlice<Derived, OtherRows_, OtherCols_>::Result operator()(TRange<OtherRows_> const& I, TRange<OtherCols_> const& J) const
    { return sub(I, J);}

    /** swap two elements: only for vectors an points. */
    void swap(int i, int  j) { std::swap(this->elt(i), this->elt(j)); }
    /** @param pos1, pos2 positions of the columns to swap */
    void swapCols(int pos1, int pos2)
    {
#ifdef STK_BOUNDS_CHECK
      if (this->beginCols() > pos1)
      { STKOUT_OF_RANGE_2ARG(ICArray::swapCols,pos1, pos2,beginCols() > pos1);}
      if (this->endCols() <= pos1)
      { STKOUT_OF_RANGE_2ARG(ICArray::swapCols,pos1, pos2,endCols() <= pos1);}
      if (this->beginCols() > pos2)
      { STKOUT_OF_RANGE_2ARG(ICArray::swapCols,pos1, pos2,beginCols() > pos2);}
      if (this->endCols() <= pos2)
      { STKOUT_OF_RANGE_2ARG(ICArray::swapCols,pos1, pos2,endCols() <= pos2);}
#endif
      // swap allocator
      allocator_.swapCols(pos1, pos2);
    }
    /** @param pos1, pos2 positions of the rows to swap */
    void swapRows(int pos1, int pos2)
    {
#ifdef STK_BOUNDS_CHECK
      if (this->beginRows() > pos1)
      { STKOUT_OF_RANGE_2ARG(ICArray::swapRows,pos1, pos2,beginRows() > pos1);}
      if (this->endRows() <= pos1)
      { STKOUT_OF_RANGE_2ARG(ICArray::swapRows,pos1, pos2,endRows() <= pos1);}
      if (this->beginRows() > pos2)
      { STKOUT_OF_RANGE_2ARG(ICArray::swapRows,pos1, pos2,beginRows() > pos2);}
      if (this->endRows() <= pos2)
      { STKOUT_OF_RANGE_2ARG(ICArray::swapRows,pos1, pos2,endRows() <= pos2);}
#endif
      // swap allocator
      allocator_.swapRows(pos1, pos2);
    }
    /** exchange this with T.
     *  @param T the container to exchange with this
     **/
    void exchange(Derived& T) { allocator_.exchange(T.allocator_);}
    /** move T to this.
     *  @param T the array to move
     **/
    void move(Derived const& T) { allocator_.move(T.allocator_);}
    /** shift the Array.
     *  @param beginRows,beginCols  first indexes of the rows and columns
     **/
    Derived& shift(int beginRows, int beginCols)
    {
      STK_STATIC_ASSERT((structure_ == (int)Arrays::array2D_)
                      ||(structure_ == (int)Arrays::lower_triangular_)
                      ||(structure_ == (int)Arrays::upper_triangular_)
                      ,YOU_CANNOT_USED_THIS_METHOD_WITH_THIS_KIND_OF_ARRAY);
      if (!hidden::CheckShift<Derived, structure_>::shift(this->asDerived(), beginRows, beginCols)) return this->asDerived();
      if (this->isRef())
      { STKRUNTIME_ERROR_2ARG(ICArray::shift,beginRows,beginCols,cannot operate on reference);}
      allocator_.shift(beginRows, beginCols);
      return this->asDerived();
    }
    /** shift the Array.
     *  @param firstIdx first index of the vector/point/diagonal/square array.
     *  @note if this method is used with arrays, upper triangular and lower triangular
     *  arrays, both indexes will be shifted.
     **/
    Derived& shift(int firstIdx)
    {
      if (!hidden::CheckShift<Derived, structure_>::shift(this->asDerived(), firstIdx)) return this->asDerived();
      if (this->isRef())
      { STKRUNTIME_ERROR_1ARG(ICArray::shift,firstIdx,cannot operate on reference);}
      allocator_.shift(firstIdx);
      return this->asDerived();
    }
    /** resize the Array.
     *  @param I, J range of the rows and columns
     **/
    Derived& resize(Range const& I, Range const& J)
    {
      STK_STATIC_ASSERT( (structure_ == (int)Arrays::array2D_)
                       ||(structure_ == (int)Arrays::square_)
                       ||(structure_ == (int)Arrays::diagonal_)
                       ||(structure_ == (int)Arrays::lower_triangular_)
                       ||(structure_ == (int)Arrays::upper_triangular_)
                       ||(structure_ == (int)Arrays::symmetric_)
                       ||(structure_ == (int)Arrays::lower_symmetric_)
                       ||(structure_ == (int)Arrays::upper_symmetric_)
                      ,YOU_CANNOT_USED_THIS_METHOD_WITH_THIS_KIND_OF_ARRAY);
      if (!hidden::CheckShift<Derived, structure_>::isAllowed(this->asDerived(), I, J))
      { STKRUNTIME_ERROR_2ARG(ICArray::resize,I,J,not permited);}
      if (!hidden::CheckShift<Derived, structure_>::resize(this->asDerived(), I, J)) return this->asDerived();
      if (this->isRef())
      { STKRUNTIME_ERROR_2ARG(ICArray::resize,I,J,cannot operate on reference);}
      allocator_.resize(I.size(), J.size()).shift(I.begin(), J.begin());
      return this->asDerived();
    }
    /** Resize the vector/point/diagonal/square array.
     *  @param I Range of the vector
     *  @note if this method is used with arrays, upper triangular and lower triangular
     *  arrays, both ranges will be resized.
     **/
    template<int Size_>
    Derived& resize(TRange<Size_> const& I)
    {
      if (!hidden::CheckShift<Derived, structure_>::resize(this->asDerived(), I)) return this->asDerived();
      if (this->isRef())
      { STKRUNTIME_ERROR_1ARG(ICArray::resize,I,cannot operate on reference);}
      allocator_.resize(I.size()).shift(I.begin());
      return this->asDerived();
    }
    /** Resize the vector/point/diagonal/square array.
     *  @param size Range of the vector/point/diagonal/square array
     *  @note if this method is used with arrays, upper triangular and lower triangular
     *  arrays, both ranges will be resized.
     **/
    Derived& resize(int size)
    {
      if (!hidden::CheckShift<Derived, structure_>::resize(this->asDerived(), size)) return this->asDerived();
      if (this->isRef())
      { STKRUNTIME_ERROR_1ARG(ICArray::resize,size,cannot operate on reference);}
      allocator_.resize(size);
      return this->asDerived();
    }
};

} // namespace STK

#endif /* STK_ICARRAY_H */
