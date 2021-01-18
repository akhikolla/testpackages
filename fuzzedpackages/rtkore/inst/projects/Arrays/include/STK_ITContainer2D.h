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
 * created on: 10 août 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ITContainer2D.h
 *  @brief In this file we define the IContainer2Dand ITContainer2D interface classes.
 **/

#ifndef STK_ITCONTAINER2D_H
#define STK_ITCONTAINER2D_H

#include <Sdk/include/STK_IRecursiveTemplate.h>
#include <Sdk/include/STK_Macros.h>

#include "STK_ArraysTraits.h"
#include "STK_Arrays_Util.h"

namespace STK
{

/** @ingroup Arrays
 *  @brief Interface base class for 2D containers.
 *
 *  The IContainer2D class is the base class for all two-dimensional containers
 *  storing data and which cannot be part of an expression.
 *  A two-dimensional container is defined by an horizontal range of index
 *  for the columns and a vertical range of index for the rows.
 *
 *  This Interface base class stores the ranges and allows to derived classes to
 *  manipulate these ranges.
 **/
template<int SizeRows_, int SizeCols_>
class IContainer2D
{
  public:
    /** Type of the Range for the rows */
    typedef TRange<SizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<SizeCols_> ColRange;

    /** Default constructor. cols_ = 1:0 and rows_ = 1:0. */
    IContainer2D(): rows_(), cols_() {}
    /** Constructor with specified ranges
     *  @param I the vertical range
     *  @param J the horizontal range
     **/
    inline IContainer2D( RowRange const& I, ColRange const& J): rows_(I), cols_(J) {}
    /** Copy constructor
     *  @param T the container to copy
     **/
    inline IContainer2D( IContainer2D const& T): rows_(T.rows_), cols_(T.cols_) {}
    /** destructor. **/
    ~IContainer2D() {}

    /** @return the columns range */
    inline ColRange const &cols() const { return cols_;}
    /** @return the index of the first column */
    inline int beginCols() const { return cols_.begin();}
    /** @return the ending index of the columns */
    inline int endCols() const { return cols_.end();}
    /** @return the number of column */
    inline int sizeCols() const { return cols_.size();}

    /** @return the range of the rows */
    inline RowRange const& rows() const { return rows_;}
    /** @return the index of the first row */
    inline int beginRows() const { return rows_.begin();}
    /** @return the ending index of rows */
    inline int endRows() const { return rows_.end();}
    /** @return the number of rows */
    inline int sizeRows() const { return rows_.size();}

    /** @return the index of the last column */
    inline int lastIdxCols() const { return cols_.lastIdx();}
    /** @return the index of the last row */
    inline int lastIdxRows() const { return rows_.lastIdx();}

    /** @return @c true if the container is empty, @c false otherwise */
    bool empty() const { return (cols_.empty() || rows_.empty());}

  protected:
    /** Set the first index of the rows and columns.
     *  @param rbeg, cbeg the first index of the rows and columns
     **/
    void shift( int rbeg, int cbeg) { rows_.shift(rbeg); cols_.shift(cbeg);}
    /** Set the ranges of the container.
     *  @param I,J the vertical and horizontal range
     **/
    void setRanges(RowRange const& I = RowRange(), ColRange const& J = ColRange())
    { rows_ = I; cols_ =J;}
    // rows
    /** Set the range of the number of rows.
     *  @param I the range of the rows number
     **/
    void setRows( RowRange const& I = RowRange()) { rows_ = I;}
    /** Set the first index of the rows.
     *  @param beg the first index of the rows
     **/
    void shiftBeginRows( int beg) { rows_.shift(beg);}
    /** Increment the range of the number of rows.
     *  @param inc the increment to apply
     **/
    void incRangeRows( int inc) { rows_.inc(inc);}
    /** Increment the first index of the number of rows.
     *  @param inc the increment to apply
     **/
    void incBeginRows( int inc) { rows_.incFirst(inc);}
    /** Decrement the first index of the number of rows.
     *  @param dec the decrement to apply
     **/
    void decBeginRows( int dec) { rows_.decFirst(dec);}
    /** Increment the end of the number of rows.
     *  @param inc the increment to apply
     **/
    void incEndRows( int inc) { rows_.incLast(inc);}
    /** Decrement the end of the number of rows.
     *  @param dec the decrement to apply
     **/
    void decEndRows( int dec) { rows_.decLast(dec);}
    // cols
    /** Set the columns range.
     * @param J the columns range
     **/
    void setCols( ColRange const& J = ColRange()) { cols_ = J;}
    /** Shift the columns first index to beg.
     *  @param beg the new first index
     **/
    void shiftBeginCols( int beg) { cols_.shift(beg);}
    /** Increment the columns range.
     *  @param inc the increment to apply the range
     **/
    void incRangeCols( int inc) { cols_.inc(inc);}
    /** increment the first index of the columns.
     *  @param inc the increment to apply
     **/
    void incBeginCols( int inc) { cols_.incFirst(inc);}
    /** Decrement the columns first index.
     *  @param dec the decrement to apply
     **/
    void decBeginCols( int dec) { cols_.decFirst(dec);}
    /** Increment the last index of the columns.
     *  @param inc the increment to apply
     **/
    void incEndCols( int inc)  { cols_.incLast(inc);}
    /** Decrement the last index of the columns.
     *  @param dec the decrement to apply
     **/
    void decEndCols( int dec) { cols_.decLast(dec);}

    /** exchange this container with T
     *  @param T the container to exchange with this
     **/
     void exchange(IContainer2D& T)
     {
       std::swap(T.rows_, this->rows_ );
       std::swap(T.cols_, this->cols_ );
     }
     /** Increment the end of the number of rows.
      *  @param inc the increment to apply
      **/
     void incLastIdxRows( int inc) { rows_.incLast(inc);}
     /** Decrement the end of the number of rows.
      *  @param dec the decrement to apply
      **/
     void decLastIdxRows( int dec) { rows_.decLast(dec);}
     /** Increment the last index of the columns.
      *  @param inc the increment to apply
      **/
     void incLastIdxCols( int inc)  { cols_.incLast(inc);}
     /** Decrement the last index of the columns.
      *  @param dec the decrement to apply
      **/
     void decLastIdxCols( int dec) { cols_.decLast(dec);}

  private:
    /** Vertical range : Range of the indexes for the rows. */
    RowRange rows_;
    /** Horizontal range : Range of the indexes for the columns. */
    ColRange cols_;
};

/** @ingroup Arrays
 *  @brief Interface base class for homogeneous 2D containers like allocators.
 *
 * The ITContainer2D class is the template base class for all
 * homogeneous two-dimensional containers containing element of type @c Type
 * where Type is note necessarily a scalar. Some methods assume that derived class
 * is not part of an expression and is not constant, so that it can be
 * #- shifted,
 * #- resized
 * #- accessed in modification.
 *
 * Implement the curious recursive template paradigm : the template
 * parameter @c Derived is the name of the class that
 * implements @c ITContainer2D. For example
 * <code>
 * template<class Type>
 * class Derived: public ITContainer2D< Derived<Type> >
 * {...}
 * </code>
 *
 * Functions used in this interface and to implement in
 * derived class if necessary have the following definitions:
 * @code
 *   Type& elt2Impl(int i, int j);
 *   ConstReturnType elt2Impl(int i, int j) const;
 *   Type& elt1Impl(int pos);
 *   ConstReturnType elt1Impl(int pos) const;
 *   Type& elt0Impl()
 *   ConstReturnType elt0Impl(int pos) const;
 *   void shift1Impl(int beg);
 *   void shift2Impl(int beginRows, int beginCols);
 *   Derived& resize1Impl(int size);
 *   resize2Impl(sizeRows, sizeCols);
 * @endcode
 *
 *
 * @sa IContainer2D, ICArray
 *
 * @note The constant getter @c elt1Impl(pos) have to return a reference as we
 * are using derived classes for storing any kind of data.
 **/
template <class Derived>
class ITContainer2D: public IContainer2D< hidden::Traits<Derived>::sizeRows_, hidden::Traits<Derived>::sizeCols_>
                   , public IRecursiveTemplate<Derived>
{
  public:
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef typename hidden::Traits<Derived>::ConstReturnType ConstReturnType;

    typedef typename hidden::Traits<Derived>::Row Row;
    typedef typename hidden::Traits<Derived>::Col Col;

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

    /** Type of the Base container */
    typedef IContainer2D<sizeRows_, sizeCols_ > Base2D;
    /** Type of the Base container */
    typedef IRecursiveTemplate<Derived> Base;

  protected:
    /** Default constructor.*/
    ITContainer2D(): Base2D(), Base() {}
    /** constructor with specified Range.
     *  @param I,J range of the rows and columns
     **/
    inline ITContainer2D( RowRange const& I, ColRange const& J): Base2D(I, J), Base() {}
    /** Copy constructor.
     *  @param T the container to copy
     **/
    inline ITContainer2D( ITContainer2D const& T): Base2D(T), Base() {}
    /** destructor. */
    ~ITContainer2D() {}

  public:
    /** @return the element (i,j) of the 2D container.
     *  @param i, j indexes of the element to get
     **/
    inline Type& elt(int i, int j)
    {
#ifdef STK_BOUNDS_CHECK
      if (this->beginRows() > i) { STKOUT_OF_RANGE_2ARG(ITContainer2D::elt, i, j, beginRows() > i);}
      if (this->endRows() <= i)  { STKOUT_OF_RANGE_2ARG(ITContainer2D::elt, i, j, endRows() <= i);}
      if (this->beginCols() > j) { STKOUT_OF_RANGE_2ARG(ITContainer2D::elt, i, j, beginCols() > j);}
      if (this->endCols() <= j)  { STKOUT_OF_RANGE_2ARG(ITContainer2D::elt, i, j, endCols() <= j);}
#endif
      return this->asDerived().elt2Impl(i,j);
    }
    /** @return a constant reference on element (i,j) of the 2D container
     *  @param i, j indexes of the element to get
     **/
    inline ConstReturnType elt(int i, int j) const
    {
#ifdef STK_BOUNDS_CHECK
      if (this->beginRows() > i) { STKOUT_OF_RANGE_2ARG(ITContainer2D::elt, i, j, beginRows() > i);}
      if (this->endRows() <= i)  { STKOUT_OF_RANGE_2ARG(ITContainer2D::elt, i, j, endRows() <= i);}
      if (this->beginCols() > j) { STKOUT_OF_RANGE_2ARG(ITContainer2D::elt, i, j, beginCols() > j);}
      if (this->endCols() <= j)  { STKOUT_OF_RANGE_2ARG(ITContainer2D::elt, i, j, endCols() <= j);}
#endif
      return this->asDerived().elt2Impl(i,j);
    }
    /** @return a reference on the ith element
     *  @param i index of the element to get
     **/
    inline Type& elt(int i) { return this->asDerived().elt1Impl(i);}
    /** @return the constant ith element
     *  @param i index of the element to get
     *  @note bounds check cannot be done there as number_ does not have
     *  begin() and end() implemented
     **/
    inline ConstReturnType elt(int i) const { return this->asDerived().elt1Impl(i);}
    /** @return a reference on the number */
    inline Type& elt() { return this->asDerived().elt0Impl();}
    /** @return a constant reference on the number */
    inline ConstReturnType elt() const  { return this->asDerived().elt0Impl();}

    /** resize the container
     *  @param sizeRows, sizeCols size of the rows and columns
     **/
    Derived& resize(int sizeRows, int sizeCols)
    { return this->asDerived().resize2Impl(sizeRows, sizeCols);}
    /** Resize 1D container
     *  @param size the size to set to the vector
     **/
    Derived& resize(int size) { return this->asDerived().resize1Impl(size);}

    /** shift the first indexes of the container
     *  @param firstRow, firstCol indexes of the first row and first column
     **/
    void shift( int firstRow, int firstCol) { this->asDerived().shift2Impl(firstRow, firstCol);}
    /** shift the first indexes of the 1D container
     *  @param beg the index of the first row or column
     **/
    void shift(int beg) { this->asDerived().shift1Impl(beg);}

};

} // namespace STK

#endif /* STK_ITCONTAINER2D_H */
