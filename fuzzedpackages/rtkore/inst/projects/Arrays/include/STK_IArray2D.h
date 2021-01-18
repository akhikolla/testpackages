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

/** @file STK_IArray2D.h
 *  @brief Interface base class for the Array1D, this is an internal header file,
 *  included by other Containers library headers.
 *
 *  You should not attempt to use it directly but rather used one of the
 *  derived class like Array2D, except if you want to create your own
 *  Container Class.
 **/

#ifndef STK_IARRAY2D_H
#define STK_IARRAY2D_H

#include "STK_IArray2DBase.h"

namespace STK
{

/** @ingroup Arrays
  * @brief Interface base class for two-dimensional arrays.
  *
  * A IArray2D is a specialized interface class for two-dimensional
  * arrays stored in columns. All derived class from @c IArray2D
  * access to a column using a @c Type* ptr.
  *
  * All memory is allocated there using the @c rangeRowsInCol method.
  *
  * The available methods for manipulating the derived (2D) arrays are
  * @code
  *     void shift( int rbeg, int cbeg);
  *     void pushBackRows( int n=1);
  *     void pushBackRows( int n=1);
  *     void pushFrontRows( int n=1);
  *     void insertRows(int pos, int n=1);
  *     void eraseRows(int pos, int n=1);
  *     void pushBackCols( int n=1);
  *     void pushFrontCols( int n=1);
  *     void insertCols(int pos, int n=1);
  *     void eraseCols(int pos, int n=1);
  *     Derived& resize( Range const& I, Range const& J);
  *     template<class Other>
  *     Derived& pushFrontRows(ExprBase<Other> const& other);
  *     template<class Other>
  *     Derived& pushBackRows(ExprBase<Other> const& other);
  *     template<class Other>
  *     Derived& pushFrontCols(ExprBase<Other> const& other);
  *     template<class Other>
  *     Derived& pushBackCols(ExprBase<Other> const& other);
  *     template<class Other>
  *     Derived& pushFrontCols(ExprBase<Other> const& other);
  *     template<class Other>
  *     Derived& pushBackCols(ITContainer1D<Other> const& other);
  * @endcode
  *
  * @tparam Derived is the name of the class implementing an @c IArray2D.
  * @sa Array2D, @sa Array2DDiagonal, @sa Array2DLowerTriangular,
  * @sa Array2DUpperTriangular, @sa Array2DPoint, @sa Array2DVector;
  * @sa Array2DSquare.
 **/
template < class  Derived  >
class IArray2D: public IArray2DBase< typename hidden::Traits<Derived>::Type*, Derived>
{
  public:
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef typename hidden::Traits<Derived>::Row  Row;
    typedef typename hidden::Traits<Derived>::Col  Col;

    enum
    {
      structure_ = hidden::Traits<Derived>::structure_,
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };

     typedef IArray2DBase< Type*, Derived> Base;

 protected:
    /** Default constructor */
    IArray2D(): Base() {}
    /** Constructor with specified ranges
     *  @param I range of the Rows
     *  @param J range of the Cols
     **/
    IArray2D( Range const& I, Range const& J): Base(I, J)
    { initializeCols(J);}
    /** Copy constructor
     *  @param T the array to copy
     *  @param ref true if we wrap T
     **/
    IArray2D( const IArray2D& T, bool ref =false): Base(T, ref)
    {
      if (!ref)
      {
        // initialize the Columns and Rows
        initializeCols(T.cols());
        for (int j=T.beginCols(); j<T.endCols(); j++)
        { copyColumnForward(T, j, j);}
      }
    }
    /** constructor by reference, ref_=1.
     *  @param T the array to copy
     *  @param I,J ranges of the rows and columns to wrap
     **/
    template<class OtherArray>
    IArray2D( IArray2D<OtherArray> const& T, Range const& I, Range const& J)
           : Base(T, I, J)
    {}
    /** Wrapper constructor the Container is a ref.
     *  @param q pointer on data
     *  @param I range of the Rows to wrap
     *  @param J range of the Columns to wrap
     **/
    IArray2D( Type** q, Range const& I, Range const& J): Base(q, I, J) {}
    /** destructor.
     *  free the vertically allocated memory (the columns). The horizontally
     *  allocated memory is handled by the Allocator class.
     **/
    ~IArray2D() { if (!this->isRef()) this->freeCols(this->cols());}

  public:
    /** set a value to the whole array */
    Derived& setValue(Type const& v)
    {
      for (int j=this->beginCols(); j<this->endCols(); j++)
      {
        Type* p(this->data(j));
        for (int i=this->rangeCols_[j].begin(); i<this->rangeCols_[j].end(); i++) p[i]= v;
      }
      return this->asDerived();
    }
    /** move T to this.
     *  @note : T is not modified but just set as a reference of the data it was responsible.
     *  @param T the array to move.
     **/
     Derived& move(Derived const& T)
     {
       if (this->asPtrDerived() == &T) return this->asDerived();
       if (!this->isRef()) { freeCols(this->cols());}
       // move Base part
       Base::move(T);
       return this->asDerived();
     }
    /** clear the object.
     *  This will free all allocated memory and reset all range to Range().
     **/
    void clear()
    {
      // Nothing to do for reference
      if (this->isRef()) return;
      // free the Rows memory
      this->freeMem();
      this->setRanges();
      // initialize if necessary
      this->mallocCols(this->cols());
      initializeCols(this->cols());
    }
    /** @brief Set new beginning indexes to the array.
     *  @param rbeg, cbeg the indexes of the first row and first column to set
     **/
    void shift( int rbeg, int cbeg)
    {
      // move begin of the columns
      this->shiftBeginCols(cbeg);
      // move begin of the rows
      this->shiftBeginRows(rbeg);
    }
    /** resize the array.
     *
     *  @note The implicit assumption made by this method is that it is easiest
     *  and faster to add column than add rows to the 2D array.
     *
     * @param I the new range for the rows of the array
     * @param J the new range for the columns of the array
     **/
    Derived& resize( Range const& I, Range const& J)
    {
      // check if there is something to do
      if ((this->rows() == I) && (this->cols() == J)) return this->asDerived();
      if (this->isRef())
      { STKRUNTIME_ERROR_2ARG(IArray2D::resize,I,J,cannot operate on reference);}
      //  translate beg
      this->shift(I.begin(), J.begin());
      // check again if there is something to do
      if ((this->rows() == I) && (this->cols() == J)) return this->asDerived();
      // just clear empty container
      if (I.size()<=0 || J.size() <= 0) { this->clear(); return this->asDerived();}
      // number of rows and columns to delete or add
      int rinc = I.end() - this->endRows();
      int cinc = J.end() - this->endCols();
      // check if we add columns
      if ((cinc >=0)) // work first on rows as we add columns
      {
        if (rinc < 0)
        { this->popBackRows(-rinc);}  // less rows
        else
        { this->pushBackRows(rinc);} // more rows
        this->pushBackCols(cinc); // add columns
        return this->asDerived();
      }
      // work first on columns as we remove column
      this->popBackCols(-cinc); // remove columns
      if (rinc < 0) this->popBackRows(-rinc); // less rows
      else          this->pushBackRows(rinc); // more rows
      return this->asDerived();
    }
    /** New first index for the object.
     *  @param beg the index of the first element to set
     **/
    void shift( int beg) { this->asDerived().shift1D(beg);}
    /** @return the resized row or column vector
     *  @param I the new range for the vector/point
     **/
    Derived& resize( Range const& I) { return this->asDerived().resize1D(I);}
    /** New first index for the Rfws of the array.
     *  @param rbeg the index of the first row to set
     **/
    void shiftBeginRows( int rbeg)
    {
      // compute increment
      int rinc = rbeg - this->beginRows();
      // if there is something to do
      if (rinc == 0) return;
      // is this structure just a pointer?
      if (this->isRef())
      { STKRUNTIME_ERROR_1ARG(IArray2D::shiftBeginRows,rbeg,cannot operate on reference);}
      // translate rows_()
      Base::shiftBeginRows(rbeg);
      // For all cols, move begin
      for (int j=this->beginCols(); j<this->endCols(); j++)
      { shiftCol(j, this->rangeCols_[j].begin()+rinc);}
    }
    /** Add n Rows to the array.
     *  @param n number of Rows to add
     **/
    void pushBackRows( int n=1)
    {
      // if n<=0 nothing to do
      if (n <= 0) return;
      // is this structure just a pointer?
      if (this->isRef())
      { STKRUNTIME_ERROR_1ARG(IArray2D::pushBackRows,n,cannot operate on reference);}
      // If the array have no rows : create its
      if (this->sizeRows() <=0)
      {
        // update the range of the array
        this->incLastIdxRows(n);
        // initialize the array
        this->initializeCols(this->cols());
      }
      else
      {
        // update the range of the rows
        this->incLastIdxRows(n);
        // allocate new Rows for each Col
        for (int j=this->beginCols(); j<this->endCols(); j++)
        {
          // compute range from the leaf
          Range range(this->asDerived().rangeRowsInCol(j));
          // if there is no column or the end is less than the array
          // end
          if ((range.size()>0)&&(range.lastIdx()>this->lastIdxRows()-n))
          {
            // if the column is empty create it
            if (this->rangeCols_[j].size()<=0)
            {
              this->initializeCol(j, range);
            }
            else
            {
              // compute position
             int pos(this->lastIdxRows()-n+1);
              // add elts
              insertRowsToCol(j, pos, range.lastIdx() - pos +1);
            }
          }
        }
      }
    }
    /** Insert n Rows in front of the array.
     *  @param n number of elements to insert (default 1)
     **/
    void pushFrontRows(int n =1)
    {
      // if n<=0 nothing to do
      if (n <= 0) return;
      // is this structure just a pointer?
      if (this->isRef())
      { STKRUNTIME_ERROR_1ARG(IArray2D::pushFrontRows,n,cannot operate on reference);}
      int pos = this->beginRows();      // update the range of the rows
      this->incLastIdxRows(n);
      // allocate new Rows for each Col
      for (int j=this->beginCols(); j<this->endCols(); j++)
      {
        // check position
        if ( (pos >= this->rangeCols_[j].begin())
           ||(pos <= this->rangeCols_[j].end())
           )
        { insertRowsToCol(j, pos, n);}
      }
    }
    /** Insert n Rows at the position pos of the array.
     *  If pos is outside the range of a column, then the
     *  method do nothing.
     *  @param pos index where to insert Rows
     *  @param n number of elements to insert (default 1)
     **/
    void insertRows(int pos, int n =1)
    {
      // if n<=0 nothing to do
      if (n <= 0) return;
      // is this structure just a pointer?
      if (this->isRef())
      { STKRUNTIME_ERROR_2ARG(IArray2D::insertRows,pos,n,cannot operate on reference);}
      if (this->beginRows() > pos)
      { STKOUT_OF_RANGE_2ARG(IArray2D::insertRows,pos,n,beginRows() > pos);}
      if (this->endRows() < pos)
      { STKOUT_OF_RANGE_2ARG(IArray2D::insertRows,pos,n,endRows() < pos);}
      // update the range of the rows
      this->incLastIdxRows(n);
      // allocate new Rows for each Col
      for (int j=this->beginCols(); j<this->endCols(); j++)
      {
        // check position
        if ( (pos >= this->rangeCols_[j].begin())
           ||(pos <= this->rangeCols_[j].end())
           )
        { insertRowsToCol(j, pos, n);}
      }
    }
    /** Delete n latest rows of the array.
     *  @param n number of rows to delete
     **/
    void popBackRows( int n = 1)
    {
      // if n<=0 nothing to do
      if (n <= 0) return;
      // is this structure just a pointer?
      if (this->isRef())
      { STKRUNTIME_ERROR_1ARG(IArray2D::popBackRows,n,cannot operate on reference);}
      if (this->sizeRows() < n)
      { STKOUT_OF_RANGE_1ARG(IArray2D::popBackRows,n,sizeRows() < n);}
      this->decLastIdxRows(n);
      // decrease range of each Col
      for (int j= this->beginCols(); j< this->endCols(); j++)
        eraseRowsToCol(j, this->endRows(), n);
    }

    /** Delete n Rows at the pos index to the array.
     *  @param pos index where to delete elements
     *  @param n number of elements to delete (default 1)
    **/
    void eraseRows(int pos, int n=1)
    {
      // if n==0 nothing to do
      if (n<=0) return;
      // is this structure just a pointer?
      if (this->isRef())
      { STKRUNTIME_ERROR_2ARG(IArray2D::eraseRows,pos,n,cannot operate on reference);}
      if (this->beginRows() > pos)
      { STKOUT_OF_RANGE_2ARG(IArray2D::eraseRows,pos,n,beginRows() > pos);}
      if (this->lastIdxRows() < pos)
      { STKOUT_OF_RANGE_2ARG(IArray2D::eraseRows,pos,n,lastIdxRows() < pos);}
      if (this->lastIdxRows() < pos+n-1)
      { STKOUT_OF_RANGE_2ARG(IArray2D::eraseRows,pos,n,lastIdxRows() < pos+n-1);}
      // update each Col
      for (int j=this->beginCols(); j<this->endCols(); j++)
        eraseRowsToCol(j, pos, n);
      // update dimensions
      this->decLastIdxRows(n);
    }
    /** Add n Columns at the end of the array.
     *  @param n the number of Columns to add
     **/
    void pushBackCols(int n = 1)
    {
      // if n<=0 nothing to do
      if (n <= 0) return;
      // is this structure just a pointer?
      if (this->isRef())
      { STKRUNTIME_ERROR_1ARG(IArray2D::pushBackCols,n,cannot operate on reference);}
      // If the array have no Columns : create its
      if (this->sizeCols() <=0)
      {
        this->incLastIdxCols(n);
        this->mallocCols( this->cols());
        initializeCols( this->cols());
      }
      else // else insert to the end of the array
      { insertCols(this->endCols(), n);}
    }
    /** Insert n Columns at the beginning of the array.
     *  @param n the number of column to insert
     **/
    void pushFrontCols(int n =1)
    {
      // if n<=0 nothing to do
      if (n <= 0) return;
      // is this structure just a pointer?
      if (this->isRef())
      { STKRUNTIME_ERROR_1ARG(IArray2D::pushFrontCols,n,cannot operate on reference);}
      // compute horizontal range of the array after insertion
      Range range_ho(this->cols());
      range_ho.incLast(n);
      // allocate, if necessary, the mem for the Cols
      if (this->availableCols() < range_ho.size()) //  not enough space
      {
        // exchange with Taux
        Derived Taux;
        this->exchange(Taux);
        // initialize columns of the array
        try
        {
          this->mallocCols(range_ho);
        }
        catch (Exception const& error)   // if an error occur
        {
          this->exchange(Taux);   // restore array
          throw error;            // and send again the exception
        }
        // set the ranges
        this->setCols(range_ho);
        this->setRows(Taux.rows());
        // translate and copy last Columns from Taux to this
        for (int k=Taux.lastIdxCols(); k>=Taux.beginCols(); k--)
          this->transferColumn(Taux, k+n, k);
      }
      else // enough space -> shift the last Cols
      {
        Range addRange(this->endCols(), n);
        // insert default capacity for the new Columns
        this->availableRows_.insert(addRange, 0);
        // insert default range for the new Columns
        this->rangeCols_.insert(addRange, Range());
        // update range_
        this->incLastIdxCols(n);
        // translate data
        for (int k=this->lastIdxCols()-n; k>=this->beginCols(); k--)
          this->transferColumn( this->asDerived(), k+n, k);
      }
      // initialize the rows for n first Columns
      this->initializeCols(Range(this->beginCols(), n));
    }
    /** Insert n Columns at the index pos to the array.
     *  @param pos the position of the inserted Columns
     *  @param n the number of column to insert
     **/
    void insertCols(int pos, int n =1)
    {
      // if n<=0 nothing to do
      if (n <= 0) return;
      // is this structure just a pointer?
      if (this->isRef())
      { STKRUNTIME_ERROR_2ARG(IArray2D::insertCols,pos,n,cannot operate on reference);}
      if (this->beginCols() > pos)
      { STKOUT_OF_RANGE_2ARG(IArray2D::insertCols,pos,n,beginCols() > pos);}
      if (this->endCols() < pos)
      { STKOUT_OF_RANGE_2ARG(IArray2D::insertCols,pos,n,endCols() < pos);}
      // compute horizontal range of the array after insertion
      Range range_ho(this->cols());
      range_ho.incLast(n);
      // allocate, if necessary, the mem for the Cols
      if (this->availableCols() < range_ho.size()) //  not enough space
      {
        // exchange with Taux
        Derived Taux;
        this->exchange(Taux);
        // initialize columns of the array
        try
        {
          this->mallocCols(range_ho);
        }
        catch (Exception & error)   // if an error occur
        {
          this->exchange(Taux);   // restore array
          throw error;        // and send again the exception
        }
        // set the ranges
        this->setCols(range_ho);
        this->setRows(Taux.rows());
        // move first Columns from Taux to this
        for (int k=Taux.beginCols(); k<pos; k++)
          this->transferColumn(Taux, k, k);
        // translate and copy last Columns from Taux to this
        for (int k=Taux.lastIdxCols(); k>=pos; k--)
          this->transferColumn(Taux, k+n, k);
      }
      else // enough space -> shift the last Cols
      {
        Range addRange(this->lastIdxCols()+1, n);
        // insert default capacity for the new Columns
        this->availableRows_.insert(addRange, 0);
        // insert default range for the new Columns
        this->rangeCols_.insert(addRange, Range());
        // update range_
        this->incLastIdxCols(n);
        // translate data
        for (int k=this->lastIdxCols()-n; k>=pos; k--)
          this->transferColumn( this->asDerived(), k+n, k);
      }
      // initialize the rows for the Cols, this->availableRows_, this->rangeCols_
      // in the range pos:pos+n-1
      this->initializeCols(Range(pos, n));
    }
    /** Delete last Columns of the array
     *  @param n the number of Columns to delete
     **/
    void popBackCols( int n =1)
    {
      // if n<=0 nothing to do
      if (n<=0) return;
      // is this structure just a pointer?
      if (this->isRef())
      { STKRUNTIME_ERROR_1ARG(IArray2D::popBackCols,n,cannot operate on reference);}
      if (this->sizeCols() < n)
      { STKOUT_OF_RANGE_1ARG(IArray2D::popBackCols,n,sizeCol() < n);}
      // delete each col
      this->freeCols(Range(this->lastIdxCols()-n+1, this->lastIdxCols(), 0));
      // update this->availableRows_
      this->availableRows_.popBack(n);
      // update this->rangeCols_
      this->rangeCols_.popBack(n);
      // update cols
      this->decLastIdxCols(n);
      // if there is no more Cols
      if (this->sizeCols() == 0) this->freeMem();
    }
    /** Delete n Columns at the specified position of the array.
     *  @param pos the position of the deleted Columns
     *  @param n the number of column to delete
     **/
    void eraseCols(int pos, int n = 1)
    {
      if (n<=0) return; // if n<=0 nothing to do
      // is this structure just a pointer?
      if (this->isRef())
      { STKRUNTIME_ERROR_2ARG(IArray2D::eraseCols,pos,n,cannot operate on reference);}
      if (this->beginCols() > pos)
      { STKOUT_OF_RANGE_2ARG(IArray2D::eraseCols,pos,n,beginCols() > pos);}
      if (this->lastIdxCols() < pos)
      { STKOUT_OF_RANGE_2ARG(IArray2D::eraseCols,pos,n,lastIdxCols() < pos);}
      if (this->lastIdxCols() < pos+n-1)
      { STKOUT_OF_RANGE_2ARG(IArray2D::eraseCols,pos,n,lastIdxCols() < pos+n-1);}
      // delete each col
      this->freeCols(Range(pos, n));
      // update cols_
      this->decLastIdxCols(n);
      // shift Cols
      for (int k=pos; k<this->endCols(); k++)
        this->data(k) = this->data(k+n);
      // update this->availableRows_
      this->availableRows_.erase(pos, n);
      // update this->rangeCols_
      this->rangeCols_.erase(pos, n);
      // if there is no more Cols
      if (this->sizeCols() == 0) this->freeMem();
    }
    /** Update the columns of the array in the specified range.
     *  @param J range of the column to update
     **/
    void update(Range const& J)
    {
      if (this->beginCols() > J.begin())
      { STKOUT_OF_RANGE_1ARG(IArray2D::update,J,beginCols() > J.begin());}
      if (this->endCols() < J.end())
      { STKOUT_OF_RANGE_1ARG(IArray2D::update,J,endCols() < J.end());}
      for ( int icol = J.begin(); icol < J.end() ; ++icol)
      { update(icol);}
    }

    /** Update the column of the array in the specified position.
     *  @param col index of the column to update
     **/
    void update(int col)
    {
      if (this->beginCols() > col)
      { STKOUT_OF_RANGE_1ARG(IArray2D::update,col,beginCols() > col);}
      if (this->lastIdxCols() < col)
      { STKOUT_OF_RANGE_1ARG(IArray2D::update,col,lastIdxCols() < col);}
      if (this->asDerived().rangeRowsInCol(col) != this->rangeCol(col))
      { resizeCol(col, this->asDerived().rangeRowsInCol(col));}
    }

    /** overwrite @c this with @c src.
     *  @note If the size match, @c this is not resized, and in this case,
     *  the method take care of the possibly of overlapping.
     *  @param src the array to copy
     *  @return a copy of src
     **/
    Derived& assign( IArray2D const& src);

    /** merge (by value) the array other with this.
     *  @param other the array to merge to this
     **/
    template<class Other>
    Derived& pushFrontCols(ExprBase<Other> const& other)
    {
      // check if the array is empty
      if (this->empty())
      {
        this->asDerived() = other.asDerived();
        return this->asDerived();
      }
      // this is not empty
      if (other.rows() != this->rows())
      { STKRUNTIME_ERROR_NO_ARG(IArray2D::pushFrontCols,range of the rows are different);}
      // if the array is not empty we add the column and copy other inside
      int size = other.sizeCols(), first = this->beginCols();
      insertCols(first,size);
      for (int j0= first, j1= other.beginCols(); j1 < other.endCols(); ++j0, ++j1)
      {
        for (int i=other.beginRows(); i<other.endRows(); i++)
          (*this)(i, j0) = other(i,j1);
      }
      // return this
      return this->asDerived();
    }
    /** merge (by value) the array other with this.
     *  @param other the array to merge to this
     **/
    template<class Other>
    Derived& pushBackCols(ExprBase<Other> const& other)
    {
      // check if the array is empty
      if (this->empty())
      {
        this->asDerived() = other.asDerived();
        return this->asDerived();
      }
      // this is not empty
      if (other.rows() != this->rows())
      { STKRUNTIME_ERROR_NO_ARG(IArray2D::pushBackCols,range of the rows are different);}
      // if the array is not empty we add the column and copy other inside
     int size = other.sizeCols(), first = this->lastIdxCols()+1;
      pushBackCols(size);
      for (int j0= first, j1= other.beginCols(); j1 < other.endCols(); ++j0, ++j1)
      {
        for (int i=other.beginRows(); i<other.endRows(); i++)
          (*this)(i, j0) = other(i,j1);
      }
      // return this
      return this->asDerived();
    }
    /** Specialization for Array1D. merge (by value) the array other with this
     *  @param other the column to add to this
     **/
    template<class Other>
    Derived& pushBackCols(IArray1D<Other> const& other)
    {
      // check if the array is empty
      if (this->empty())
      {
        resize(other.rows(),1);
        for (int i=other.begin(); i<other.end(); i++)
          (*this)(i, baseIdx) = other[i];
        return this->asDerived();
      }
      // not empty
      if (other.rows() != this->rows())
      { STKRUNTIME_ERROR_NO_ARG(IArray2D::pushBackCols(other),other.rows() != rows());}
      // if the array is not empty we add the column and copy other inside
      int size = other.sizeCols(), first = this->lastIdxCols()+1;
      pushBackCols(size);
      for (int i=other.begin(); i<other.end(); i++)
        (*this)(i, first) = other[i];
      // return this
      return this->asDerived();
    }
    /** set other at the beginning of this (concatenate). Perform a copy of the
     *  values stored in other to this.
     *  @param other the array to add
     *  @note the size and the type have to match
     **/
    template<class Other>
    Derived& pushFrontRows(ExprBase<Other> const& other)
    {
      // check if the array is empty
      if (this->empty())
      {
        this->asDerived() = other.asDerived();
        return this->asDerived();
      }
      // not empty
      if (other.cols() != this->cols())
        STKRUNTIME_ERROR_NO_ARG(IArray2D::pushBackRows,range of the columns are different);
      // add nbRow to existing rows
      int nbRow = other.sizeRows();
      pushFrontRows(nbRow);
      for (int j=this->beginCols(); j< this->endCols(); ++j)
      {
        for (int i=this->beginRows(), iOther= other.beginRows(); iOther<other.endRows(); ++i, ++iOther)
        { this->elt(i,j) = other.elt(iOther,j);}
      }
      // return this
      return this->asDerived();
    }
    /** set other at the end of this (concatenate). Perform a copy of the
     *  values stored in other to this.
     *  @param other the array to add back
     *  @note the size and the type have to match
     **/
    template<class Other>
    Derived& pushBackRows(ExprBase<Other> const& other)
    {
      // check if the array is empty
      if (this->empty())
      {
        this->asDerived() = other.asDerived();
        return this->asDerived();
      }
      // not empty
      if (other.cols() != this->cols())
        STKRUNTIME_ERROR_NO_ARG(IArray2D::pushBackRows,range of the columns are different);
      // add nbRow to existing rows
      int nbRow = other.sizeRows();
      pushBackRows(nbRow);
      for (int j=this->beginCols(); j< this->endCols(); ++j)
      {
        // start from the end in order to avoid
        for (int i=this->lastIdxRows(), iOther= other.lastIdxRows(); iOther>=other.beginRows(); --i, --iOther)
        { this->elt(i,j) = other.elt(iOther,j);}
      }
      // return this
      return this->asDerived();
    }
    /** @brief function for reserving memory in all the columns
     *  @param sizeRows the size to reserve for the rows
     *  @param sizeCols the size to reserve for the columns
     **/
    void reserve(int sizeRows, int sizeCols)
    {
      this->reserveCols(sizeCols);
      reserveRows(sizeRows);
    }
    /** @brief function for reserving memory in all the columns
     *  @param size the size to reserve
     **/
    void reserveRows(int size) { reserveRows(this->cols(), size);}

    // for one dimension containers
    /** STL compatibility: insert element @c v in the range @c I of the Array.
     *  @param v,I value and range of indexes
     **/
    void insert( Range const& I, Type const& v)
    { // insert defined for vector, point and diagonal arrays
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);
      this->asDerived().insertElt(I.begin(), I.size());
      for (int i=I.begin(); i<I.end(); i++) this->elt(i) = v;
    }
    /** STL compatibility: push front an element.
     *  @param v value to push front
     **/
    void push_front(Type const& v)
    { // push_front defined for vector, point and diagonal arrays
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);
      insert(Range(this->begin(), 1), v);
    }
    /** STL compatibility: append an element v.
     *  @param v value to append back
     **/
    void push_back(Type const& v)
    { // push_back defined for vector, point and diagonal arrays
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);
      this->asDerived().pushBack();
      this->back() = v;
    }

  protected:
    /** copy forward the column @c srcCol of @c src in the column @c dstCol of this. */
    void copyColumnForward(IArray2D const& src, int jDst, int jSrc)
    {
      Type *dp =this->data(jDst), *sp =src.data(jSrc);
      const int tfirst(this->rangeCols_[jDst].begin());
      const int sfirst(src.rangeCols_[jSrc].begin());
      const int send (src.rangeCols_[jSrc].end());
      for ( int it=tfirst, is=sfirst; is<send; it++, is++) dp[it] = sp[is];
    }
    /** copy backward the column @c jSrc of @c src in the column @c jDst of this. */
    void copyColumnBackward(IArray2D const& src, int jDst, int jSrc)
    {
      Type *dp =this->data(jDst), *sp =src.data(jSrc);
      const int tlast (this->rangeCols_[jDst].lastIdx())
              , sfirst(src.rangeCols_[jSrc].begin())
              , slast (src.rangeCols_[jSrc].lastIdx());

      for ( int it=tlast, is=slast; is>=sfirst; it--, is--) dp[it] = sp[is];
    }
    /** Memory deallocation.
     *  This method clear all allocated memory. The range of the Cols
     *  is set to (beginHo_:beginHo_-1). The range of the Rows remain
     *  unmodified.
     **/
    void freeMem()
    {
      // Nothing to do for reference
      if (this->isRef()) return;
      // free the Rows memory
      this->freeCols(this->cols());
      // liberate horizontally
      this->freeRows();
    }

  private:
    /** @brief internal function for reserving memory in a range of columns.
     *  @param J range of the Columns to initialize
     *  @param size the size to reserve
     **/
    void reserveRows(Range const& J, int size)
    {
      for (int j=J.begin(); j<J.end(); j++)
      {
        try
        {
          this->reserveRowsToCol(j, size);
        }
        catch(Exception const& error)// if an error occur just stop iterations
        {
#ifdef STK_DEBUG
          stk_cout << STKERROR_2ARG2(IArray2D::reserveRows,J, size,memory allocation failed.);
#endif
          // and throw an exception
          STKRUNTIME_ERROR_2ARG(IArray2D::reserveRows,J, size,memory allocation failed.);
        }
      }
    }
    /** @brief main function for memory allocation and initialization of the columns.
     *  The capacity for the Rows have to be set before calling this method.
     *  @param J vertical range of the Columns to initialize
     **/
    void initializeCols(Range const& J)
    {
      for (int j=J.begin(); j<J.end(); j++)
      {
        try
        {
          // initialize the column with the range specific to the array
          this->initializeCol(j, this->rangeRowsInCol(j));
        }
        catch (Exception const& error)   // if an error occur
        {
          // free each column allocated
          for (int k=J.begin(); k<j; k++) this->freeCol(k);
          // put default for the other Cols
          for (int k=j; k<J.end(); k++) this->data(k) = 0;
          // and throw an exception
          throw error;
        }
      }
    }
    /** @brief internal method for initializing a column.
     *
     *  Method for the the allocation of memory of the col
     *  pos with the given range.
     *  @param pos the index of the column to initialize
     *  @param I   range of the Col
     **/
    void initializeCol(int pos, Range const& I)
    {
      if (I.size() <=0)
      {
        // set default for ptr
        this->data(pos) = 0;
        // set default value for this->availableRows_[pos]
        this->availableRows_[pos] = 0;
        // set default value for this->rangeCols_[pos]
        this->rangeCols_[pos] = I;
        // return
        return;
      }
      // compute the size necessary (cannot be 0)
      int size = Arrays::evalSizeCapacity(I.size());
      // try to allocate memory
      try
      {
        this->data(pos) = new Type[size];
      }
      catch (std::bad_alloc & error)  // if an alloc error occur
      {
        // set default for ptr
        this->data(pos) = 0;
        // set default value for this->availableRows_[pos]
        this->availableRows_[pos] = 0;
        // set default value for this->rangeCols_[pos]
        this->rangeCols_[pos] = Range();
        // and throw an exception
        STKRUNTIME_ERROR_2ARG(IArray2D::initializeCol,pos, I,memory allocation failed.);
      }
      // increment ptr of the column
      this->data(pos) -= I.begin();
      // set size for this->availableRows_[pos]
      this->availableRows_[pos] = size;
      // set value for this->rangeCols_[pos]
      this->rangeCols_[pos] = I;
    }
    /** vertical memory deallocation.
     *  @param J range of the columns to liberate.
     **/
    void freeCols(Range const& J)
    { for (int j=J.begin(); j<J.end(); j++) { freeCol(j);}}
    /** @brief Internal method for memory deallocation.
     *  @param col the number of the column to free
     **/
    void freeCol(int col)
    {
      if (this->data(col)) // if there is a column at this position
      {
        // increment the ptr
        this->data(col) += this->rangeCols_[col].begin();
        // delete allocated mem for the column col
        delete [] this->data(col);
        // set default value for ptr
        this->data(col) =0;
        // set default value for this->availableRows_[col]
        this->availableRows_[col] = 0;
        // set default value for this->rangeCols_[col]
        this->rangeCols_[col] = Range();
      }
    }
    /** @brief internal method for translating a column.
     *
     *  Method for the the allocation of memory of the column
     *  pos with the given range.
     *  @param col the index of the column to translate
     *  @param beg new begin of the column
     **/
    void shiftCol( int col, int beg)
    {
      // check if there is data
      if (this->data(col))
      { this->data(col) -= (beg - this->rangeCols_[col].begin());}
      // translate this->rangeCols_
      this->rangeCols_[col].shift(beg);
    }
    /** @brief Internal method for reserving rows to a specified column.
     *
     *  reserve @c size memory place to the column @c col of the array.
     *  @param col index of the column
     *  @param size to reserve
     **/
    void reserveRowsToCol( int col, int size)
    {
      // nothing to do
      if (this->availableRows_[col] > size) return;
      // wrap old Col
      Type* p_oldCol(this->data(col));
      // create new Col
      try
      {
        this->data(col) = new Type[size];
        this->availableRows_[col] = size;
      }
      catch (std::bad_alloc & error)  // if an alloc error occur
      {
        this->data(col) = p_oldCol;
        STKRUNTIME_ERROR_2ARG(IArray2D::reserveRowsToCol,col,size,memory allocation failed.);
      }
      // increment ptr of the column
      this->data(col) -= this->rangeCols_[col].begin();
      // ger ptr on the new col
      Type* p_newCol = this->data(col);
      // copy Elts
      for (int i=this->rangeCols_[col].begin(); i<this->rangeCols_[col].end(); i++)
        p_newCol[i] = p_oldCol[i];
      // if there is allocated memory, liberate it
      if (p_oldCol)
      {
        p_oldCol += this->rangeCols_[col].begin();
        delete [] p_oldCol;
      }
    }
    /** @brief Internal method for resizing a column with a specified range.
     *
     *  This method resize the column @c col to the desired range using:
     * - @c shiftCol
     * - either @c popBackRowsToCol or @c pushBackRowsToCol if needed.
     *  @param col index of the column
     *  @param I range to set to the column
    **/
    void resizeCol( int col, Range const& I)
    {
      // check if there is something to do
      if (this->rangeCol(col) == I) return;
      // shift to the desired first index
      shiftCol(col, I.begin());
      // compute difference of size
      int inc = this->rangeCol(col).size() - I.size();
      // nothing to do
      if (inc == 0) return;
      // add row
      if (inc < 0) { pushBackRowsToCol(col, -inc);}
      else         { popBackRowsToCol(col, inc);}
    }
    /** @brief Internal method for inserting rows to a specified column.
     *
     *  Insert n Rows at the position pos to the column column of the
     *  array. No check is done about the index.
     *  @param col column index
     *  @param pos index where to insert Rows
     *  @param n number of elements to insert (default 1)
     **/
    void insertRowsToCol( int col, int pos, int n =1)
    {
      // wrap old Column
      Type* p_oldCol(this->data(col));
      // get vertical range of the Column
      Range oldRange(this->rangeCols_[col]);
      // update range
      this->rangeCols_[col].incLast(n);
      // allocate if necessary the Col
      if (this->availableRows_[col] < this->rangeCols_[col].size())
      {
        // create new Col
        this->initializeCol(col, this->rangeCols_[col]);
        // if there was data, copy and liberate
        if (p_oldCol)
        {
          // get ptr on the new col
          Type* p_newCol(this->data(col));
          // copy first Elts
          for (int k=oldRange.begin(); k<pos; k++) p_newCol[k] = p_oldCol[k];
          // translate and copy last Elts
          for (int k=oldRange.lastIdx(); k>=pos; k--) p_newCol[k+n] = p_oldCol[k];
          // increment ptr_col
          p_oldCol += oldRange.begin();
          // and free old col
          delete [] p_oldCol;
        }
      }
      else // enough space
      {
        // translate last Elts
        for (int k=oldRange.lastIdx(); k>=pos; k--)
          p_oldCol[k+n] = p_oldCol[k];
      }
    }
    /** @brief Internal method for appending rows to a specified column.
     *
     *  Push back n Rows at the end of the column @c col of the array.
     *  @param col column index
     *  @param n number of elements to append (default 1)
     **/
    void pushBackRowsToCol( int col, int n =1)
    {
      // wrap old Col
      Type* p_oldCol(this->data(col));
      // get vertical range of the Col
      Range oldRange(this->rangeCols_[col]);
      // compute vertical range of the Col after insertion
      this->rangeCols_[col].incLast(n);
      // allocate if necessary the Col
      if (this->availableRows_[col] < this->rangeCols_[col].size())
      {
        // create new Col
        this->initializeCol(col, this->rangeCols_[col]);
        // ger ptr on the new col
        Type* p_newCol(this->data(col));
        // copy Elts
        for (int k=oldRange.begin(); k<oldRange.end(); k++) p_newCol[k] = p_oldCol[k];
        // if there is allocated memory, liberate it
        if (p_oldCol)
        {
          p_oldCol += oldRange.begin();
          delete [] p_oldCol;
        }
      }
    }
    /** @brief Internal method for deleting rows from a specified column.
     *
     *  Delete n Rows at the position @c pos to the column @c col of the array.
     *  No check is done about indexes. It is possible to remove data
     *  outside the range of the column. In this case it is assumed
     *  that the data are known and there was no necessity to store
     *  them inside the array (think to a triangular matrix).
     *
     *  @param col index of the Column
     *  @param pos index where to delete elements
     *  @param n number of elements to delete (default 1)
    **/
    void eraseRowsToCol( int col, int pos, int n=1)
    {
      // check trivial cases
      if (this->rangeCols_[col].lastIdx() < pos) return;
      if (this->rangeCols_[col].begin()> pos+n-1)
      { shiftCol( col, this->rangeCols_[col].begin() - n); return;}
      // find the exisiting rows to delete
      Range rangeDel(pos, n);
      rangeDel.inf(this->rangeCols_[col]);
      if (rangeDel == this->rangeCols_[col]) { freeCol(col); return;}
      // shift data, rangeDel is inside the rang of the column
      Type* p_col(this->data(col));
      for ( int k=rangeDel.begin(), k1=rangeDel.end(); k1<this->rangeCols_[col].end(); k++, k1++)
      {  p_col[k]   = p_col[k1];}
      // update size of the range
      this->rangeCols_[col].decLast(rangeDel.size());
      // and shift if necessary
      if (pos < rangeDel.begin())
      { shiftCol( col, this->rangeCols_[col].begin() - (n-rangeDel.size()));}
    }
    /** @brief Internal method for deleting last rows to a specified column.
     *
     *  Delete the  n latest Rows to the array.
     *
     *  @param col index of the Column
     *  @param n number of elements to delete (default is 1)
    **/
    void popBackRowsToCol( int col, int n=1)
    {
      // check if there is something to do
      if (n <= 0) return;
      // update range
      this->rangeCols_[col].decLast(n);
      // free mem if necessary
      if (this->rangeCols_[col].size()==0) freeCol(col);
    }
};


/** overwrite @c this with @c src.
 *  @note If the size match, @c this is not resized, and in this case,
 *  the method take care of the possibly of overlapping.
 *  @param src the array to copy
 **/
template < class  Derived  >
Derived& IArray2D<Derived>::assign( IArray2D const& src)
{
  // Resize if necessary.
  if ( (this->sizeRows() != src.sizeRows()) ||(this->sizeCols() != src.sizeCols()) )
  { this->resize(src.rows(), src.cols());}
  // Copy without overlapping
  if (src.beginRows()>=this->beginRows())
  {
    if (src.beginCols()>this->beginCols())
    {
      for ( int jSrc=src.beginCols(), jDst=this->beginCols(); jSrc<src.endCols(); jDst++, jSrc++)
      { this->copyColumnForward(src, jDst, jSrc);}
      return this->asDerived();
    }
    for ( int jSrc=src.lastIdxCols(), jDst=this->lastIdxCols(); jSrc>=src.beginCols(); jDst--, jSrc--)
    { this->copyColumnForward(src, jDst, jSrc);}
    return this->asDerived();
  }
  // src.beginRows()<this->beginRows()
  if (src.beginCols()>=this->beginCols())
  {
    for ( int jSrc=src.beginCols(), jDst=this->beginCols(); jSrc<src.endCols(); jDst++, jSrc++)
    { this->copyColumnBackward(src, jDst, jSrc);}
    return this->asDerived();
  }
  // src.beginCols()<this->beginCols()
  for ( int jSrc=src.lastIdxCols(), jDst=this->lastIdxCols(); jSrc>=src.beginCols(); jDst--, jSrc--)
  { this->copyColumnBackward(src, jDst, jSrc);}

  return this->asDerived();
}
} // namespace STK

#endif
// STK_IARRAY2D_H
