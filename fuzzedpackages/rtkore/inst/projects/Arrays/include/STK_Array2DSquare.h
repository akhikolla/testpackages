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
 * Purpose:  Define the Array2DSquare class.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_Array2DSquare.h
  * @brief In this file, we define Array2DSquare class.
 **/

#ifndef STK_ARRAY2DSQUARE_H
#define STK_ARRAY2DSQUARE_H

#include "STK_IArray2D.h"

namespace STK
{
// forward declaration
template< typename Type> class Array2DSquare;
template< typename Type> class Array2D;
template< typename Type> class Array2DPoint;
template< typename Type> class Array2DVector;

typedef Array2DSquare<Real>   ArraySquareX;
typedef Array2DSquare<double> ArraySquareXd;
typedef Array2DSquare<int>    ArraySquareXi;

namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of the Traits class for the Array2DSquare class.
 **/
template<class Type_>
struct Traits<Array2DSquare<Type_> >
{
  private:
    class Void {};
  public:
    typedef Array2DPoint<Type_>  Row;
    typedef Array2DVector<Type_> Col;
    typedef Array2DPoint<Type_>  SubRow;
    typedef Array2DVector<Type_> SubCol;
    typedef Array2D<Type_>       SubArray;
    typedef Void                 SubVector;

    typedef Type_                Type;
    typedef typename RemoveConst<Type>::Type const& ConstReturnType;

    enum
    {
      structure_ = Arrays::square_,
      orient_    = Arrays::by_col_,
      sizeRows_  = UnknownSize,
      sizeCols_  = UnknownSize,
      size_      = UnknownSize,
      storage_   = Arrays::dense_ // always dense
    };
};

} // namespace hidden

/** @ingroup Arrays
  * @brief Derivation of the Array2DSquare class for square arrays of
  * Real.
  *
  * A Array2DSquare is a column oriented two dimensional
  * container of Real with the same number of rows and columns.
  *
  * The range of the rows and the columns is the same.
  **/
template<class Type_>
class Array2DSquare: public IArray2D< Array2DSquare<Type_> >
{
  public:
    /** Type for the Interface Class. */
    typedef IArray2D< Array2DSquare<Type_> > Base;
    typedef ArrayBase< Array2DSquare<Type_> > LowBase;

    typedef typename hidden::Traits< Array2DSquare<Type_> >::Row Row;
    typedef typename hidden::Traits< Array2DSquare<Type_> >::Col Col;
    typedef typename hidden::Traits< Array2DSquare<Type_> >::SubRow SubRow;
    typedef typename hidden::Traits< Array2DSquare<Type_> >::SubCol SubCol;
    typedef typename hidden::Traits< Array2DSquare<Type_> >::SubVector SubVector;
    typedef typename hidden::Traits< Array2DSquare<Type_> >::SubArray SubArray;

    typedef typename hidden::Traits< Array2DSquare<Type_> >::Type Type;
    typedef typename hidden::Traits< Array2DSquare<Type_> >::ConstReturnType ConstReturnType;

    enum
    {
      structure_ = hidden::Traits< Array2DSquare<Type_> >::structure_,
      orient_    = hidden::Traits< Array2DSquare<Type_> >::orient_,
      sizeCols_  = hidden::Traits< Array2DSquare<Type_> >::sizeCols_,
      sizeRows_  = hidden::Traits< Array2DSquare<Type_> >::sizeRows_,
      size_      = hidden::Traits< Array2DSquare<Type_> >::size_,
      storage_   = hidden::Traits< Array2DSquare<Type_> >::storage_
    };

    /** Default constructor with cols_=(1:0) and rows_=(1:0). */
    Array2DSquare(): Base() {}
    /** Default constructor with specified Range
     *  @param I range of the Rows and Cols
     **/
    Array2DSquare( Range const& I): Base(I, I) {}
    /** constructor with cols_and rows_ givens,
     *  initialization with a constant.
     *  @param I range of the Rows and Cols
     *  @param v initial value of the container
     **/
    Array2DSquare( Range const& I, Type const& v): Base(I, I)
    { LowBase::setValue(v);}
    /** Copy constructor.
     *  @param T the container to copy
     *  @param ref true if T is wrapped
     **/
    Array2DSquare( Array2DSquare const&T, bool ref=false): Base(T, ref) {}
    /** constructor by reference, ref_=1 in the range given by I.
     *  @param T the Container to wrap
     *  @param I range of the container to wrap
     **/
    Array2DSquare( Array2DSquare const& T, Range const& I): Base(T, I, I) {}
    /** Copy constructor using an expression.
     *  @param T the container to wrap
     **/
    template<class OtherDerived>
    Array2DSquare( ExprBase<OtherDerived> const& T): Base()
    { LowBase::operator=(T);}
    /** destructor. */
    ~Array2DSquare() {}
    /** New beginning index for the object.
     *  @param beg first index of the container
     **/
    void shift1D(int beg)
    { Base::shift(beg, beg);}
    /** New size for the container.
     *  @param I range of the columns and rows of the container
     **/
    Array2DSquare& resize1D( Range const& I)
    { Base::resize(I, I); return *this;}
    /** Insert n rows and column at the given position to the container.
     *  @param pos,n position and number of rows and columns to insert
     **/
    void insert( int pos, int n =1)
    {
      Base::insertRows(pos, n);
      Base::insertCols(pos, n);
    }
    /** Delete n rows and columns at the specified position to
     *  the container.
     *  @param pos,n position and number of rows and columns to erase
     **/
    void erase( int pos, int n=1)
    {
      Base::eraseCols(pos, n);
      Base::eraseRows(pos, n);
    }
    /** Add n rows and columns to the container.
     *  @param n number of Rows and Cols to add
     **/
    void pushBack(int n=1)
    {
      Base::pushBackRows(n);
      Base::pushBackCols(n);
    }
    /** Delete n rows and columns at the end of the container.
     *  @param n number of Rows and Cols to delete
     **/
    void popBack(int n=1)
    {
      Base::popBackCols(n);
      Base::popBackRows(n);
    }
    /** operator = : overwrite the CArray with the Right hand side T.
     *  @param T the container to copy
     **/
    template<class Rhs>
    Array2DSquare& operator=(ExprBase<Rhs> const& T) { return LowBase::operator=(T);}
    /** overwrite the Array2D with T.
     *  @param T the container to copy
     **/
    Array2DSquare& operator=( Array2DSquare const& T) { return LowBase::assign(T);}
    /** set the container to a constant value.
     *  @param v the value to set
     **/
    Array2DSquare& operator=(Type const& v) { return LowBase::setValue(v);}
};

} // namespace STK


#endif
// STK_ARRAY2DSQUARE_H
