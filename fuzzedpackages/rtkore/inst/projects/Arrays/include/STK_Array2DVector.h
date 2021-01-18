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
 * Project: stkpp::Arrays
 * Purpose:  Define the Array2DVector class.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_Array2DVector.h
  * @brief A Array2DVector is a one dimensional horizontal container
  *
  * An Array2DVector is an implementation of the interface IArray2D.
  * It's a one dimensional horizontal container.
 **/

#ifndef STK_ARRAY2DVECTOR_H
#define STK_ARRAY2DVECTOR_H

#include "STK_IArray2D.h"

namespace STK
{

template<typename> class Array2DPoint;
template<typename> class Array2DVector;


/** @ingroup Arrays
  * @brief final class for a Real vertical container.
  *
  * A Vector is a column oriented 1D container of Real.
  **/
typedef Array2DVector<Real>   Vector; // [DEPRECATED]
typedef Array2DVector<Real>   VectorX;
typedef Array2DVector<double> VectorXd;
typedef Array2DVector<int>    VectorXi;

namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of the Traits class for the Array2DVector class.
 **/
template<class Type_>
struct Traits< Array2DVector<Type_> >
{
  typedef Array2DPoint<Type_>  Row;
  typedef Array2DVector<Type_> Col;
  typedef Array2DPoint<Type_>  SubRow;
  typedef Array2DVector<Type_> SubCol;
  typedef Array2DVector<Type_> SubArray;
  typedef Array2DVector<Type_> SubVector;

  typedef Type_                Type;
  typedef typename RemoveConst<Type>::Type const& ConstReturnType;

  enum
  {
    structure_ = Arrays::vector_,
    orient_    = Arrays::by_col_,
    sizeCols_  = 1,
    sizeRows_  = UnknownSize,
    size_      = UnknownSize,
    storage_ = Arrays::dense_ // always dense
  };
};

} // namespace hidden

/** @ingroup Arrays
 *  @brief template one dimensional horizontal Array.
 *
 *  An Array2DVector is a Vertical container of a single column.
 *
 *  By default the index of the first element is 0 but this can be
 *  modified using the appropriate constructor or using the method @c shift.
 *
 *  @sa Array2DPoint
 **/
template<class Type_>
class Array2DVector: public IArray2D< Array2DVector<Type_> >
{
  public:
    typedef IArray2D< Array2DVector<Type_> > Base;
    typedef ArrayBase< Array2DVector<Type_> > LowBase;

    typedef typename hidden::Traits<Array2DVector<Type_> >::Row Row;
    typedef typename hidden::Traits<Array2DVector<Type_> >::Col Col;
    typedef typename hidden::Traits<Array2DVector<Type_> >::SubRow SubRow;
    typedef typename hidden::Traits<Array2DVector<Type_> >::SubCol SubCol;
    typedef typename hidden::Traits<Array2DVector<Type_> >::SubVector SubVector;
    typedef typename hidden::Traits<Array2DVector<Type_> >::SubArray SubArray;

    typedef typename hidden::Traits<Array2DVector<Type_> >::Type Type;
    typedef typename hidden::Traits<Array2DVector<Type_> >::ConstReturnType ConstReturnType;

    enum
    {
      structure_ = hidden::Traits< Array2DVector<Type_> >::structure_,
      orient_    = hidden::Traits< Array2DVector<Type_> >::orient_,
      sizeCols_  = hidden::Traits< Array2DVector<Type_> >::sizeCols_,
      sizeRows_  = hidden::Traits< Array2DVector<Type_> >::sizeRows_,
      size_      = hidden::Traits< Array2DVector<Type_> >::size_,
      storage_   = hidden::Traits< Array2DVector<Type_> >::storage_
    };

    /** Default constructor */
    Array2DVector(): Base( Range(), Range(1)) {}
    /** constructor with specified range.
     *  @param I range of the container
     **/
    Array2DVector( Range const& I) :Base(I, Range(1)) {}
    /** constructor with specified range, initialization with a constant.
     *  @param I range of the container
     *  @param v initial value of the container
     **/
    Array2DVector( Range const& I, Type const& v): Base(I, Range(1))
    { LowBase::setValue(v);}
    /** Copy constructor
     *  @param T the container to copy
     *  @param ref true if this is a wrapper of T
     **/
    Array2DVector( const Array2DVector &T, bool ref =false)
                : Base(T, ref) {}
    /** constructor by reference, ref_=1.
     *  @param T the container to wrap
     *  @param I the columns range to wrap
     **/
    Array2DVector( const Array2DVector& T, Range const& I)
                : Base(T, I, T.cols())
    {}
    /** constructor by reference, ref_=1.
     *  @param T the container to wrap
     *  @param I the range of the data to wrap
     *  @param col the index of the column to wrap
     **/
    template<class OtherArray>
    Array2DVector( IArray2D<OtherArray> const& T, Range const& I, int col)
               : Base(T, I, Range(col, 1))
    {}
    /** Copy constructor using an expression.
     *  @param T the container to wrap
     **/
    template<class OtherDerived>
    Array2DVector( ExprBase<OtherDerived> const& T): Base( Range(), Range(1))
    { LowBase::operator=(T);}
    /** constructor by reference, ref_=1.
     *  @param p_data a pointer on the data to wrap
     *  @param I the range of the data to wrap
     *  @param col the index of the column to wrap
     **/
     Array2DVector( Type** p_data, Range const& I, int col)
                 : Base(p_data, I, Range(col, 1))
    {}
    /** destructor. */
    ~Array2DVector() {}
    /** @return a constant reference on the ith element
     *  @param i index of the element (const)
     **/
    Type const & elt1Impl( int i) const { return this->data(this->beginCols())[i];}
    /** @return a reference on the ith element
     *  @param i index of the element
     **/
    inline Type& elt1Impl( int i) { return this->data(this->beginCols())[i];}
    /** New first index for the object.
     *  @param rbeg the index of the first row to set
     **/
    void shift1D( int rbeg) { Base::shift(rbeg, this->beginCols());}
    /**  Resize the container.
     *  @param I the range to set to the container
     **/
    Array2DVector<Type>& resize1D( Range const& I)
    { Base::resize(I, this->cols()); return *this;}
    /** Add n elements to the container.
     *  @param n number of elements to add
     **/
    void pushBack( int n=1) { Base::pushBackRows(n);}
    /** Delete n elements at the pos index to the container.
     *  @param pos index where to delete elements
     *  @param n number of elements to delete (default 1)
    **/
    void erase( int pos, int const& n=1) { Base::eraseRows(pos, n);}
    /** Insert n elements at the position pos of the container. The bound
     *  end_ should be modified at the very end of the insertion as pos
     *  can be a reference to it.
     *  @param pos index where to insert elements
     *  @param n number of elements to insert (default 1)
     **/
    void insertElt(int pos, int const& n =1)
    { Base::insertRows(pos, n);}
    /** operator = : overwrite the CArray with the Right hand side T.
     *  @param T the container to copy
     **/
    template<class Rhs>
    Array2DVector& operator=(ExprBase<Rhs> const& T) { return LowBase::operator=(T);}
    /** operator = : overwrite the CArray with the Right hand side T.
     *  @param T the container to copy
     **/
    Array2DVector& operator=(Array2DVector const& T) { return LowBase::assign(T);}
    /** set the container to a constant value.
     *  @param v the value to set
     **/
    Array2DVector& operator=(Type const& v) { return LowBase::setValue(v);}
};

} // namespace STK

#endif // STK_ARRAY2DVECTOR_H
