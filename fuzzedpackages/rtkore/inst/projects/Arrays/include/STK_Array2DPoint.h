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
 * Purpose:  Define the Array2DPoint class.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_Array2DPoint.h
  * @brief A Array2DPoint is a one dimensional horizontal container
 **/

#ifndef STK_ARRAY2DPOINT_H
#define STK_ARRAY2DPOINT_H

#include "STK_IArray2D.h"

namespace STK
{

template<typename> class Array2DPoint;
template<typename> class Array2DVector;

/** @ingroup Arrays
  * @brief final class for a Real horizontal container.
  * A Point is a row oriented 1D container of Real.
  */
typedef Array2DPoint<Real>   PointX;
typedef Array2DPoint<double> PointXd;
typedef Array2DPoint<int>    PointXi;

namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of the Traits class for Array2DPoint class.
 **/
template<class Type_>
struct Traits< Array2DPoint<Type_> >
{
  typedef Array2DPoint<Type_>  Row;
  typedef Array2DVector<Type_> Col;
  typedef Array2DPoint<Type_>  SubRow;
  typedef Array2DVector<Type_> SubCol;
  typedef Array2DPoint<Type_>  SubArray;
  typedef Array2DPoint<Type_>  SubVector;

  typedef Type_                Type;
  typedef typename RemoveConst<Type>::Type const& ConstReturnType;

  enum
  {
    structure_ = Arrays::point_,
    orient_    = Arrays::by_col_,
    sizeRows_  = 1,
    sizeCols_  = UnknownSize,
    size_      = UnknownSize,
    storage_   = Arrays::dense_ // always dense
  };
};

} // namespace hidden

/** @ingroup Arrays
 *  @brief template one dimensional horizontal Array.
 *
 * An Array2DPoint is an implementation of the interface IArray2D.
 * It's a one dimensional row-vector and is refered as a point.
 *
 *  By default the index of the first element is 0 but this can be
 *  modified using the appropriate constructor or using the method @c shift.
 **/
template<class Type_>
class Array2DPoint: public IArray2D< Array2DPoint<Type_> >
{
  public:
    typedef IArray2D< Array2DPoint<Type_> > Base;
    typedef ArrayBase< Array2DPoint<Type_> > LowBase;

    typedef typename hidden::Traits< Array2DPoint<Type_> >::Row Row;
    typedef typename hidden::Traits< Array2DPoint<Type_> >::Col Col;
    typedef typename hidden::Traits< Array2DPoint<Type_> >::SubRow SubRow;
    typedef typename hidden::Traits< Array2DPoint<Type_> >::SubCol SubCol;
    typedef typename hidden::Traits< Array2DPoint<Type_> >::SubVector SubVector;
    typedef typename hidden::Traits< Array2DPoint<Type_> >::SubArray SubArray;

    typedef typename hidden::Traits< Array2DPoint<Type_> >::Type Type;
    typedef typename hidden::Traits< Array2DPoint<Type_> >::ConstReturnType ConstReturnType;

    enum
    {
      structure_ = hidden::Traits< Array2DPoint<Type_> >::structure_,
      orient_    = hidden::Traits< Array2DPoint<Type_> >::orient_,
      sizeCols_  = hidden::Traits< Array2DPoint<Type_> >::sizeCols_,
      sizeRows_  = hidden::Traits< Array2DPoint<Type_> >::sizeRows_,
      size_      = hidden::Traits< Array2DPoint<Type_> >::size_,
      storage_   = hidden::Traits< Array2DPoint<Type_> >::storage_
    };

    /** Default constructor */
    Array2DPoint(): Base(Range(1), Range()) {}
    /** constructor with specified range.
     *  @param J range of the container
     **/
    Array2DPoint( Range const& J): Base(Range(1), J) {}
    /** constructor with specified range, initialization with a constant.
     *  @param J range of the container
     *  @param v initial value of the container
     **/
    Array2DPoint( Range const& J, Type const& v): Base(Range(1), J)
    { LowBase::setValue(v);}
    /** Copy constructor
     *  @param T the container to copy
     *  @param ref true if this is a wrapper of T
     **/
    Array2DPoint( Array2DPoint const& T, bool ref =false)
               : Base(T, ref)
    {}
    /** constructor by reference, ref_=1.
     *  @param T the container to wrap
     *  @param J the columns range to wrap
     **/
    Array2DPoint( Array2DPoint const& T, Range const& J)
               : Base(T, T.rows(), J) {}
    /** constructor by reference, ref_=1.
     *  @param T the container to wrap
     *  @param J the range of the data to wrap
     *  @param row the index of the row to wrap
     **/
    template<class OtherArray>
    Array2DPoint( IArray2D<OtherArray> const& T, Range const& J, int row)
               : Base(T, Range(row, 1), J)
    {}
    /** Copy constructor using an expression.
     *  @param T the container to wrap
     **/
    template<class OtherDerived>
    Array2DPoint( ExprBase<OtherDerived> const& T): Base(Range(1), Range())
    { LowBase::operator=(T);}
    /** constructor by reference, ref_=1.
     *  @param p_data a pointer on the data to wrap
     *  @param J the range of the data to wrap
     *  @param row the index of the row to wrap
     **/
     Array2DPoint( Type** p_data, Range const& J, int row)
                : Base(p_data, Range(row, 1), J) {}
    /** destructor. */
    ~Array2DPoint() {}
    /** @return a constant reference on the jth element
     *  @param j index of the element (const)
     **/
    Type const & elt1Impl(int const& j) const { return this->data(j)[this->beginRows()];}
    /** @return a reference on the jth element
     *  @param j index of the element
     **/
    inline Type& elt1Impl(int const& j) { return this->data(j)[this->beginRows()];}
    /** New first indexes for the object.
     *  @param cbeg the index of the first column to set
     **/
    void shift1D(int const& cbeg) { Base::shift(this->beginRows(), cbeg);}
    /**  Resize the container.
     *  @param J the range to set to the container
     **/
    Array2DPoint<Type>& resize1D(Range const& J)
    { Base::resize(this->rows(), J); return *this;}
    /** Add n elements to the container.
     *  @param n number of elements to add
     **/
    void pushBack( int n=1) { Base::pushBackCols(n);}
    /** Delete n elts at the pos index to the container.
     *  @param pos, n position and number of elements to delete
    **/
    void erase( int pos, int n=1)
    { Base::eraseCols(pos, n);}
    /** Insert n elts at the position pos of the container. The bound
     *  end_ should be modified at the very end of the insertion as pos
     *  can be a reference to it.
     *  @param pos, n position and number of elements to insert
     **/
    void insertElt( int pos, int n =1)
    { Base::insertCols(pos, n);}
    /** operator = : overwrite the CArray with the Right hand side T.
     *  @param T the container to copy
     **/
    template<class Rhs>
    Array2DPoint& operator=(ExprBase<Rhs> const& T) { return LowBase::operator=(T);}
    /** operator = : overwrite the CArray with the Right hand side T.
     *  @param T the container to copy
     **/
    Array2DPoint& operator=(const Array2DPoint &T) { return LowBase::assign(T);}
    /** set the container to a constant value.
     *  @param v the value to set
     **/
    Array2DPoint& operator=(Type const& v) { return LowBase::setValue(v);}
};

} // namespace STK

#endif // STK_ARRAY2DPOINT_H
