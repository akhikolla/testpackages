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
 * Purpose:  Define the Array2DNumber class.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_Array2DNumber.h
  * @brief A Array2DNumber is a one dimensional horizontal container
 **/

#ifndef STK_ARRAY2DNUMBER_H
#define STK_ARRAY2DNUMBER_H

#include "STK_IArray2D.h"

namespace STK
{

template<typename> class Array2DNumber;

/** @ingroup Arrays
  * @brief final class for a Real horizontal container.
  * A Number is a row oriented 1D container of Real.
  */
typedef Array2DNumber<Real>   Number;
typedef Array2DNumber<Real>   NumberX;
typedef Array2DNumber<double> NumberXd;
typedef Array2DNumber<int>    NumberXi;

namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of the Traits class for Array2DNumber class.
 **/
template<class Type_>
struct Traits< Array2DNumber<Type_> >
{
  typedef Array2DNumber<Type_>  Row;
  typedef Array2DNumber<Type_>  Col;
  typedef Array2DNumber<Type_>  SubRow;
  typedef Array2DNumber<Type_>  SubCol;
  typedef Array2DNumber<Type_>  SubArray;
  typedef Array2DNumber<Type_>  SubVector;

  typedef Type_                Type;
  typedef typename RemoveConst<Type>::Type const& ConstReturnType;

  enum
  {
    structure_ = Arrays::number_,
    orient_    = Arrays::by_col_,
    sizeRows_  = 1,
    sizeCols_  = 1,
    size_      = 1,
    storage_   = Arrays::dense_ // always dense
  };
};

} // namespace hidden

/** @ingroup Arrays
 *  @brief template number Array.
 *
 * An Array2DNumber is an implementation of the interface IArray2D.
 * It's a number container and is referred as a number.
 *
 *  By default the index of the first element is 1 but this can be
 *  modified using the appropriate constructor or using the method @c shift.
 **/
template<class Type_>
class Array2DNumber: public IArray2D< Array2DNumber<Type_> >
{
  public:
    typedef IArray2D< Array2DNumber<Type_> > Base;
    typedef ArrayBase< Array2DNumber<Type_> > LowBase;

    typedef typename hidden::Traits< Array2DNumber<Type_> >::Row Row;
    typedef typename hidden::Traits< Array2DNumber<Type_> >::Col Col;
    typedef typename hidden::Traits< Array2DNumber<Type_> >::SubRow SubRow;
    typedef typename hidden::Traits< Array2DNumber<Type_> >::SubCol SubCol;
    typedef typename hidden::Traits< Array2DNumber<Type_> >::SubVector SubVector;
    typedef typename hidden::Traits< Array2DNumber<Type_> >::SubArray SubArray;

    typedef typename hidden::Traits< Array2DNumber<Type_> >::Type Type;
    typedef typename hidden::Traits< Array2DNumber<Type_> >::ConstReturnType ConstReturnType;

    enum
    {
       structure_ = hidden::Traits< Array2DNumber >::structure_,
       orient_    = hidden::Traits< Array2DNumber >::orient_,
       sizeRows_  = hidden::Traits< Array2DNumber >::sizeRows_,
       sizeCols_  = hidden::Traits< Array2DNumber >::sizeCols_,
       size_      = hidden::Traits< Array2DNumber >::size_,
       storage_   = hidden::Traits< Array2DNumber >::storage_
    };

    /** Default constructor */
    Array2DNumber(): Base(Range(1), Range(1)) {}
    /** constructor with specified range, initialization with a constant.
     *  @param v initial value of the container
     **/
    Array2DNumber( Type const& v): Base(Range(1), Range(1))
    { LowBase::setValue(v);}
    /** Copy constructor
     *  @param T the container to copy
     *  @param ref true if this is a wrapper of T
     **/
    Array2DNumber( Array2DNumber const& T, bool ref =false): Base(T, ref) {}
    /** constructor by reference, ref_=1.
     *  @param T the container to wrap
     *  @param row, col row and column indexes to wrap
     **/
    template<class OtherArray>
    Array2DNumber( IArray2D<OtherArray> const& T, int row, int col)
               : Base(T, Range(row, 1), Range(col, 1))
    {}
    /** Copy constructor using an expression.
     *  @param T the container to wrap
     **/
    template<class OtherDerived>
    Array2DNumber( ExprBase<OtherDerived> const& T): Base(Range(1), Range(1))
    { LowBase::operator=(T);}
    /** constructor by reference, ref_=1.
     *  @param p_data a pointer on the data to wrap
     *  @param row the index of the row to wrap
     *  @param col the index of the column to wrap
     **/
     Array2DNumber( Type** p_data, int row, int col)
                 : Base(p_data, Range(row, 1), Range(col, 1)) {}
    /** destructor. */
    ~Array2DNumber() {}
    /** @return a constant reference on the jth element
     **/
    inline Type const& elt0Impl() const { return this->data(this->beginCols())[this->beginRows()];}
    /** @return a reference on the jth element
     **/
    inline Type& elt0Impl() { return this->data(this->beginCols())[this->beginRows()];}
    /** New first indexes for the object.
     *  @param cbeg the index of the first column to set
     **/
    void shift1D(int const& cbeg) { Base::shift(this->beginRows(), cbeg);}
    /** operator = : overwrite the CArray with the Right hand side T.
     *  @param T the container to copy
     **/
    template<class Rhs>
    Array2DNumber& operator=(ExprBase<Rhs> const& T) { return LowBase::operator=(T);}
    /** operator = : overwrite the CArray with the Right hand side T.
     *  @param T the container to copy
     **/
    Array2DNumber& operator=(const Array2DNumber &T) { return LowBase::assign(T);}
    /** set the container to a constant value.
     *  @param v the value to set
     **/
    Array2DNumber& operator=(Type const& v) { return LowBase::setValue(v);}
};

} // namespace STK

#endif // STK_ARRAY2DNUMBER_H
