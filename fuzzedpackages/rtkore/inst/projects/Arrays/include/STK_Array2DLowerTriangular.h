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
 * Purpose:  Define the Array2DLowerTriangular class.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_Array2DLowerTriangular.h
  * @brief In this file we define the Array2DLowerTriangular class
 **/

#ifndef STK_MATRIXLOWERTRIANGULAR_H
#define STK_MATRIXLOWERTRIANGULAR_H

#include "STK_IArray2D.h"

namespace STK
{
// forward declaration
template<typename> class Array2DLowerTriangular;
template<typename> class Array2DPoint;
template<typename> class Array2DVector;

/** @ingroup Arrays
  * @brief Specialization of the Array2D class for lower triangular
  * matrices.
  *
  * An Array2DLowerTriangular is a column oriented 2D container of Real
  * whcih is lower triangular.
 **/
typedef Array2DLowerTriangular<Real>   ArrayLowerTriangularXX;
typedef Array2DLowerTriangular<double> ArrayLowerTriangularXXd;
typedef Array2DLowerTriangular<int>    ArrayLowerTriangularXXi;

namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of the Traits class for the Array2DLowerTriangular class.
 **/
template<class Type_>
struct Traits< Array2DLowerTriangular<Type_> >
{  private:
    class Void {};
  public:
    typedef Array2DPoint<Type_>           Row;
    typedef Array2DVector<Type_>          Col;
    typedef Array2DPoint<Type_>           SubRow;
    typedef Array2DVector<Type_>          SubCol;
    typedef Array2DLowerTriangular<Type_> SubArray;
    typedef Void                          SubVector;

    typedef Type_                Type;
    typedef typename RemoveConst<Type>::Type const& ConstReturnType;

   enum
   {
     structure_ = Arrays::lower_triangular_,
     orient_    = Arrays::by_col_,
     sizeRows_  = UnknownSize,
     sizeCols_  = UnknownSize,
     storage_ = Arrays::dense_ // always dense
   };
};

}

/** @ingroup Arrays
  * @brief Declaration of the lower triangular matrix class.
  *
  * A Array2DLowerTriangular is a column oriented 2D lower triangular
  * container of element. It is possible to add/remove rows and columns
  * but in this case the container will no more be triangular.
  * The container can be set lower triangular again using the method
  * IArray2D::update().
  *
  * One can give any range for the rows/columns in the constructor but as the
  * diagonal of the matrix is defined as @em i=j then the lower part of the matrix
  * is the set of index [i,j], with i<j.
 **/
template<class Type_ >
class Array2DLowerTriangular: public IArray2D< Array2DLowerTriangular<Type_> >
{
  public:
    /** Type for the Interface Class.*/
    typedef IArray2D< Array2DLowerTriangular<Type_> > Base;
    typedef ArrayBase < Array2DLowerTriangular<Type_> > LowBase;

    typedef typename hidden::Traits<Array2DLowerTriangular<Type_> >::Row Row;
    typedef typename hidden::Traits<Array2DLowerTriangular<Type_> >::Col Col;
    typedef typename hidden::Traits<Array2DLowerTriangular<Type_> >::SubRow SubRow;
    typedef typename hidden::Traits<Array2DLowerTriangular<Type_> >::SubCol SubCol;
    typedef typename hidden::Traits<Array2DLowerTriangular<Type_> >::SubVector SubVector;
    typedef typename hidden::Traits<Array2DLowerTriangular<Type_> >::SubArray SubArray;

    typedef typename hidden::Traits<Array2DLowerTriangular<Type_> >::Type Type;
    typedef typename hidden::Traits<Array2DLowerTriangular<Type_> >::ConstReturnType ConstReturnType;

   enum
   {
     structure_ = Arrays::lower_triangular_,
     orient_    = Arrays::by_col_,
     sizeRows_  = UnknownSize,
     sizeCols_  = UnknownSize,
     storage_ = Arrays::dense_ // always dense
   };

    /** Default constructor */
    Array2DLowerTriangular(): Base() {}
    /** Constructor with specified ranges
     *  @param I range of the Rows
     *  @param J range of the Cols
     **/
    Array2DLowerTriangular( Range const& I, Range const& J): Base(I, J) {}
    /** constructor with rows_ and rageHo_ specified, initialization with a
     *  specified value.
     *  @param I range of the Rows
     *  @param J range of the Cols
     *  @param v initial value in the container
     **/
    Array2DLowerTriangular( Range const& I, Range const& J, Type const& v)
                         : Base(I, J)
    { LowBase::setValue(v);}
    /** Copy constructor
     *  @param T the container to copy
     *  @param ref true if T is wrapped
     **/
    Array2DLowerTriangular( const Array2DLowerTriangular &T, bool ref=false)
                         : Base(T, ref) {}
    /** constructor by reference, ref_=1.
     *  @param T the container to wrap
     *  @param I range of the Rows to wrap
     *  @param J range of the Cols to wrap
     **/
    Array2DLowerTriangular( const Base& T, Range const& I, Range const& J)
                         : Base(T, I, J) {}
    /** Wrapper constructor Contruct a reference container.
     *  @param q pointer on the data
     *  @param I range of the  Rows to wrap
     *  @param J range of the Cols to wrap
     **/
    Array2DLowerTriangular( Type** q, Range const& I, Range const& J)
                         : Base(q, I, J) {}
    /** Copy constructor using an expression.
     *  @param T the container to wrap
     **/
    template<class OtherDerived>
    Array2DLowerTriangular( ExprBase<OtherDerived> const& T): Base()
    { LowBase::operator=(T);}
    /** destructor : use destructor of IArray2D. */
    ~Array2DLowerTriangular() {}
    /** operator = : overwrite the CArray with the Right hand side T.
     *  @param T the container to copy
     **/
    template<class Rhs>
    Array2DLowerTriangular& operator=(ExprBase<Rhs> const& T)
    { return LowBase::operator=(T);}
    /** operator = : overwrite the Array2DLowerTriangular with T.
     *  @param T the container to copy
     **/
    Array2DLowerTriangular& operator=( Array2DLowerTriangular const& T) { return LowBase::assign(T);}
    /** operator= : set the container to a constant value.
     *  @param v the value to set
     **/
    Array2DLowerTriangular& operator=(Type const& v){ return LowBase::setValue(v);}
    /** New beginning index for the object.
     *  @param beg first index of the container
     **/
    void shift1D(int beg)
    { Base::shift(beg, beg);}
    /** New size for the container.
     *  @param I range of the columns and rows of the container
     **/
    Array2DLowerTriangular& resize1D( Range const& I)
    { Base::resize(I, I); return *this;}
};

} // namespace STK

#endif
// STK_MATRIXLOWERTRIANGULAR_H
