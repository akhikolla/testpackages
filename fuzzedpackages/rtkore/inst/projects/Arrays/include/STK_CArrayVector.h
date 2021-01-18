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
 * created on: 25 nov. 2011
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_CArrayVector.h
 *  @brief In this file we implement the final class CArrayVector.
 **/

#ifndef STK_CARRAYVECTOR_H
#define STK_CARRAYVECTOR_H

#include "STK_ICArray.h"

namespace STK
{
// forward declaration
template< typename Type, int SizeRows_=UnknownSize, bool Orient_ = Arrays::by_col_> class CArrayVector;

template< typename Type, int SizeRows_, int SizeCols_, bool Orient_> class CArray;
template< typename Type, int Size_, bool Orient_> class CArraySquare;
template< typename Type, int SizeCols_, bool Orient_> class CArrayPoint;
template< typename Type, bool Orient_> class CArrayNumber;

// useful typedef
typedef CArrayVector<Real, UnknownSize, Arrays::by_col_>   CVectorX;
typedef CArrayVector<Real, 2, Arrays::by_col_>             CVector2;
typedef CArrayVector<Real, 3, Arrays::by_col_>             CVector3;
typedef CArrayVector<double, UnknownSize, Arrays::by_col_> CVectorXd;
typedef CArrayVector<double, 2, Arrays::by_col_>           CVector2d;
typedef CArrayVector<double, 3, Arrays::by_col_>           CVector3d;
typedef CArrayVector<int, UnknownSize, Arrays::by_col_>    CVectorXi;
typedef CArrayVector<int, 2, Arrays::by_col_>              CVector2i;
typedef CArrayVector<int, 3, Arrays::by_col_>              CVector3i;

typedef CArrayVector<Real, UnknownSize, Arrays::by_row_>   CVectorByRowX;
typedef CArrayVector<Real, 2, Arrays::by_row_>             CVectorByRow2;
typedef CArrayVector<Real, 3, Arrays::by_row_>             CVectorByRow3;
typedef CArrayVector<double, UnknownSize, Arrays::by_row_> CVectorByRowXd;
typedef CArrayVector<double, 2, Arrays::by_row_>           CVectorByRow2d;
typedef CArrayVector<double, 3, Arrays::by_row_>           CVectorByRow3d;
typedef CArrayVector<int, UnknownSize, Arrays::by_row_>    CVectorByRowXi;
typedef CArrayVector<int, 2, Arrays::by_row_>              CVectorByRow2i;
typedef CArrayVector<int, 3, Arrays::by_row_>              CVectorByRow3i;

namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of the Traits class for CArrayVector class.
 */
template<typename Type_, int SizeRows_, bool Orient_>
struct Traits< CArrayVector<Type_, SizeRows_, Orient_> >
{
    typedef CArrayNumber<Type_, Orient_> Row;
    typedef CArrayVector<Type_, SizeRows_, Orient_> Col;

    // The CAllocator have to have the same structure than the CArray
    typedef CAllocator<Type_, SizeRows_, 1, Orient_> Allocator;

    typedef Type_                Type;
    typedef typename RemoveConst<Type>::Type const& ConstReturnType;

    enum
    {
      structure_ = Arrays::vector_,
      orient_    = Orient_,
      sizeRows_  = SizeRows_,
      sizeCols_  = 1,
      size_      = SizeRows_,
      storage_   = Arrays::dense_
    };
};


} // namespace hidden

/** @ingroup Arrays
 * @brief specialization for the vector case.
 */
template <typename Type_, int SizeRows_, bool Orient_>
class CArrayVector: public ICArray < CArrayVector<Type_, SizeRows_, Orient_> >
{
  public:
    typedef ICArray < CArrayVector<Type_, SizeRows_, Orient_> > Base;
    typedef ArrayBase < CArrayVector<Type_, SizeRows_, Orient_> > LowBase;

    typedef typename hidden::Traits< CArrayVector<Type_, SizeRows_, Orient_> >::Row Row;
    typedef typename hidden::Traits< CArrayVector<Type_, SizeRows_, Orient_> >::Col Col;
    typedef typename hidden::Traits< CArrayVector<Type_, SizeRows_, Orient_> >::Type Type;
    typedef typename hidden::Traits< CArrayVector<Type_, SizeRows_, Orient_> >::ConstReturnType ConstReturnType;

    enum
    {
      structure_ = hidden::Traits< CArrayVector<Type_, SizeRows_, Orient_> >::structure_,
      orient_    = hidden::Traits< CArrayVector<Type_, SizeRows_, Orient_> >::orient_,
      sizeRows_  = hidden::Traits< CArrayVector<Type_, SizeRows_, Orient_> >::sizeRows_,
      sizeCols_  = hidden::Traits< CArrayVector<Type_, SizeRows_, Orient_> >::sizeCols_,
      size_      = hidden::Traits< CArrayVector<Type_, SizeRows_, Orient_> >::size_,
      storage_   = hidden::Traits< CArrayVector<Type_, SizeRows_, Orient_> >::storage_
    };

    /** Default constructor. */
    CArrayVector(): Base() {}
    /** constructor with specified dimension.
     *  @param sizeRows size of the rows
     **/
    CArrayVector( int sizeRows): Base(sizeRows, 1) {}
    /** constructor with specified ranges.
     *  @param range range of the rows
     **/
    CArrayVector( Range range): Base(range.size(), 1)
    { this->shift(range.begin());}
    /** constructor with sizeRows specified, initialization with a constant.
     *  @param sizeRows size of the rows
     *  @param v initial value of the container
     **/
    CArrayVector( int sizeRows, Type const& v): Base(sizeRows, 1, v) {}
    /** constructor with range specified, initialization with a constant.
     *  @param range range of the rows
     *  @param v initial value of the container
     **/
    CArrayVector( Range range, Type const& v): Base(range.size(), 1, v)
    { this->shift(range.begin());}
    /** Copy constructor
     *  @param T the container to copy
     *  @param ref true if T is wrapped
     **/
    CArrayVector( const CArrayVector &T, bool ref=false): Base(T, ref) {}
    /** constructor by reference, ref_=1.
     *  @param T the container to wrap
     *  @param I the columns range to wrap
     **/
    CArrayVector( CArrayVector const& T, Range const& I)
               : Base(T.allocator(), I, T.cols())
    {}
    /** constructor by reference, ref_=1.
     *  @param T the container to wrap
     *  @param I the range of the data to wrap
     *  @param col the index of the col to wrap
     **/
    template<class OtherArray>
    CArrayVector( ICArray<OtherArray> const& T, Range const& I, int col)
               : Base(T.allocator(), I, Range(col, 1))
    {}
    /** wrapper constructor for 0 based C-Array.
     *  @param q pointer on the array
     *  @param nbRow number of rows
     **/
    CArrayVector( Type* const& q, int nbRow): Base(q, nbRow, 1) {}
    /** constructor by reference.
     *  @param allocator the allocator to wrap
     **/
    template<class OtherAllocator>
    CArrayVector( ITContainer2D<OtherAllocator> const& allocator): Base(allocator.asDerived()) {}
    /** Copy constructor using an expression.
     *  @param T the container to wrap
     **/
    template<class OtherDerived>
    CArrayVector( ExprBase<OtherDerived> const& T): Base(T.size(), 1) { LowBase::operator=(T);}
    /** destructor. */
    ~CArrayVector() {}
    /** operator= : set the container to a constant value.
     *  @param v the value to set
     **/
    CArrayVector& operator=(Type const& v) { return LowBase::setValue(v);}
    /** operator = : overwrite the CArrayPoint with the Right hand side T.
     *  @param T the container to copy
     **/
    template<class Rhs>
    CArrayVector& operator=( ExprBase<Rhs> const& T) { return LowBase::assign(T);}
    /** operator = : overwrite the CArray with the Right hand side rhs.
     *  @param rhs the container to copy
     **/
    CArrayVector& operator=(CArrayVector const& rhs) { return LowBase::assign(rhs);}
};

} // namespace STK


#endif /* STK_CARRAYVECTOR_H */
