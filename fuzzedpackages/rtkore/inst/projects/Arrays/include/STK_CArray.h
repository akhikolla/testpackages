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

/** @file STK_CArray.h
 *  @brief In this file we implement the final class CArray.
 **/

#ifndef STK_CARRAY_H
#define STK_CARRAY_H

#include "STK_ICArray.h"

namespace STK
{
// forward declarations
template< typename Type, int SizeRows_ = UnknownSize, int SizeCols_ = UnknownSize, bool Orient_ = Arrays::by_col_> class CArray;
template< typename Type, int Size_, bool Orient_> class CArraySquare;
template< typename Type, int SizeCols_, bool Orient_> class CArrayPoint;
template< typename Type, int SizeRows_, bool Orient_> class CArrayVector;
template< typename Type, bool Orient_> class CArrayNumber;

// useful typedef
typedef CArray<Real, UnknownSize, UnknownSize, Arrays::by_col_> CArrayXX;
typedef CArray<Real, UnknownSize, 2, Arrays::by_col_>           CArrayX2;
typedef CArray<Real, UnknownSize, 3, Arrays::by_col_>           CArrayX3;
typedef CArray<Real, 2, UnknownSize, Arrays::by_col_>           CArray2X;
typedef CArray<Real, 3, UnknownSize, Arrays::by_col_>           CArray3X;
typedef CArray<Real, 2, 2, Arrays::by_col_>                     CArray22;
typedef CArray<Real, 3, 3, Arrays::by_col_>                     CArray33;
typedef CArray<double, UnknownSize, UnknownSize, Arrays::by_col_>CArrayXXd;
typedef CArray<double, UnknownSize, 2, Arrays::by_col_>         CArrayX2d;
typedef CArray<double, UnknownSize, 3, Arrays::by_col_>         CArrayX3d;
typedef CArray<double, 2, UnknownSize, Arrays::by_col_>         CArray2Xd;
typedef CArray<double, 3, UnknownSize, Arrays::by_col_>         CArray3Xd;
typedef CArray<double, 2, 2, Arrays::by_col_>                   CArray22d;
typedef CArray<double, 3, 3, Arrays::by_col_>                   CArray33d;
typedef CArray<int, UnknownSize, UnknownSize, Arrays::by_col_>  CArrayXXi;
typedef CArray<int, UnknownSize, 2, Arrays::by_col_>            CArrayX2i;
typedef CArray<int, UnknownSize, 3, Arrays::by_col_>            CArrayX3i;
typedef CArray<int, 2, UnknownSize, Arrays::by_col_>            CArray2Xi;
typedef CArray<int, 3, UnknownSize, Arrays::by_col_>            CArray3Xi;
typedef CArray<int, 2, 2, Arrays::by_col_>                      CArray22i;
typedef CArray<int, 3, 3, Arrays::by_col_>                      CArray33i;

typedef CArray<Real, UnknownSize, UnknownSize, Arrays::by_row_> CArrayByRowXX;
typedef CArray<Real, UnknownSize, 2, Arrays::by_row_>           CArrayByRowX2;
typedef CArray<Real, UnknownSize, 3, Arrays::by_row_>           CArrayByRowX3;
typedef CArray<Real, 2, UnknownSize, Arrays::by_row_>           CArrayByRow2X;
typedef CArray<Real, 3, UnknownSize, Arrays::by_row_>           CArrayByRow3X;
typedef CArray<Real, 2, 2, Arrays::by_row_>                     CArrayByRow22;
typedef CArray<Real, 3, 3, Arrays::by_row_>                     CArrayByRow33;
typedef CArray<double, UnknownSize, UnknownSize, Arrays::by_row_> CArrayByRowXXd;
typedef CArray<double, UnknownSize, 2, Arrays::by_row_>         CArrayByRowX2d;
typedef CArray<double, UnknownSize, 3, Arrays::by_row_>         CArrayByRowX3d;
typedef CArray<double, 2, UnknownSize, Arrays::by_row_>         CArrayByRow2Xd;
typedef CArray<double, 3, UnknownSize, Arrays::by_row_>         CArrayByRow3Xd;
typedef CArray<double, 2, 2, Arrays::by_row_>                   CArrayByRow22d;
typedef CArray<double, 3, 3, Arrays::by_row_>                   CArrayByRow33d;
typedef CArray<int, UnknownSize, UnknownSize, Arrays::by_row_>  CArrayByRowXXi;
typedef CArray<int, UnknownSize, 2, Arrays::by_row_>            CArrayByRowX2i;
typedef CArray<int, UnknownSize, 3, Arrays::by_row_>            CArrayByRowX3i;
typedef CArray<int, 2, UnknownSize, Arrays::by_row_>            CArrayByRow2Xi;
typedef CArray<int, 3, UnknownSize, Arrays::by_row_>            CArrayByRow3Xi;
typedef CArray<int, 2, 2, Arrays::by_row_>                      CArrayByRow22i;
typedef CArray<int, 3, 3, Arrays::by_row_>                      CArrayByRow33i;

namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of the Traits class for CArray class.
 */
template<typename Type_, int SizeRows_, int SizeCols_, bool Orient_>
struct Traits< CArray<Type_, SizeRows_, SizeCols_, Orient_> >
{
    typedef CArrayPoint<Type_, SizeCols_, Orient_>  Row;
    typedef CArrayVector<Type_, SizeRows_, Orient_> Col;

    // The CAllocator have to have the same structure than the CArray
    typedef CAllocator<Type_, SizeRows_, SizeCols_, Orient_> Allocator;

    typedef Type_                Type;
    typedef typename RemoveConst<Type>::Type const& ConstReturnType;

    enum
    {
      structure_ = Arrays::array2D_,
      orient_    = Orient_,
      sizeRows_  = SizeRows_,
      sizeCols_  = SizeCols_,
      storage_   = Arrays::dense_
    };
};


} // namespace hidden

/** @ingroup Arrays
 * @brief A CArray is a two dimensional array with a continuous storage like
 * a C-array.
 */
template <typename Type_, int SizeRows_, int SizeCols_, bool Orient_>
class CArray: public ICArray < CArray<Type_, SizeRows_, SizeCols_, Orient_> >
{
  public:
    typedef ICArray < CArray<Type_, SizeRows_, SizeCols_, Orient_> > Base;
    typedef ArrayBase < CArray<Type_, SizeRows_, SizeCols_, Orient_> > LowBase;

    typedef typename hidden::Traits< CArray<Type_, SizeRows_, SizeCols_, Orient_> >::Row Row;
    typedef typename hidden::Traits< CArray<Type_, SizeRows_, SizeCols_, Orient_> >::Col Col;
    typedef typename hidden::Traits< CArray<Type_, SizeRows_, SizeCols_, Orient_> >::Type Type;
    typedef typename hidden::Traits< CArray<Type_, SizeRows_, SizeCols_, Orient_> >::ConstReturnType ConstReturnType;

    enum
    {
      structure_ = hidden::Traits< CArray<Type_, SizeRows_, SizeCols_, Orient_> >::structure_,
      orient_    = hidden::Traits< CArray<Type_, SizeRows_, SizeCols_, Orient_> >::orient_,
      sizeRows_  = hidden::Traits< CArray<Type_, SizeRows_, SizeCols_, Orient_> >::sizeRows_,
      sizeCols_  = hidden::Traits< CArray<Type_, SizeRows_, SizeCols_, Orient_> >::sizeCols_,
      storage_   = hidden::Traits< CArray<Type_, SizeRows_, SizeCols_, Orient_> >::storage_
    };

    /** Default constructor. */
    CArray(): Base() {}
    /** constructor with specified dimension.
     *  @param sizeRows, sizeCols size of the rows and columns
     **/
    CArray( int sizeRows, int sizeCols): Base(sizeRows, sizeCols) {}
    /** constructor with specified ranges.
     *  @param rows, cols range of the rows and columns
     **/
    CArray( Range rows, Range cols): Base(rows.size(), cols.size())
    { this->shift(rows.begin(), cols.begin());}
    /** constructor with specified size, initialization with a constant.
     *  @param sizeRows, sizeCols size of the rows and columns
     *  @param v initial value of the container
     **/
    CArray( int sizeRows, int sizeCols, Type const& v): Base(sizeRows, sizeCols, v) {}
    /** constructor with specified ranges, initialization with a constant.
     *  @param rows, cols range of the rows and columns
     *  @param v initial value of the container
     **/
    CArray( Range rows, Range cols, Type const& v): Base(rows.size(), cols.size(), v)
    { this->shift(rows.begin(), cols.begin());}
    /** Copy constructor
     *  @param T the container to copy
     *  @param ref true if T is wrapped
     **/
    CArray( CArray const& T, bool ref=false): Base(T, ref) {}
    /** Copy constructor by reference, ref_=1.
     *  @param T the container to wrap
     *  @param I,J range of the rows and columns to wrap
     **/
    template<class OtherArray>
    CArray( ICArray<OtherArray> const& T, Range const& I, Range const& J)
          : Base(T.allocator(), I, J) {}
    /** wrapper constructor for 0 based C-Array.
     *  @param q pointer on the array
     *  @param sizeRows, sizeCols size of the rows and columns
     **/
    CArray( Type* const& q, int sizeRows, int sizeCols): Base(q, sizeRows, sizeCols) {}
    /** constructor by reference.
     *  @param allocator the allocator to wrap
     **/
    template<class OtherAllocator>
    CArray( ITContainer2D<OtherAllocator> const& allocator): Base(allocator.asDerived()) {}
    /** Copy constructor using an expression.
     *  @param T the container to wrap
     **/
    template<class OtherDerived>
    CArray( ExprBase<OtherDerived> const& T): Base(T.sizeRows(), T.sizeCols()) { LowBase::operator=(T);}
    /** destructor. */
    ~CArray() {}
    /** operator= : set the container to a constant value.
     *  @param v the value to set
     **/
    CArray& operator=(Type const& v) { return LowBase::setValue(v);}
    /** operator = : overwrite this with the Right hand side rhs.
     *  @param rhs the container to copy
     **/
    template<class Rhs>
    CArray& operator=(ExprBase<Rhs> const& rhs) { return LowBase::assign(rhs.asDerived());}
    /** operator = : overwrite the CArray with the Right hand side rhs.
     *  @param rhs the container to copy
     **/
    CArray& operator=(CArray const& rhs) { return LowBase::assign(rhs);}
};

/** @ingroup Arrays
 *  output stream for CAllocator.
 *  @param s the output stream
 *  @param V the CArray to write
 **/
template <typename Type, int SizeRows_, int SizeCols_, bool Orient_>
ostream& operator<<(ostream& s, const CAllocator<Type, SizeRows_, SizeCols_, Orient_>& V)
{
  CArray<Type, SizeRows_, SizeCols_, Orient_> wrap(V);
  return out2D(s,wrap);
}


} // namespace STK


#endif /* STK_CARRAY_H */
