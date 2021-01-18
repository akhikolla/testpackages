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

/** @file STK_CArraySquare.h
 *  @brief In this file we implement the final class CArraySquare.
 **/

#ifndef STK_CARRAYSQUARE_H
#define STK_CARRAYSQUARE_H

#include "STK_ICArray.h"

namespace STK
{
// forward declarations
template< typename Type_, int Size_ = UnknownSize, bool Orient_ = Arrays::by_col_> class CArraySquare;
template< typename Type, int SizeRows_, int SizeCols_, bool Orient_> class CArray;
template< typename Type, int SizeCols_, bool Orient_> class CArrayPoint;
template< typename Type, int SizeRows_, bool Orient_> class CArrayVector;
template< typename Type, bool Orient_> class CArrayNumber;

// useful typedef
typedef CArraySquare<Real, UnknownSize, Arrays::by_col_>   CSquareX;
typedef CArraySquare<Real, 2, Arrays::by_col_>             CSquare2;
typedef CArraySquare<Real, 3, Arrays::by_col_>             CSquare3;
typedef CArraySquare<double, UnknownSize, Arrays::by_col_> CSquareXd;
typedef CArraySquare<double, 2, Arrays::by_col_>           CSquare2d;
typedef CArraySquare<double, 3, Arrays::by_col_>           CSquare3d;
typedef CArraySquare<int, UnknownSize, Arrays::by_col_>    CSquareXi;
typedef CArraySquare<int, 2, Arrays::by_col_>              CSquare2i;
typedef CArraySquare<int, 3, Arrays::by_col_>              CSquare3i;

typedef CArraySquare<Real, UnknownSize, Arrays::by_row_>   CSquareByRowX;
typedef CArraySquare<Real, 2, Arrays::by_row_>             CSquareByRow2;
typedef CArraySquare<Real, 3, Arrays::by_row_>             CSquareByRow3;
typedef CArraySquare<double, UnknownSize, Arrays::by_row_> CSquareByRowXd;
typedef CArraySquare<double, 2, Arrays::by_row_>           CSquareByRow2d;
typedef CArraySquare<double, 3, Arrays::by_row_>           CSquareByRow3d;
typedef CArraySquare<int, UnknownSize, Arrays::by_row_>    CSquareByRowXi;
typedef CArraySquare<int, 2, Arrays::by_row_>              CSquareByRow2i;
typedef CArraySquare<int, 3, Arrays::by_row_>              CSquareByRow3i;

namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of the Traits class for CArray class.
 */
template<typename Type_, int Size_, bool Orient_>
struct Traits< CArraySquare<Type_, Size_, Orient_> >
{
    typedef CArrayPoint<Type_, Size_, Orient_> Row;
    typedef CArrayVector<Type_, Size_, Orient_> Col;

    // The CAllocator have to have the same structure than the CArray
    typedef CAllocator<Type_, Size_, Size_, Orient_> Allocator;

    typedef Type_                Type;
    typedef typename RemoveConst<Type>::Type const& ConstReturnType;

    enum
    {
      structure_ = Arrays::square_,
      orient_    = Orient_,
      sizeRows_  = Size_,
      sizeCols_  = Size_,
      size_      = Size_,
      storage_   = Arrays::dense_
    };
};


} // namespace hidden

/** @ingroup Arrays
 *  @brief specialization for the square case.
 */
template <typename Type_, int Size_, bool Orient_>
class CArraySquare: public ICArray < CArraySquare<Type_, Size_, Orient_> >
{
  public:
    typedef ICArray < CArraySquare<Type_, Size_, Orient_> > Base;
    typedef ArrayBase < CArraySquare<Type_, Size_, Orient_> > LowBase;

    typedef typename hidden::Traits< CArraySquare <Type_, Size_> >::Row Row;
    typedef typename hidden::Traits< CArraySquare <Type_, Size_> >::Col Col;

    typedef typename hidden::Traits< CArraySquare <Type_, Size_> >::Type Type;
    typedef typename hidden::Traits< CArraySquare <Type_, Size_> >::ConstReturnType ConstReturnType;

    enum
    {
      structure_ = hidden::Traits< CArraySquare <Type_, Size_, Orient_> >::structure_,
      orient_    = hidden::Traits< CArraySquare <Type_, Size_, Orient_> >::orient_,
      sizeRows_  = hidden::Traits< CArraySquare <Type_, Size_, Orient_> >::sizeRows_,
      sizeCols_  = hidden::Traits< CArraySquare <Type_, Size_, Orient_> >::sizeCols_,
      size_      = hidden::Traits< CArraySquare <Type_, Size_, Orient_> >::size_,
      storage_   = hidden::Traits< CArraySquare <Type_, Size_, Orient_> >::storage_
    };
    /** Default constructor. */
    CArraySquare(): Base() {}
    /** constructor with specified dimension.
     *  @param size range of the columns
     **/
    CArraySquare( int const& size): Base(size, size) {}
    /** constructor with specified ranges.
     *  @param range range of the rows and columns
     **/
    CArraySquare( Range range): Base(range.size(), range.size())
    { this->shift(range.begin());}
    /** constructor with specified dimension, initialization with a constant.
     *  @param size range of the columns
     *  @param v initial value of the container
     **/
    CArraySquare( int size, Type const& v): Base(size, size, v) {}
    /** constructor with specified ranges, initialization with a constant.
     *  @param range range of the rows and columns
     *  @param v initial value of the container
     **/
    CArraySquare( Range range, Type const& v): Base(range.size(), range.size(), v)
    { this->shift(range.begin());}
    /** Copy constructor
     *  @param T the container to copy
     *  @param ref true if T is wrapped
     **/
    CArraySquare( CArraySquare const& T, bool ref=false): Base(T, ref) {}
    /** wrapper constructor for 0 based C-Array.
     *  @param q pointer on the array
     *  @param size number of rows/columns
     **/
    CArraySquare( Type* const& q, int size): Base(q, size, size) {}
    /** constructor by reference.
     *  @param allocator the allocator to wrap
     **/
    template<class OtherAllocator>
    CArraySquare( ITContainer2D<OtherAllocator> const& allocator): Base(allocator.asDerived()) {}
    /** Copy constructor using an expression.
     *  @param T the expression to copy
     **/
    template<class OtherDerived>
    CArraySquare( ExprBase<OtherDerived> const& T): Base(T.sizeRows(), T.sizeCols())
    { LowBase::operator=(T);}
    /** destructor. */
    ~CArraySquare() {}
    /** operator= : set the container to a constant value.
     *  @param v the value to set
     **/
    CArraySquare& operator=(Type const& v) { return LowBase::setValue(v);}
    /** operator = : overwrite the CArray with the Right hand side T.
     *  @param T the container to copy
     **/
    template<class Rhs>
    CArraySquare& operator=( ExprBase<Rhs> const& T) { return LowBase::assign(T);}
    /** operator = : overwrite the CArray with the Right hand side rhs.
     *  @param rhs the container to copy
     **/
    CArraySquare& operator=(CArraySquare const& rhs) { return LowBase::assign(rhs);}
};

} // namespace STK


#endif /* STK_CARRAYSQUARE_H */
