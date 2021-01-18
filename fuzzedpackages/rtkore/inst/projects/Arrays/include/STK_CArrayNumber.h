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

/** @file STK_CArrayNumber.h
 *  @brief In this file we implement the final class CArrayNumber.
 **/

#ifndef STK_CARRAYNUMBER_H
#define STK_CARRAYNUMBER_H

#include "STK_ICArray.h"

namespace STK
{
// forward declaration
template< typename Type, bool Orient_ = Arrays::by_col_> class CArrayNumber;
template< typename Type, int SizeRows_, int SizeCols_, bool Orient_> class CArray;
template< typename Type, int Size_, bool Orient_> class CArraySquare;
template< typename Type, int SizeCols_, bool Orient_> class CArrayPoint;
template< typename Type, int SizeRows_, bool Orient_ > class CArrayVector;

// typedef for CArrayVector Real is by default double, but can be float
typedef CArrayNumber<Real, Arrays::by_col_>   CNumber;
typedef CArrayNumber<double, Arrays::by_col_> CNumberd;
typedef CArrayNumber<int, Arrays::by_col_>    CNumberi;

typedef CArrayNumber<Real, Arrays::by_row_>   CNumberByRow;
typedef CArrayNumber<double, Arrays::by_row_> CNumberByRowd;
typedef CArrayNumber<int, Arrays::by_row_>    CNumberByRowi;

namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of the Traits class for CArrayNumber class.
 */
template<typename Type_, bool Orient_>
struct Traits< CArrayNumber<Type_, Orient_> >
{
    typedef CArrayNumber<Type_, Orient_> Row;
    typedef CArrayNumber<Type_, Orient_> Col;

    // The CAllocator have to have the same structure than the CArrayNumber
    typedef CAllocator<Type_, 1, 1, Orient_> Allocator;

    typedef Type_                Type;
    typedef typename RemoveConst<Type>::Type const& ConstReturnType;

    enum
    {
      structure_ = Arrays::number_,
      orient_    = Orient_,
      sizeRows_  = 1,
      sizeCols_  = 1,
      size_      = 1,
      storage_   = Arrays::dense_
    };
};


} // namespace hidden


/** @ingroup Arrays
 * @brief specialization for the number case.
 */
template <typename Type_, bool Orient_>
class CArrayNumber: public ICArray < CArrayNumber<Type_, Orient_> >
{
  public:
    typedef ICArray < CArrayNumber<Type_, Orient_> > Base;
    typedef ArrayBase < CArrayNumber<Type_, Orient_> > LowBase;

    typedef typename hidden::Traits< CArrayNumber<Type_, Orient_> >::Row Row;
    typedef typename hidden::Traits< CArrayNumber<Type_, Orient_> >::Col Col;
    typedef typename hidden::Traits< CArrayNumber<Type_, Orient_> >::Type Type;
    typedef typename hidden::Traits< CArrayNumber<Type_, Orient_> >::ConstReturnType ConstReturnType;

    enum
    {
      structure_ = hidden::Traits< CArrayNumber<Type_, Orient_> >::structure_,
      orient_    = hidden::Traits< CArrayNumber<Type_, Orient_> >::orient_,
      sizeRows_  = hidden::Traits< CArrayNumber<Type_, Orient_> >::sizeRows_,
      sizeCols_  = hidden::Traits< CArrayNumber<Type_, Orient_> >::sizeCols_,
      size_      = hidden::Traits< CArrayNumber<Type_, Orient_> >::size_,
      storage_   = hidden::Traits< CArrayNumber<Type_, Orient_> >::storage_
    };
    /** Default constructor. */
    CArrayNumber(): Base() {}
    /** constructor with an initial value, initialization with a constant.
     *  @param v initial value of the container
     **/
    CArrayNumber( Type const& v): Base(1, 1, v) {}
    /** Copy constructor
     *  @param T the container to copy
     *  @param ref true if T is wrapped
     **/
    CArrayNumber( const CArrayNumber &T, bool ref=false): Base(T, ref) {}
    /** wrapper constructor for 0 based C-Array.
     *  @param q pointer on the array
     **/
    CArrayNumber( Type* const& q): Base(q, 1, 1) {}
    /** constructor by reference.
     *  @param allocator the allocator to wrap
     **/
    template<class OtherAllocator>
    CArrayNumber( OtherAllocator const& allocator): Base(allocator) {}
    /** Copy constructor using an expression.
     *  @param T the container to wrap
     **/
    template<class OtherDerived>
    CArrayNumber( ExprBase<OtherDerived> const& T): Base(1, 1) { LowBase::operator=(T);}
    /** destructor. */
    ~CArrayNumber() {}
    /** operator= : set the container to a constant value.
     *  @param v the value to set
     **/
    CArrayNumber& operator=(Type const& v) { return LowBase::setValue(v);}
    /** operator = : overwrite the CArrayNumber with the Right hand side T.
     *  @param T the container to copy
     **/
    template<class Rhs>
    CArrayNumber& operator=(ExprBase<Rhs> const& T) { return LowBase::assign(T);}
    /** operator = : overwrite the CArrayNumber with the Right hand side rhs.
     *  @param rhs the container to copy
     **/
    CArrayNumber& operator=(CArrayNumber const& rhs) { return LowBase::assign(rhs);}
};

} // namespace STK


#endif /* STK_CARRAYNUMBER_H */
