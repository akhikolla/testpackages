/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2017  Serge Iovleff, Université Lille 1, Inria

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
 * created on: 10 mars 2017
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_SAllocator.h
 *  @brief In this file we implement the sparse allocators.
 **/

#ifndef STK_SALLOCATOR_H
#define STK_SALLOCATOR_H

#include "STK_MemSAllocator.h"

namespace STK
{
/** @ingroup Arrays
 *  @brief Allocator for sparse Array classes.
 *  The data are stored either in the compressed sparse row (CSR) format
 *  or compressed sparse column (CSC) format.
 */
template<typename Type, int SizeRows_, int SizeCols_, bool Orient_>
class SAllocator;

namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of the Traits class for CAllocator.
 */
template< typename Type_, int SizeRows_, int SizeCols_, bool Orient_>
struct Traits< SAllocator<Type_, SizeRows_, SizeCols_, Orient_> >
{
  typedef Type_  Type;
  typedef typename RemoveConst<Type_>::Type const& ConstReturnType;

  enum
  {
    structure_ = (SizeRows_==1) ? (SizeCols_== 1) ? Arrays::number_ : Arrays::point_
                                : (SizeCols_== 1) ? Arrays::vector_
                                                  : (SizeRows_ == SizeCols_ && SizeRows_!= UnknownSize) ?
                                                    Arrays::square_ : Arrays::array2D_,
    orient_    = Orient_,
    sizeRows_  = SizeRows_,
    sizeCols_  = SizeCols_,
    storage_   = Arrays::sparse_,
    sizePtr_   = (Orient_ == Arrays::by_col_) ? SizeCols_+1: SizeRows_+1,
    sizeProd_  = ProductSizeRowsBySizeCols<SizeRows_, SizeCols_>::prod_, // theoretical size
    nzMax_     = (sizeProd_ > maxFixedSize ? UnknownSize : sizeProd_),
  };

  typedef SAllocator<Type_, 1, SizeCols_, Orient_> Row;
  typedef SAllocator<Type_, SizeRows_, 1, Orient_> Col;

  typedef MemSAllocator<Type_, sizePtr_, nzMax_> Allocator;

};

} // namespace hidden


// forward declaration
template<class Derived, bool Orient_> class OrientedSAllocator;


/** Specialization for the row oriented storage scheme */
template<class Derived>
class OrientedSAllocator<Derived, Arrays::by_col_>: public ITContainer2D<Derived>
{
  public:
    typedef SAllocatorBase<Derived> Base;
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef typename hidden::Traits<Derived>::Row Row;
    typedef typename hidden::Traits<Derived>::Col Col;

    enum
    {
      structure_ = hidden::Traits<Derived>::structure_,
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;
    typedef Array1D<int>  VectorIdx;

    /** constructor with specified ranges */
    OrientedSAllocator( Range const& I, Range const& J): Base(I, J) {}
    /** copy constructor */
    OrientedSAllocator( OrientedSAllocator const& A, bool ref): Base(A)
    { if (!ref) allocator_.assign(A.allocator_);}
    /** Reference constructor */
    template<class OtherDerived>
    inline OrientedSAllocator( OrientedCAllocator<OtherDerived, Arrays::by_row_> const& A
                             , Range const& I, Range const& J)
                             : Base(I, J), ldx_(A.ldx()), allocator_(A.allocator(), true)
    {}


};

template<class Derived>
class OrientedSAllocator<Derived, Arrays::by_row_>: public ITContainer2D<Derived>
{
  public:
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef typename hidden::Traits<Derived>::Row Row;
    typedef typename hidden::Traits<Derived>::Col Col;

    enum
    {
      structure_ = hidden::Traits<Derived>::structure_,
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;
    typedef Array1D<int>  VectorIdx;

};
/** @brief
 *
 */
template<typename Type, int SizeRows_, int SizeCols_, bool Orient_>
class SAllocator: public OrientedSAllocator<SAllocator<Type, SizeRows_, SizeCols_, Orient_>, Orient_>
{
  public:
    SAllocator ( int n_rows, int n_cols, int nz_max = 0) {}
    ~SAllocator();
};

} // namespace STK

#endif /* STK_SALLOCATOR_H */
