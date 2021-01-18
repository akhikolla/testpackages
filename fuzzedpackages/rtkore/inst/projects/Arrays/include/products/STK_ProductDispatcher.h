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
 * created on: 25 déc. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ProductDispatcher.h
 *  @brief In this file we select the product method to use.
 **/

#ifndef STK_PRODUCTDISPATCHER_H
#define STK_PRODUCTDISPATCHER_H

#include "STK_ProductRaw.h"
#include "STK_ArrayByVectorProduct.h"
#include "STK_ArrayByArrayProduct.h"

namespace STK
{

// forward declarations (needed because there is no CAllocator with this structure)
template<typename> class Array2DUpperTriangular;
template<typename> class Array2DLowerTriangular;

namespace hidden
{

//------------------------------------------------------------
// general lhs case. result is general except if rhs is vector
/** specialization for lhs is array2D */
template<typename Lhs, typename Rhs, int RStructure_>
struct ProductTraits<Lhs, Rhs, Arrays::array2D_, RStructure_>
{
  enum
  {
    structure_ = Arrays::array2D_,
    sizeRows_  = Traits<Lhs>::sizeRows_,
    sizeCols_  = Traits<Rhs>::sizeCols_,
    // use rhs orientation by default, otherwise orientation of the panel size
    orient_    = ( Traits<Lhs>::sizeRows_ == UnknownSize || Traits<Rhs>::sizeCols_ == UnknownSize)
                 ? Traits<Rhs>::orient_
                 : int(Traits<Lhs>::sizeRows_) > int(Traits<Rhs>::sizeCols_)
                   ? int(Traits<Lhs>::orient_) : int(Traits<Rhs>::orient_),
    storage_  = ( Traits<Lhs>::storage_ == int(Arrays::dense_))||(Traits<Rhs>::storage_ == int(Arrays::dense_))
                ? int(Arrays::dense_) : int(Arrays::sparse_)
    , useForRows_ = Arrays::useLhsSize_
    , useForCols_ = Arrays::useRhsSize_
  };
  typedef typename hidden::Promote< typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type const& ConstReturnType;
  typedef CAllocator<Type, sizeRows_, sizeCols_, orient_> Allocator;
};

//------------------------------------------------------------
// general lhs case. result is general except if rhs is vector
/** specialization for lhs is array2D */
template<typename Lhs, typename Rhs>
struct ProductTraits<Lhs, Rhs, Arrays::array2D_, Arrays::square_>
{
  enum
  {
    structure_ = Arrays::array2D_,
    sizeRows_  = Traits<Lhs>::sizeRows_,
    sizeCols_  = Traits<Lhs>::sizeCols_ != UnknownSize ?
                 int(Traits<Lhs>::sizeCols_) : int(Traits<Rhs>::sizeCols_),
    // use rhs orientation by default, otherwise orientation of the panel size
    orient_    = ( Traits<Lhs>::sizeRows_ == UnknownSize || Traits<Rhs>::sizeCols_ == UnknownSize)
                 ? Traits<Rhs>::orient_
                 : int(Traits<Lhs>::sizeRows_) > int(Traits<Rhs>::sizeCols_)
                   ? int(Traits<Lhs>::orient_) : int(Traits<Rhs>::orient_),
    storage_  = ( Traits<Lhs>::storage_ == int(Arrays::dense_)) || (Traits<Rhs>::storage_ == int(Arrays::dense_))
                ? int(Arrays::dense_) : int(Arrays::sparse_)
    , useForRows_ = Arrays::useLhsSize_
    , useForCols_ = Traits<Lhs>::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
  };

  typedef typename hidden::Promote< typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type const& ConstReturnType;
  typedef CAllocator<Type, sizeRows_, sizeCols_, orient_> Allocator;
};

//----------------------------------
// square lhs case. result is general except if rhs square
template<typename Lhs, typename Rhs, int RStructure_>
struct ProductTraits<Lhs, Rhs, Arrays::square_, RStructure_>
{
  enum
  {
    structure_ = Arrays::array2D_,
    sizeRows_ = int(Traits<Lhs>::sizeRows_) < int(Traits<Rhs>::sizeRows_)
              ? int(Traits<Lhs>::sizeRows_) : int(Traits<Rhs>::sizeRows_),

    sizeCols_ = Traits<Rhs>::sizeCols_,
    // use rhs orientation by default, otherwise orientation of the panel size
    orient_    = ( Traits<Lhs>::sizeRows_ == UnknownSize || Traits<Rhs>::sizeCols_ == UnknownSize)
                 ? Traits<Rhs>::orient_
                 : int(Traits<Lhs>::sizeRows_) > int(Traits<Rhs>::sizeCols_)
                   ? int(Traits<Lhs>::orient_) : int(Traits<Rhs>::orient_),
    storage_  = ( Traits<Lhs>::storage_ == int(Arrays::dense_)) || (Traits<Rhs>::storage_ == int(Arrays::dense_))
                ? int(Arrays::dense_) : int(Arrays::sparse_)
    , useForRows_ = Traits<Lhs>::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
    , useForCols_ = Arrays::useRhsSize_
  };
  typedef typename hidden::Promote< typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type const& ConstReturnType;

  typedef CAllocator<Type, sizeRows_,sizeCols_, orient_> Allocator;
};

//----------------------------------
// square rhs case. result is general except if rhs square
template<typename Lhs, typename Rhs, int LStructure_>
struct ProductTraits<Lhs, Rhs, LStructure_, Arrays::square_>
{
  enum
  {
    structure_ = Arrays::array2D_,
    sizeRows_  = Traits<Lhs>::sizeRows_,

    sizeCols_  = int(Traits<Lhs>::sizeCols_) < int(Traits<Rhs>::sizeCols_)
               ? int(Traits<Lhs>::sizeCols_) : int(Traits<Rhs>::sizeCols_),

    // use rhs orientation by default, otherwise orientation of the panel size
    orient_    = ( Traits<Lhs>::sizeRows_ == UnknownSize || Traits<Rhs>::sizeCols_ == UnknownSize)
                 ? Traits<Rhs>::orient_
                 : int(Traits<Lhs>::sizeRows_) > int(Traits<Rhs>::sizeCols_)
                   ? int(Traits<Lhs>::orient_) : int(Traits<Rhs>::orient_),
    storage_  = ( Traits<Lhs>::storage_ == int(Arrays::dense_)) || (Traits<Rhs>::storage_ == int(Arrays::dense_))
                ? int(Arrays::dense_) : int(Arrays::sparse_)
   , useForRows_ = Arrays::useLhsSize_
   , useForCols_ = Traits<Lhs>::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
  };
  typedef typename hidden::Promote< typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type const& ConstReturnType;

  typedef CAllocator<Type, sizeRows_,sizeCols_, orient_> Allocator;
};

template<typename Lhs, typename Rhs>
struct ProductTraits<Lhs, Rhs, Arrays::square_, Arrays::square_>
{
  enum
  {
    structure_ = Arrays::square_,
    sizeRows_ = (int(Traits<Lhs>::sizeRows_) < int(Traits<Rhs>::sizeRows_))
                ? int(Traits<Lhs>::sizeRows_) : int(Traits<Rhs>::sizeRows_),
    sizeCols_ = (int(Traits<Lhs>::sizeCols_) < int(Traits<Rhs>::sizeCols_))
                ? int(Traits<Lhs>::sizeCols_) : int(Traits<Rhs>::sizeCols_),
    // use rhs orientation by default, otherwise orientation of the panel size
    orient_    = ( Traits<Lhs>::sizeRows_ == UnknownSize || Traits<Rhs>::sizeCols_ == UnknownSize)
                 ? Traits<Rhs>::orient_
                 : int(Traits<Lhs>::sizeRows_) > int(Traits<Rhs>::sizeCols_)
                   ? int(Traits<Lhs>::orient_) : int(Traits<Rhs>::orient_),
    storage_  = ( Traits<Lhs>::storage_ == int(Arrays::dense_)) || (Traits<Rhs>::storage_ == int(Arrays::dense_))
                ? int(Arrays::dense_) : int(Arrays::sparse_)
    , useForRows_ = Traits<Lhs>::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
    , useForCols_ = Traits<Rhs>::sizeCols_ != UnknownSize ? Arrays::useRhsSize_ : Arrays::useLhsSize_
  };
  typedef typename hidden::Promote< typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type const& ConstReturnType;

  typedef CAllocator<Type, sizeRows_, sizeCols_, orient_> Allocator;
};

//----------------------------------
// lower triangular case
template<typename Lhs, typename Rhs, int RStructure_>
struct ProductTraits<Lhs, Rhs, Arrays::lower_triangular_, RStructure_>
{
  enum
  {
    structure_ = Arrays::array2D_,
    sizeRows_ = Traits<Lhs>::sizeRows_,
    sizeCols_ = Traits<Rhs>::sizeCols_,
    // use rhs orientation by default, otherwise orientation of the panel size
    orient_    = ( Traits<Lhs>::sizeRows_ == UnknownSize || Traits<Rhs>::sizeCols_ == UnknownSize)
                 ? Traits<Rhs>::orient_
                 : int(Traits<Lhs>::sizeRows_) > int(Traits<Rhs>::sizeCols_)
                   ? int(Traits<Lhs>::orient_) : int(Traits<Rhs>::orient_),
    storage_  = ( Traits<Lhs>::storage_ == int(Arrays::dense_)) || (Traits<Rhs>::storage_ == int(Arrays::dense_))
                ? int(Arrays::dense_) : int(Arrays::sparse_)
    , useForRows_ = Arrays::useLhsSize_
    , useForCols_ = Arrays::useRhsSize_
  };
  typedef typename hidden::Promote< typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type const& ConstReturnType;

  typedef CAllocator<Type, sizeRows_, sizeCols_, orient_> Allocator;
};

template<typename Lhs, typename Rhs>
struct ProductTraits<Lhs, Rhs, Arrays::lower_triangular_, Arrays::square_>
{
  enum
  {
    structure_ = Arrays::array2D_,
    sizeRows_ = Traits<Lhs>::sizeRows_,
    sizeCols_ = (int(Traits<Lhs>::sizeCols_) < int(Traits<Rhs>::sizeCols_))
                ? int(Traits<Lhs>::sizeCols_) : int(Traits<Rhs>::sizeCols_),
    // use rhs orientation by default, otherwise orientation of the panel size
    orient_    =  ( Traits<Lhs>::sizeRows_ == UnknownSize || Traits<Rhs>::sizeCols_ == UnknownSize)
                 ? Traits<Rhs>::orient_
                 : int(Traits<Lhs>::sizeRows_) > int(Traits<Rhs>::sizeCols_)
                   ? int(Traits<Lhs>::orient_) : int(Traits<Rhs>::orient_),
    storage_  = ( Traits<Lhs>::storage_ == int(Arrays::dense_)) || (Traits<Rhs>::storage_ == int(Arrays::dense_))
                ? int(Arrays::dense_) : int(Arrays::sparse_)
    , useForRows_ = Arrays::useLhsSize_
    , useForCols_ = Traits<Lhs>::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
  };
  typedef typename hidden::Promote< typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type const& ConstReturnType;

  typedef CAllocator<Type, sizeRows_, sizeCols_, orient_> Allocator;
};

template<typename Lhs, typename Rhs>
struct ProductTraits<Lhs, Rhs, Arrays::lower_triangular_, Arrays::lower_triangular_>
{
  enum
  {
    structure_ = Arrays::lower_triangular_,
    sizeRows_  = Traits<Lhs>::sizeRows_,
    sizeCols_  = Traits<Rhs>::sizeCols_,
    // use rhs orientation by default, otherwise orientation of the panel size
    orient_    =  ( Traits<Lhs>::sizeRows_ == UnknownSize || Traits<Rhs>::sizeCols_ == UnknownSize)
                 ? Traits<Rhs>::orient_
                 : int(Traits<Lhs>::sizeRows_) > int(Traits<Rhs>::sizeCols_)
                   ? int(Traits<Lhs>::orient_) : int(Traits<Rhs>::orient_),
    storage_  = ( Traits<Lhs>::storage_ == int(Arrays::dense_)) || (Traits<Rhs>::storage_ == int(Arrays::dense_))
                ? int(Arrays::dense_) : int(Arrays::sparse_)
   , useForRows_ = Arrays::useLhsSize_
   , useForCols_ = Arrays::useRhsSize_
};
  typedef typename hidden::Promote< typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type const& ConstReturnType;

  typedef Array2DLowerTriangular<Type> Allocator; // no CAllocator
};

//----------------------------------
// upper triangular case
template<typename Lhs, typename Rhs, int RStructure_>
struct ProductTraits<Lhs, Rhs, Arrays::upper_triangular_, RStructure_>
{
  enum
  {
    structure_ = Arrays::array2D_,
    sizeRows_ = Traits<Lhs>::sizeRows_,
    sizeCols_ = Traits<Rhs>::sizeCols_,
    // use rhs orientation by default, otherwise orientation of the panel size
    orient_    =  ( Traits<Lhs>::sizeRows_ == UnknownSize || Traits<Rhs>::sizeCols_ == UnknownSize)
                 ? Traits<Rhs>::orient_
                 : int(Traits<Lhs>::sizeRows_) > int(Traits<Rhs>::sizeCols_)
                   ? int(Traits<Lhs>::orient_) : int(Traits<Rhs>::orient_),
    storage_  = ( Traits<Lhs>::storage_ == int(Arrays::dense_)) || (Traits<Rhs>::storage_ == int(Arrays::dense_))
                ? int(Arrays::dense_) : int(Arrays::sparse_)
    , useForRows_ = Arrays::useLhsSize_
    , useForCols_ = Arrays::useRhsSize_
  };
  typedef typename hidden::Promote< typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type const& ConstReturnType;

  typedef CAllocator<Type, sizeRows_ , sizeCols_, orient_> Allocator;
};

template<typename Lhs, typename Rhs>
struct ProductTraits<Lhs, Rhs, Arrays::upper_triangular_, Arrays::square_>
{
  enum
  {
    structure_ = Arrays::array2D_,
    sizeRows_ = Traits<Lhs>::sizeRows_,
    sizeCols_ = (int(Traits<Lhs>::sizeCols_) < int(Traits<Rhs>::sizeCols_))
                ? int(Traits<Lhs>::sizeCols_) : int(Traits<Rhs>::sizeCols_),
    // use rhs orientation by default, otherwise orientation of the panel size
    orient_    = ( Traits<Lhs>::sizeRows_ == UnknownSize || Traits<Rhs>::sizeCols_ == UnknownSize)
                 ? Traits<Rhs>::orient_
                 : int(Traits<Lhs>::sizeRows_) > int(Traits<Rhs>::sizeCols_)
                   ? int(Traits<Lhs>::orient_) : int(Traits<Rhs>::orient_),
    storage_  = ( Traits<Lhs>::storage_ == int(Arrays::dense_)) || (Traits<Rhs>::storage_ == int(Arrays::dense_))
                ? int(Arrays::dense_) : int(Arrays::sparse_)
   , useForRows_ = Arrays::useLhsSize_
   , useForCols_ = Traits<Lhs>::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
  };
  typedef typename hidden::Promote< typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type const& ConstReturnType;

  typedef CAllocator<Type, sizeRows_, sizeCols_, orient_> Allocator;
};

template<typename Lhs, typename Rhs>
struct ProductTraits<Lhs, Rhs, Arrays::upper_triangular_, Arrays::upper_triangular_>
{
  enum
  {
    structure_ = Arrays::upper_triangular_,
    sizeRows_  = Traits<Lhs>::sizeRows_,
    sizeCols_  = Traits<Rhs>::sizeCols_,
    // use rhs orientation by default, otherwise orientation of the panel size
    orient_    = ( Traits<Lhs>::sizeRows_ == UnknownSize || Traits<Rhs>::sizeCols_ == UnknownSize)
                 ? Traits<Rhs>::orient_
                 : int(Traits<Lhs>::sizeRows_) > int(Traits<Rhs>::sizeCols_)
                   ? int(Traits<Lhs>::orient_) : int(Traits<Rhs>::orient_),
    storage_  = ( Traits<Lhs>::storage_ == int(Arrays::dense_)) || (Traits<Rhs>::storage_ == int(Arrays::dense_))
                ? int(Arrays::dense_) : int(Arrays::sparse_)
   , useForRows_ = Arrays::useLhsSize_
   , useForCols_ = Arrays::useRhsSize_
  };
  typedef typename hidden::Promote< typename Lhs::Type, typename Rhs::Type>::result_type Type;
  typedef typename RemoveConst<Type>::Type const& ConstReturnType;

  typedef Array2DUpperTriangular<Type> Allocator;  // no CAllocator
};

/** @ingroup hidden
 *  Dispatcher allowing to choose the way to multiply two expressions.
 *
 *  @note In the some cases, e.g. when some of the structures are triangular for
 *  examples, we use the usual matrix multiplication formula. This can be enhanced
 *  in future versions.
 **/
template < class Lhs, class Rhs, class Result
         , int lhsStructure_ = hidden::Traits<Lhs>::structure_
         , int RhsStructure_ = hidden::Traits<Rhs>::structure_ >
struct ProductDispatcher
{
  enum
  {
    // structure_ = Traits<Result>::structure_,
    sizeRows_  = Traits<Result>::sizeRows_,
    sizeCols_  = Traits<Result>::sizeCols_,
    orient_    = Traits<Result>::orient_,
    storage_   = Traits<Result>::storage_
  };
  typedef MultCoefImpl<Lhs, Rhs, Result> MultCoeff;

  /** loop over the columns of rhs first */
  static void run(Lhs const& lhs, Rhs const& rhs, Result& res )
  {
    if (lhs.sizeRows()<rhs.sizeCols())
    {
      for (int j=rhs.beginCols(); j< rhs.endCols(); j++)
      {
        Range const r = res.rangeRowsInCol(j);
//         int const begin = res.rangeRowsInCol(j).begin(), end = res.rangeRowsInCol(j).end();
        for (int i=r.begin(); i< r.end(); i++)
        { MultCoeff::dot(lhs, rhs, res, i, j);}
      }
    }
    else
    {
      for (int i=lhs.beginRows(); i< lhs.endRows(); i++)
      {
        Range const r = res.rangeColsInRow(i);
        for (int j=r.begin(); j< r.end(); j++)
        { MultCoeff::dot(lhs, rhs, res, i, j);}
 //       stk_cout << "(begin,end)=("<<begin<<","<<end<<")\n";
      }
 //     stk_cout << "product done\n";
    }
  }
};

/** Specialization for the array2d by array2D case. */
template <class Lhs, class Rhs, class Result>
struct ProductDispatcher<Lhs, Rhs, Result, Arrays::array2D_, Arrays::array2D_>
{
  enum
  {
    structure_ = Arrays::array2D_,
    sizeRows_  = Traits<Result>::sizeRows_,
    sizeCols_  = Traits<Result>::sizeCols_,
    orient_    = Traits<Result>::orient_,
    storage_   = Traits<Result>::storage_
  };
  static void run(Lhs const& lhs, Rhs const& rhs, Result& res )
  {
    if (MultCoefImpl<Lhs, Rhs, Result>::multDispatcher(lhs, rhs, res)) return;
    (lhs.sizeRows()<rhs.sizeCols()) ? BlockByPanel<Lhs,Rhs,Result>::run(lhs, rhs, res)
                                    : PanelByBlock<Lhs,Rhs,Result>::run(lhs, rhs, res);
  }
};

template <class Lhs, class Rhs, class Result>
struct ProductDispatcher<Lhs, Rhs, Result, Arrays::array2D_, Arrays::square_>
{
  enum
  {
    structure_ = Arrays::array2D_,
    sizeRows_  = Traits<Result>::sizeRows_,
    sizeCols_  = Traits<Result>::sizeCols_,
    orient_    = Traits<Result>::orient_,
    storage_   = Traits<Result>::storage_
  };
  static void run(Lhs const& lhs, Rhs const& rhs, Result& res )
  {
    if (MultCoefImpl<Lhs, Rhs, Result>::multDispatcher(lhs, rhs, res)) return;
    (lhs.sizeRows()<rhs.sizeCols()) ? BlockByPanel<Lhs,Rhs,Result>::run(lhs, rhs, res)
                                    : PanelByBlock<Lhs,Rhs,Result>::run(lhs, rhs, res);
  }
};

template <class Lhs, class Rhs, class Result>
struct ProductDispatcher<Lhs, Rhs, Result, Arrays::square_, Arrays::square_>
{
  enum
  {
    structure_ = Arrays::square_,
    sizeRows_  = Traits<Result>::sizeRows_,
    sizeCols_  = Traits<Result>::sizeCols_,
    orient_    = Traits<Result>::orient_,
    storage_   = Traits<Result>::storage_
  };
  static void run(Lhs const& lhs, Rhs const& rhs, Result& res )
  {
    if (MultCoefImpl<Lhs, Rhs, Result>::multDispatcher(lhs, rhs, res)) return;
    (lhs.sizeRows()<rhs.sizeCols()) ? BlockByPanel<Lhs,Rhs,Result>::run(lhs, rhs, res)
                                    : PanelByBlock<Lhs,Rhs,Result>::run(lhs, rhs, res);
  }
};

template <class Lhs, class Rhs, class Result>
struct ProductDispatcher<Lhs, Rhs, Result, Arrays::square_, Arrays::array2D_>
{
  enum
  {
    structure_ = Arrays::array2D_,
    sizeRows_  = Traits<Result>::sizeRows_,
    sizeCols_  = Traits<Result>::sizeCols_,
    orient_    = Traits<Result>::orient_,
    storage_   = Traits<Result>::storage_
  };
  static void run(Lhs const& lhs, Rhs const& rhs, Result& res )
  {
    if (MultCoefImpl<Lhs, Rhs, Result>::multDispatcher(lhs, rhs, res)) return;
    (lhs.sizeRows()<rhs.sizeCols()) ? BlockByPanel<Lhs,Rhs,Result>::run(lhs, rhs, res)
                                    : PanelByBlock<Lhs,Rhs,Result>::run(lhs, rhs, res);
  }
};

template <class Lhs, class Rhs, class Result>
struct ProductDispatcher<Lhs, Rhs, Result, Arrays::array2D_, Arrays::vector_>
{
  enum
  {
    structure_ = Arrays::vector_,
    sizeRows_  = Traits<Result>::sizeRows_,
    sizeCols_  = Traits<Result>::sizeCols_,
    orient_    = Traits<Result>::orient_,
    storage_   = Traits<Result>::storage_
  };
  static void run(Lhs const& lhs, Rhs const& rhs, Result& res )
  { bv<Lhs,Rhs,Result>::run(lhs, rhs, res);}
};

template <class Lhs, class Rhs, class Result>
struct ProductDispatcher<Lhs, Rhs, Result, Arrays::square_, Arrays::vector_>
{
  enum
  {
    structure_ = Arrays::vector_,
    sizeRows_  = Traits<Result>::sizeRows_,
    sizeCols_  = Traits<Result>::sizeCols_,
    orient_    = Traits<Result>::orient_,
    storage_   = Traits<Result>::storage_
  };
  static void run(Lhs const& lhs, Rhs const& rhs, Result& res )
  { bv<Lhs,Rhs,Result>::run(lhs, rhs, res);}
};

template <class Lhs, class Rhs, class Result, int lhsStructure_>
struct ProductDispatcher<Lhs, Rhs, Result, lhsStructure_, Arrays::vector_>
{
  enum
  {
    structure_ = Arrays::vector_,
    sizeRows_  = Traits<Result>::sizeRows_,
    sizeCols_  = Traits<Result>::sizeCols_,
    orient_    = Traits<Result>::orient_,
    storage_   = Traits<Result>::storage_
  };
  typedef MultCoefImpl<Lhs, Rhs, Result> MultCoeff;
  static void run(Lhs const& lhs, Rhs const& rhs, Result& res )
  {
    for (int i=lhs.beginRows(); i< lhs.endRows(); i++)
    { MultCoeff::dot(lhs, rhs, res, i);}
  }
};


template <class Lhs, class Rhs, class Result, int RhsStructure_>
struct ProductDispatcher<Lhs, Rhs, Result, Arrays::point_, RhsStructure_>
{
  enum
  {
    structure_ = Arrays::point_,
    sizeRows_  = Traits<Result>::sizeRows_,
    sizeCols_  = Traits<Result>::sizeCols_,
    orient_    = Traits<Result>::orient_,
    storage_   = Traits<Result>::storage_
  };
  typedef MultCoefImpl<Lhs, Rhs, Result> MultCoeff;
  static void run(Lhs const& lhs, Rhs const& rhs, Result& res )
  {
    for (int j=rhs.beginCols(); j< rhs.endCols(); j++)
    { MultCoeff::dot(lhs, rhs, res, j);}
  }
};

template <class Lhs, class Rhs, class Result>
struct ProductDispatcher<Lhs, Rhs, Result, Arrays::point_, Arrays::array2D_>
{
  enum
  {
    structure_ = Arrays::point_,
    sizeRows_  = Traits<Result>::sizeRows_,
    sizeCols_  = Traits<Result>::sizeCols_,
    orient_    = Traits<Result>::orient_,
    storage_   = Traits<Result>::storage_
  };
  static void run(Lhs const& l, Rhs const& r, Result& res )
  { MultPointArray<Lhs, Rhs, Result>::run(l, r, res);}
};

template <class Lhs, class Rhs, class Result>
struct ProductDispatcher<Lhs, Rhs, Result, Arrays::point_, Arrays::square_>
{
  enum
  {
    structure_ = Arrays::point_,
    sizeRows_  = Traits<Result>::sizeRows_,
    sizeCols_  = Traits<Result>::sizeCols_,
    orient_    = Traits<Result>::orient_,
    storage_   = Traits<Result>::storage_
  };
  static void run(Lhs const& l, Rhs const& r, Result& res )
  { MultPointArray<Lhs, Rhs, Result>::run(l, r, res);}
};

} // namespace hidden

} // namespace STK

#endif /* STK_PRODUCTDISPATCHER_H */
