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

    Contact: S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  stkpp::Arrays
 * created on: 26 nov. 2007
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Array1D.h
  * @brief In this file we define and implement the final class Array1D
 **/

#ifndef STK_ARRAY1D_H
#define STK_ARRAY1D_H

#include "STK_IArray1D.h"

namespace STK
{
template<class Type_, int Size_ = UnknownSize > class Array1D;

namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of the Traits class for Array1D class.
 **/
template<class Type_, int Size_>
struct Traits< Array1D<Type_, Size_> >
{
  typedef Array1D<Type_, 1>     Row;
  typedef Array1D<Type_, Size_> Col;
  typedef Array1D<Type_, UnknownSize> SubVector;

  typedef Type_ Type;
  typedef typename RemoveConst<Type>::Type const& ConstReturnType;

  enum
  {
    structure_ = Arrays::vector_,
    orient_    = Arrays::by_col_,
    sizeCols_  = 1,
    sizeRows_  = Size_,
    size_      = Size_,
    storage_   = Arrays::dense_ // always dense
  };

  typedef TRange<size_> RowRange;
  typedef TRange<1>     ColRange;

  typedef RandomIterator1D<Array1D<Type_, Size_> > Iterator;
  typedef ConstRandomIterator1D<Array1D<Type_, Size_> > ConstIterator;

  typedef std::reverse_iterator<Iterator> ReverseIterator;
  typedef std::reverse_iterator<ConstIterator> ConstReverseIterator;
};

} // namespace hidden

}

//#include "iterators/STK_RandomIterator1D.h"

namespace STK
{
/** @ingroup Arrays
 *  @brief template one dimensional Arrays.
 *
 * An Array1D is a template non-oriented container (even if the @c Traits
 * struct define it as column oriented) implementing the interface
 * class @c IArray1D.
 *
 * @note It is a final class for the curious recursive paradigm.
 *
 * @tparam Type of the objects stored in the @c Array1D
 * @tparam Size_ size of the vector if it known. Default is @c UnknownSize
 **/
template<class Type_, int Size_ >
class Array1D: public IArray1D< Array1D<Type_, Size_> >
{
  public:
    typedef IArray1D< Array1D<Type_, Size_> > Base;

    typedef typename hidden::Traits< Array1D<Type_, Size_> >::Row Row;
    typedef typename hidden::Traits< Array1D<Type_, Size_> >::Col Col;
    typedef typename hidden::Traits< Array1D<Type_, Size_> >::SubVector SubVector;

    enum
    {
      structure_ = hidden::Traits< Array1D<Type_, Size_> >::structure_,
      orient_    = hidden::Traits< Array1D<Type_, Size_> >::orient_,
      sizeCols_  = hidden::Traits< Array1D<Type_, Size_> >::sizeCols_,
      sizeRows_  = hidden::Traits< Array1D<Type_, Size_> >::sizeRows_,
      size_      = hidden::Traits< Array1D<Type_, Size_> >::size_,
      storage_   = hidden::Traits< Array1D<Type_, Size_> >::storage_
    };

    typedef typename hidden::Traits< Array1D <Type_, Size_> >::Type Type;
    typedef typename hidden::Traits< Array1D <Type_, Size_> >::ConstReturnType ConstReturnType;

    typedef typename hidden::Traits< Array1D<Type_, Size_> >::Iterator Iterator;
    typedef typename hidden::Traits< Array1D<Type_, Size_> >::ConstIterator ConstIterator;
    typedef typename hidden::Traits< Array1D<Type_, Size_> >::ReverseIterator ReverseIterator;
    typedef typename hidden::Traits< Array1D<Type_, Size_> >::ConstReverseIterator ConstReverseIterator;

    typedef TRange<size_> RowRange;
    typedef TRange<1>     ColRange;

    // Compatibility naming scheme with STL
    typedef Iterator iterator;
    typedef ConstIterator const_iterator;
    typedef ReverseIterator reverse_iterator;
    typedef ConstReverseIterator const_reverse_iterator;

    /** Default constructor. */
    Array1D(): Base(){}
    /** constructor with a specified Range
     *  @param I range of the container
     **/
    Array1D( Range const& I): Base(I) {}
    /** Misc constructor with beg and end, initialization with a constant.
     *  @param I,v range and initial value of the container
     **/
    Array1D( Range const& I, Type const& v): Base(I, v) {}
    /** Copy constructor
     *  @param T the container to copy
     *  @param ref @c true if T is wrapped
     **/
    Array1D( Array1D const&T, bool ref =false): Base(T, ref) {}
    /** constructor by reference, ref_=1.
     *  @param T,I the container and the range of data to wrap
     **/
    Array1D( Array1D const& T, Range const& I): Base(T, I) {}
    /** constructor by reference, ref_=1.
     *  @param T,I the container and the range of data to wrap
     **/
    template<int OtherSize>
    Array1D( Array1D<Type, OtherSize> const& T, Range const& I): Base(T, I) {}
    /** Copy an other type of array/expression in an Array1D.
     *  @param T the array/expression to copy
     **/
    template<class OtherArray>
    Array1D( ExprBase<OtherArray> const& T): Base(T.asDerived()) {}
    /** Wrapper constructor: the container is a reference of a C-Array.
     *  @param q, I pointer and range of data
     **/
    Array1D( Type* q, Range const& I): Base(q, I) {}
    /** destructor: allocated memory is liberated by MemAllocator base class. */
    ~Array1D() {}
    /** operator = : overwrite the Array1D with T.
     *  @param T the container to copy
     **/
    Array1D& operator=(Array1D const& T)
    { return this->assign(T);}
    /** Copy an other type of 1D array in an Array1D.
     *  @param T the array to copy
     **/
    template<class OtherArray>
    Array1D& operator=(ITContainer1D<OtherArray> const& T)
    {
      // check size
      if (this->range()!=T.range()) this->resize(T.range());
      for (int i=this->begin(); i<this->end(); i++) this->elt(i) = T.elt(i);
      return *this;
    }
    /** Copy an other type of array/expression in an Array1D.
     *  @param T the array/expression to copy
     **/
    template<class OtherArray>
    Array1D& operator=(ExprBase<OtherArray> const& T)
    {
      // check size
      if (this->size()!=T.size()) this->resize(T.range());
      for (int i=this->begin(); i<this->end(); i++) this->elt(i)= T.elt(i);
      return *this;
    }
    /** operator= : set the container to a constant value.
     *  @param v the value to set
     **/
    Array1D& operator=(Type const& v)
    {
      for (int i=this->begin(); i<this->end(); i++) this->elt(i)= v;
      return *this;
    }
};

} // namespace STK

#endif // STK_ARRAY1D_H
