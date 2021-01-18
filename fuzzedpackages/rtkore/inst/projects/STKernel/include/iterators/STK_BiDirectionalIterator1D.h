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
 * created on: 4 déc. 2015
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_BiDirectionalIterator1D.h
 *  @brief In this file we implement the BiDirectionalIterator1D class.
 **/

#ifndef STK_BIDIRECTIONALITERATOR1D_H
#define STK_BIDIRECTIONALITERATOR1D_H

#include "STK_IteratorBase.h"

namespace STK
{

// forward declaration
template<class Array> struct BiDirectionalIterator1D;
template<class Array> struct ConstBiDirectionalIterator1D;

namespace hidden
{
  /** @ingroup hidden
   *  @brief Specialization for the RandomIterator1D iterator class
   **/
  template<class Array>
  struct IteratorTraits<BiDirectionalIterator1D<Array> >
  {
    typedef std::bidirectional_iterator_tag iterator_category;
    typedef typename Array::Type value_type;
    typedef int difference_type;
    typedef value_type* pointer;
    typedef value_type& reference;
};

  /** @ingroup hidden
   *  @brief Specialization for the RandomIterator1D iterator class
   **/
  template<class Array>
  struct IteratorTraits<ConstBiDirectionalIterator1D<Array> >
  {
    typedef std::bidirectional_iterator_tag iterator_category;
    typedef typename Array::Type value_type;
    typedef int difference_type;
    typedef value_type const* pointer;
    typedef value_type const& reference;
};

} // namespace hidden


/** @ingroup Arrays
 *  @brief BiDirectionalIterator1D allows to loop over the element of
 *  one dimensional list containers.
 *
 *  @sa STK::List1D
 **/
template<class Array>
struct BiDirectionalIterator1D: public IteratorBase< BiDirectionalIterator1D<Array> >
{
    typedef  IteratorBase< BiDirectionalIterator1D<Array> > Base;
    typedef typename Base::iterator_category iterator_category;
    typedef typename Base::value_type value_type;
    typedef typename Base::reference reference;
    typedef typename Base::pointer pointer;

    using Base::pos_;

  public:
    // creating
    BiDirectionalIterator1D( Array& list1D, int pos)
                          : Base(pos), list1D_(list1D) {}
    BiDirectionalIterator1D( BiDirectionalIterator1D const& it)
                          : Base(it), list1D_(it.list1D_) {}
    ~BiDirectionalIterator1D() {}
    BiDirectionalIterator1D& operator=(BiDirectionalIterator1D const& it)
    { list1D_ = it.list1D_; pos_= it.pos_; return *this;}
    // getting
    reference operator*() { return list1D_[pos_]; }
    pointer operator->()  { return &(list1D_[pos_]); }
    // misc
    friend void swap(BiDirectionalIterator1D& lhs, BiDirectionalIterator1D& rhs)
    {
      Base::swap(lhs, rhs);
      std::swap(lhs.list1D_, rhs.list1D_);
    }

  private:
    Array& list1D_;
};

/** @ingroup Arrays
 *  @brief ConstBiDirectionalIterator1D allows to loop over the element of
 *  one dimensional list containers.
 *
 *  @sa STK::List1D
 **/
template<class Array>
struct ConstBiDirectionalIterator1D:  public IteratorBase< ConstBiDirectionalIterator1D<Array> >
{
    typedef  IteratorBase< ConstBiDirectionalIterator1D<Array> > Base;
    typedef typename Base::iterator_category iterator_category;
    typedef typename Base::value_type value_type;
    typedef typename Base::reference reference;
    typedef typename Base::pointer pointer;

    using Base::pos_;

  public:
    // creating
    ConstBiDirectionalIterator1D( Array const& list1D, int pos)
                               : Base(pos), list1D_(list1D) {}
    ConstBiDirectionalIterator1D( ConstBiDirectionalIterator1D const& it)
                               : Base(it), list1D_(it.list1D_) {}
    ~ConstBiDirectionalIterator1D() {}
    ConstBiDirectionalIterator1D& operator=(ConstBiDirectionalIterator1D const& it)
    { list1D_ = it.list1D_; pos_= it.pos_; return *this;}
    // getting
    reference operator*() { return list1D_[pos_]; }
    pointer operator->()  { return &(list1D_[pos_]); }
    // misc
    friend void swap(ConstBiDirectionalIterator1D& lhs, ConstBiDirectionalIterator1D& rhs)
    {
      Base::swap(lhs, rhs);
      std::swap(lhs.list1D_, rhs.list1D_);
    }

  private:
    Array const& list1D_;
};

} // namespace STK

#endif /* STK_BIDIRECTIONALITERATOR1D_H */
