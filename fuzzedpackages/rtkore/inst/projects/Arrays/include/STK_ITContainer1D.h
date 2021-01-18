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
 * Purpose:  Define the Interface 1D template Container class.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_ITContainer1D.h
 *  @brief in this file we define the interface class ITContainer1D.
 **/

#ifndef STK_ITCONTAINER1D_H
#define STK_ITCONTAINER1D_H


#include <Sdk/include/STK_IRecursiveTemplate.h>
#include <Sdk/include/STK_Macros.h>
#include <STKernel/include/iterators/STK_BiDirectionalIterator1D.h>
#include <STKernel/include/iterators/STK_RandomIterator1D.h>

#include "STK_ArraysTraits.h"
#include "STK_Arrays_Util.h"


namespace STK
{
/** @ingroup Arrays
 *  @brief Interface base class for homogeneous 1D containers.
 *
 * The ITContainer1D class is the template base class for all
 * homogeneous one-dimensional containers containing element of type @c Type
 * where Type is not necessarily a scalar. It assumes that the derived class
 * cannot be part of an expression and is not constant, so that it can be
 * #- shifted,
 * #- resized
 * #- accessed in modification.
 *
 * Implement the curious recursive template paradigm: the template parameter
 * c Derived is the name of the class that implements @c ITContainer1D.
 * For example
 * <code>
 * template<class Type>
 * class Derived: public ITContainer1D< Derived<Type> >
 * {...}
 * </code>
 *
 * The pseudo virtual function defined in this interface and to be implemented
 * by derived classes have the following definitions
 * @code
 *   Type& elt1Impl(int const& pos);
 *   Type const& elt1Impl(int const& pos) const;
 *   void shiftImpl(int const& beg);
 *   Derived& resizeImpl(Range const& I);
 * @endcode
 * It is required that all these methods have to be implemented in Derived classes.
 *
 * @sa List1D, Array1D
 **/
template <class Derived>
class ITContainer1D: public IRecursiveTemplate<Derived>
{
  public:
    enum
    {
      size_ = hidden::Traits<Derived>::size_
    };

    typedef typename hidden::Traits<Derived>::Type Type;
    typedef typename hidden::Traits<Derived>::ConstReturnType ConstReturnType;

    typedef typename hidden::Traits<Derived>::Row Row;
    typedef typename hidden::Traits<Derived>::Col Col;
    typedef typename hidden::Traits<Derived>::SubVector SubVector;

    typedef typename hidden::Traits<Derived>::Iterator Iterator;
    typedef typename hidden::Traits<Derived>::ConstIterator ConstIterator;
    typedef typename hidden::Traits<Derived>::ReverseIterator ReverseIterator;
    typedef typename hidden::Traits<Derived>::ConstReverseIterator ConstReverseIterator;

    typedef TRange<size_> Range1D;

  protected:
    /** Default constructor */
    ITContainer1D(): range_() {}
    /** constructor with a specified range.
     *  @param I : the range of the container
     **/
    ITContainer1D( Range const& I): range_(I) {}
    /** Copy constructor
     *  @param T the container to copy
     **/
    ITContainer1D( ITContainer1D const& T): range_(T.range_) {}
    /** destructor. */
    ~ITContainer1D() {}

  public:
    /**  @return the range of the container */
    inline Range1D const& range() const  { return range_;}
    /** @return the index of the first element */
    inline int begin() const { return range_.begin();}
    /**  @return the ending index of the elements */
    inline int end() const { return range_.end();}
    /**  @return the size of the container */
    inline int size() const { return range_.size();}

    /**  @return the index of the last element */
    inline int lastIdx() const { return range_.lastIdx();}

    /** @return an iterator on the beginning of the container */
    Iterator beginIterator() { return Iterator(this->asDerived(), begin());}
    /**  @return an iterator on the end of the container */
    Iterator endIterator() { return Iterator(this->asDerived(), end());}

    /** @return an iterator on the beginning of the container */
    ConstIterator beginConstIterator() const { return ConstIterator(this->asDerived(), begin());}
    /**  @return an iterator on the end of the container */
    ConstIterator endConstIterator() const { return ConstIterator(this->asDerived(), end());}

    /** @return a reversed iterator on the beginning of the reversed container */
    ReverseIterator rbeginIterator() { return ReverseIterator(endIterator());}
    /**  @return a reversed iterator on the end of the reversed container */
    ReverseIterator rendIterator() { return ReverseIterator(beginIterator());}

    /** @return a reversed iterator on the beginning of the reversed container */
    ConstReverseIterator rbeginConstIterator() const { return ConstReverseIterator(endConstIterator());}
    /**  @return a reversed iterator on the end of the reversed container */
    ConstReverseIterator rendConstIterator() const { return ConstReverseIterator(beginConstIterator());}

    /** Is there some data ?
     *  @return @c true if the container is empty, @c false otherwise
     **/
    bool empty() const { return range_.empty();}

    /** @return the ith element for vector_, point_ and diagonal_ containers
     *  @param i index of the element to get
     **/
    inline Type& elt(int i)
    {
#ifdef STK_BOUNDS_CHECK
      if (begin() > i) { STKOUT_OF_RANGE_1ARG(ITContainer1D::elt, i, begin() > i);}
      if (end() <= i)  { STKOUT_OF_RANGE_1ARG(ITContainer1D::elt, i, end() <= i);}
#endif
      return this->asDerived().elt1Impl(i);
    }
    /** @return a constant reference on the ith element for vector_, point_ and diagonal_ containers
     *  @param i index of the element to get
     **/
    inline ConstReturnType elt(int i) const
    {
#ifdef STK_BOUNDS_CHECK
      if (begin() > i) { STKOUT_OF_RANGE_1ARG(ITContainer1D::elt, i, begin() > i);}
      if (end() <= i)  { STKOUT_OF_RANGE_1ARG(ITContainer1D::elt, i, end() <= i);}
#endif
      return this->asDerived().elt1Impl(i);
    }
    /** @return the element ith element
     *  @param i index of the element to get
     **/
    inline Type& operator[](int i)
    {
#ifdef STK_BOUNDS_CHECK
      if (begin() > i) { STKOUT_OF_RANGE_1ARG(ITContainer1D::operator[], i, begin() > i);}
      if (end() <= i)  { STKOUT_OF_RANGE_1ARG(ITContainer1D::operator[], i, end() <= i);}
#endif
      return elt(i);
    }
    /** @return a constant reference on the ith  element
     *  @param i index of the element to get
     **/
    inline ConstReturnType operator[](int i) const
    {
#ifdef STK_BOUNDS_CHECK
      if (begin() > i) { STKOUT_OF_RANGE_1ARG(ITContainer1D::operator[], i, begin() > i);}
      if (end() <= i)  { STKOUT_OF_RANGE_1ARG(ITContainer1D::operator[], i, end() <= i);}
#endif
      return elt(i);
    }
    /** @return safely the jth element
     *  @param i index of the element
     **/
    inline Type& at(int i)
    {
      if (begin() > i) { STKOUT_OF_RANGE_1ARG(ITContainer1D::at, i, begin() > i);}
      if (end() <= i)  { STKOUT_OF_RANGE_1ARG(ITContainer1D::at, i, end() <= i);}
      return elt(i);
    }
    /** @return safely the constant jth element
     *  @param i index of the element
     **/
    ConstReturnType at(int i) const
    {
      if (begin() > i)
      { STKOUT_OF_RANGE_1ARG(ITContainer1D::at, i, begin() > i);}
      if (end() <= i)
      { STKOUT_OF_RANGE_1ARG(ITContainer1D::at, i, end() <= i);}
      return elt(i);
    }
    /** Access to many elements.
     *  @param I the range of the elements
     *  @return a reference container with the elements of this in the range I
     **/
    inline SubVector sub(Range const& I) const
    {
#ifdef STK_BOUNDS_CHECK
      if ((I.begin()<begin()))
      { STKOUT_OF_RANGE_1ARG(ITContainer1D::sub,I,I.begin()<begin());}
      if ((I.end()>end()))
      { STKOUT_OF_RANGE_1ARG(ITContainer1D::sub,I,I.end()>end());}
#endif
      return SubVector(this->asDerived(),I);
    }

    /** @return a reference on the first element. */
    inline Type& front() { return elt(begin());}
    /** @return a constant reference on the first element */
    inline ConstReturnType front() const { return elt(begin());}
    /** @return a reference on the last element */
    inline Type& back() { return elt(lastIdx());}
    /** @return a constant reference on the last element */
    inline ConstReturnType back() const { return elt(lastIdx());}

    /**  @param beg the index of the first column to set */
    void shift(int beg)
    {
      this->asDerived().shiftImpl(beg);
      range_.shift(beg);
    }
    /** @return the resized container.
     *  @param I the range of the container
     **/
    Derived& resize(Range const& I) { return this->asDerived().resizeImpl(I);}

  protected:
    /** exchange this container with T
     * @param T the container to exchange with T
     **/
     void exchange(ITContainer1D& T) { std::swap(T.range_, range_ );}
    /** Set range of the rows of the container.
     *  @param I the range to set (default empty)
     **/
    void setRange(Range const& I = Range1D()) { range_ = I;}
    /** increment the range of the container (can be negative).
     *  @param inc increment to apply to the range
     **/
    void incRange(int inc =1) { range_.inc(inc);}
    /** increment the beginning of the container (can be negative).
     *  @param inc the increment to apply to the beginning of the range
     **/
    void incFirst(int inc =1) { range_.incFirst(inc);}
    /** decrement the beginning of the container.
     *  @param inc the decrement to apply to the beginning of the range
     **/
    void decFirst(int inc =1) { range_.decFirst(inc);}
    /** increment the end of the container (can be negative).
     *  @param inc the increment to apply to the end of the range
     **/
    void incLast(int inc =1) { range_.incLast(inc);}
    /** decrement the end of the container.
     *  @param inc the decrement to apply to the end of the range
     **/
    void decLast(int inc =1) { range_.decLast(inc);}

  private:
    /** range of the array. */
    Range1D range_;
};

} // namespace STK

#endif // STK_ITCONTAINER1D_H
