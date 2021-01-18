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
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_IArray1D.h
 *  @brief Interface base class for the Array1D, this is an internal header file,
 *  included by other Containers library headers.
 *
 *  You should not attempt to use it directly but rather used one of the
 *  derived class like Array1D, except if you want to create your own
 *  Container Class.
 **/

#ifndef STK_IARRAY1D_H
#define STK_IARRAY1D_H

#include "STK_ExprBase.h"
#include "allocators/STK_MemAllocator.h"
#include "STK_ITContainer1D.h"

namespace STK
{

/** @ingroup Arrays
 *  @brief template one dimensional Array.
 *
 * An IArray1D is a template one column container implementing the interface
 * base class ITContainer1D.
 **/
template<class Derived >
class IArray1D: public ITContainer1D<Derived>
{
  public:
    enum
    {
      size_ = hidden::Traits<Derived>::size_,
      isFixedSize_ = (size_ != UnknownSize)
    };
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef typename hidden::Traits<Derived>::ConstReturnType ConstReturnType;

    typedef typename hidden::Traits<Derived>::RowRange RowRange;
    typedef typename hidden::Traits<Derived>::ColRange ColRange;

    typedef MemAllocator<Type, size_> Allocator;
    typedef ITContainer1D< Derived > Base;

    using Base::range;
    using Base::begin;
    using Base::end;
    using Base::size;
    using Base::elt;
    using Base::setRange;

  protected:

    /** Default constructor. */
    IArray1D();
    /** constructor with a specified Range.
      *  @param I range of the container
     **/
    IArray1D( Range const& I);
    /** Misc constructor with first and last, initialization with a constant.
     *  @param I range of the container
     *  @param v initial value of the container
     **/
    IArray1D( Range const& I, Type const& v);
    /** Copy constructor
     *  @param T the container to copy
     *  @param ref true if T is wrapped
     **/
    IArray1D( const IArray1D &T, bool ref =false);
    /** Copy constructor
     *  @param T the container to copy
     **/
    template<class OtherDerived>
    IArray1D( ExprBase<OtherDerived> const& T);
    /** constructor by reference, ref_=1.
     *  @param T, I the container and the range of data to wrap
     **/
    IArray1D( IArray1D const& T, Range const& I);
    /** constructor by reference, ref_=1.
     *  @param T,I the container and the range of data to wrap
     **/
    template<class OtherDerived>
    IArray1D( IArray1D<OtherDerived> const& T, Range const& I);

    /** destructor: allocated memory is liberated by MemAllocator base class.*/
    ~IArray1D() {}

  public:
    /** @return @c true if *this is reference container, @c false otherwise */
    inline bool isRef() const { return allocator_.isRef();}
    /** Modify the state of the container: this become a reference (if ref is
     *  @c true) or the owner of the data (if ref is @c false).
     *  @note To use with care in order to avoid memory leak
     *  @param ref : has top be false if this own its own data
     **/
    inline void setRef(bool ref) const { allocator_.setRef(ref);}

    /** @return a constant pointer on the data set*/
    inline Type* const& p_data() const { return allocator_.p_data();}
    /** @return a constant reference on the allocator */
    Allocator const& allocator() const { return allocator_;}

    /**  @return the range of the rows of the container */
    inline RowRange const& rows() const  { return range();}
     /** @return the index of the first element */
    inline int beginRows() const { return begin();}
    /**  @return the ending index of the elements */
    inline int endRows() const { return end();}
    /**  @return the size of the container */
    inline int sizeRows() const  { return size();}

    /** @return the Horizontal range (1 column) */
    inline ColRange cols() const { return ColRange(1);}
    /** @return the index of the first column */
    inline int beginCols() const { return baseIdx;}
    /**  @return the index of the ending column */
    inline int endCols() const  { return baseIdx+1;}
    /** @return the number of columns */
    inline int sizeCols() const  { return 1;};

    /**  @return the index of the last element */
    inline int lastIdxRows() const  { return this->lastIdx();}
    /**  @return the index of the last element */
    inline int lastIdxCols() const  { return baseIdx;}

    /** @return the maximum possible number of elements without reallocation*/
    int capacity() const { return isRef() ? 0 : allocator_.size();}

    /** access to an element
     *  @param pos index of the element
     *  @return a reference on the element to modify
     **/
    inline Type& elt1Impl(int pos) { return allocator_.data(pos);}
    /** access to a constant element
     *  @param pos index of the const element
     *  @return a constant reference on the element
     **/
    inline ConstReturnType elt1Impl(int pos) const { return allocator_.data(pos);}

    /** New beginning index for the object.
     *  @param beg the index of the first column to set
     **/
    void shiftImpl(int beg = baseIdx);
    /**  Resize the container.
     * - call @c shift
     * - call @c pushBack if there will be more elements
     * - call @c popBack if three will be less elements
     * @param I the range to set to the Array1D
     **/
    Derived& resizeImpl(Range const& I);
    /** reserve internal memory for at least size elements.
     *  @param size number of elements to reserve
     **/
    void reserve(int size);
    /** Clear the object. Memory is liberated and the
     *  range of the Container is set to 0:-1 or 1:0 (@see baseIdx).
     **/
    void clear();
    /** move T to this.
     *  @note : T is not modified but just set as a reference of the data it was owner.
     *  @param T the container to move to this.
     **/
    void move(Derived const& T);
    /** Add n Elements to the end of the container.
     *  @param n number of elements to add
     **/
    Derived& pushBack( int n=1);
    /** Delete last elts of the container.
     *  @param n number of elts to delete
     **/
    Derived& popBack(int n = 1);
    /** Delete n elements at the pos index to the container.
     *  @param pos index where to delete elements
     *  @param n number of elements to delete (default 1)
     **/
    Derived& erase(int pos, int n=1);
    /** Insert n elements at the position pos of the container.
     *  @param pos,n index where to insert the @c n elements (default is 1)
     **/
    Derived& insertElt( int pos, int n =1);
    /** STL compatibility: Insert element @c v at position @c pos of the Array.
     *  @param pos position to insert elements
     *  @param v value to insert
     **/
    Derived& insert( int pos, Type const& v);
    /** STL compatibility: Insert element @c v in the range @c I of the Array.
     *  @param I range of the index where to insert elements
     *  @param v value to insert
     **/
    Derived& insert( Range const& I, Type const& v);
    /** STL compatibility: push front an element.
     *  @param v value to append
     **/
    Derived& push_front(Type const& v);
    /** STL compatibility: append an element v.
     *  @param v value to append
     **/
    Derived& push_back(Type const& v);
    /** Swapping the pos1 elt and the pos2 elt.
     *  @param pos1,pos2 positions of the elements to swap
     **/
    void swap(int pos1, int pos2);
    /** exchange this Container with T.
     *  @param T the Array to exchange with this
     **/
    void exchange(IArray1D &T);
    /** overwrite @c this with @c src.
     *  @note If the size match, @c this is not resized, and in this case,
     *  the method take care of the possibility of overlapping.
     *  @param src the container to assign
     **/
    Derived& assign( IArray1D const& src);

    /** set a value to this container.
     *  @param value the value to set
     **/
    Derived& setValue(Type const& value);

  protected:
    /** function for memory allocation and initialization.
     *  This method will free all allocated memory owned by this container before initialization.
     *  @param I range of the container
     **/
    void initialize(RowRange const& I);
    /** function for memory allocation and initialization.
     *  The range is not set in this method. If a bad_alloc occur, we set the
     *  range of the container to default before throwing it.
     *  @param I range of the container
     **/
    void allocate(RowRange const& I);
    /** Method for memory deallocation. Memory is liberated and the
     *  range of the Container is set to begin:begin-1.
     **/
    void freeMem();

  private:
    Allocator allocator_;
};

template<class Derived >
IArray1D<Derived>::IArray1D(): Base(), allocator_(Arrays::evalRangeCapacity(range())) {}

template<class Derived >
IArray1D<Derived>::IArray1D( Range const& I)
                           : Base(I), allocator_(Arrays::evalRangeCapacity(range())) {}

template<class Derived >
IArray1D<Derived>::IArray1D( Range const& I, Type const& v)
                           : Base(I)
                           , allocator_(Arrays::evalRangeCapacity(range()))
{ for(int i=begin(); i<end(); i++) allocator_.data(i) = v;}

/* Copy constructor
 *  @param T the container to copy
 *  @param ref
 **/
template<class Derived >
IArray1D<Derived>::IArray1D( const IArray1D &T, bool ref)
                           : Base(T)
                           , allocator_(T.allocator_, ref)
{
  if (!ref) // copy data
  {
    allocate(T.range());
    allocator_.memcpy(begin(), T.allocator_, range());
  }
}

/* Copy constructor
 *  @param T the container to copy
 **/
template<class Derived >
template<class OtherDerived>
IArray1D<Derived>::IArray1D( ExprBase<OtherDerived> const& T)
                           : Base(T.range())
                           , allocator_(Arrays::evalRangeCapacity(T.range()))
{
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(OtherDerived);
  for (int i=begin(); i<end(); i++) elt(i)= T.elt(i);
}

template<class Derived >
IArray1D<Derived>::IArray1D( IArray1D const& T, Range const& I)
                          : Base(I)
                          , allocator_(T.allocator_, true)
{}
/* constructor by reference, ref_=1.
 *  @param T,I the container and the range of data to wrap
 **/
template<class Derived >
template<class OtherDerived>
IArray1D<Derived>::IArray1D( IArray1D<OtherDerived> const& T, Range const& I)
                           : Base(I), allocator_(T.allocator(), true) {}

template<class Derived >
void IArray1D<Derived>::shiftImpl(int beg)
{
  // compute increment
  int inc = beg - begin();
  if (inc == 0) return;
  if (isRef())  // is this structure just a pointer ?
  { STKRUNTIME_ERROR_1ARG(IArray1D::shiftImpl,beg,cannot operate on references);}
  // translate range and data
  this->incRange(inc);
  allocator_.shift(beg);
}

template<class Derived >
Derived& IArray1D<Derived>::resizeImpl(Range const& I)
{
  // check if there is something to do
  if ( this->range() == I) return this->asDerived();
  if (isRef()) { STKRUNTIME_ERROR_1ARG(IArray1D::resizeImpl,I,cannot operate on references);}
  // translate
  shiftImpl(I.begin());
  // compute number of elements to delete or add
  const int inc = I.end() - end();
  // adjust size of the container
  if (inc > 0) pushBack(inc);  // more elements
  else         popBack(-inc);  // less elements
  return this->asDerived();
}

template<class Derived >
void IArray1D<Derived>::reserve(int size)
{
  // nothing to do
  if (size < this->capacity() || isFixedSize_) return;
  // is this structure a ptr ?
  if (isRef()) { STKRUNTIME_ERROR_1ARG(IArray1D::reserve,size,cannot operate on references);}
  allocator_.realloc(Range(begin(), size));
}

template<class Derived >
void IArray1D<Derived>::clear()
{
  if (isRef()) return;   // Nothing to do for ref
  freeMem();  // Free Mem
  setRange(); // Set dimension to default
}

template<class Derived >
void IArray1D<Derived>::move(Derived const& T)
{
  if (this->asPtrDerived() == &T) return; // avoid move on itself
  if (!isRef()) { freeMem();}
  allocator_.move(T.allocator_);  // move Allocator part
  setRange(T.range());            // Set ITContainer1D part
}

template<class Derived >
Derived& IArray1D<Derived>::pushBack( int n)
{
#ifdef STK_ARRAYS_VERY_VERBOSE
    stk_cout << _T("Entering IArray1D<Derived>::pushBack(") << n << _T(")\n");
#endif
  // checks
  if (n <= 0) return this->asDerived();
  if (isRef()) { STKRUNTIME_ERROR_1ARG(IArray1D::pushBack,n,cannot operate on references);}
  if (this->empty())  // If container is empty : create it
  { initialize(RowRange(begin(), n));}
  else // otherwise add element from end() position
  { insertElt(end(), n);}
  return this->asDerived();
}

template<class Derived >
Derived& IArray1D<Derived>::popBack(int n)
{
  // checks
  if (n <= 0 || isFixedSize_) return this->asDerived();
  if (isRef())
  { STKRUNTIME_ERROR_1ARG(IArray1D::popBack,n,cannot operate on reference);}
  this->decLast(n);
  if (this->size() <= 0) this->freeMem(); // release mem if there's no more elts
  return this->asDerived();
}

template<class Derived >
Derived& IArray1D<Derived>::erase(int pos, int n)
{
  // checks
  if (n<=0) return this->asDerived();
  if (isRef()) { STKRUNTIME_ERROR_2ARG(IArray1D::erase,pos, n,cannot operate on reference);}
#ifdef STK_BOUNDS_CHECK
  if (begin() > pos)
  { STKOUT_OF_RANGE_2ARG(IArray1D::erase,pos, n,begin() > pos);}
  if (end() <= pos)
  { STKOUT_OF_RANGE_2ARG(IArray1D::erase,pos, n,end() <= pos);}
  if (end() < pos+n)
  { STKOUT_OF_RANGE_2ARG(IArray1D::erase,pos, n,end() < pos+n);}
#endif
  // translate remaining elements and update dimensions
  allocator_.memmove(pos, _R(pos+n, end()-1));
  this->decLast(n);
  if (size() <= 0) this->freeMem(); // if empty release allocated memory
  return this->asDerived();
}

template<class Derived >
Derived& IArray1D<Derived>::insertElt( int pos, int n)
{
  // checks
  if (n <= 0 ) return this->asDerived();
  if (isRef()) { STKRUNTIME_ERROR_2ARG(IArray1D::insertElt,pos,n,cannot operate on references);}
#ifdef STK_BOUNDS_CHECK
  if (begin() > pos)
  { STKOUT_OF_RANGE_2ARG(IArray1D::insertElt,pos, n,begin() > pos);}
  if (end() < pos)
  { STKOUT_OF_RANGE_2ARG(IArray1D::insertElt,pos, n,end() < pos);}
#endif

  // allocate, if necessary, memory for the elements
  if ( (capacity() < this->size()+n) && !isFixedSize_ )
  {
    // temporary container
    IArray1D Taux;
    exchange(Taux);
    // compute range of the container after insertion
    RowRange range(Taux.range());
    range.incLast(n);
    // allocate
    try { allocate(range);}
    catch (Exception const& error)   // if an error occur
    {
      exchange(Taux); // restore this
      throw error;    // and send again the Exception
    }
    // set range
    setRange(Taux.range());
    // copy original elements
    allocator_.memcpy(Taux.begin(), Taux.allocator_, Range(Taux.begin(), pos - begin()) );
    allocator_.memcpy(pos+n, Taux.allocator_, Range(pos, end()-pos) );
  }
  else // enough space -> shift the last elements
  {
    if (!isFixedSize_)
    { allocator_.memmove(pos+n, Range(pos, end() - pos));}
    else
    { allocator_.memmove(pos+n, Range(pos, end() - pos - n));}
  }
  this->incLast(n);
  return this->asDerived();
}

/* STL compatibility: Insert element @c v at position @c pos of the Array.
 *  @param pos position to insert elements
 *  @param v value to insert
 **/
template<class Derived >
Derived& IArray1D<Derived>::insert( int pos, Type const& v)
{
  insertElt(pos, 1);
  elt(pos) = v;
  return this->asDerived();
}

template<class Derived >
Derived&  IArray1D<Derived>::insert( Range const& I, Type const& v)
{
  insertElt(I.begin(), I.size());
  for (int i=I.begin(); i<I.end(); i++) elt(i) = v;
  return this->asDerived();
}

template<class Derived >
Derived& IArray1D<Derived>::push_front(Type const& v)
{
  insert(Range(begin(), 1), v);
  return this->asDerived();
}

template<class Derived >
Derived& IArray1D<Derived>::push_back(Type const& v)
{
  pushBack();
  this->back() = v;
  return this->asDerived();
}

template<class Derived >
void IArray1D<Derived>::swap(int pos1, int pos2)
{
#ifdef STK_BOUNDS_CHECK
  if (begin() > pos1) { STKOUT_OF_RANGE_2ARG(IArray1D::swap,pos1,pos2,begin()>pos1);}
  if (end() <= pos1)  { STKOUT_OF_RANGE_2ARG(IArray1D::swap,pos1,pos2,end()<=pos1);}
  if (begin() > pos2) { STKOUT_OF_RANGE_2ARG(IArray1D::swap,pos1,pos2,begin()>pos2);}
  if (end() <= pos2)  { STKOUT_OF_RANGE_2ARG(IArray1D::swap,pos1,pos2,end()<=pos2);}
#endif
  std::swap(elt(pos1), elt(pos2));
}

template<class Derived >
void IArray1D<Derived>::exchange(IArray1D &T)
{
  allocator_.exchange(T.allocator_);
  Base::exchange(T);
}

template<class Derived >
Derived& IArray1D<Derived>::assign( IArray1D const& src)
{
  if (p_data() == src.p_data()) { return this->asDerived();}
  // Resize if necessary.
  if ( size() != src.size() ) { this->resize(src.range());}
  allocator_.memcpy(begin(), src.allocator_, src.range());
  return this->asDerived();
}

/* set a value to this container.
 *  @param value the value to set
 **/
template<class Derived >
Derived& IArray1D<Derived>::setValue(Type const& v)
{
  for(int i=begin(); i<end(); i++) elt(i) = v;
  return this->asDerived();
}

template<class Derived >
void IArray1D<Derived>::initialize(RowRange const& I)
{
  allocate(I);
  allocator_.setRef(false);
  setRange(I);
}

template<class Derived >
void IArray1D<Derived>::allocate(RowRange const& I)
{
  try
  {
    allocator_.malloc(Arrays::evalRangeCapacity(I));
  }
  catch (Exception const& error)
  {
    setRange(); // if an error occur set default range
    throw error;
  }
}

template<class Derived >
void IArray1D<Derived>::freeMem()
{
  if (isRef()) return;  // Nothing to do for ref
  allocator_.free();
  setRange(Range(begin(),0)); // set range to default (Base)
}

} // namespace STK

#endif // STK_IARRAY1D_H
