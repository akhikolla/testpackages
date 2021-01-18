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
 * Purpose:  Define the Base Interface for the Array classes.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_MemAllocator.h
  * @brief In this file we define the class MemAllocator
 **/

#ifndef STK_MEMALLOCATOR_H
#define STK_MEMALLOCATOR_H

#include <cstring>

#include <STKernel/include/STK_Range.h>
#include <Sdk/include/STK_Macros.h>
#include <Sdk/include/STK_MetaTemplate.h>
#include "../STK_IContainerRef.h"

namespace STK
{

namespace hidden
{

/** @ingroup hidden Helper class allowing to use std::memcpy or std::memmove for
 * fundamental types
 **/
template<int, class Type_> struct memChooser;

/** @ingroup Specialization for fundamental types
 *  copy or move memory using standard C copy and move functions  */
template<class Type_> struct memChooser<1, Type_>
{
  static Type_* memcpy(Type_* p, int pos, Type_* q, Range const& range)
  { return static_cast<Type_*>(std::memcpy( p+pos, q+range.begin(), sizeof(Type_)*range.size()));}
  static Type_* memmove(Type_* p, int pos, Range const& range)
  { return static_cast<Type_*>(std::memmove( p + pos, p+ range.begin(), sizeof(Type_)*range.size()));}
};

/** @ingroup Specialization for other types using loop and operator= */
template<class Type_> struct memChooser<0, Type_>
{
    /** copy mem using loop */
  static Type_* memcpy(Type_* p, int pos, Type_* q, Range const& range)
  {
    Type_* ptr = p+pos-range.begin();
    for (int k=range.begin(); k<range.end(); k++) { ptr[k] = q[k];}
    return p;
  }
  /** move mem without overlapping using loop */
  static Type_* memmove(Type_* p, int pos, Range const& range)
  {
    if (range.begin() < pos)
    {
      for(int iSrc=range.lastIdx(), iDst=pos+range.size()-1; iSrc>=range.begin(); iSrc--, iDst--)
      { p[iDst] = p[iSrc];}
    }
    else
    {
      for(int iSrc=range.begin(), iDst=pos; iSrc<range.end(); iSrc++, iDst++)
      { p[iDst] = p[iSrc];}
    }
   return p;
  }
};

} // namespace hidden
/** @ingroup Arrays
 *  @brief template base class for all Allocator classes.
 *
 *  The MemAllocator class is the base class of all memory allocators. This
 *  class manages the main pointer on the data. It derives from the IContainerRef
 *  class as an array stored in memory can always be wrapped in some way or be
 *  a wrapper of some data stored in memory.
 *
 *  This class can also be used as a concrete class.
 *  @tparam Type_ can be any type of data that can be stored in memory.
 *  @tparam Size_ size of the data if it is known at compile time
 **/
template<typename Type_, int Size_>
struct MemAllocator: public IContainerRef
{
  enum { isNumeric_ = hidden::IsArithmetic<Type_>::yes};
  typedef Type_  Type;
  typedef typename hidden::RemoveConst<Type>::Type const& ConstReturnType;

  typedef TRange<Size_> AllocatorRange;
  using IContainerRef::isRef;
  using IContainerRef::setRef;

  /** Default constructor. */
  MemAllocator();
  /** constructor with specified Range
   *  @param I range of the data
   **/
  MemAllocator( Range const& I);
  /** Copy constructor. We don't know, how the user classes want to copy the
   *  data if this is not a reference.
   *  @param T the array to copy or reference
   *  @param ref is this a wrapper of T ?
   *  @warning if ref is @c false, the derived class is responsible of the data copy
   **/
  MemAllocator( MemAllocator const& T, bool ref = false);
  /** Copy constructor. We don't know, how the user classes want to copy the
   *  data if this is not a reference.
   *  @param T the array to copy or reference
   *  @param ref is this a wrapper of T ?
   *  @warning if ref is @c false, the derived class is responsible of the data copy
   **/
  template<int OtherSize_>
  inline MemAllocator( MemAllocator<Type, OtherSize_> const& T, bool ref = false);
  /** constructor by reference.
   *  @param T,I the allocator and range to reference
   **/
  MemAllocator( MemAllocator const& T, Range const& I);
  /** constructor by reference.
   *  @param T,I the allocator and range to reference
   **/
  template<int OtherSize_>
  inline MemAllocator( MemAllocator<Type, OtherSize_> const& T, Range const& I);
  /** @brief Wrapper or copy constructor.
   *  @param q,I ptr and range of the data to wrap
   *  @param ref is this a wrapper ? If @c true data will not be freed when this
   *  object is released
   **/
  MemAllocator( Type* const& q, Range const& I, bool ref);
  /** Wrapper or copy constructor: second form. This constructor assumes the
   *  data as a C-like array. Thus first index is 0.
   *  @param q, size ptr and size of the data to wrap
   *  @param ref is this a wrapper ? If @c true data will not be freed when this
   *  object is released
   **/
  MemAllocator( Type* const& q, int size, bool ref);
  /** destructor. */
  ~MemAllocator();

  /** @return the range of the data*/
  inline AllocatorRange const& range() const { return range_;}
  /** @return the first index of the data. */
  inline int begin() const { return range_.begin();}
  /**@return the ending index of the data */
  inline int end() const { return range_.end();}
  /** @return the size of the data */
  inline int size() const { return range_.size();}

  /** Get the a constant pointer on p_data_ */
  inline Type* const& p_data() const { return p_data_;}
  /** Get the pointer on p_data_ */
  inline Type*& p_data() { return p_data_;}

  /** Get the const element number pos.
   *  @param pos the position of the element we get
   **/
  inline ConstReturnType data( int pos) const
  {
#ifdef STK_BOUNDS_CHECK
    if (pos < begin())
    { STKOUT_OF_RANGE_1ARG(MemAllocator::data,pos,MemAllocator::begin() > pos);}
    if (pos >= end())
    { STKOUT_OF_RANGE_1ARG(MemAllocator::data,pos,MemAllocator::end() <= pos);}
#endif
    return p_data_[pos];
  }
  /** Get the element number pos.
   *  @param pos the position of the element we get
   **/
  inline Type& data(int pos)
  {
#ifdef STK_BOUNDS_CHECK
    if (pos < begin())
    { STKOUT_OF_RANGE_1ARG(MemAllocator::data,pos,MemAllocator::begin() > pos);}
    if (pos >= end())
    { STKOUT_OF_RANGE_1ARG(MemAllocator::data,pos,MemAllocator::end() <= pos);}
#endif
    return p_data_[pos];
  }
  /** swap two elements of the Allocator.
   *  @param pos1, pos2 the positions of the first and second element
   **/
  void swap(int pos1, int pos2)
  { std::swap(p_data_[pos1], p_data_[pos2]);}
  /** @brief main ptr memory allocation.
   *  @param I range of the data allocated
   **/
  template<int OtherSize>
  void malloc( TRange<OtherSize> const& I);
  /** @brief function for main ptr memory reallocation.
   *
   *  If the size requested is greater than the allocated size,
   *  the Type stored are saved and copied using the operator=. the Type
   *  class have to provide this operator.
   *
   *  If the size requested is lesser than the allocated size, only
   *  the first elements fitting in the container are copied.
   *  @param I range of the data to reserve
   **/
  template<int OtherSize>
  void realloc( TRange<OtherSize> const& I);
  /** function for main ptr memory deallocation. */
  void free();
  /** function copying a part of allocator T in this.
   *  @param pos position where will be copied data
   *  @param T,range the array of data and the range of the data to copy  */
  template<int OtherSize__>
  void memcpy(int pos, MemAllocator<Type, OtherSize__> const& T, Range const& range);

  /** function copying a part of allocator T2 in allocator T2.
   *  @param T1,pos,T2,range the arrays the position and the range of the data to copy
   **/
  template<int OtherSize_1_, int OtherSize_2_>
  static void memcpy( MemAllocator<Type, OtherSize_1_>& T1, int pos
                    , MemAllocator<Type, OtherSize_2_> const& T2, Range const& range)
  { T1.memcpy(pos, T2, range);}

  /** function moving a part of allocator T.
   *  @param pos,range position and range in form [begin,end) to move */
  void memmove(int pos, Range const& range);
  /** function moving a part of allocator T.
   *  @param T,pos,range data, position and range of data to move */
  template<int OtherSize__>
  static void memmove( MemAllocator<Type, OtherSize__>& T, int pos, Range const& range)
  {
#ifdef STK_BOUNDS_CHECK
    if (pos < T.begin()){ STKOUT_OF_RANGE_1ARG(MemAllocator::memmove,pos,T.begin() > pos);}
    if (pos >= T.end()) { STKOUT_OF_RANGE_1ARG(MemAllocator::memmove,pos,T.end() <= pos);}
    if (!T.range().isContaining(range))
    { STKOUT_OF_RANGE_1ARG(MemAllocator::memmove,range,range not in T.range());}
#endif
    std::memmove( T.p_data_ + pos, T.p_data_+ range.begin(), sizeof(Type)*range.size());
  }

  /** exchange this with T.
   *  @param T the container to exchange with T
   **/
  MemAllocator& exchange(MemAllocator &T);
  /** @brief copy the Allocator T by value.
   *  The memory is free and the Allocator T is physically copied in this.
   *  @param T the allocator to copy by value
   *  @return a copy of T
   **/
  MemAllocator& assign( MemAllocator const& T);
  /** @brief move the Allocator T to this.
   *  The memory of this is freed and T becomes a reference of this. This
   *  method allow to move the data of T to this without using physical copy.
   *
   *  @param T the allocator to move to this
   *  @return this object.
   *  @note the data member ref_ is mutable so that T can be passed as a
   *  constant reference.
   **/
  MemAllocator& move( MemAllocator const& T);
  /** shift the first index of the data to first.
   *  @param first the index of the first data to set
   **/
  void shift(int const& first);
  /** @brief Set address and range of allocated data.
   *  This method is to be used when the memory have been allocated outside.
   *  If allocator wrap allocated memory, it is not freed.
   *  @param p_data the address to set
   *  @param range range of the data
   *  @param ref is p_data_ a wrapper ?
   **/
  void setPtr( Type* p_data,  Range const& range, bool ref)
  { p_data_ = p_data; range_ = range; setRef(ref);}

  protected:
    /** Main pointer on the data. */
    Type* p_data_;

  private:
    /** Set array members to default values. */
    void setDefault()
    {
      p_data_ = 0;
      range_ = AllocatorRange();
      setRef(false);
    }
    /** Increment the address of the data.
     *  @param inc the increment to apply
     **/
    void incPtr( int inc)
    {
      if (p_data_) { p_data_ += inc;}
      range_.dec(inc);
    }
    /** Decrement the address of the data.
     *  @param dec the increment to apply
     **/
    void decPtr( int dec)
    {
      if (p_data_) { p_data_ -= dec;}
      range_.inc(dec);
    }
    /** Range of the data */
    AllocatorRange range_;
};

/* Default constructor. */
template<typename Type, int Size_>
MemAllocator<Type,Size_>::MemAllocator(): IContainerRef(false), p_data_(0), range_() {}
/* constructor with specified Range
 *  @param I range of the data
 **/
template<typename Type, int Size_>
MemAllocator<Type,Size_>::MemAllocator( Range const& I): IContainerRef(false), p_data_(0), range_(I)
{ malloc(I);}
/* Copy constructor. We don't know, how the user classes want to copy the
 *  data if this is not a reference.
 *  @param T the array to copy or reference
 *  @param ref is this a wrapper of T ?
 *  @note if ref is @c false, the derived class is responsible of the data copy
 **/
template<typename Type, int Size_>
MemAllocator<Type,Size_>::MemAllocator( MemAllocator const& T, bool ref)
             : IContainerRef(ref)
             , p_data_(ref ? T.p_data_: 0)
             , range_(T.range_)
{/* derived class has to copy the data if ref==false */}
/* Copy constructor. We don't know, how the user classes want to copy the
 *  data if this is not a reference.
 *  @param T the array to copy or reference
 *  @param ref is this a wrapper of T ?
 *  @note if ref is @c false, the derived class is responsible of the data copy
 **/
template<typename Type, int Size_>
template<int OtherSize_>
inline MemAllocator<Type,Size_>::MemAllocator( MemAllocator<Type, OtherSize_> const& T, bool ref)
                    : IContainerRef(ref)
                    , p_data_(ref ? T.p_data(): 0)
                    , range_(T.range())
{ /* derived class has to copy the data if ref==false */}
/* constructor by reference.
 *  @param T,I the allocator and range to reference
 **/
template<typename Type, int Size_>
MemAllocator<Type,Size_>::MemAllocator( MemAllocator const& T, Range const& I)
             : IContainerRef(true), p_data_(T.p_data_), range_(I)
{}
/* constructor by reference.
 *  @param T,I the allocator and range to reference
 **/
template<typename Type, int Size_>
template<int OtherSize_>
inline MemAllocator<Type,Size_>::MemAllocator( MemAllocator<Type, OtherSize_> const& T, Range const& I)
                                             : IContainerRef(true), p_data_(T.p_data()), range_(I)
{}
/* @brief Wrapper or copy constructor.
 *  @param q,I ptr and range of the data to wrap
 *  @param ref is this a wrapper ? If @c true data will not be freed when this
 *  object is released
 **/
template<typename Type, int Size_>
MemAllocator<Type,Size_>::MemAllocator( Type* const& q, Range const& I, bool ref)
             : IContainerRef(ref), p_data_(q), range_(I)
{ /* derived class have to copy the data if ref==false */}
/* Wrapper or copy constructor: second form. This constructor assumes the
 *  data as a C-like array. Thus first index is 0.
 *  @param q, size ptr and size of the data to wrap
 *  @param ref is this a wrapper ? If @c true data will not be freed when this
 *  object is released
 **/
template<typename Type, int Size_>
MemAllocator<Type,Size_>::MemAllocator( Type* const& q, int size, bool ref)
             : IContainerRef(ref), p_data_(q), range_(AllocatorRange(0,size))
{ /* derived class have to copy the data if ref==false */}
/** destructor. */
template<typename Type, int Size_>
MemAllocator<Type,Size_>::~MemAllocator() { free();}

/* exchange this with T.
 *  @param T the container to exchange with T
 **/
template<typename Type, int Size_>
MemAllocator<Type,Size_>& MemAllocator<Type,Size_>::exchange(MemAllocator<Type,Size_> &T)
{
  std::swap(p_data_, T.p_data_);
  std::swap(range_, T.range_);
  IContainerRef::exchange(T);
  return *this;
}
/* @brief copy the Allocator T by value.
 *  The memory is free and the Allocator T is physically copied in this.
 *  @param T the allocator to copy by value
 *  @return a copy of this
 **/
template<typename Type, int Size_>
MemAllocator<Type,Size_>& MemAllocator<Type,Size_>::assign( MemAllocator<Type,Size_> const& T)
{
  // allocate memory if necessary
  malloc(T.range_);
  // copy values
  for (int pos= T.begin(); pos < T.end(); ++pos)
  { p_data_[pos] = T.p_data_[pos];}
  return *this;
}
/* @brief move the Allocator T to this.
 *  The memory is free and T become a reference of this. This method allow
 *  to steal the data of T without physical copy.
 *
 *  @return this object.
 *  @note the data member ref_ is mutable so that T can be passed as a
 *  constant reference.
 *  @param T the allocator to move to this
 **/
template<typename Type, int Size_>
MemAllocator<Type, Size_>& MemAllocator<Type,Size_>::move( MemAllocator<Type, Size_> const& T)
{
  if (this == &T) return *this;
  free();
  setPtr(T.p_data_, T.range_, T.isRef());
  T.setRef(true); // T become a reference of the data it own
  return *this;
}
template<typename Type, int Size_>
void MemAllocator<Type,Size_>::shift(int const& first)
{
  // check if there is something to do
  if (first == begin()) return;
  // check for reference
  if (isRef())
    STKRUNTIME_ERROR_1ARG(MemAllocator::shift, first, cannot operate on reference.);
  // compute increment
  int inc = first - begin();
  // translate data
  decPtr(inc);
}
template<typename Type, int Size_>
template<int OtherSize>
void MemAllocator<Type,Size_>::malloc( TRange<OtherSize> const& I)
{
  // there is no necessity to allocate if data is already allocated, range_
  // is the same and the data is not owned by an other allocator
  if ((range_ == I)&&(p_data_)&&(!isRef())) return;
  // free any existing data
  free();
  if (I.size() <= 0)  // check size
  {
    setPtr(0, I, false);
    return;
  }
  // allocate memory
  try
  { // set (p_data_, range_, ref_)
    setPtr(new Type[I.size()], Range(0, I.size()), false);
  }
  catch (std::bad_alloc const& error)
  {
    setDefault();
    STKRUNTIME_ERROR_1ARG(MemAllocator::malloc, I, memory allocation failed);
  }
  decPtr(I.begin());
}

template<typename Type, int Size_>
template<int OtherSize>
void MemAllocator<Type,Size_>::realloc( TRange<OtherSize> const& I)
{
  // there is no necessity to allocate if data is already allocated, range_
  // is the same and the data is not owned by an other allocator
  if ((range_ == I)&&(p_data_)&&(!isRef())) return;
  if (I.size() <= 0)
  {
    free();
    setPtr(0, I, false);
    return;
  }
  try
  {
    // allocate memory and apply increment
    Type* p  = new Type[I.size()];
    p -= I.begin();
    // copy data and liberate old memory and set (p_data_, range_, ref_)
    range_ = inf(range_, I);
    for (int i = range_.begin(); i<range_.end(); ++i) { p[i] = p_data_[i];}
    free();
    setPtr(p, I, false);
  }
  catch (std::bad_alloc const& error)
  { STKRUNTIME_ERROR_1ARG(MemAllocator::realloc, I, memory allocation failed);}
}
/* function for main ptr memory deallocation. */
template<typename Type, int Size_>
void MemAllocator<Type,Size_>::free()
{
  // nothing to do for reference
  if (isRef()) return;
  // if there is elts
  if (p_data_)
  {
    incPtr(begin());   // translate
    delete [] p_data_; // erase
    setDefault();      // set default values
  }
}

template<typename Type, int Size_>
template<int OtherSize__>
void MemAllocator<Type,Size_>::memcpy(int pos, MemAllocator<Type, OtherSize__> const& T, Range const& range)
{
  if (range.size() <= 0) return;
#ifdef STK_BOUNDS_CHECK
  if (pos < begin()) { STKOUT_OF_RANGE_1ARG(MemAllocator::memcpy,pos,begin() > pos);}
  if (pos >= end())  { STKOUT_OF_RANGE_1ARG(MemAllocator::memcpy,pos,end() <= pos);}
  if (!T.range().isContaining(range))
  { STKOUT_OF_RANGE_1ARG(MemAllocator::memcopy,range,range not in T.range());}
#endif
  hidden::memChooser<isNumeric_, Type>::memcpy(p_data_, pos, T.p_data_, range);
}

template<typename Type, int Size_>
void MemAllocator<Type,Size_>::memmove(int pos, Range const& range)
{
  if (range.size() <= 0) return;
#ifdef STK_BOUNDS_CHECK
  if (pos < begin())
  { STKOUT_OF_RANGE_1ARG(MemAllocator::memmove,pos,MemAllocator::begin() > pos);}
  if (pos >= end())
  { STKOUT_OF_RANGE_1ARG(MemAllocator::memmove,pos,MemAllocator::end() <= pos);}
  if (!range_.isContaining(range))
  { STKOUT_OF_RANGE_1ARG(MemAllocator::memmove,range,range not in range_);}
#endif
  hidden::memChooser<isNumeric_, Type>::memmove( p_data_, pos, range);
}


} // namespace STK

#endif /* STK_MEMALLOCATOR_H */
