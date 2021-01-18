

#ifndef __TYPELISTNORMAL_H__
#define __TYPELISTNORMAL_H__

class NullType {};

struct EmptyType {};

template < class T, class U >
class Typelist
{
 public:
  typedef T Head;
  typedef U Tail;
};


// calculate length
template <class TList> class Length {};
template <> class Length<NullType>
{
 public:
  enum { Value = 0 };
};
template <class T, class U>
class Length< Typelist<T,U> >
{
 public:
  enum { Value = 1 + Length<U>::Value };
};


// define typelists of different lengths (upto twenty for now - easy to extend)
//note : don't remove the space before the closing bracket( >) , removal of this space will cause 
//a build error on the windows platforms  
#define TYPELIST_1(T1) Typelist<T1,NullType>
#define TYPELIST_2(T1,T2) Typelist<T1, TYPELIST_1(T2) >
#define TYPELIST_3(T1,T2,T3) Typelist<T1,TYPELIST_2(T2,T3) >
#define TYPELIST_4(T1,T2,T3,T4) Typelist<T1,TYPELIST_3(T2,T3,T4) >
#define TYPELIST_5(T1,T2,T3,T4,T5) Typelist<T1,TYPELIST_4(T2,T3,T4,T5) >
#define TYPELIST_6(T1,T2,T3,T4,T5,T6) Typelist<T1,TYPELIST_5(T2,T3,T4,T5,T6) >
#define TYPELIST_7(T1,T2,T3,T4,T5,T6,T7) Typelist<T1,TYPELIST_6(T2,T3,T4,T5,T6,T7) >
#define TYPELIST_8(T1,T2,T3,T4,T5,T6,T7,T8) Typelist<T1,TYPELIST_7(T2,T3,T4,T5,T6,T7,T8) >
#define TYPELIST_9(T1,T2,T3,T4,T5,T6,T7,T8,T9) Typelist<T1,TYPELIST_8(T2,T3,T4,T5,T6,T7,T8,T9) >
#define TYPELIST_10(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10) Typelist<T1,TYPELIST_9(T2,T3,T4,T5,T6,T7,T8,T9,T10) >
#define TYPELIST_11(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11) Typelist<T1,TYPELIST_10(T2,T3,T4,T5,T6,T7,T8,T9,T10,T11) >
#define TYPELIST_12(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12) Typelist<T1,TYPELIST_11(T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12) >
#define TYPELIST_13(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13) Typelist<T1,TYPELIST_12(T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13) >
#define TYPELIST_14(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14) Typelist<T1,TYPELIST_13(T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14) >
#define TYPELIST_15(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15) Typelist<T1,TYPELIST_14(T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15) >
#define TYPELIST_16(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16) Typelist<T1,TYPELIST_15(T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16) >
#define TYPELIST_17(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17) Typelist<T1,TYPELIST_16(T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17) >
#define TYPELIST_18(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17,T18) Typelist<T1,TYPELIST_17(T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17,T18) >
#define TYPELIST_19(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17,T18,T19) Typelist<T1,TYPELIST_18(T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17,T18,T19) >
#define TYPELIST_20(T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17,T18,T19,T20) Typelist<T1,TYPELIST_19(T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17,T18,T19,T20) >

template <class TList, unsigned int nIndex> 
class TypeAt;

template <class Head, class Tail>
class TypeAt<Typelist<Head,Tail>, 0>
{
  public:
  typedef Head Result;
};


template <class Head, class Tail, unsigned int i>
class TypeAt<Typelist<Head,Tail>, i>
{
  public:
  typedef typename TypeAt<Tail,i-1>::Result Result;
};


#ifdef SUN


namespace Private
{
  template <bool>
    class  SelectImpl
    {
    public:
      template <class T, class U>
	class In
	{
	public:
	  typedef T Result;
	};
    };
  
  template <>
    class SelectImpl<false>
    {
    public:
      template <class T, class U>
	class In
	{
	public:
	  typedef U Result;
	};
    };
  
}	// end of namespace private

template <bool flag, typename T, typename U>
    class Select
    {
    public:
      typedef typename Private::SelectImpl<flag>::template In<T, U>::Result Result;
    };

template <class TList, unsigned int i, class DefType = NullType> class TypeAtNonStrict;

namespace Private
{
  // if TList is not NullType, check if Index is 0.
  // if Index is 0, the result is TList::Head
  // if Index is > 0, the result is the result of appliying TypeAtNonStrict
  // to the list's and Index-1
  template <class TList>
	class TypeAtNonStrictImpl
	{
	public:
	  template <class DefType, unsigned int Index>
	    class In
	    {
	    public:
	      typedef typename Select
		<
		Index == 0,				// The condition
		typename TList::Head,	// true-case
		typename TypeAtNonStrict<typename TList::Tail, Index-1, DefType>::Result
		>::Result Result;
	    };
	};
  
  // if TList is NullType the result is *always* the specified DefaultType.
  template <>
    class TypeAtNonStrictImpl<NullType>
    {
    public:
      template <class DefType, unsigned int Index>
	class In
	{
	public:
	  typedef DefType Result;
	};
    };
  
}	// end of namespace Private

template <class TList, unsigned int i, class DefType>
class TypeAtNonStrict
{
 public:
  typedef typename
    Private::TypeAtNonStrictImpl<TList>::template In<DefType, i>::Result Result;
};

#else


template <class TList, unsigned int nIndex, typename DefaultType = NullType> 
class TypeAtNonStrict
{
  public:
  typedef DefaultType Result;
};

template <class Head, class Tail, typename DefaultType>
class TypeAtNonStrict<Typelist<Head,Tail>, 0, DefaultType>
{
  public:
  typedef Head Result;
};

template <class Head, class Tail, unsigned int i, typename DefaultType>
class TypeAtNonStrict<Typelist<Head,Tail>, i, DefaultType>
{
  public:
  typedef typename TypeAtNonStrict<Tail,i-1,DefaultType>::Result Result;
};



#endif








template <class TList, class T>  class Erase;

template <class T>
class Erase<NullType,T>
{
 public:
  typedef NullType Result;
};

template <class T, class Tail>
class Erase<Typelist<T, Tail>, T>
{
 public:
  typedef Tail Result;
};

template <class Head, class Tail, class T>
class Erase<Typelist<Head, Tail>, T>
{
 public:
  typedef Typelist<Head,typename Erase<Tail,T>::Result> Result;
};


template <class TList> class NoDuplicates;

template<>
class NoDuplicates<NullType>
{
 public:
  typedef NullType Result;
};


template<class Head, class Tail>
class NoDuplicates<Typelist<Head,Tail> >
{
 private:
  typedef typename NoDuplicates<Tail>::Result L1;
  typedef typename Erase<L1,Head>::Result L2;

 public:
  typedef Typelist<Head,L2> Result;
};


template <class TList, class T>  class IndexOf;

template <class T>
class IndexOf<NullType,T>
{
 public:
  enum {Value = -1};
};


template <class T, class Tail>
class IndexOf<Typelist<T,Tail>,T>
{
 public:
  enum {Value = 0};
};

template <class Head, class Tail, class T>
class IndexOf<Typelist<Head,Tail>,T>
{
 private:
  enum {Temp = IndexOf<Tail,T>::Value};
 public:
  enum {Value = Temp == -1 ? -1 : 1 + Temp};
};

template <typename T>
class Type2Type
{
 public:
  typedef T OriginalType;
};

template <int v>
class Int2Type
{
 public:
  enum { value = v };
};

#endif


