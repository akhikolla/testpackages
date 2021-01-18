
#ifndef ___CPPAD_EIGEN_H___
#define ___CPPAD_EIGEN_H___

#include <Eigen/Dense>
#include <cppad/example/cppad_eigen.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>

#include "generics.h"


using CppAD::AD;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::Array;


// scalar types
typedef double Scalar_0; 
typedef AD<Scalar_0> Scalar_1; 
typedef AD<Scalar_1> Scalar_2; 
typedef AD<Scalar_2> Scalar_3; 
typedef AD<Scalar_3> Scalar_4; 
typedef AD<Scalar_4> Scalar_5; 
typedef AD<Scalar_5> Scalar_6; 
typedef AD<Scalar_6> Scalar_7; 
typedef AD<Scalar_7> Scalar_8; 
typedef AD<Scalar_8> Scalar_9; 

typedef TYPELIST_10(Scalar_0,Scalar_1,Scalar_2,Scalar_3,Scalar_4,Scalar_5,Scalar_6,Scalar_7,Scalar_8,Scalar_9) ScalarTypes;

// matrix types
typedef Matrix<Scalar_0,Dynamic,Dynamic> Mat_0;
typedef Matrix<Scalar_1,Dynamic,Dynamic> Mat_1;
typedef Matrix<Scalar_2,Dynamic,Dynamic> Mat_2;
typedef Matrix<Scalar_3,Dynamic,Dynamic> Mat_3;
typedef Matrix<Scalar_4,Dynamic,Dynamic> Mat_4;
typedef Matrix<Scalar_5,Dynamic,Dynamic> Mat_5;
typedef Matrix<Scalar_6,Dynamic,Dynamic> Mat_6;
typedef Matrix<Scalar_7,Dynamic,Dynamic> Mat_7;
typedef Matrix<Scalar_8,Dynamic,Dynamic> Mat_8;
typedef Matrix<Scalar_9,Dynamic,Dynamic> Mat_9;



typedef TYPELIST_10(Mat_0,Mat_1,Mat_2,Mat_3,Mat_4,Mat_5,Mat_6,Mat_7,Mat_8,Mat_9) MatTypes;

#define MATRIX(T) Matrix<T,Dynamic,Dynamic>
#define MULTIARG(T) std::vector<MATRIX(T) >
#define MULTIMULTIARG(T) std::vector<MULTIARG(T)>
#define TENSOR(T) Array<MATRIX(T),Dynamic,Dynamic>
#define TRIPLE(T) boost::tuple<MULTIARG(T),MATRIX(T),TENSOR(T) >
#define MULTITRIPLE(T) std::vector<TRIPLE(T)>
#define FUNCTION(T) boost::function<MULTIARG(T) (const MULTIARG(T)&)>
#define AFUNCTION(T) boost::function<TRIPLE(T) (const MULTIARG(T)&)>
#define MULTIAFUNCTION(T) std::vector<AFUNCTION(T)>





template <bool b>
struct StaticAssert {};

// template specialized on true
template <>
struct StaticAssert<true>
{
    static void Assert() {}
};


template <class T1, class T2, class TList, int N = -1>
class asADN_impl
{
public:
  static T1 doit(const T2& in)
  {
    enum{NextLevel = IndexOf<TList,T1>::Value + 1};
    enum{Difference = N-1};
    return Value(asADN_impl<typename TypeAt<TList,NextLevel>::Result ,T2,TList,Difference>::doit(in));
  } 
};

template <class T1, class T2, class TList>
class asADN_impl<T1,T2,TList,1>
{
public:
  static T1 doit(const T2& in)
  {
    return Value(in);
  } 
};

template <class T1, class T2, class TList>
class asADN_impl<T1,T2,TList,0>
{
public:
  static T1 doit(const T2& in)
  {
 return in;
  } 
};



template <class T1,class T2>
T1 asADN(const T2& in)
{
enum{ Difference = IndexOf<ScalarTypes,T2>::Value - IndexOf<ScalarTypes,T1>::Value};
StaticAssert<Difference >= 0>::Assert();
return asADN_impl<T1,T2,ScalarTypes,Difference>::doit(in);

}


template <class Tout,class Tin>
Tout Demoter(const Tin& in)
{
unsigned int nrows = in.rows();
unsigned int ncols = in.cols();
Tout out(nrows,ncols);
for(unsigned int i = 0; i < nrows; i++)
  {
    for(unsigned int j = 0; j < ncols; j++)
     {
       out(i,j) = asADN<typename TypeAt<ScalarTypes,IndexOf<MatTypes,Tout>::Value>::Result>(in(i,j));
     }    
  }
return out;
}

template <class Tout,class Tin>
Tout Promoter(const Tin& in)
{
unsigned int nrows = in.rows();
unsigned int ncols = in.cols();
Tout out(nrows,ncols);
for(unsigned int i = 0; i < nrows; i++)
  {
    for(unsigned int j = 0; j < ncols; j++)
     {
       out(i,j) = in(i,j);
     }    
  }
return out;

}

template <class Tout,class Tin>
Tin Retainer(const Tin& in)
{
return in;
}



template <class TList,template <class,class> class Unit> class PairWithSelf;

template <class T1,class T2,class T3,template <class,class> class Unit> 
class PairWithSelf<Typelist<T1,Typelist<T2,T3> >,Unit>  : public PairWithSelf<Typelist<T2,T3>,Unit> , public Unit<T1,T1> {};

template <class T1,class T2,template <class,class> class Unit> 
class PairWithSelf<Typelist<T1,Typelist<T2,NullType> >,Unit>  : public Unit<T1,T1> , public Unit<T2,T2> {};


template <class TList,template <class,class> class Unit> class PairRootWithAllAndBranch;
                                                               
template <class TList,template <class,class> class Unit> class PairRootWithAll;


template <class T1,class T2,class T3,template <class,class> class Unit> 
class PairRootWithAllAndBranch<Typelist<T1,Typelist<T2,T3> >,Unit> : public PairRootWithAll<Typelist<T1,T3>,Unit>, 
								     public PairRootWithAllAndBranch<Typelist<T2,T3>,Unit>, 
								     public Unit<T1,T2>,
								     public Unit<T2,T1>
{
  public:
  PairRootWithAllAndBranch() : Unit<T1,T2>(boost::function<T1 (const T2&)>(Demoter<T1,T2>))  ,
			       Unit<T2,T1>(boost::function<T2 (const T1&)>(Promoter<T2,T1>)) {}
};


template <class T1,class T2,template <class,class> class Unit> 
class PairRootWithAllAndBranch<Typelist<T1,Typelist<T2,NullType> >,Unit> :  public Unit<T1,T2> , 
									    public Unit<T2,T1>
{
public:
  PairRootWithAllAndBranch() : Unit<T1,T2>(boost::function<T1 (const T2&)>(Demoter<T1,T2>)) ,
			       Unit<T2,T1>(boost::function<T2 (const T1&)>(Promoter<T2,T1>)) {}
};


template <class T1,class T2,class T3,template <class,class> class Unit> 
class PairRootWithAll<Typelist<T1,Typelist<T2,T3> >,Unit> : public PairRootWithAll<Typelist<T1,T3>,Unit>, 
							    public Unit<T1,T2>,
							    public Unit<T2,T1>
{
  public:
  PairRootWithAll() : Unit<T1,T2>(boost::function<T1 (const T2&)>(Demoter<T1,T2>))  ,
                          Unit<T2,T1>(boost::function<T2 (const T1&)>(Promoter<T2,T1>)) {}
};


template <class T1,class T2,template <class,class> class Unit> 
class PairRootWithAll<Typelist<T1,Typelist<T2,NullType> >,Unit> :  public Unit<T1,T2> , 
								   public Unit<T2,T1>
{
public:
  PairRootWithAll() : Unit<T1,T2>(boost::function<T1 (const T2&)>(Demoter<T1,T2>)) ,
		      Unit<T2,T1>(boost::function<T2 (const T1&)>(Promoter<T2,T1>)) {}
};



template <class TList,template <class,class> class Unit> 
class Converter : public PairRootWithAllAndBranch<TList,Unit> , public PairWithSelf<TList,Unit> {};


template <class T1, class T2>
class ConverterImpl
{
private:
  boost::function<T1 (const T2&)> m_method;  
public:
  T1 Method(const T2& in) { return m_method(in); }
  ConverterImpl()   
  {
    m_method = boost::function<T1 (const T2&)>(Retainer<T1,T2>);
  }
  ConverterImpl(boost::function<T1 (const T2&)> method) 
  {
    m_method = method;
  }

};

template <class out_type,class in_type>
out_type Convert(const in_type& in)
{
Converter<MatTypes,ConverterImpl>& c = CGenericSingleton<Converter<MatTypes,ConverterImpl> >::Instance();
return c.ConverterImpl<out_type,in_type>::Method(in);
}


#endif

