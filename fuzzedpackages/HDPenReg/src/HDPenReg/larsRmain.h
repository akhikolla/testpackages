#ifndef _LARSRMAIN_H
#define _LARSRMAIN_H

#include <RTKpp.h>

#include "lars/PathState.h"
#include "lars/functions.h"
#include "lars/Path.h"
#include "lars/Lars.h"
#include "lars/Cvlars.h"
#include "lars/Fusion.h"


#if defined(LARS_DEBUG) || defined(CVLARS_DEBUG)|| defined(FUSION_DEBUG)

template< class Type>
void printHead(STK::CArray<Type> const& A, std::string  const& name)
{
  stk_cout << "print: " << name << _T("\n";);
  stk_cout << name << _T(".isRef() =")                 << A.isRef()  << _T("\n");
  stk_cout << name << _T(".allocator().isRef() =")     << A.allocator().isRef() << _T("\n");
  stk_cout << name << _T(".rows() =")                  << A.rows()  << _T("\n");
  stk_cout << name << _T(".cols() =")                  << A.cols()  << _T("\n");
}

template< class Type>
void printHead(STK::Array2DUpperTriangular<Type> const& A, std::string  const& name)
{
  stk_cout << "print: " << name << _T("\n";);
  stk_cout << name << _T(".isRef() =")                 << A.isRef()  << _T("\n");
  stk_cout << name << _T(".allocator().isRef() =")     << A.allocator().isRef() << _T("\n");
  stk_cout << name << _T(".rangeCols().isRef() =")     << A.rangeCols().isRef() << _T("\n");
  stk_cout << name << _T(".rows() =")                  << A.rows()  << _T("\n");
  stk_cout << name << _T(".cols() =")                  << A.cols()  << _T("\n");
  stk_cout << name << _T(".allocator().range() =")     << A.allocator().range()  << _T("\n");
  stk_cout << name << _T(".rangeCols().range() =")     << A.rangeCols().range() << _T("\n");
  stk_cout << name << _T(".availableCols() =")         << A.availableCols()  << _T("\n");
  stk_cout << name << _T(".rangeCols() = ")            << A.rangeCols();
//  stk_cout << name << _T(".allocator() = ") << Array1D< Array1D<Type>* >(A.allocator(),A.cols(), true)  << _T("\n");
}
template< class Type>
void printHead(STK::Array2D<Type> const& A, std::string  const& name)
{
  stk_cout << "print: " << name << _T("\n";);
  stk_cout << name << _T(".isRef() =")                 << A.isRef()  << _T("\n");
  stk_cout << name << _T(".allocator().isRef() =")     << A.allocator().isRef() << _T("\n");
  stk_cout << name << _T(".rangeCols().isRef() =")     << A.rangeCols().isRef() << _T("\n");
  stk_cout << name << _T(".rows() =")                  << A.rows()  << _T("\n");
  stk_cout << name << _T(".cols() =")                  << A.cols()  << _T("\n");
  stk_cout << name << _T(".allocator().range() =")     << A.allocator().range()  << _T("\n");
  stk_cout << name << _T(".rangeCols().range() =")     << A.rangeCols().range() << _T("\n");
  stk_cout << name << _T(".availableCols() =")         << A.availableCols()  << _T("\n");
//  stk_cout << name << _T(".rangeCols() = ")            << A.rangeCols();
//  stk_cout << name << _T(".allocator() = ") << Array1D< Array1D<Type>* >(A.allocator(),A.cols(), true)  << _T("\n");
}

template< class TYPE>
void printHead(STK::Array2DVector<TYPE> const& A, std::string const& name)
{
  stk_cerr << "print: " << name << _T("\n";);
  stk_cerr << name << _T(".isRef() =")                 << A.isRef()  << _T("\n");
  stk_cerr << name << _T(".allocator().isRef() =")     << A.allocator().isRef() << _T("\n");
  stk_cerr << name << _T(".rangeCols().isRef() =")     << A.rangeCols().isRef() << _T("\n");
  stk_cerr << name << _T(".cols() =")                  << A.cols()  << _T("\n");
  stk_cerr << name << _T(".rows() =")                  << A.rows()  << _T("\n");
  stk_cerr << name << _T(".allocator().range() =")     << A.allocator().range()  << _T("\n");
  stk_cerr << name << _T(".rangeCols().range() =")     << A.rangeCols().range() << _T("\n");
  stk_cerr << name << _T(".availableCols() =")         << A.availableCols()  << _T("\n");
  stk_cerr << name << _T(".rangeCols() = ")            << A.rangeCols() << _T("\n");
//  stk_cerr << name << _T("=\n")                        << A.transpose() << _T("\n\n");
}
template< class TYPE>
void printHead(STK::CArrayVector<TYPE> const& A, std::string const& name)
{
  stk_cerr << "print: " << name << _T("\n";);
  stk_cerr << name << _T(".isRef() =")                 << A.isRef()  << _T("\n");
  stk_cerr << name << _T(".allocator().isRef() =")     << A.allocator().isRef() << _T("\n");
  stk_cerr << name << _T(".cols() =")                  << A.cols()  << _T("\n");
  stk_cerr << name << _T(".rows() =")                  << A.rows()  << _T("\n");
//  stk_cerr << name << _T("=\n")                        << A.transpose() << _T("\n\n");
}
template< class TYPE>
void printHead(STK::Array1D<TYPE> const& A, std::string const& name)
{
  stk_cerr << "print: " << name << _T("\n";);
  stk_cerr << name << _T(".isRef() =")                 << A.isRef()  << _T("\n");
  stk_cerr << name << _T(".cols() =")                  << A.cols()  << _T("\n");
  stk_cerr << name << _T(".rows() =")                  << A.rows()  << _T("\n");
}

template< class Derived>
void print( STK::ICArray<Derived> const& A, std::string const& name, std::string const& message)
{
  stk_cerr << _T("==>") << message<<_T("\n");
  printHead(A.asDerived(), name);
}
template< class Derived>
void print(STK::IArray2D<Derived> const& A, std::string const& name, std::string const& message)
{
  stk_cerr << _T("==>") << message<<_T("\n");
  printHead(A.asDerived(), name);
//  stk_cout << name << _T("=\n") << A << _T("\n\n");
}
template< class Type>
void print(STK::Array1D<Type> const& A, std::string const& name, std::string const& message)
{
  stk_cerr << _T("==>") << message<<_T("\n");
  printHead(A, name);
//  stk_cout << name << _T("=\n") << A << _T("\n\n");
}
#endif


RcppExport SEXP larsmain(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP maxStep, SEXP intercept, SEXP eps);
RcppExport SEXP fusionmain(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP maxStep, SEXP intercept, SEXP eps);
RcppExport SEXP cvlarsmain(SEXP data, SEXP response, SEXP nbIndiv, SEXP nbVar, SEXP maxStep, SEXP intercept, SEXP eps, SEXP nbFold, SEXP partition, SEXP index, SEXP mode);

#endif
