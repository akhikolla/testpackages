/*
  This is the global head file. Every file in ROPTLIB will include this file.

  ---- WH
*/

#ifndef DEF_H
#define DEF_H

// macros for going between fortran and C++, assumes a 
// m x n matrix being accesse
#define ROWM(x, y, n) (x*n + y) //row-major
#define COLM(x, y, m) (y*m + x) //column-major

//#define MATLAB_MEX_FILE//For debug---

// macros for viewing complier macros
#define XSTR(x) STR(x)
#define STR(x) #x


/* 
   If all the test files are included in a project, then only uncomment one of them to specify which test problem is run.
*/

//#define TESTEUCFRECHETMEAN
//#define TESTEUCQUADRATIC
//#define TESTPRODUCT
//#define TESTSPHERERAYQUO
//#define TESTSTIEBROCKETT
//#define TESTGRASSRQ
//#define TESTCSO
//#define TESTSTIESOFTICA
//#define TESTTESTSPARSEPCA
//#define TESTWEIGHTEDLOWRANK
//#define TESTELASTICCURVESRO
//#define TESTMYMATRIX
//#define TESTSPDMEAN
//#define TESTPRESHAPEPATHSTRAIGHTEN
//#define TESTSPDTENSORDL

//#define TESTSIMPLEEXAMPLE
//#define TESTPRODUCTEXAMPLE

#include <cmath>
/*If ROPTLIB is not compiled in Matlab, then the following wrapper functions of blas and lapack
  are included.*/
// #ifndef MATLAB_MEX_FILE
#if !defined(MATLAB_MEX_FILE) && !defined(R_BUILD)

// blas and lapack related
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <cstdio>

#define integer int
#define doublecomplex Rcomplex

#define zdotc_ zdotc
#define zgegs_ zgegs
#define ztgsyl_ ztgsyl
#define zgees_ zgees
#define ztrtrs_ ztrtrs
#define zgemm_ zgemm
#define zscal_ zscal
#define zgeqp3_ zgeqp3
#define zunmqr_ zunmqr


#include <dgemm.h>
#include <dgetrf.h>
#include <dgetrs.h>
#include <dgemv.h>
#include <dcopy.h>
#include <ddot.h>
#include <dscal.h>
#include <daxpy.h>
#include <dger.h>
#include <dgeqp3.h>
#include <dorgqr.h>
#include <dormqr.h>
#include <dtrsm.h>
#include <dlarfx.h>
#include <ddot.h>
#include <dgesdd.h>
#include <dgesvd.h>
#include <dsymv.h>
#include <dgetri.h>
#include <dlapmt.h>
#include <dgees.h>
#include <dnrm2.h>
#include <dgesv.h>
#include <dsyevx.h>
#include <dlamch.h>
#include <dpotrf.h>
#include <dtrtrs.h>
#include <dsyevd.h>
#include <dsyevr.h>
#include <dsyev.h>

#include <zdotc.h>
#include <zgegs.h>
#include <ztgsyl.h>
#include <zgees.h>
#include <ztrtrs.h>
#include <zgemm.h>
#include <zscal.h>
#include <zgeqp3.h>
#include <zunmqr.h>
#endif // end of ifndef MATLAB_MEX_FILE


#if defined(R_BUILD)

// Much of the ROPTLIB code needs BLAS and Lapack, so let them all use the
// R implementation.
//
// Many ROPTLIB files use std::cout or printf. CRAN people do not like this, so
// here we redefine printf to Rprintf, and provide a generic-looking OUTSTREAM
// which really uses Rcpp::Rcout. Instead of including ostream, R.h, and
// Rstreambuf.h, we could simply include Rcpp.h; however, this causes very slow
// compilation and a bloated .so file at the end. This seems to be caused by
// lots of symbol information being packed into the .o files. (This can be
// verified by using the strip command on the .so file).
//
// Even without including Rcpp.h from each source file, the .so at the end is
// quite large. The CRAN people do not seem to have a good solution for this
// currently:
// http://r.789695.n4.nabble.com/compile-c-code-in-an-R-package-without-g-td4725700.html

#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <ostream>
#include <R.h>
#include <Rcpp/iostream/Rstreambuf.h>

#define integer int
#define doublecomplex Rcomplex
#define OUTSTREAM Rcpp::Rcout
//#define printf Rprintf

#endif

#ifdef _WIN64 // The following code is compiled only when this library is compiled in Windows (64-bit only)
/*If the code is compile under DEBUG mode, then test whether there is memory leakage or not*/
#ifdef _DEBUG
#define DEBUG_CLIENTBLOCK   new( _CLIENT_BLOCK, __FILE__, __LINE__)
/*Use my code to help checking memory leakage. One has to define a global variable:
  std::map<integer *, integer> *CheckMemoryDeleted;
  before running the code.
*/
//#define CHECKMEMORYDELETED
#else
#define DEBUG_CLIENTBLOCK
#endif

/*This is used for checking the memory leakage in windows system*/
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>

/*This is used for checking the memory leakage in windows system if the code is run in DEBUG mode*/
#ifdef _DEBUG
#define new DEBUG_CLIENTBLOCK
#endif

#elif _WIN32 // The following code is compiled only when this library is compiled in Windows (both 32-bit and 64-bit only)
//define something for Windows (32-bit and 64-bit, this part is common)
#elif __APPLE__ // The following code is compiled only when this library is compiled in MAC
#include "TargetConditionals.h"
#if TARGET_IPHONE_SIMULATOR
// iOS Simulator
#elif TARGET_OS_IPHONE
// iOS device
#elif TARGET_OS_MAC
// Other kinds of Mac OS
#else
// Unsupported platform
#endif
#elif __linux// The following code is compiled only when this library is compiled in Linux system
// linux
#ifdef __GNUC__

#endif // end of __GNUC__


#elif __unix // all unices not caught above
// Unix
#elif __posix
// POSIX
#endif // end of checking platforms


#if __linux || __unix || __posix

#if __cplusplus < 201103L

#if __clang__ 

#if __has_feature(cxx_nullptr)
#else
#define NO_NULL_PTR
#endif

#else
#define NO_NULL_PTR
#endif

#endif


#ifdef NO_NULL_PTR
const class {
 public:
  template<class T> // convertible to any type
    operator T*(void) const // of null non-member
    {
      return 0;
    } // pointer...
  template<class C, class T> // or any type of null
    operator T C::*(void) const // member pointer...
  {
    return 0;
  }
 private:
  void operator&(void) const; // whose address can't be taken
} nullptr = {};
#endif

#endif

/*
// check compiler for c++11 compliance.  If this compliance doesn't exist, then emulate nullptr per Scott Meyers c++ book
#if __cplusplus < 201103L

// if using clang then nullptr exists for version 3 and greater
#if __clang__ && __clang_major__ >= 3
#pragma message "Using clang version 3 or higher, assuming nullptr defined."
#pragma message "The value of cplusplus: " XSTR(__cplusplus)
#pragma message "The version of clang: " XSTR(__clang_version__)

#else

const class {
 public:
  template<class T> // convertible to any type
    operator T*(void) const // of null non-member
    {
      return 0;
    } // pointer...
  template<class C, class T> // or any type of null
    operator T C::*(void) const // member pointer...
  {
    return 0;
  }
 private:
  void operator&(void) const; // whose address can't be taken
} nullptr = {};
#endif
#endif
#endif
*/


/*If ROPTLIB is compiled in Matlab, then removing the underscore to make the wrappers consistant.*/
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#include "blas.h"
#include "lapack.h"
#define integer ptrdiff_t
#define dgemm_ dgemm
#define dgetrf_ dgetrf
#define dgetrs_ dgetrs
#define dgemv_ dgemv
#define dcopy_ dcopy
#define ddot_ ddot
#define dscal_ dscal
#define daxpy_ daxpy
#define dger_ dger
#define dgeqp3_ dgeqp3
#define dorgqr_ dorgqr
#define dormqr_ dormqr
#define dtrsm_ dtrsm
#define dlarfx_ dlarfx
#define dgesdd_ dgesdd
#define dgesvd_ dgesvd
#define dsymv_ dsymv
#define dgetri_ dgetri
#define dgees_ dgees
#define dnrm2_ dnrm2
#define dgesv_ dgesv
#define dsyevx_ dsyevx
#define dlamch_ dlamch
#define dpotrf_ dpotrf
#define dtrtrs_ dtrtrs
#define dsyevd_ dsyevd
#define dsyevr_ dsyevr
#define dsyev_ dsyev

#define zdotc_ zdotc
#define zgegs_ zgegs
#define ztgsyl_ ztgsyl
#define zgees_ zgees
#define ztrtrs_ ztrtrs
#define zgemm_ zgemm
#define zscal_ zscal
#define zgeqp3_ zgeqp3
#define zunmqr_ zunmqr
#endif // end of ifdef MATLAB_MEX_FILE

/*Help to debug the code*/
// #include "ForDebug.h"

/*For obtaining the lower bound, upper bound of numbers of double precision*/
#include <climits>
#include <limits>
/*For obtain the computational time*/
#include "Timer.h"

#include <map>
#include <string>

typedef std::map<std::string, double> PARAMSMAP;

/*Define the number PI*/
#ifndef PI
#define PI 3.14159265358979323846264
#endif

#endif // end of DEF_H
