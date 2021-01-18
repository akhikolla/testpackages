#ifndef _basics_BASICS_H
#define _basics_BASICS_H

#include <Rcpp.h>
#include <vector>

// extern "C" tells the compiler that this function should be accessible from anywhere
// input an SEXP, return an int
// These can be found only from WITHIN C++

extern "C" float maxval(SEXP maxent01sub) ;
extern "C" float minval(SEXP maxent01sub) ;
extern "C" bool all_ints_equal(std::vector<int> myvector1, std::vector<int> myvector2) ;
extern "C" bool any_ints_equal(std::vector<int> myvector1, std::vector<int> myvector2) ;
extern "C" bool all_ints_found(std::vector<int> myvector1, std::vector<int> myvector2) ;

// This causes this CRAN check warning in CRAN:
// 'merge_int_vectors' has C-linkage specified, but returns 
// user-defined type 'std::vector<int>' which is incompatible with C
//extern "C" std::vector<int> merge_int_vectors(std::vector<int> myvector1, std::vector<int> myvector2) ;
std::vector<int> merge_int_vectors(std::vector<int> myvector1, std::vector<int> myvector2) ;


extern "C" void moncombn(int* combmat, int* n, int* m);
extern "C" void moncombn_zerostart(int* combmat, int* n, int* m);
extern "C" int nChoosek( int n, int k );

extern "C" int get_missing_int(std::vector<int> myvector1, std::vector<int> myvector2);
extern "C" void printvec (std::vector<int> myvector1);
extern "C" void printBoolVec (std::vector<bool> myvector1);

//RcppExport int maxval(SEXP maxent01sub) ;
//RcppExport int minval(SEXP maxent01sub) ;



// RcppExport, on the other hand, can be found from R
/*
 * note : RcppExport is an alias to `extern "C"` defined by Rcpp.
 *
 * It gives C calling convention to the convolve function so that 
 * it can be called from .Call in R. Otherwise, the C++ compiler mangles the 
 * name of the function and .Call can't find it.
 *
 * It is only useful to use RcppExport when the function is intended to be called
 * by .Call. See the thread http://thread.gmane.org/gmane.comp.lang.r.rcpp/649/focus=672
 * on Rcpp-devel for a misuse of RcppExport
 */
RcppExport SEXP mult2probvect(SEXP leftprobs, SEXP rightprobs) ;
RcppExport SEXP convolve3cpp(SEXP a, SEXP b) ;

#endif
