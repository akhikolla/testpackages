#ifndef _cpp_calc_anclikes_sp_cpp_calc_anclikes_sp_H
#define _cpp_calc_anclikes_sp_cpp_calc_anclikes_sp_H

#include <Rcpp.h>
#include "basics.h"

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
RcppExport SEXP cpp_calc_anclikes_sp() ;
RcppExport SEXP cpp_calc_anclikes_sp_COOprobs() ;
RcppExport SEXP cpp_calc_anclikes_sp_COOweights_faster() ;

RcppExport SEXP cpp_calc_anclikes_sp_rowsums() ;

RcppExport SEXP cpp_areas_list_to_states_list() ;

RcppExport SEXP moncombn();

RcppExport SEXP cpp_combn_zerostart();

RcppExport SEXP cpp_states_list_to_DEmat();
RcppExport SEXP cpp_states_list_to_DEmat_COO();

RcppExport SEXP cpp_calc_rowsums_for_COOweights_columnar();
RcppExport SEXP cpp_calc_anclikes_sp_using_COOprobs();
RcppExport SEXP cpp_calc_splitlikes_using_COOweights_columnar();

RcppExport SEXP cpp_calc_loglike_sp_fast();

#endif
