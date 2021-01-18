/*
 * Author: Andreas Alfons
 *         Erasmus University Rotterdam
 */

#ifndef _ccaPP_CCA_H
#define _ccaPP_CCA_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include "cor.h"

using namespace Rcpp;

// functions to export to R
RcppExport SEXP R_fastCor(SEXP R_x, SEXP R_y, SEXP R_method, SEXP R_control);	// for testing
RcppExport SEXP R_ccaPP(SEXP R_x, SEXP R_y, SEXP R_k, SEXP R_method, SEXP R_corControl,
  	SEXP R_algorithm, SEXP R_ppControl, SEXP R_standardize, SEXP R_fallback);

#endif
