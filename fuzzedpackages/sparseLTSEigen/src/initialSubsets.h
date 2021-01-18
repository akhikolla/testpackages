/*
 * Author: Andreas Alfons
 *         Erasmus University Rotterdam
 */

#ifndef _sparseLTSEigen_INITIALSUBSETS_H
#define _sparseLTSEigen_INITIALSUBSETS_H

#define EIGEN_NO_DEBUG

#include <sparseLTSEigen.h>
#include "fastLasso.h"
#include "utils.h"

// functions to export to R
RcppExport SEXP R_sparseSubsets(SEXP R_x, SEXP R_y, SEXP R_lambda, SEXP R_h, 
    SEXP R_subsets, SEXP R_normalize, SEXP R_intercept, SEXP R_eps, 
    SEXP R_useGram);

#endif
