/*
 * Author: Andreas Alfons
 *         KU Leuven
 */

#ifndef _ccaPP_COR_H
#define _ccaPP_COR_H

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include "utils.h"
#include "fastCorKendall.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// functions to export to R
RcppExport SEXP R_corPearson(SEXP R_x, SEXP R_y);
RcppExport SEXP R_corSpearman(SEXP R_x, SEXP R_y, SEXP R_consistent);
RcppExport SEXP R_corKendall(SEXP R_x, SEXP R_y, SEXP R_consistent);
RcppExport SEXP R_corQuadrant(SEXP R_x, SEXP R_y, SEXP R_consistent);
RcppExport SEXP R_corM(SEXP R_x, SEXP R_y, SEXP R_prob,
		SEXP initial, SEXP R_tol);

// functions to be used within C++
double corPearson(const vec& x, const vec& y);
double corSpearman(const vec& x, const vec& y, const bool& consistent);
double corKendall(const vec& x, const vec& y, const bool& consistent);
double corQuadrant(const vec& x, const vec& y, const bool& consistent);
double corM(const vec& x, const vec& y, const double& prob,
		const string& initial, const double& tol);

#endif
