#pragma once
#include <Rcpp.h>

using namespace Rcpp;

List tfCox_onelambda(int ord, double alpha, double lambda, IntegerVector discrete, int ndis,
                     double stepSize, NumericMatrix X, NumericMatrix init_theta, 
                     NumericVector Time, IntegerVector status, IntegerVector indx_time, 
                     IntegerVector tie,  int n, int p, IntegerMatrix Perm, 
                     IntegerMatrix Rank, IntegerVector thin, NumericVector vec_xtol, 
                     double tol, int niter, int backtracking);

void tf_main(double *x, double *y, double *w, int n, int ord, double lambda, 
             int thinning, double x_tol, double *beta, double *weight, int *m); 
  
NumericVector gradient(NumericMatrix theta, int n, int p,
                       IntegerVector status, IntegerVector indx_time, 
                       IntegerVector tie);

double negloglik(NumericMatrix theta, int n, int p, 
                 IntegerVector status, IntegerVector indx_time, 
                 IntegerVector tie);
