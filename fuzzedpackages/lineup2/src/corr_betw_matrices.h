#ifndef CORR_BETW_MATRICES_H
#define CORR_BETW_MATRICES_H

// correlations between columns of two matrices

#include <Rcpp.h>


// calculate correlation between columns of x and corresponding columns of y
Rcpp::NumericVector corr_betw_matrices_paired(const Rcpp::NumericMatrix& x,
                                              const Rcpp::NumericMatrix& y);

// for each column of left matrix, find the column in the right matrix
// with the highest correlation
Rcpp::List corr_betw_matrices_unpaired_bestright(const Rcpp::NumericMatrix& x,
                                                 const Rcpp::NumericMatrix& y);

// return correlations between column of left matrix and column of right matrix
// that exceed corr_threshold
Rcpp::List corr_betw_matrices_unpaired_bestpairs(const Rcpp::NumericMatrix& x,
                                                 const Rcpp::NumericMatrix& y,
                                                 const double corr_threshold);

// calculate full set of correlations between columns of x and columns of y
Rcpp::NumericMatrix corr_betw_matrices_unpaired_all(const Rcpp::NumericMatrix& x,
                                                    const Rcpp::NumericMatrix& y);

#endif // CORR_BETW_MATRICES_H
