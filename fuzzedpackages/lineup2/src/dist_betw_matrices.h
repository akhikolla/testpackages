#ifndef DIST_BETW_MATRICES_H
#define DIST_BETW_MATRICES_H

#include <Rcpp.h>

// calculate the distance between columns of matrix x and columns of matrix y
// d(i,j) = sqrt{ sum_k [x(k,i) - y(k,j)]^2 }
Rcpp::NumericMatrix rmsd_betw_matrices(const Rcpp::NumericMatrix& x,
                                       const Rcpp::NumericMatrix& y);


// like rmsd_betw_matrices but using the average absolute difference
Rcpp::NumericMatrix mad_betw_matrices(const Rcpp::NumericMatrix& x,
                                      const Rcpp::NumericMatrix& y);

#endif // DIST_BETW_MATRICES_H
