// standardize the columns of a matrix
#ifndef FSCALE_H
#define FSCALE_H

#include <Rcpp.h>

// fscale: standardize the vector x using only the values where
//         both x and y are not missing
Rcpp::NumericVector fscale(const Rcpp::NumericVector& x,
                           const Rcpp::NumericVector& y);

#endif // FSCALE_H
