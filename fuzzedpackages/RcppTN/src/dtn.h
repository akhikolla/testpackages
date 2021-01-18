// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

/// Vectorized Truncated Normal Density Function

# include <Rcpp.h>

# ifndef DTN_HPP
# define DTN_HPP

void dtn(Rcpp::NumericVector &X,
         Rcpp::NumericVector &Mean,
         Rcpp::NumericVector &Sd,
         Rcpp::NumericVector &Low,
         Rcpp::NumericVector &High,
         Rcpp::NumericVector &Dens
         ) ;

# endif
