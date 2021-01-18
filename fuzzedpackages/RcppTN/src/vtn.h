# ifndef VTN_H
# define VTN_H

# include <Rcpp.h>

void vtn(Rcpp::NumericVector &Mean,
         Rcpp::NumericVector &Sd,
         Rcpp::NumericVector &Low,
         Rcpp::NumericVector &High,
         Rcpp::NumericVector &Vars
    ) ;

# endif
