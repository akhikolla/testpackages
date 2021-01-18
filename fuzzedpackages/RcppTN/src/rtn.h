# ifndef RTN_H
# define RTN_H

# include <Rcpp.h>

void rtn(Rcpp::NumericVector &Mean,
         Rcpp::NumericVector &Sd,
         Rcpp::NumericVector &Low,
         Rcpp::NumericVector &High,
         Rcpp::NumericVector &Draws
    ) ;

# endif
