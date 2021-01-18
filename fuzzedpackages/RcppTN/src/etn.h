# ifndef ETN_H
# define ETN_H

# include <Rcpp.h>

void etn(Rcpp::NumericVector &Mean,
         Rcpp::NumericVector &Sd,
         Rcpp::NumericVector &Low,
         Rcpp::NumericVector &High,
         Rcpp::NumericVector &Exps
    ) ;

# endif
