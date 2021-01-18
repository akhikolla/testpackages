# ifndef ENTTN_HPP
# define ENTTN_HPP

# include <Rcpp.h>

void enttn(Rcpp::NumericVector &Mean,
           Rcpp::NumericVector &Sd,
           Rcpp::NumericVector &Low,
           Rcpp::NumericVector &High,
           Rcpp::NumericVector &Ents
           ) ;

# endif
