#ifndef DEBUG_H_
#define DEBUG_H_

//#include<Rcpp.h>
#include <RcppArmadillo.h>

namespace debug {
  void print(const arma::colvec& input)
  {
    Rcpp::Rcout << input << std::endl;
  }
  
  void print(const arma::uvec& input)
  {
    Rcpp::Rcout << input << std::endl;
  }
  
  void print(const arma::mat& input)
  {
    Rcpp::Rcout << input << std::endl;
  }
}

#endif // DEBUG_H_
