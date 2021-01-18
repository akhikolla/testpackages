#include <RcppArmadillo.h>
#include <iostream>
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::cube read_cube(std::string x){
  std::ifstream file(x);

  arma::cube C;
  C.load(file, arma::arma_binary);

  file.close();

  return C;
}

// This function is to load the raw_sam processed externally by C++
