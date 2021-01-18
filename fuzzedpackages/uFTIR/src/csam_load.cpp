#include <RcppArmadillo.h>
#include <fstream>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::cube csam_load(char const* filename)
{
  std::string fname(filename);
  arma::cube A;
  A.load(fname, arma::arma_binary);
  
  return A;
}
