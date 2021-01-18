#include <RcppArmadillo.h>
#include <fstream>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
int mosaic_sam_write(arma::cube A, char const* filename)
{
  A.save(filename, arma::arma_binary);
  return 0;
}

