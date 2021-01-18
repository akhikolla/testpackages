#ifndef RUNCONT_h
#define RUNCONT_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

arma::vec RunCont(arma::mat& xc, arma::mat& xg, arma::vec& y, double lamb1, double lamb2, 
					arma::vec bc0, arma::vec bg0, double r, arma::mat& a, int p, int pc, char method);
#endif
