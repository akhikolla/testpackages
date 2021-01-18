#ifndef SURVCD_h
#define SURVCD_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

arma::vec RunNetSurv(arma::mat& xc, arma::mat& xg, arma::vec& y, double lamb1, double lamb2, arma::vec bc0, arma::vec bg0, double r, 
                        arma::mat& a, int p, int pc, bool robust);
arma::vec RunMCPSurv(arma::mat xc, arma::mat xg, arma::vec y, double lamb1, arma::vec bc, arma::vec bg, double r, int p, int pc, bool robust);
arma::vec RunLassoSurv(arma::mat xc, arma::mat xg, arma::vec y, double lamb1, arma::vec bc, arma::vec bg, int p, int pc, bool robust);

#endif
