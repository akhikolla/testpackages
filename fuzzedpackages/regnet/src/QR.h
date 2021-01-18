#ifndef QR_h
#define QR_h

#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

arma::vec QRWMR(const arma::mat& x, const arma::vec& y, arma::vec b);

#endif
