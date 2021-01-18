#ifndef HMB_CPP_HMB_H
#define HMB_CPP_HMB_H

Rcpp::List cpp_hmb(
  const arma::vec& Y_S,
  const arma::mat& X_S,
  const arma::mat& X_Sa,
  const arma::mat& Z_Sa,
  const arma::mat& Z_U);

#endif
