#ifndef HMB_CPP_GHMB_H
#define HMB_CPP_GHMB_H

Rcpp::List cpp_ghmb(
  const arma::vec& y_S,
  const arma::mat& X_S,
  const arma::mat& X_Sa,
  const arma::mat& Z_Sa,
  const arma::mat& Z_U,
  const arma::mat& Omega_S,
  const arma::mat& Omega_Sa,
  const arma::mat& Sigma_Sa);

#endif
