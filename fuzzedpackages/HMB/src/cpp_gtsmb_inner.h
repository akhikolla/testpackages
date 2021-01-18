#ifndef HMB_CPP_GTSMB_INNER_H
#define HMB_CPP_GTSMB_INNER_H

Rcpp::List cpp_gtsmb_inner(
  const arma::vec& y_S,
  const arma::mat& X_S,
  const arma::mat& X_Sa,
  const arma::mat& Z_Sa,
  const arma::mat& Z_U,
  const arma::mat& Omega_S,
  const arma::cube& Phi_Sa
);

#endif
