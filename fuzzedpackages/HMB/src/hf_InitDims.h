#ifndef HMB_HF_INITDIMS
#define HMB_HF_INITDIMS

hmbdims InitDims (
  const arma::vec& Y_S,
  const arma::mat& X_S,
  const arma::mat& X_Sa,
  const arma::mat& Z_Sa,
  const arma::mat& Z_U,
  const bool hasIntercepts = TRUE
);

hmbdims InitDims_Z_S (
  const arma::vec& Y_S,
  const arma::mat& X_S,
  const arma::mat& X_Sa,
  const arma::mat& Z_S,
  const arma::mat& Z_Sa,
  const arma::mat& Z_U,
  const bool hasIntercepts = TRUE
);

#endif
