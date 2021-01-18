#ifndef HMB_HF_MUVAR
#define HMB_HF_MUVAR

void MuVar(
  long double* return_var,
  const arma::mat& Predictors,
  const arma::mat& muCov,
  const arma::mat& muVarCov
);

#endif
