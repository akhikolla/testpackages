likelihood <- function(omega, nugget, M.mx, X, w, eta_tilde, T_val, corr, orthogonalGP){

  sigma_hat <- omega[1]
  sigma2 <- sigma_hat^2
  rho_hat <- omega[2:length(omega)]
  R.mx <- corr_matrix(X, rho_hat, corr) + diag(nugget/sigma2, nrow(X))
  if(orthogonalGP) R.mx <- R.mx - orthogonalize(X, rho_hat, corr)

  out <- likelihood_fun(M.mx, X, w, eta_tilde, R.mx, sigma2, T_val)

  return(out)
}
