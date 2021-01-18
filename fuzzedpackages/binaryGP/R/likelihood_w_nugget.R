likelihood_w_nugget <- function(omega, M.mx, X, w, eta_tilde, T_val, corr, orthogonalGP){

  sigma_hat <- omega[1]
  sigma2_hat <- sigma_hat^2
  rho_hat <- omega[2:(length(omega)-1)]
  nugget_hat <- omega[length(omega)]
  R.mx <- corr_matrix(X, rho_hat, corr) + diag(nugget_hat/sigma2_hat, nrow(X))
  if(orthogonalGP) R.mx <- R.mx - orthogonalize(X, rho_hat, corr)

  out <- likelihood_fun(M.mx, X, w, eta_tilde, R.mx, sigma2_hat, T_val)

  return(out)
}
