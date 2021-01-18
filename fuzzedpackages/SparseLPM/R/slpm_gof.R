slpm_gof <- function(var_pars)
{
  M <- nrow(var_pars$alpha_u_tilde)
  N <- nrow(var_pars$alpha_v_tilde)
  K <- ncol(var_pars$alpha_u_tilde)
  theta <- matrix(0,M,N)
  for (i in 1:M) for (j in 1:N)
  {
    k <- which.max(var_pars$lambda_tilde[i,j,])
    theta[i,j] = 1 / (var_pars$alpha_u_tilde[i,k] - var_pars$alpha_v_tilde[j,k])^2
  }
  theta
}
