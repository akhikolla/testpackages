slpm_gen <- function(M, N, K, hyper_pars = NULL)
{
  delta <- a <- b <- rep(1,K)
  if (!is.null(hyper_pars)) 
  {
    delta = hyper_pars$delta
    a = hyper_pars$a_gamma
    b = hyper_pars$b_gamma
  }
  gamma <- rep(NA,K)
  for (k in 1:K) gamma[k] = rgamma(n = 1, shape = a[k], rate = b[k])
  U <- matrix(NA, M, K)
  V <- matrix(NA, N, K)
  for (k in 1:K)
  {
    U[,k] = rnorm(n = M, mean = 0, sd = sqrt(1/gamma[k]))
    V[,k] = rnorm(n = N, mean = 0, sd = sqrt(1/gamma[k]))
  }
  lambda <- as.numeric(rdirichlet(n = 1, alpha = delta))
  adj <- matrix(0,M,N)
  for (i in 1:M) for (j in 1:N) for (k in 1:K)
  {
    d <- (U[i,k] - V[j,k])^2
    adj[i,j] = adj[i,j] + lambda[k] / d
  }
  list(adj = adj, U = U, V = V, lambda = lambda, gamma = gamma)
}

