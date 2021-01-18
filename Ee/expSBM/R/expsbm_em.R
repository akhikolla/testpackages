expSBM_EM <- function(N, edgelist, Z, lambda, mu, nu, directed = F, trunc = T, tol = 0.001, n_iter_max = 100, verbose = F)
{
  cpp_expSBM_EM(N, cbind(edgelist[,1]-1,edgelist[,2]-1,edgelist[,3],edgelist[,4]), Z, lambda, mu, nu, directed, trunc, tol, n_iter_max, verbose)
}
