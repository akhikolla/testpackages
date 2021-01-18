expSBM_ELBO <- function(N, edgelist, Z, lambda, mu, nu, directed = F, trunc = T, verbose = F)
{
  cpp_expSBM_ELBO(N, cbind(edgelist[,1]-1,edgelist[,2]-1,edgelist[,3],edgelist[,4]), Z, lambda, mu, nu, directed, trunc, verbose)
}

