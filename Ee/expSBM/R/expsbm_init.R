blur <- function(mat, value = 0.5)
{
  K <- ncol(mat)
  v <- 1/K + (1 - 1/K) * value# the actual value is somewhere between a completely random and a hard clustering
  res <- mat
  if (ncol(mat) > 1) for (i in 1:nrow(mat)) for (j in 1:ncol(mat)) 
  {
    if (mat[i,j] == 1) res[i,j] = v
    else res[i,j] = (1 - v) / (K-1)
  }
  res
}

expSBM_init <- function(edgelist, K, method = c("random", "SBM_binary", "SBM_gaussian", "spectral"), sbm_density = 0.5, blur_value = 1)
{
  method_arg = match.arg(method)
  N <- max( c(edgelist[,1],edgelist[,2]) )
  res <- NULL
  if (method_arg == "random") res = sample(1:K, N, T)
  if (method_arg == "SBM_binary") res = expSBM_init_sbm_binary(edgelist, K, sbm_density)
  if (method_arg == "SBM_gaussian") res = expSBM_init_sbm_gaussian(edgelist, K)
  if (method_arg == "spectral") res = expSBM_init_spectral(edgelist, K)
  blur(mat = mclust::unmap(res), value = blur_value)
}

aggregate_values <- function(edgelist)
{
  N <- max(edgelist[,1:2])
  aggr_mu <- aggr_nu <- matrix(0, N, N)
  for (index in 1:nrow(edgelist))
  {
    i <- edgelist[index,1]
    j <- edgelist[index,2]
    if (edgelist[index,3] == 1) aggr_mu[i,j] = aggr_mu[i,j] + edgelist[index,4]
    if (edgelist[index,3] == 0) aggr_nu[i,j] = aggr_nu[i,j] + edgelist[index,4]
  }
  list(aggr_value_A0 = aggr_nu, aggr_value_A1 = aggr_mu,
       log_aggr_value_A0 = ifelse(aggr_nu != 0, log(aggr_nu), 0), 
       log_aggr_value_A1 = ifelse(aggr_mu != 0, log(aggr_mu), 0))
}

expSBM_init_sbm_binary <- function(edgelist, K, sbm_density = 0.5)
{
  values <- aggregate_values(edgelist)$aggr_value_A1
  adj <- ifelse(values > quantile(values, probs = 1 - sbm_density), 1, 0)
  res <- BM_bernoulli(membership_type = "SBM", adj = adj, verbosity = 0, plotting = "", explore_min = K, explore_max = K)
  res$estimate()
  soft_allocations <- res$memberships[[which.max(res$ICL)]]$Z
  hard_allocations <- apply(soft_allocations, 1, which.max)
  hard_allocations
}

expSBM_init_sbm_gaussian <- function(edgelist, K)
{
  values <- aggregate_values(edgelist)
  res <- BM_gaussian_multivariate(membership_type = "SBM", adj = list(values$log_aggr_value_A0, values$log_aggr_value_A1), verbosity = 0, plotting = "", explore_min = K, explore_max = K)
  res$estimate()
  soft_allocations <- res$memberships[[which.max(res$ICL)]]$Z
  hard_allocations <- apply(soft_allocations, 1, which.max)
  hard_allocations
}

expSBM_init_spectral <- function(edgelist, K)
{
  N <- max( c(edgelist[,1],edgelist[,2]) )
  matu <- matrix(0, N, N)
  matv <- matrix(0, N, N)
  for (index in 1:nrow(edgelist))
  {
    i <- edgelist[index,1]
    j <- edgelist[index,2]
    A <- edgelist[index,3]
    value <- edgelist[index,4]
    if (A == 0) 
    {
      matv[i,j] = matv[i,j] + value
      matv[j,i] = matv[j,i] + value
    }
    if (A == 1) 
    {
      matu[i,j] = matu[i,j] + value
      matu[j,i] = matu[j,i] + value
    }
  }
  mat <- matu / (matu + matv)
  diag(mat) = 0
  d <- apply(mat, 1, sum) + 0.00001
  L <- diag(d) - mat
  L <- diag(d^-0.5) %*% L %*% diag(d^-0.5)
  ev <- eigen(L, symmetric = TRUE)
  xx <- ev$vectors[,(ncol(ev$vectors)-(K-1)):ncol(ev$vectors)]
  fit <- Mclust(data = xx, G = K, verbose = F)
  fit$classification
}


