expSBM_select <- function(K_max, N, edgelist, method = "SBM_gaussian", directed = F, trunc = T, tol = 0.001, n_iter_max = 100, init_blur_value = 1, verbose = F)
{
  res <- list()
  res$fitted_models <- list()
  elbo_values <- icl_values <- rep(NA,K_max)
  for (k in 1:K_max)
  {
    lambda_init <- as.numeric(rdirichlet(1, rep(1,k)))
    Z_init_collapsed <- expSBM_init(edgelist = edgelist, K = k, method = method, blur_value = init_blur_value)
    K_collapsed <- sum(colSums(Z_init_collapsed) > 0)
    Z_init <- matrix(0,N,k)
    if (k > 1) Z_init[,1:K_collapsed] = Z_init_collapsed else Z_init = Z_init_collapsed
    mu_init <- matrix(rgamma(k*k,1,1),k,k)
    if (!directed) for (g in 1:k) for (h in 1:k) if (g > h) mu_init[g,h] = mu_init[h,g]
    nu_init <- matrix(rgamma(k*k,1,1),k,k)
    if (!directed) for (g in 1:k) for (h in 1:k) if (g > h) nu_init[g,h] = nu_init[h,g]
    res$fitted_models[[k]] <- output <- expSBM_EM(N, edgelist, Z_init, lambda_init, mu_init, nu_init, directed, trunc, tol, n_iter_max, verbose)
    elbo_values[k] = output$elbo_values[length(output$elbo_values)]
    Z_map_vec <- apply(output$Z_star, 1, which.max)
    Z_map <- matrix(0,N,k)
    for (i in 1:N) Z_map[i,Z_map_vec[i]] = 1
    marginal_likelihood_contribution <- expSBM_ELBO(N, edgelist, Z_map, output$lambda_star, output$mu_star, output$nu_star, directed, trunc, verbose)$elbo_value
    prior_contribution <- sum(Z_map %*% log(ifelse(output$lambda_star > 0, output$lambda_star, 1)))
    if (directed) icl_values[k] = marginal_likelihood_contribution + prior_contribution - k^2 * log(nrow(edgelist)) - 0.5 * (k-1) * log(N)
    if (!directed) icl_values[k] = marginal_likelihood_contribution + prior_contribution - (k*(k+1)/2) * log(nrow(edgelist)) - 0.5 * (k-1) * log(N)
  }
  res$icl_values <- icl_values
  res$K_star <- which.max(icl_values)
  res$best_model <- res$fitted_models[[res$K_star]]
  res
}
