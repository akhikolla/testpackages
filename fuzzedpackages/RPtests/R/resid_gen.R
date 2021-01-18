# Generate lasso residuals
# Returns n by B matrix of scaled residuals
# Optional error_mat supplies n by B matrix of errors
# If signal is specified, beta_out need not be
# Returns an n by B matrix of scaled residuals
resid_gen_lasso <- function(x, beta_out, sigma_est, B=249L, rand_gen = rnorm, mc.cores=1L,
                            error_mat, lam0=NULL, signal, ...) {
  n <- nrow(x)
  if (missing(signal)) signal <- as.numeric(x %*% beta_out$beta)

  if (missing(error_mat)) error_mat <- sigma_est * matrix(rand_gen(n*B), n, B)

  rep_fun <- function(b) {
    y_star <- as.numeric(signal + error_mat[, b])
    return(resid_lasso(x, y_star, lam0=lam0, ...))
  }

  output <- parallel::mclapply(1:B, rep_fun, mc.cores=mc.cores)
  output <- matrix(unlist(output), nrow=n, ncol=B)
  return(output)
}

# Returns scaled residuals
# If beta_out is supplied it is used to calculate the residauls
# beta_out can be either a beta vector or a list with beta component
resid_lasso <- function(x, y, beta_out, lam0=NULL, rescale=TRUE, ...) {
  y <- y - mean(y)
  if (missing(beta_out)) {
    beta_hat <- sqrt_lasso(x, y, lam0 = lam0, ...)
  } else {
    if (is.list(beta_out)) {
      beta_hat <- beta_out$beta
    } else {
      beta_hat <- beta_out
    }
  }
  if (rescale) {
    return(Scale(y - x %*% beta_hat))
  } else {
    return(as.numeric(y - x %*% beta_hat))
  }
}

comp_sigma_est <- function(x, y, beta_star) {
  # This will typically be an overestimate of sigma when b_star is biased
  sqrt(mean((y - mean(y) - x %*% beta_star)^2))
}

# Generate OLS residuals
# Returns n by B matrix of scaled residuals
# proj (I-P) is optional
# An intercept term is added
resid_gen_ols <- function(x, proj, B=250, rand_gen = rnorm, error_mat = NULL) {
  n <- nrow(x)
  if (missing(proj)) {
    x <- cbind(1, x) # add intercept
    u <- svd(x, nv=0)$u
    proj <- diag(rep(1, n)) - u%*%t(u) # I - P
  }
  n <- nrow(proj)
  error_mat <- if(is.null(error_mat)) matrix(rand_gen(B*n), ncol=B) else error_mat
  resid_sim <- scale(proj %*% error_mat)*sqrt(n/(n-1))
  return(resid_sim)
}

# OLS residuals
# returns a list with components
# resid - scaled residuals, proj - (I-P)
# An intercept term is added
# y can potentially be a matrix
resid_ols <- function(x, y, incl_proj = FALSE) {
  n <- nrow(x)
  x <- cbind(1, x) # add intercept
  u <- svd(x, nv=0)$u
  proj <- diag(rep(1, n)) - u%*%t(u) # I - P
  if (incl_proj) {
    if (is.matrix(y)) {
      return(list("resid"=scale(proj %*% y) * sqrt((n-1)/2), "proj"=proj))
    } else {
      return(list("resid"=Scale(as.numeric(proj %*% y)), "proj"=proj))
    }
  } else {
    if (is.matrix(y)) {
      return("resid"=scale(proj %*% y) * sqrt((n-1)/2))
    } else {
      return(Scale(as.numeric(proj %*% y)))
    }
  }
}

Scale <- function(resid) {
  resid <- as.numeric(resid - mean(resid))
  resid / sqrt(mean((resid)^2))
}
