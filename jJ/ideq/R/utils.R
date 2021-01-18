# Utility functions for dstm.R

# ifelse for arguments
`%else%` <- function(a, b) if ( is.null(a) || all(is.na(a)) ) b else a

# check for NAs
check.na <- function(x, x_name=deparse(substitute(x))) {
  if (any(is.na(x)))
    stop(paste(x_name, "may not contain NAs"))
}

# check dims
check.dim <- function(x, P, x_name=deparse(substitute(x)), dim_name="P") {
  if (any(dim(x) != P))
    stop(paste(x_name, "must be of dimension", dim_name, "by", dim_name))
}

# positive numeric
check.numeric.scalar <- function(x, x_name=deparse(substitute(x)), x_min=0,
                                 x_min_name=deparse(substitute(x_min)),
                                 strict_inequality=TRUE) {
  if (strict_inequality) {
    if (!is.numeric(x) || x <= x_min)
      stop(paste(x_name, "must be numeric >", x_min))
  } else {
    if (!is.numeric(x) || x < x_min)
      stop(paste(x_name, "must be numeric >=", x_min))
  }
}

# numeric vector
check.numeric.vector <- function(x, P, x_name=deparse(substitute(x)),
                                 dim_name="P") {
  if (!is.numeric(x) || !is.vector(x))
    stop(paste(x_name, "must be numeric vector"))
  if (length(x) != P)
    stop(paste(x_name, "must be of length", dim_name))
  check.na(x, x_name)
}

# numeric matrix
check.numeric.matrix <- function(x, P, x_name=deparse(substitute(x)),
                                 dim_name="P") {
  if (!is.numeric(x) || !is.matrix(x))
    stop(paste(x_name, "must be numeric matrix"))
  check.dim(x, P, x_name, dim_name)
  check.na(x, x_name)
}

# covariance matrix
is.cov.matrix <- function(x) {
  matrixcalc::is.positive.semi.definite(x) &&
    matrixcalc::is.symmetric.matrix(x)
}
check.cov.matrix <- function(x, P, x_name=deparse(substitute(x)), 
                             dim_name="P") {
  if (!is.cov.matrix(x))
    stop(paste(x_name, "must be valid variance-covariance matrix"))
  check.dim(x, P, x_name, dim_name)
  check.na(x, x_name)
}

# Cholesky inverse
chol_inv <- function(x) chol2inv(chol(x))

# Processing parameters common to dstm_eof() and dstm_ide()
process_common_params <- function(params, proc_error, P, sample_sigma2) {
  # Set m_0
  m_0 <- params[["m_0"]] %else% rep(0, P)
  check.numeric.vector(m_0, P)

  # Set C_0
  C_0 <- params[["C_0"]] %else% diag(P)
  check.cov.matrix(C_0, P)
  
  # Observation Error; creates alpha_sigma2, beta_sigma2, sigma2
  alpha_sigma2 <- beta_sigma2 <- sigma2 <- NA
  if (sample_sigma2) {
    alpha_sigma2 <- params[["alpha_sigma2"]] %else% 5
    beta_sigma2  <- params[["beta_sigma2"]]  %else% 4
    check.numeric.scalar(alpha_sigma2)
    check.numeric.scalar(beta_sigma2)
  } else {
    sigma2 <- params[["sigma2"]] %else% 1
    check.numeric.scalar(sigma2)
  }
  
  alpha_lambda <- beta_lambda <- df_W <- NA
  scale_W <- matrix()
  
  if (proc_error == "IW") {
    df_W <- params[["df_W"]] %else% (2*P)
    scale_W <- params[["scale_W"]] %else% diag(df_W, P)
    check.cov.matrix(scale_W, P)
    check.numeric.scalar(df_W, x_min=P-1)
  } else if (proc_error == "Discount") {
    alpha_lambda <- params[["alpha_lambda"]] %else% 100
    beta_lambda <- params[["beta_lambda"]] %else% 1
    check.numeric.scalar(alpha_lambda)
    check.numeric.scalar(beta_lambda)
    # NOTE: default prior for lambda corresponds to lambda=1/99
    # Or, in terms of the discount factor, delta=0.99
    
  } else {
    stop("proc_error must be \"IW\" or \"Discount\"")
  }
  
  list(
    m_0=m_0, C_0=C_0,
    alpha_sigma2=alpha_sigma2, beta_sigma2=beta_sigma2, sigma2=sigma2,
    scale_W=scale_W, df_W=df_W, alpha_lambda=alpha_lambda, beta_lambda=beta_lambda
  )
}

################################################################################
# Functions for dstm_ide()
scale_col <- function(x, L, x_range) {
  if (is.na(x_range)) x_range <- diff(range(x))
  x <- L/x_range * x
  x <- x + (L/2 - max(x))
  x
}

scale_all <- function(x, L, ranges=rep(NA, ncol(x))) {
  for ( i in seq(ncol(x)) ) x[, i] <- scale_col(x[, i], L, ranges[i])
  x
}

gen_grid <- function(kernel_locs, L) {
  increment <- L / (kernel_locs)
  extreme <- (L - increment) / 2
  x <- seq(-extreme, extreme, increment)
  g <- expand.grid(x, x)
  g <- as.matrix(g)
  colnames(g) <- c("x", "y")
  g
}

################################################################################
# Functions for predict.dstm()
update_C <- function(C_prev, G) {
  C_new <- array(-1, dim=dim(C_prev))
  n_samples <- dim(C_new)[3]
  for (i in seq(n_samples)) C_new[,,i] <- G[,,i] %*% C_prev[,,i] %*% t(G[,,i])
  C_new
}

calc_W <- function(lambda, C_T) {
  W <- array(-1, dim=dim(C_T))
  n_samples <- dim(C_T)[3]
  for (i in seq(n_samples)) W[,,i] <- lambda[i] * C_T[,,i]
  W
}

next_thetas <- function(thetas_prev, G, W) {
  n_samples <- ncol(thetas_prev)
  E_t <- sapply(seq(n_samples), function(i) G[,,i] %*% thetas_prev[,i])
  sapply(seq(n_samples), function(i) mvtnorm::rmvnorm(1, E_t[,i], W[,,i]))
}
