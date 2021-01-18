#' Test significance of single predictors
#'
#' Compute p-values for the significance of each variable in \code{x}.
#'
#' @param x Input matrix with \code{nobs} rows, each an observation vector.
#' @param y Response variable; shoud be a numeric vector.
#' @param x_alt Optional: a matrix with jth column the sparse projection of the
#'   jth column of x on all its other columns i.e. the output of
#'   \code{\link{sparse_proj}}. If not supplied this is computed by the
#'   function.
#' @param B Number of bootstrap samples. If set to 0, the asymptotic ditribution
#'   is used for calibration.
#' @param rand_gen A function to generate the simulated errors up to an unknown
#'   scale factor. It must permit calling as \code{rand_gen(nobs*B)}. Determines
#'   the form of errors in each of the null models, though the results are
#'   broadly insensitive to this choice. The default \code{rnorm} equates to
#'   null hypotheses of (sparse) Gaussian linear models. Setting
#'   \code{rand_gen=NULL} resamples residuals to generate simulated errors and
#'   approximates nulls of i.i.d. errors with unknown distributions.
#' @param mc.cores Number of cores to use.
#' @return A vector of p-values for each variable.
#' @references Shah, R. D., Buhlmann, P. (2016) \emph{Goodness of fit tests for
#'   high-dimensional linear models} \url{http://arxiv.org/abs/1511.03334}
#' @seealso \code{\link{RPtest}} and \code{\link{sparse_proj}}
#' @examples
#' x <- scale(matrix(rnorm(50*100), 50, 100))
#' x <- scale(x)
#' y <- as.numeric(x[, 1:5] %*% rep(1, 5) + rnorm(nrow(x)))
#' out <- RPtest_single(x=x, y=y, B=25)
#' @export
RPtest_single <- function(x, y, x_alt, B=100L, rand_gen=rnorm, mc.cores=1L) {

  # Params not currently used
  normal_pvals=TRUE
  sigma_est_fac=1

  # Check x, y, x_alt and standardise
  if (!is.matrix(x) | any(is.na(x))) stop("x should be a matrix with at least two columns
                                          and no missing values")
  np <- dim(x)
  if (is.null(np) | (np[2] <= 1L))
    stop("x should be a matrix with at least two columns")
  n <- as.integer(np[1])
  p <- as.integer(np[2])
  y <- as.numeric(y)
  if (length(y) != n) stop("y must have nrow(x) components")
  if (!missing(x_alt)) {
    if (!is.matrix(x_alt) | nrow(x_alt) != n | ncol(x_alt) != p | any(is.na(x_alt)))
      stop("x_alt must have the same dimensions as x and have no missing values")
  } else {
    x_alt <- sparse_proj(x)
  }

  x <- scale(x)
  y <- y - mean(y)

  # Estimate beta(s) and generate true residuals
  beta_out <- comp_beta_star_all(x, y)
  sigma_est <- sigma_est_fac * comp_sigma_est(x, y, beta_out$beta_star_mat[, "main"])

  # Only use main beta to estimate sigma
  resid <- resid_lasso_all(x, y, beta_out=beta_out)

  # Asymptotic distribution calibration
  if (B==0) {
    x_alt <- n*scale(x_alt)/sqrt(n-1)
    test <- test_all_cov(x_alt, resid, test_func=cov_test_twoside)
    pvals_out <- pvals_all_normal_noboot(test)
    return(pvals_out)
  }

  # Nonparametrix bootstrap residuals (only use main beta residuals)
  if (is.null(rand_gen)) {
    rand_gen <- function(n_samp) {
      sample(as.numeric(resid$resid_mat[, 1]), n_samp, replace=TRUE)
    }
  }

  # Bootstrap residuals
  resid_sim <- resid_gen_lasso_all(x, beta_out, sigma_est, B=B, rand_gen=rand_gen, mc.cores=mc.cores)

  # Calculate test statistics
  test_sim <- test_all_cov_sim(x_alt, resid_sim, resid$nonzero, test_func=cov_test_twoside)
  test <- test_all_cov(x_alt, resid, test_func=cov_test_twoside)

  # Compute vector of -values
  pvals_out <- pvals_all_normal(test_sim, test)

  return(pvals_out)
}


# First compute a full lasso fit beta_main
# Then for each variable that is non zero, we compute a separate lasso fit
# excluding that variable
# We return a list with 1st comp a matrix of fits, 2nd comp a logical vector of the nonzero vars
# 3rd comp a matrix describing the folds used
# non_zero cand is a p-logical vector with T for if we need to re-estimate
# when that variable is non-zero
comp_beta_star_all <- function(x, y, non_zero_cand=rep(TRUE, ncol(x)), ...) {
  beta_main <- sqrt_lasso(x, y, ...)

  nonzero <- (beta_main != 0) & non_zero_cand
  beta_star_mat <- matrix(nrow = length(beta_main), ncol = sum(nonzero)+1)
  colnames(beta_star_mat) <- c("main", which(nonzero))
  beta_star_mat[, 1] <- beta_main
  rm(beta_main)

  i <- 2
  w_nonzero <- which(nonzero)
  for(var in w_nonzero) {
    beta_star_mat[, i] <-sqrt_lasso(x, y, exclude=var, ...)
    i <- i+1
  }
  return(list("beta_star_mat"=beta_star_mat, "nonzero"=nonzero))
}

# Returns a list with components
# resid_mat: Scaled residuals after regression ommitting each variable.
# Cols of this matrix correspond to cols of beta_star_mat
resid_lasso_all <- function(x, y, beta_out, non_zero_cand=rep(TRUE, ncol(x)), ...) {
  n <- nrow(x)
  if (missing(beta_out)) {
    beta_out <- comp_beta_star_all(x, y, non_zero_cand=non_zero_cand, ...)
  }
  resid_mat <- scale(y - x %*%  beta_out$beta_star_mat, center = FALSE) * sqrt(n/(n-1))
  return(list("resid_mat"=resid_mat, "nonzero"=beta_out$nonzero))
}

# Equivalent of resid_gen_lasso, but for residuals omitting each variable in
# turn. Output takes form of a list
# First comp main is a list of B lists. Each of the B lists is the output of
# resid_lasso_all
# Remaining components give a matrix of residuals omitting each variable
# selected in the beta_out$beta_star_mat[, 1]
resid_gen_lasso_all <- function(x, beta_out, sigma_est, B=50L,
                                rand_gen = rnorm, mc.cores=1L, ...) {
  n <- nrow(x)
  # Prepare output
  n_beta_star <- ncol(beta_out$beta_star_mat)
  output <- vector(mode = "list", length = n_beta_star)
  names(output) <- c("main", which(beta_out$nonzero))

  not_nonzero <- !beta_out$nonzero

  error_mat <- sigma_est*matrix(rand_gen(n*B), n, B)

  # main first
  signal <- as.numeric(x %*% beta_out$beta_star_mat[, "main"])
  rep_fun <- function(b) {
    y_star <- as.numeric(signal + error_mat[, b])
    return(resid_lasso_all(x, y_star, non_zero_cand=not_nonzero, ...))
  }
  output$main <- parallel::mclapply(1:B, rep_fun, mc.cores=mc.cores)

  # all others
  w_nonzero <- which(beta_out$nonzero)
  i <- 2
  for (var in w_nonzero) {
    output[[i]] <- resid_gen_lasso(x, list("beta"=beta_out$beta_star_mat[, i]),
                                   sigma_est, B=B, rand_gen = rand_gen, mc.cores=mc.cores,
                                   error_mat=error_mat, exclude=var, ...)
    i <- i + 1
  }
  return(output)
}

# B by p matrix of absolute values of correlations
# x could be a residualised version of the original design matrix
test_all_cov_sim <- function(x, resid_list_all, nonzero, test_func=cov_test_twoside()) {
  p <- ncol(x); n <- nrow(x)
  B <- length(resid_list_all$main)
  output <- matrix(nrow=B, ncol=p)
  for (var in seq_len(p)) {
    if (nonzero[var]) {
      output[, var] <- test_func(x[, var], resid_list_all[[as.character(var)]])
    } else {
      for (b in seq_len(B)) {
        if (resid_list_all$main[[b]]$nonzero[var]) {
          output[b, var] <-
            test_func(x[, var], resid_list_all$main[[b]]$resid_mat[, as.character(var)])
        } else {
          output[b, var] <- test_func(x[, var], resid_list_all$main[[b]]$resid_mat[, "main"])
        }
      }
    }
  }
  return(output)
}

# p-vector of absolute values of correlations
test_all_cov <- function(x, resid_all, test_func=cov_test_twoside()) {
  p <- ncol(x); n <- nrow(x)
  output <- numeric(p)
  for (var in seq_len(p)) {
    if (resid_all$nonzero[var]) {
      output[var] <- test_func(x[, var], resid_all$resid_mat[, as.character(var)])
    } else {
      output[var] <- test_func(x[, var], resid_all$resid_mat[, "main"])
    }
  }
  return(output)
}

# p-vector of p vals
pvals_all_normal <- function(test_stat_sim, test_stat) {
  p <- length(test_stat)
  B <- nrow(test_stat_sim)
  output <- numeric(p)
  for (var in seq_len(p)) {
    output[var] <- 2 * pnorm(abs(test_stat[var]-mean(test_stat_sim[, var])),
                             sd=sd(test_stat_sim[, var]), lower.tail = FALSE)
  }
  return(output)
}

pvals_all_normal_noboot <- function(test_stat) {
  p <- length(test_stat)
  output <- numeric(p)
  for (var in seq_len(p)) {
    output[var] <- 2 * pnorm(abs(test_stat[var]),
                             sd=1, lower.tail = FALSE)
  }
  return(output)
}
