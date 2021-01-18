# ========================= find_lambda_one_d_rcpp =========================

#' Selecting Box-Cox parameter lambda in the one-dimensional case using
#' C++ via Rcpp
#'
#' Finds a value of the Box-Cox transformation parameter lambda for which
#' the (positive univariate) random variable with log-density logf
#' has a density closer to that of a Gaussian random variable.
#' Works by estimating a set of quantiles of the distribution implied
#' by logf and treating those quantiles as data in a standard Box-Cox
#' analysis.  In the following we use theta to denote the argument of
#' logf on the original scale and phi on the Box-Cox transformed scale.
#'
#' @param logf A pointer to a compiled C++ function returning the log
#'   of the target density \eqn{f}.
#' @param ... further arguments to be passed to \code{logf} and related
#'   functions.
#' @param ep_bc A (positive) numeric scalar. Smallest possible value of phi
#'   to consider.  Used to avoid negative values of phi.
#' @param min_phi,max_phi  Numeric scalars.  Smallest and largest values
#'   of phi at which to evaluate logf, i.e. the range of values of phi over
#'   which to evaluate logf. Any components in min_phi that are not positive
#'   are set to ep_bc.
#' @param num A numeric scalar. Number of values at which to evaluate logf.
#' @param xdiv A numeric scalar.  Only values of phi at which the density f is
#'   greater than the (maximum of f) / \code{xdiv} are used.
#' @param probs A numeric scalar. Probabilities at which to estimate the
#'   quantiles of that will be used as data to find lambda.
#' @param lambda_range A numeric vector of length 2.  Range of lambda over
#'   which to optimise.
#' @param phi_to_theta A pointer to a compiled C++ function returning
#'   (the inverse) of the transformation from theta to phi used to ensure
#'   positivity of phi prior to Box-Cox transformation.  The argument is
#'   phi and the returned value is theta.  If \code{phi_to_theta}
#'   is undefined at the input value then the function should return NA.
#' @param user_args A list of numeric components providing arguments to
#'   the user-supplied functions \code{phi_to_theta} and \code{log_j}.
#' @param log_j A pointer to a compiled C++ function returning the log of
#'  the Jacobian of the transformation from theta to phi, i.e. based on
#'  derivatives of phi with respect to theta. Takes theta as its argument.
#'  If this is not supplied then a constant Jacobian is used.
#' @details The general idea is to estimate quantiles of f corresponding to a
#'   set of equally-spaced probabilities in \code{probs} and to use these
#'   estimated quantiles as data in a standard estimation of the Box-Cox
#'   transformation parameter \code{lambda}.
#'
#'   The density f is first evaluated at \code{num} points equally spaced over
#'   the interval (\code{min_phi}, \code{max_phi}).  The continuous density f
#'   is approximated by attaching trapezium-rule estimates of probabilities
#'   to the midpoints of the intervals between the points.  After standardizing
#'   to account for the fact that f may not be normalized,
#'   (\code{min_phi}, \code{max_phi}) is reset so that values with small
#'   estimated probability (determined by \code{xdiv}) are excluded and the
#'   procedure is repeated on this new range.  Then the required quantiles are
#'   estimated by inferring them from a weighted empirical distribution
#'   function based on treating the midpoints as data and the estimated
#'   probabilities at the midpoints as weights.
#' @return A list containing the following components
#'   \item{lambda}{A numeric scalar.  The value of \code{lambda}.}
#'   \item{gm}{A numeric scalar.  Box-cox scaling parameter, estimated by the
#'     geometric mean of the quantiles used in the optimisation to find the
#'     value of lambda.}
#'   \item{init_psi}{A numeric scalar.  An initial estimate of the mode of the
#'     Box-Cox transformed density}
#'   \item{sd_psi}{A numeric scalar.  Estimates of the marginal standard
#'     deviations of the Box-Cox transformed variables.}
#'  \item{phi_to_theta}{as detailed above (only if \code{phi_to_theta} is
#'    supplied)}
#'  \item{log_j}{as detailed above (only if \code{log_j} is supplied)}
#'  \item{user_args}{as detailed above (only if \code{user_args} is supplied)}
#'
#' @references Box, G. and Cox, D. R. (1964) An Analysis of Transformations.
#'  Journal of the Royal Statistical Society. Series B (Methodological), 26(2),
#'  211-252.
#' @references Andrews, D. F. and Gnanadesikan, R. and Warner, J. L. (1971)
#'  Transformations of Multivariate Data, Biometrics, 27(4).
#' @references Eddelbuettel, D. and Francois, R. (2011). Rcpp: Seamless
#'  R and C++ Integration. \emph{Journal of Statistical Software},
#'  \strong{40}(8), 1-18.
#'  \url{https://www.jstatsoft.org/v40/i08/}.
#' @references Eddelbuettel, D. (2013). \emph{Seamless R and C++ Integration
#'  with Rcpp}, Springer, New York. ISBN 978-1-4614-6867-7.
#' @examples
#'
#' # Log-normal density ===================
#'
#' # Note: the default value of max_phi = 10 is OK here but this will not
#' # always be the case.
#'
#' ptr_lnorm <- create_xptr("logdlnorm")
#' mu <- 0
#' sigma <- 1
#' lambda <- find_lambda_one_d_rcpp(logf = ptr_lnorm, mu = mu, sigma = sigma)
#' lambda
#' x <- ru_rcpp(logf = ptr_lnorm, mu = mu, sigma = sigma, log = TRUE, d = 1,
#'              n = 1000, trans = "BC", lambda = lambda)
#'
#' # Gamma density ===================
#'
#' alpha <- 1
#' # Choose a sensible value of max_phi
#' max_phi <- qgamma(0.999, shape = alpha)
#' # [I appreciate that typically the quantile function won't be available.
#' # In practice the value of lambda chosen is quite insensitive to the choice
#' # of max_phi, provided that max_phi is not far too large or far too small.]
#'
#' ptr_gam <- create_xptr("logdgamma")
#' lambda <- find_lambda_one_d_rcpp(logf = ptr_gam, alpha = alpha,
#'                                  max_phi = max_phi)
#' lambda
#' x <- ru_rcpp(logf = ptr_gam, alpha = alpha, d = 1, n = 1000, trans = "BC",
#'              lambda = lambda)
#'
#' alpha <- 0.1
#' # NB. for alpha < 1 the gamma(alpha, beta) density is not bounded
#' # So the ratio-of-uniforms emthod can't be used but it may work after a
#' # Box-Cox transformation.
#' # find_lambda_one_d() works much better than find_lambda() here.
#'
#' max_phi <- qgamma(0.999, shape = alpha)
#' lambda <- find_lambda_one_d_rcpp(logf = ptr_gam, alpha = alpha,
#'                                  max_phi = max_phi)
#' lambda
#' x <- ru_rcpp(logf = ptr_gam, alpha = alpha, d = 1, n = 1000, trans = "BC",
#'              lambda = lambda)
#' \donttest{
#' plot(x)
#' plot(x, ru_scale = TRUE)
#' }
#' @seealso \code{\link{ru_rcpp}} to perform ratio-of-uniforms sampling.
#' @seealso \code{\link{find_lambda_rcpp}} to produce (somewhat) automatically
#'   a list for the argument \code{lambda} of \code{ru} for any value of
#'   \code{d}.
#'
#' @export
find_lambda_one_d_rcpp <- function(logf, ..., ep_bc = 1e-4, min_phi = ep_bc,
                              max_phi = 10, num = 1001L, xdiv = 100,
                              probs = seq(0.01, 0.99, by = 0.01),
                              lambda_range = c(-3, 3), phi_to_theta = NULL,
                              log_j = NULL, user_args = list()) {
  # Check that logf is an external pointer.
  is_pointer <- (class(logf) == "externalptr")
  if (!is_pointer) {
    stop("logf must be an external pointer to a function")
  }
  # Check if phi_to_theta and log_j are supplied then they are external
  # pointers.
  is_pointer <- (class(phi_to_theta) == "externalptr")
  if (!is_pointer & !is.null(phi_to_theta)) {
    stop("phi_to_theta must be an external pointer to a function or NULL")
  }
  is_pointer <- (class(log_j) == "externalptr")
  if (!is_pointer & !is.null(log_j)) {
    stop("log_j must be an external pointer to a function or NULL")
  }
  # Check that max_phi > min_phi in all cases
  if (any(max_phi-min_phi <= 0)) {
    stop("max_phi must be larger than min_phi elementwise.")
  }
  if (any(max_phi <= 0)) {
    stop("all components of max_phi must be positive")
  }
  # Set to ep_bc any non-positive elements in min_phi.
  min_phi <- ifelse(min_phi > 0, min_phi, ep_bc)
  # If phi_to_theta is supplied save them to return them later
  trans_list <- list()
  if (!is.null(phi_to_theta)) {
    trans_list$phi_to_theta <- phi_to_theta
    trans_list$log_j <- log_j
    if (is.null(log_j)) {
      log_j = create_log_jac_xptr("log_none_jac")
      trans_list$log_j <- log_j
    }
  } else {
    phi_to_theta = null_phi_to_theta_xptr("no_trans")
    log_j = create_log_jac_xptr("log_none_jac")
  }
  trans_list$user_args <- user_args
  # Set num equally-spaced values of x in [min_phi, max_phi]
  x <- seq(min_phi, max_phi, len = num)
  # Set parameters for passing to cpp_log_rho_4().
  pars <- list(...)
  # Function to determine how deep a list is, i.e. how many layers
  # of listing it has.
  list_depth <- function(x) {
    ifelse(is.list(x), 1L + max(sapply(x, list_depth)), 0L)
  }
  # Find the depth of pars.
  if (length(pars) > 0) {
    pars_depth <- list_depth(pars)
  } else {
    pars_depth <- 0
  }
  # If the user has supplied a list rather than individual components
  # then remove the extra layer of list and retrieve the original
  # variable names.
  if (pars_depth > 1) {
    par_names <- names(pars)
    pars <- unlist(pars, recursive = FALSE)
    # Remove the user's list name, if they gave it one.
    if (!is.null(par_names)) {
      keep_name <- nchar(par_names) + 2
      names(pars) <- substring(names(pars), keep_name)
    }
  }
  tpars <- list()
  ptpfun <- create_psi_to_phi_xptr("no_trans")
  d <- 1
  logf_args <- list(psi_mode = rep(0, d), rot_mat = diag(d), hscale = 0,
                    logf = logf, pars = pars, tpars = tpars, ptpfun = ptpfun,
                    phi_to_theta = phi_to_theta, log_j = log_j,
                    user_args = user_args)
  # Calculate the density (weights) at these values
  log_w <- do.call(rcpp_apply, c(list(x = matrix(x)), logf_args))
  # Shift log_w so that it has a maximum of 0, to try to avoid underflow.
  log_w <- log_w - max(log_w, na.rm = TRUE)
  # Evaluate the density values.
  w <- exp(log_w)
  # Standardise, so that the weights sum to 1
  w <- w / sum(w)
  n <- length(w)
  # Calculate the mean value on each interval between xs
  wbar <- (w[1:(n-1)] + w[2:n]) / 2
  # Calculate the widths of the intervals
  xdiff <- x[2]-x[1]
  # Find the midpoints of the intervals
  xmid <- (x[1:(n-1)] + x[2:n]) / 2
  # Approximate the probability on each interval (by area of trapezium)
  areas <- wbar * xdiff
  # Which areas are greater than the maximum area / xdiv?
  r <- range(xmid[areas > max(areas) / xdiv])
  # Reset the xs over this new range
  x <- seq(r[1], r[2], len = num)
  # Recalculate the weights and the areas and midpoints
  log_w <- do.call(rcpp_apply, c(list(x = matrix(x)), logf_args))
  log_w <- log_w - max(log_w, na.rm = TRUE)
  w <- exp(log_w)
  w <- w / sum(w)
  n <- length(w)
  wbar <- (w[1:(n-1)] + w[2:n]) / 2
  xdiff <- x[2]-x[1]
  xmid <- (x[1:(n-1)] + x[2:n]) / 2
  areas <- wbar * xdiff
  w <- areas / sum(areas)
  # Estimate the 100*probs% quantiles of the density
  qs <- stats::quantile(wecdf(xmid, w), probs = probs)
  qs <- matrix(qs, ncol = 1)
  # Use the quantiles qs as a sample of data to estimate lambda
  # Set to 1 all the input weights associated with the quantiles
  w_ones <- rep(1, nrow(qs))
  temp <- optim_box_cox(x = qs, w = w_ones, lambda_range = lambda_range)
  # Calculate an estimate of the mode of the transformed density
  fy <- w * (xmid / temp$gm) ^ (1 - temp$lambda)
  mode_y <- which.max(fy)
  temp$init_psi <- box_cox(xmid[mode_y], lambda = temp$lambda, gm = temp$gm)
  # Calculate an estimate of the standard deviation of the transformed density
  low_q <- qs[which.min(abs(probs-stats::pnorm(-1 / 2)))]
  up_q <- qs[which.min(abs(probs-stats::pnorm(1 / 2)))]
  temp$sd_psi <- up_q - low_q
  return(c(temp, trans_list))
}

# ========================= find_lambda_rcpp =========================

#' Selecting Box-Cox parameter lambda for general d using C++ via Rcpp.
#'
#' Finds a value of the Box-Cox transformation parameter lambda for which
#' the (positive) random variable with log-density logf has a density
#' closer to that of a Gaussian random variable.
#' In the following we use theta to denote the argument of
#' logf on the original scale and phi on the Box-Cox transformed scale.
#'
#' @param logf A pointer to a compiled C++ function returning the log
#'   of the target density \eqn{f}.
#' @param ... further arguments to be passed to \code{logf} and related
#'   functions.
#' @param d A numeric scalar. Dimension of f.
#' @param n_grid A numeric scalar.  Number of ordinates for each variable in
#'   phi.  If this is not supplied a default value of
#'   ceiling(2501 ^ (1 / d)) is used.
#' @param ep_bc A (positive) numeric scalar. Smallest possible value of phi
#'   to consider.  Used to avoid negative values of phi.
#' @param min_phi,max_phi  Numeric vectors.  Smallest and largest values
#'   of phi at which to evaluate logf, i.e. the range of values of phi over
#'   which to evaluate logf. Any components in min_phi that are not positive
#'   are set to ep_bc.
#' @param which_lam A numeric vector.  Contains the indices of the components
#'   of phi that ARE to be Box-Cox transformed.
#' @param lambda_range A numeric vector of length 2.  Range of lambda over
#'   which to optimise.
#' @param init_lambda A numeric vector of length 1 or d.  Initial value of
#'   lambda used in the search for the best lambda.  If \code{init_lambda}
#'   is a scalar then \code{rep(init_lambda, d)} is used.
#' @param phi_to_theta A pointer to a compiled C++ function returning
#'   (the inverse) of the transformation from theta to phi used to ensure
#'   positivity of phi prior to Box-Cox transformation.  The argument is
#'   phi and the returned value is theta.  If \code{phi_to_theta}
#'   is undefined at the input value then the function should return NA.
#' @param user_args A list of numeric components providing arguments to
#'   the user-supplied functions \code{phi_to_theta} and \code{log_j}.
#' @param log_j A pointer to a compiled C++ function returning the log of
#'  the Jacobian of the transformation from theta to phi, i.e. based on
#'  derivatives of phi with respect to theta. Takes theta as its argument.
#' @details The general idea is to evaluate the density f on a d-dimensional
#'  grid, with \code{n_grid} ordinates for each of the \code{d} variables.
#'  We treat each combination of the variables in the grid as a data point
#'  and perform an estimation of the Box-Cox transformation parameter
#'  \code{lambda}, in which each data point is weighted by the density
#'  at that point.  The vectors \code{min_phi} and \code{max_phi} define the
#'  limits of the grid and \code{which_lam} can be used to specify that only
#'  certain components of phi are to be transformed.
#' @return A list containing the following components
#'   \item{lambda}{A numeric vector.  The value of \code{lambda}.}
#'   \item{gm}{A numeric vector.  Box-cox scaling parameter, estimated by the
#'     geometric mean of the values of phi used in the optimisation to find
#'     the value of lambda, weighted by the values of f evaluated at phi.}
#'   \item{init_psi}{A numeric vector.  An initial estimate of the mode of the
#'     Box-Cox transformed density}
#'   \item{sd_psi}{A numeric vector.  Estimates of the marginal standard
#'     deviations of the Box-Cox transformed variables.}
#'  \item{phi_to_theta}{as detailed above (only if \code{phi_to_theta} is
#'    supplied)}
#'  \item{log_j}{as detailed above (only if \code{log_j} is supplied)}
#'  \item{user_args}{as detailed above (only if \code{user_args} is supplied)}
#' @references Box, G. and Cox, D. R. (1964) An Analysis of Transformations.
#'  Journal of the Royal Statistical Society. Series B (Methodological), 26(2),
#'  211-252.
#' @references Andrews, D. F. and Gnanadesikan, R. and Warner, J. L. (1971)
#'  Transformations of Multivariate Data, Biometrics, 27(4).
#' @references Eddelbuettel, D. and Francois, R. (2011). Rcpp: Seamless
#'  R and C++ Integration. \emph{Journal of Statistical Software},
#'  \strong{40}(8), 1-18.
#'  \url{https://www.jstatsoft.org/v40/i08/}.
#' @references Eddelbuettel, D. (2013). \emph{Seamless R and C++ Integration
#'  with Rcpp}, Springer, New York. ISBN 978-1-4614-6867-7.
#' @examples
#'
#' # Log-normal density ===================
#' # Note: the default value max_phi = 10 is OK here but this will not always
#' # be the case
#' ptr_lnorm <- create_xptr("logdlnorm")
#' mu <- 0
#' sigma <- 1
#' lambda <- find_lambda_rcpp(logf = ptr_lnorm, mu = mu, sigma = sigma)
#' lambda
#' x <- ru_rcpp(logf = ptr_lnorm, mu = mu, sigma = sigma, d = 1, n = 1000,
#'              trans = "BC", lambda = lambda)
#'
#' # Gamma density ===================
#' alpha <- 1
#' #  Choose a sensible value of max_phi
#' max_phi <- qgamma(0.999, shape = alpha)
#' # [Of course, typically the quantile function won't be available.  However,
#' # In practice the value of lambda chosen is quite insensitive to the choice
#' # of max_phi, provided that max_phi is not far too large or far too small.]
#'
#' ptr_gam <- create_xptr("logdgamma")
#' lambda <- find_lambda_rcpp(logf = ptr_gam, alpha = alpha, max_phi = max_phi)
#' lambda
#' x <- ru_rcpp(logf = ptr_gam, alpha = alpha, d = 1, n = 1000, trans = "BC",
#'              lambda = lambda)
#'
#' \donttest{
#' # Generalized Pareto posterior distribution ===================
#'
#' n <- 1000
#' # Sample data from a GP(sigma, xi) distribution
#' gpd_data <- rgpd(m = 100, xi = -0.5, sigma = 1)
#' # Calculate summary statistics for use in the log-likelihood
#' ss <- gpd_sum_stats(gpd_data)
#' # Calculate an initial estimate
#' init <- c(mean(gpd_data), 0)
#'
#' n <- 1000
#' # Sample on original scale, with no rotation ----------------
#' ptr_gp <- create_xptr("loggp")
#' for_ru_rcpp <- c(list(logf = ptr_gp, init = init, d = 2, n = n,
#'                      lower = c(0, -Inf)), ss, rotate = FALSE)
#' x1 <- do.call(ru_rcpp, for_ru_rcpp)
#' plot(x1, xlab = "sigma", ylab = "xi")
#' # Parameter constraint line xi > -sigma/max(data)
#' # [This may not appear if the sample is far from the constraint.]
#' abline(a = 0, b = -1 / ss$xm)
#' summary(x1)
#'
#' # Sample on original scale, with rotation ----------------
#' for_ru_rcpp <- c(list(logf = ptr_gp, init = init, d = 2, n = n,
#'                       lower = c(0, -Inf)), ss)
#' x2 <- do.call(ru_rcpp, for_ru_rcpp)
#' plot(x2, xlab = "sigma", ylab = "xi")
#' abline(a = 0, b = -1 / ss$xm)
#' summary(x2)
#'
#' # Sample on Box-Cox transformed scale ----------------
#'
#' # Find initial estimates for phi = (phi1, phi2),
#' # where phi1 = sigma
#' #   and phi2 = xi + sigma / max(x),
#' # and ranges of phi1 and phi2 over over which to evaluate
#' # the posterior to find a suitable value of lambda.
#' temp <- do.call(gpd_init, ss)
#' min_phi <- pmax(0, temp$init_phi - 2 * temp$se_phi)
#' max_phi <- pmax(0, temp$init_phi + 2 * temp$se_phi)
#'
#' # Set phi_to_theta() that ensures positivity of phi
#' # We use phi1 = sigma and phi2 = xi + sigma / max(data)
#'
#' # Create an external pointer to this C++ function
#' ptr_phi_to_theta_gp <- create_phi_to_theta_xptr("gp")
#' # Note: log_j is set to zero by default inside find_lambda_rcpp()
#' lambda <- find_lambda_rcpp(logf = ptr_gp, ss = ss, d = 2, min_phi = min_phi,
#'                            max_phi = max_phi, user_args = list(xm = ss$xm),
#'                            phi_to_theta = ptr_phi_to_theta_gp)
#' lambda
#'
#' # Sample on Box-Cox transformed, without rotation
#' x3 <- ru_rcpp(logf = ptr_gp, ss = ss, d = 2, n = n, trans = "BC",
#'               lambda = lambda, rotate = FALSE)
#' plot(x3, xlab = "sigma", ylab = "xi")
#' abline(a = 0, b = -1 / ss$xm)
#' summary(x3)
#'
#' # Sample on Box-Cox transformed, with rotation
#' x4 <- ru_rcpp(logf = ptr_gp, ss = ss, d = 2, n = n, trans = "BC",
#'               lambda = lambda)
#' plot(x4, xlab = "sigma", ylab = "xi")
#' abline(a = 0, b = -1 / ss$xm)
#' summary(x4)
#'
#' def_par <- graphics::par(no.readonly = TRUE)
#' par(mfrow = c(2,2), mar = c(4, 4, 1.5, 1))
#' plot(x1, xlab = "sigma", ylab = "xi", ru_scale = TRUE,
#'   main = "mode relocation")
#' plot(x2, xlab = "sigma", ylab = "xi", ru_scale = TRUE,
#'   main = "mode relocation and rotation")
#' plot(x3, xlab = "sigma", ylab = "xi", ru_scale = TRUE,
#'   main = "Box-Cox and mode relocation")
#' plot(x4, xlab = "sigma", ylab = "xi", ru_scale = TRUE,
#'   main = "Box-Cox, mode relocation and rotation")
#' graphics::par(def_par)
#' }
#' @seealso \code{\link{ru_rcpp}} to perform ratio-of-uniforms sampling.
#' @seealso \code{\link{find_lambda_one_d_rcpp}} to produce (somewhat)
#'   automatically a list for the argument \code{lambda} of \code{ru} for the
#'   \code{d} = 1 case.
#' @export
find_lambda_rcpp <- function(logf, ..., d = 1, n_grid = NULL, ep_bc = 1e-4,
                        min_phi = rep(ep_bc, d), max_phi = rep(10, d),
                        which_lam = 1:d, lambda_range = c(-3,3),
                        init_lambda = NULL, phi_to_theta = NULL,
                        log_j = NULL, user_args = list()) {
  # Check that logf is an external pointer.
  is_pointer <- (class(logf) == "externalptr")
  if (!is_pointer) {
    stop("logf must be an external pointer to a function")
  }
  # Check if phi_to_theta and log_j are supplied then they are external
  # pointers.
  is_pointer <- (class(phi_to_theta) == "externalptr")
  if (!is_pointer & !is.null(phi_to_theta)) {
    stop("phi_to_theta must be an external pointer to a function or NULL")
  }
  is_pointer <- (class(log_j) == "externalptr")
  if (!is_pointer & !is.null(log_j)) {
    stop("log_j must be an external pointer to a function or NULL")
  }
  #
  if (!is.null(init_lambda)) {
    if (!is.vector(init_lambda)) {
      stop("init_lambda must be a vector.")
    }
    if (length(init_lambda) == 1) {
      init_lambda <- rep(init_lambda, d)
    }
    if (length(init_lambda) != d) {
      stop("init_lambda must be a vector of length 1 or d.")
    }
  }
  # Check that max_phi > min_phi in all cases
  if (any(max_phi-min_phi <= 0)) {
    stop("max_phi must be larger than min_phi elementwise.")
  }
  # Set to ep_bc any non-positive elements in min_phi and max_phi.
  min_phi <- ifelse(min_phi > 0, min_phi, ep_bc)
  max_phi <- ifelse(max_phi > 0, max_phi, ep_bc)
  # If phi_to_theta is supplied save them to return them later
  trans_list <- list()
  if (!is.null(phi_to_theta)) {
    trans_list$phi_to_theta <- phi_to_theta
    trans_list$log_j <- log_j
    if (is.null(log_j)) {
      log_j = create_log_jac_xptr("log_none_jac")
      trans_list$log_j <- log_j
    }
  } else {
    phi_to_theta = null_phi_to_theta_xptr("no_trans")
    log_j = create_log_jac_xptr("log_none_jac")
  }
  trans_list$user_args <- user_args
  # Evaluate the target density (up to a multiplicative constant) over a
  # grid of values that contains most of the probability.
  #
  # If n_grid has not been specified then set a default value
  if (is.null(n_grid)) {
    n_grid <- n_grid_fn(d)
  }
  # Set the minimum and maximum for each variable
  low_phi <- min_phi
  up_phi <- max_phi
  # Create matrix with each column giving the values of each variable
  phi <- mapply(seq, from = low_phi, to = up_phi, len = n_grid)
  # Make this matrix into a list with an entry for each column
  phi <- lapply(seq_len(ncol(phi)), function(i) phi[, i])
  # Expand into a matrix containing the grid of combinations (one
  # combination in each row).
  phi <- expand.grid(phi)
  # Set parameters for passing to cpp_log_rho_4().
  pars <- list(...)
  # Function to determine how deep a list is, i.e. how many layers
  # of listing it has.
  list_depth <- function(x) {
    ifelse(is.list(x), 1L + max(sapply(x, list_depth)), 0L)
  }
  # Find the depth of pars.
  if (length(pars) > 0) {
    pars_depth <- list_depth(pars)
  } else {
    pars_depth <- 0
  }
  # If the user has supplied a list rather than individual components
  # then remove the extra layer of list and retrieve the original
  # variable names.
  if (pars_depth > 1) {
    par_names <- names(pars)
    pars <- unlist(pars, recursive = FALSE)
    # Remove the user's list name, if they gave it one.
    if (!is.null(par_names)) {
      keep_name <- nchar(par_names) + 2
      names(pars) <- substring(names(pars), keep_name)
    }
  }
  tpars <- list()
  ptpfun <- create_psi_to_phi_xptr("no_trans")
  logf_args <- list(psi_mode = rep(0, d), rot_mat = diag(d), hscale = 0,
                    logf = logf, pars = pars, tpars = tpars, ptpfun = ptpfun,
                    phi_to_theta = phi_to_theta, log_j = log_j,
                    user_args = user_args)
  # Evaluate the target log-density at each combination in the grid.
  phi <- data.matrix(phi)
  log_w <- do.call(rcpp_apply, c(list(x = phi), logf_args))
  # Shift log_w so that it has a maximum of 0, to try to avoid underflow.
  log_w <- log_w - max(log_w, na.rm = TRUE)
  # Evaluate the density values.
  w <- exp(log_w)
  #
  # Seek marginal Box-Cox transformations such that the joint posterior is
  # closer to being bivariate normal
  #
  if (any(phi[, which_lam] <= 0)) {
    stop("Attempt to use Box-Cox transformation on non-positive value(s).")
  }
  res_bc <- optim_box_cox(x = phi, w = w, which_lam = which_lam,
                          lambda_range = lambda_range, start = init_lambda)
  lambda <- res_bc$lambda
  gm <- res_bc$gm
  phi_to_psi <- function(phi)  {
    for (j in 1:ncol(phi)) {
      phi_in <- phi[, j]
      phi[, j] <- box_cox(x = phi_in, lambda = lambda[j], gm = gm[j])
    }
    return(phi)
  }
  temp <- init_psi_calc(phi_to_psi = phi_to_psi, phi = phi, lambda = lambda,
                        gm = gm, w = w, which_lam = which_lam)
  temp <- c(list(lambda = lambda, gm = gm), temp)
  return(c(temp, trans_list))
}
