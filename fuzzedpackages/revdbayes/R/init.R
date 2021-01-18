diag_pos <- function(x) {
  # Diagonal elements of a matrix, returning NA for any negative values
  y <- diag(x)
  y[y < 0] <- NA
  return(y)
}

# ======================= change_pp_pars ========================

change_pp_pars <- function(pars, in_noy, out_noy) {
  #
  # Converts PP parameters pars = (mu, sigma, xi) from a parameterisation
  # in which the number of years (blocks) is in_noy to one in which the
  # number of years is out_noy.
  #
  # Arguments:
  #   pars    : A numeric vector. PP parameters (mu, sigma, xi).
  #   in_noy  : A numeric scalar. Input noy.
  #   out_noy : A numeric scalar. Output noy.
  #
  # Returns:
  #   A numeric vector.  Output PP parameter vector.
  #
  const <- in_noy / out_noy
  pars[1] <- pars[1] + pars[2] * box_cox(const, pars[3])
  pars[2] <- pars[2] * const ^ pars[3]
  return(pars)
}

# =========================== pp_init ===========================

pp_init <- function(data, n_exc = length(data[!is.na(data)]), thresh,
                    xm = max(data, na.rm = TRUE), noy,
                    sum_pp = sum(data, na.rm = TRUE),
                    m, xi_eq_zero = FALSE, init_ests = NULL) {
  #
  # Initial estimates for point process model parameters
  #
  # Calculates initial estimates and estimated standard errors (SEs) for
  # the point process model parameters mu, sigma and xi based on an
  # assumed realisation from this process.  Also calculates initial
  # estimates and estimated standard errors for
  # phi1 = mu,
  # phi2 = (sigma + xi * (thresh - mu) ) / sr
  # and
  # phi3 = (sigma + xi * (xm - mu) ) / sr
  # where thresh is the threshold, xm is the sample maximum and
  # sr = (xm - thresh) ^ (1/2).
  #
  # Arguments:
  #   data       : A numeric vector containing threshold exceedances.
  #   n_exc      : A numeric scalar.  The number of threshold exceedances,
  #   thresh     : A numeric scalar. The threshold. i.e. the length of data.
  #   xm         : A numeric scalar. The sample maximum.
  #   noy        : A numeric scalar.  The number of years (blocks) of data.
  #   sum_pp     : A numeric scalar. The sum of the sample threshold excesses.
  #   m          : A numeric scalar.  The number of non-missing observations
  #                in the raw data.
  #   xi_eq_zero : A logical scalar.  If TRUE assume that the shape
  #                parameter xi = 0.
  #   init_ests  : A numeric vector.  Initial estimate of
  #                theta = (mu, sigma, xi).  If supplied \code{pp_init()} just
  #                returns the corresponding initial estimate of
  #                phi = (phi1, phi2, phi3).
  # Returns:
  #   A list with components
  #     init     : A numeric vector. Initial estimates of (mu, sigma, xi).
  #     se       : A numeric vector. Estimated standard errors of
  #                (mu, sigma, xi).
  #     init_phi : A numeric vector. Initial estimates of (phi1, phi2, phi3).
  #     se_phi   : A numeric vector. Estimated standard errors of
  #                (phi1, phi2, phi3).
  #   If init_ests is supplied then only the numeric vector init_phi is
  #     returned.
  #
  sr <- sqrt(xm - thresh)
  theta_to_phi <- function(theta) {
    phi1 <- theta[1]
    phi2 <- (theta[2] + theta[3] * (thresh - theta[1])) / sr
    phi3 <- (theta[2] + theta[3] * (xm - theta[1])) / sr
    c(phi1, phi2, phi3)
  }
  #
  if (!is.null(init_ests)) {
    return(theta_to_phi(init_ests))
  }
  #
  # Estimates of PP parameters -----
  #
  # Find intial estimates of the GP parameters sigma_u and xi based on
  # excesses of threshold thresh (temp_gp$init) and their estimated
  # covariance matrix.
  #
  temp_gp <- gp_init(data = data - thresh, m = n_exc, xm = xm - thresh,
                     sum_gp = sum_pp - thresh, xi_eq_zero = xi_eq_zero)
  #
  # Estimates of GP parameters.
  gp_init <- temp_gp$init
  #
  # Transform from (thresh, sigma_u xi), which are estimates of the PP
  # parameters in which the number of blocks is equal to n_exc, to PP
  # parameters in which the number of blocks is equal to noy.
  # If noy = n_exc then const = 1 and change_pp_pars() has no effect.
  #
  init <- c(thresh, gp_init)
  init <- change_pp_pars(init, in_noy = n_exc, out_noy = noy)
  #
  # Standard errors of PP parameters -----
  #
  # Quantile corresponding to the threshold thresh
  q <- 1 - n_exc / m
  # Estimate of the variance of the threshold, when viewed as an estimator
  # of the quantile q.
  f_thresh <- (1-q) / gp_init[1]
  var_thresh <- q * (1 - q) / (m * f_thresh ^ 2)
  cov_gp <- temp_gp$cov_mtx
  cov_nexc <- matrix(0, 3, 3)
  cov_nexc[1, 1] <- var_thresh
  cov_nexc[2:3, 2:3] <- cov_gp
  const <- n_exc / noy
  m11 <- 1
  m12 <- box_cox(x = const, lambda = gp_init[2])
  m13 <- gp_init[1] * box_cox_deriv(x = const, lambda = gp_init[2])
  row1 <- c(m11, m12, m13)
  m21 <- 0
  m22 <- const ^ gp_init[2]
  m23 <- gp_init[1] * const ^ gp_init[2] * log(const)
  row2 <- c(m21, m22, m23)
  row3 <- c(0, 0, 1)
  mat <- rbind(row1, row2, row3)
  cov_pp <- mat %*% cov_nexc %*% t(mat)
  se <- sqrt(diag(cov_pp))
  #
  init_phi <- theta_to_phi(init)
  row1 <- c(1, 0, 0)
  row2 <- c(-init[3], 1, thresh - init[1]) / sr
  row3 <- c(-init[3], 1, xm - init[1]) / sr
  mat <- rbind(row1, row2, row3)
  var_phi <- mat %*% cov_pp %*% t(mat)
  se_phi <- sqrt(diag(var_phi))
  return(list(init = init, se = se, init_phi = init_phi, se_phi = se_phi))
}

# =========================== gp_init ===========================

gp_init <- function(data, m = length(data[!is.na(data)]),
                    xm = max(data, na.rm = TRUE),
                    sum_gp = sum(data, na.rm = TRUE), xi_eq_zero = FALSE,
                    init_ests = NULL, w = NULL, sumw = NULL) {
  #
  # Initial estimates for Generalized Pareto distribution parameters
  #
  # Calculates initial estimates and estimated standard errors (SEs) for the
  # generalized Pareto parameters sigma and xi based on an
  # assumed random sample from this distribution.  Also calculates
  # initial estimates and estimated standard errors for
  # phi1 = sigma and phi2 = xi + sigma / xm, where xm is the sample maximum.
  #
  # @param data A numeric vector containing positive sample values.
  # @param m A numeric scalar.  The sample size, i.e. the length of data.
  # @param xm A numeric scalar. The sample maximum.
  # @param sum_gp A numeric scalar. The sum of the sample values.
  # @param xi_eq_zero A logical scalar.  If TRUE assume that the shape
  #   parameter xi = 0.
  # @param init_ests A numeric vector.  Initial estimate of
  #   theta = (sigma, xi).  If supplied \code{gp_init()} just
  #   returns the corresponding initial estimate of phi = (phi1, phi2).
  # @details The main aim is to calculate an admissible estimate of
  #   theta = (sigma, xi), i.e. one at which the log-likelihood is finite
  #   (necessary for the posterior log-density to be finite) at the estimate,
  #   and associated estimated SEs. These are converted into estimates and SEs
  #   for phi.  The latter can be used to set values of \code{min_phi} and
  #   \code{max_phi} for input to \code{find_lambda}.
  #
  #   In the default setting (\code{xi_eq_zero = FALSE} and
  #   \code{init_ests = NULL}) the methods tried are Maximum Likelihood
  #   Estimation (MLE) (Grimshaw, 1993), Probability-Weighted Moments (PWM)
  #   (Hosking and Wallis, 1987) and Linear Combinations of Ratios of Spacings
  #   (LRS) (Reiss and Thomas, 2007, page 134) in that order.
  #
  #   For xi < -1 the likelihood is unbounded, MLE may fail when xi is not
  #   greater than -0.5 and the observed Fisher information for (sigma, xi) has
  #   finite variance only if xi > -0.25.  We use the ML estimate provided that
  #   the estimate of xi returned from \code{gp_mle} is greater than -1. We only
  #   use the SE if the MLE of xi is greater than -0.25.
  #
  #   If either the MLE or the SE are not OK then we try PWM.  We use the PWM
  #   estimate only if is admissible, and the MLE was not OK.  We use the PWM
  #   SE, but this will be \code{c(NA, NA)} is the PWM estimate of xi is > 1/2.
  #   If the estimate is still not OK then we try LRS.  As a last resort, which
  #   will tend to occur only when xi is strongly negative, we set xi = -1 and
  #   estimate sigma conditional on this.
  # @return If \code{init_ests} is not supplied by the user, a list is returned
  #   with components
  #     \item{init}{A numeric vector. Initial estimates of sigma
  #      and xi.}
  #     \item{se}{A numeric vector. Estimated standard errors of
  #      sigma and xi.}
  #     \item{cov_mat}{A numeric 2 by 2 matrix. Estimated covariance matrix
  #       associated with the estimates of sigma and xi.}
  #     \item{init_phi}{A numeric vector. Initial estimates of
  #      phi1 = sigma and phi2 = xi + sigma / xm,
  #      where xm is the maximum of \code{data}.}
  #     \item{se_phi}{A numeric vector. Estimated standard errors of
  #      phi1 and phi2.}
  #   If \code{init_ests} is supplied then only the numeric vector
  #   \code{init_phi} is returned.
  # @references Grimshaw, S. D. (1993) Computing Maximum Likelihood Estimates
  #   for the Generalized Pareto Distribution.  Technometrics, 35(2), 185-191.
  #   \url{https://doi.org/10.1080/00401706.1993.10485040}.
  # @references Hosking, J. R. M. and Wallis, J. R. (1987) Parameter and
  #   Quantile Estimation for the Generalized Pareto Distribution.
  #   Technometrics, 29(3), 339-349. \url{https://doi.org/10.2307/1269343}.
  # @references Reiss, R.-D., Thomas, M. (2007) Statistical Analysis of
  #   Extreme Values with Applications to Insurance, Finance, Hydrology and
  #   Other Fields.Birkhauser.
  #   \url{https://doi.org/10.1007/978-3-7643-7399-3}.
  # @seealso \code{\link{gev_init}} to calculate initial estimates for
  #   Generalized extreme value distribution parameters.
  #
  theta_to_phi <- function(theta) {
    c(theta[1], theta[2] + theta[1] / xm)
  }
  #
  if (!is.null(init_ests)) {
    return(theta_to_phi(init_ests))
  }
  if (xi_eq_zero) {
    s_hat <- mean(data)
    init <- c(s_hat, 0)
    init_phi <- theta_to_phi(init)
    cov_mtx <- matrix(c(2*s_hat^2, -s_hat, -s_hat, 1), 2, 2) / m
    se <- sqrt(diag(cov_mtx))
    mat <- matrix(c(1, 0, 1 / xm, 1), 2, 2, byrow = TRUE)
    var_phi <- mat %*% cov_mtx %*% t(mat)
    se_phi <- sqrt(diag(var_phi))
    return(list(init = init, se = se, init_phi = init_phi, se_phi = se_phi))
  }
  #
  ests_ok <- ses_ok <- FALSE
  # First try MLE
  temp <- gp_mle(data)
  init <- temp$mle
  # Check that the MLE is OK
  if (!is.na(init[1]) & init[2] > -1 & init[2] > -init[1]/xm &
      !is.infinite(temp$nllh)){
    ests_ok <- TRUE
    # Check whether or not we should use the SE
    if (init[2] > -0.25){
      cov_mtx <- solve(gp_obs_info(init, data))
      se <- sqrt(diag_pos(cov_mtx))
      mat <- matrix(c(1, 0, 1 / xm, 1), 2, 2, byrow = TRUE)
      var_phi <- mat %*% cov_mtx %*% t(mat)
      if (all(diag(var_phi) > 0)){
        se_phi <- sqrt(diag(var_phi))
        ses_ok <- TRUE
      }
    }
  }
  if (!ests_ok | !ses_ok){
    # Try PWM
    pwm <- gp_pwm(data)
    se <- pwm$se
    cov_mtx <- pwm$cov
    mat <- matrix(c(1, 0, 1 / xm, 1), 2, 2, byrow = TRUE)
    var_phi <- mat %*% cov_mtx %*% t(mat)
    se_phi <- sqrt(diag(var_phi))
    # Note: se and se_phi will be NA if pwm$est[2] > 1/2
    check <- gp_loglik(pars = pwm$est, data = data, m = m, xm = xm,
                        sum_gp = sum_gp)
    # If MLE wasn't OK and PWM estimate is OK then use PWM estimate
    if (!ests_ok & init[2] > -1 & !is.infinite(check)) {
      init <- pwm$est
      ests_ok <- TRUE
    }
  }
  # If estimate is not OK then try LRS
  if (!ests_ok){
    init <- gp_lrs(data)
    check <- gp_loglik(pars = init, data = data, m = m, xm = xm,
                        sum_gp = sum_gp)
    if (init[2] > -1 & !is.infinite(check)) {
      ests_ok <- TRUE
    }
  }
  # If we get here without ests_ok = TRUE then the posterior mode is
  # probably close to xi = -1.  We want to avoid xi < -1 so we set
  # xi = -1 and use the (bias-adjusted) mle of sigma conditional on xi = -1.
  if (!ests_ok) {
    init <- c(xm * (m + 1) / m, -1)
  }
  init_phi <- theta_to_phi(init)
  return(list(init = init, se = se, cov_mtx = cov_mtx, init_phi = init_phi,
              se_phi = se_phi))
}

# ============================= gev_init ==================================

gev_init <- function(data, m = length(data[!is.na(data)]),
                     x1 = min(data, na.rm = TRUE),
                     xm = max(data, na.rm = TRUE),
                     sum_gev = sum(data, na.rm = TRUE), xi_eq_zero = FALSE,
                     init_ests = NULL){
  # Initial estimates for Generalized extreme value distribution parameters
  #
  # Calculates initial estimates and estimated standard errors (SEs) for the
  # generalized extreme value parameters mu, sigma and xi based on an
  # assumed random sample from this distribution.  Also calculates
  # initial estimates and estimated standard errors for
  #   phi1 = mu,
  #   phi2 = (sigma + xi * (x1 - mu) ) / sr
  #   and
  #   phi3 = (sigma + xi * (xm - mu) ) / sr
  #   where x1 is the sample minimum, xm is the sample maximum and
  #   sr = (xm - x1) ^ (1/2).
  #
  # @param data A numeric vector containing sample maxima.
  # @param m A numeric scalar.  The sample size, i.e. the length of data.
  # @param x1 A numeric scalar. The sample minimum.
  # @param xm A numeric scalar. The sample maximum.
  # @param sum_gev A numeric scalar. The sum of the sample values.
  # @param xi_eq_zero A logical scalar.  If TRUE assume that the shape
  #   parameter xi = 0.
  # @param init_ests A numeric vector.  Initial estimate of
  #   theta = (mu, sigma, xi).  If supplied \code{gev_init()} just
  #   returns the corresponding initial estimate of phi = (phi1, phi2, phi3).
  # @details The main aim is to calculate an admissible estimate of
  #   theta = (mu, sigma, xi), i.e. one at which the log-likelihood is finite
  #   (necessary for the posterior log-density to be finite) at the estimate,
  #   and associated estimated SEs. These are converted into estimates and SEs
  #   for phi.  The latter can be used to set values of \code{min_phi} and
  #   \code{max_phi} for input to \code{find_lambda}.
  #
  #   In the default setting (\code{xi_eq_zero = FALSE} and
  #   \code{init_ests = NULL}) the methods tried are Maximum Likelihood
  #   Estimation (MLE) (Smith, 1985), Probability-Weighted Moments (PWM)
  #   (Hosking and Wallis, 1987) and Linear Combinations of Ratios of Spacings
  #   (LRS) (Reiss and Thomas, 2007, page 111) in that order.
  #
  #   For xi < -1 the likelihood is unbounded, MLE may fail when xi is not
  #   greater than -0.5 and the observed Fisher information for (mu, sigma, xi) has
  #   finite variance only if xi > -0.25.  We use the ML estimate provided that
  #   the estimate of xi returned from \code{gev_mle} is greater than -1. We only
  #   use the SE if the MLE of xi is greater than -0.25.
  #
  #   If either the MLE or the SE are not OK then we try PWM.  We use the PWM
  #   estimate only if is admissible, and the MLE was not OK.  We use the PWM SE,
  #   but this will be \code{c(NA, NA, NA)} is the PWM estimate of xi is > 1/2.
  #   If the estimate is still not OK then we try LRS.  As a last resort, which
  #   will tend to occur only when xi is strongly negative, we set xi = -1 and
  #   estimate sigma conditional on this.
  # @return If \code{init_ests} is not supplied by the user, a list is returned
  #   with components
  #     \item{init}{A numeric vector. Initial estimates of mu, sigma
  #      and xi.}
  #     \item{se}{A numeric vector. Estimated standard errors of
  #      mu, sigma and xi.}
  #     \item{init_phi}{A numeric vector. Initial estimates of
  #     phi1 = mu,
  #     phi2 = (sigma + xi * (x1 - mu) ) / sr
  #     and
  #     phi3 = (sigma + xi * (xm - mu) ) / sr
  #     where x1 is the minium of \code{data}, xm is the maximum of \code{data}
  #     and sr = (xm - x1) ^ (1/2).}
  #     \item{se_phi}{A numeric vector. Estimated standard errors of
  #      phi1, phi2 and phi3.}
  #   If \code{init_ests} is supplied then only the numeric vector
  #   \code{init_phi} is returned.
  # @references Smith, R. L. (1985) Maximum likelihood estimation in a class of
  #   nonregular cases.  Biometrika, 72(1), 67-90.
  #   \url{https://doi.org/10.1093/biomet/72.1.67}.
  # @references Hosking, J. R. M. and Wallis, J. R. and Wood, E. F. (1985)
  #   Estimation of the Generalized Extreme-Value Distribution by the Method of
  #   Probability-Weighted Moments. Technometrics, 27(3), 251-261.
  #   \url{https://doi.org/10.1080/00401706.1985.10488049}.
  # @references Reiss, R.-D., Thomas, M. (2007) Statistical Analysis of Extreme
  #   Values with Applications to Insurance, Finance, Hydrology and Other
  #   Fields. Birkhauser. \url{https://doi.org/10.1007/978-3-7643-7399-3}.
  # @seealso \code{\link{gp_init}} to calculate initial estimates for
  #   Generalized Pareto distribution parameters.
  #
  sr <- sqrt(xm - x1)
  theta_to_phi <- function(theta) {
    phi1 <- theta[1]
    phi2 <- (theta[2] + theta[3] * (x1 - theta[1])) / sr
    phi3 <- (theta[2] + theta[3] * (xm - theta[1])) / sr
    c(phi1, phi2, phi3)
  }
  #
  if (!is.null(init_ests)) {
    return(theta_to_phi(init_ests))
  }
  #
  # Initial estimates based on the Gumbel (xi = 0) case.
  sigma_init <- stats::sd(data) * sqrt(6) / pi
  mu_init <- mean(data) - sigma_init * 0.5772156649015323
  gum_init <- c(mu_init, sigma_init, 0)
  if (xi_eq_zero) {
    init <- gum_init
    init_phi <- theta_to_phi(init)
    fish <- gev_fish(gum_init)
    cov_mtx <- solve(fish) / m
    se <- sqrt(diag(cov_mtx))
    row1 <- c(1, 0, 0)
    row2 <- c(-gum_init[3], 1, x1 - gum_init[1]) / sr
    row3 <- c(-gum_init[3], 1, xm - gum_init[1]) / sr
    mat <- rbind(row1, row2, row3)
    var_phi <- mat %*% cov_mtx %*% t(mat)
    se_phi <- sqrt(diag(var_phi))
    return(list(init = init, se = se, init_phi = init_phi, se_phi = se_phi))
  }
  #
  ests_ok <- ses_ok <- FALSE
  # First try MLE
  temp <- gev_mle(gum_init, data = data, m = m, sum_gev = sum_gev)
  init <- temp$mle
  # Check that the MLE is OK
  # If mle for xi is > -1 and likelihood is non-zero at init
  if (init[3] > -1 & !is.infinite(temp$nllh)){
    ests_ok <- TRUE
    # Check whether or not we should use the SE
    if (init[3] > -0.25 & !is.null(temp$cov)){
      cov_mtx <- temp$cov
      se <- sqrt(diag_pos(cov_mtx))
      row1 <- c(1, 0, 0)
      row2 <- c(-init[3], 1, x1 - init[1]) / sr
      row3 <- c(-init[3], 1, xm - init[1]) / sr
      mat <- rbind(row1, row2, row3)
      var_phi <- mat %*% cov_mtx %*% t(mat)
      if (all(diag(var_phi) > 0)){
        se_phi <- sqrt(diag(var_phi))
        ses_ok <- TRUE
      }
    }
  }
  if (!ests_ok | !ses_ok) {
    # Try PWM
    pwm <- gev_pwm(data)
    se <- pwm$se
    cov_mtx <- pwm$cov
    row1 <- c(1, 0, 0)
    row2 <- c(-pwm$est[3], 1, x1 - pwm$est[1]) / sr
    row3 <- c(-pwm$est[3], 1, xm - pwm$est[1]) / sr
    mat <- rbind(row1, row2, row3)
    var_phi <- mat %*% cov_mtx %*% t(mat)
    se_phi <- sqrt(diag(var_phi))
    # Note: se and se_phi will be NA if pwm$est[3] > 1/2
    check <- gev_loglik(pars = pwm$est, data = data, m = m, x1 = x1, xm = xm,
                        sum_gev = sum_gev)
    # If MLE wasn't OK and PWM estimate is OK then use PWM estimate
    if (!ests_ok & init[3] > -1 & !is.infinite(check)) {
      init <- pwm$est
      ests_ok <- TRUE
    }
  }
  # If estimate is not OK then try LRS
  if (!ests_ok) {
    init <- gev_lrs(data)
    check <- gev_loglik(pars = pwm$est, data = data, m = m, x1 = x1, xm = xm,
                        sum_gev = sum_gev)
    if (init[3] > -1 & !is.infinite(check)) {
      ests_ok <- TRUE
    }
  }
  # If we get here without ests_ok = TRUE then the posterior mode is
  # probably very close to xi = -1.  Therefore, use xi = -1 and
  # estimates of mu and sigma conditional on xi = -1.
  if (!ests_ok) {
    xbar <- sum_gev / m
    mu <- ((m + 1) * xm / m + xbar) / 2
    sigma <- ((m + 1) * xm / m - xbar) / 2
    init <- c(mu, sigma, -1)
  }
  init_phi <- theta_to_phi(init)
  return(list(init = init, se = se, init_phi = init_phi, se_phi = se_phi))
}

os_init <- function(data, min_data = apply(data, 1, min, na.rm = TRUE),
                    max_data = apply(data, 1, max, na.rm = TRUE),
                    nos = sum(!is.na(data)), x1 = min(data, na.rm = TRUE),
                    xm = max(data, na.rm = TRUE),
                    sum_os = sum(data, na.rm = TRUE), xi_eq_zero = FALSE,
                    init_ests = NULL){
  sr <- sqrt(xm - x1)
  theta_to_phi <- function(theta) {
    phi1 <- theta[1]
    phi2 <- (theta[2] + theta[3] * (x1 - theta[1])) / sr
    phi3 <- (theta[2] + theta[3] * (xm - theta[1])) / sr
    c(phi1, phi2, phi3)
  }
  #
  if (!is.null(init_ests)) {
    return(theta_to_phi(init_ests))
  }
  #
  # Initial estimates based on the Gumbel (xi = 0) case, with max_data,
  # that is, the block maxima used as data.
  sigma_init <- stats::sd(max_data) * sqrt(6) / pi
  mu_init <- mean(max_data) - sigma_init * 0.5772156649015323
  gum_init <- c(mu_init, sigma_init, 0)
  if (xi_eq_zero) {
    init <- gum_init
    # Find the MLE for the xi =0 (Gumbel case).
    temp <- os_mle(gum_init[1:2], data = data, min_data = min_data, nos = nos,
                   sum_os = sum_os, gumbel = TRUE)
    init <- temp$mle
    init_phi <- theta_to_phi(init)
    cov_mtx <- temp$cov
    se <- sqrt(diag_pos(cov_mtx))
    row1 <- c(1, 0, 0)
    row2 <- c(-gum_init[3], 1, x1 - gum_init[1]) / sr
    row3 <- c(-gum_init[3], 1, xm - gum_init[1]) / sr
    mat <- rbind(row1, row2, row3)
    var_phi <- mat %*% cov_mtx %*% t(mat)
    se_phi <- sqrt(diag(var_phi))
    return(list(init = init, se = se, init_phi = init_phi, se_phi = se_phi))
  }
  #
  ests_ok <- ses_ok <- FALSE
  # Find the MLE.
  temp <- os_mle(gum_init, data = data, min_data = min_data, nos = nos,
                 sum_os = sum_os)
  init <- temp$mle
  # Check that the MLE is OK
  # If mle for xi is > -1 and likelihood is non-zero at init
  if (init[3] > -1 & !is.infinite(temp$nllh)){
    ests_ok <- TRUE
    # Check whether or not we should use the SE
    if (init[3] > -0.5 & !is.null(temp$cov)){
      cov_mtx <- temp$cov
      se <- sqrt(diag_pos(cov_mtx))
      row1 <- c(1, 0, 0)
      row2 <- c(-init[3], 1, x1 - init[1]) / sr
      row3 <- c(-init[3], 1, xm - init[1]) / sr
      mat <- rbind(row1, row2, row3)
      var_phi <- mat %*% cov_mtx %*% t(mat)
      if (all(diag(var_phi) > 0)){
        se_phi <- sqrt(diag(var_phi))
        ses_ok <- TRUE
      }
    }
  }
  init_phi <- theta_to_phi(init)
  return(list(init = init, se = se, init_phi = init_phi, se_phi = se_phi))
}

