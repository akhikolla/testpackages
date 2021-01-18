#' Internal rust functions
#'
#' Internal rust functions
#' @details
#' These functions are not intended to be called by the user.
#' @name rust-internal
#' @keywords internal
NULL

# ==================================== box_cox ================================

#' @keywords internal
#' @rdname rust-internal
box_cox <- function (x, lambda = 1, gm = 1, lambda_tol = 1e-6) {
  #
  # Computes the Box-Cox transformation of a vector.
  #
  # Args:
  #   x          : A numeric vector. (Positive) values to be Box-Cox
  #                transformed.
  #   lambda     : A numeric scalar.  Transformation parameter.
  #   gm         : A numeric scalar.  Optional scaling parameter.
  #   lambda_tol : A numeric scalar.  For abs(lambda) < lambda_tol use
  #                a Taylor series expansion.
  #
  # Returns:
  #   A numeric vector.  The transformed value
  #     (x^lambda - 1) / (lambda * gm ^ (lambda - 1))
  #
  if (abs(lambda) > lambda_tol) {
    retval <- (x ^ lambda - 1) / lambda / gm ^ (lambda - 1)
  } else {
    i <- 0:3
    retval <- sum(log(x) ^ (i + 1) * lambda ^ i / factorial(i + 1))
    retval <- retval / gm ^ (lambda - 1)
  }
  retval
}

# ================================ box_cox_vec ================================

# Version of box_cox vectorized for lambda and gm.

#' @keywords internal
#' @rdname rust-internal
box_cox_vec <- Vectorize(box_cox, vectorize.args = c("lambda", "gm"))

# =============================== optim_box_cox ===============================

#' @keywords internal
#' @rdname rust-internal
optim_box_cox <- function(x, w, lambda_range = c(-3,3), start = NULL,
                          which_lam = 1:ncol(x)) {
  #
  # Finds the optimal value of the Box-Cox transformation parameter lambda,
  # based values of a probability density function.
  #
  # Args:
  #   x            : A numeric matrix. Values at which the density is evaluated.
  #                  each row contains a combination of the ncol(x) variables.
  #                  column numbers in which_lam must contain positive values.
  #   w            : A numeric vector. Density values corresponding to each
  #                  row of x (up to proportionality).
  #   lambda_range : A numeric vector (of length 2).  Range of lambda values
  #                  over which to search.
  #   start        : A numeric vector.  Optional starting value for lambda.
  #   which_lam    : A numeric vector.  Indicates which variables to Box-Cox
  #                  transform.
  #
  # Returns: a list containing
  #   lambda       : A numeric vector.  The optimal value of lambda.
  #   gm           : A numeric vector.  Geometric mean of x weighted by w.
  #
  n_var <- ncol(x)
  neg_loglik <- function(lambda_in, x, gm, w, which_lam) {
    lambda <- rep(1, n_var)
    lambda[which_lam] <- lambda_in
    for (j in 1:n_var) {
      x[, j] <- box_cox(x = x[, j], lambda = lambda[j], gm = gm[j])
    }
    (nrow(x) / 2) * log(det(stats::cov.wt(x, w, method = "ML")$cov))
  }
  w_n <- w / mean(w)
  n_lam <- length(which_lam)
  #
  if (any(x[, which_lam] <= 0)) {
    stop("All values must be > 0")
  }
  gm_w <- rep(1, n_var)
  x_mat <- x[, which_lam, drop = FALSE]
  gm_w[which_lam] <- apply(x_mat, 2, function(x) exp(mean(w_n * log(x))))
  #
  skew_wt <- function(x, w) {
    vp <- function(p) mean(w^p)
    v1 <- vp(1)
    xbar <- mean(w * x) / v1
    mp <- function(p) mean(w * (x - xbar) ^ p) / v1
    mp(3) / mp(2)^(3 / 2)
  }
  if (is.null(start)) {
    start <- 1 - apply(x, 2, skew_wt, w = w) / 2
  }
  if (n_lam == 1L) {
    method <- "Brent"
  } else {
    method <- "L-BFGS-B"
  }
  lower <- lambda_range[1]
  upper <- lambda_range[2]
  start <- start[which_lam]
  ret <- stats::optim(start, neg_loglik, method = method, x = x, gm = gm_w,
                      w = w_n, lower = lower, upper = upper,
                      which_lam = which_lam)
  lambda <- rep(1L, n_var)
  lambda[which_lam] <- ret$par
  ret$par <- lambda
  if (ret$convergence != 0) {
    warning(paste("Convergence failure: return code =", ret$convergence))
  }
  list(lambda = ret$par, gm = gm_w)
}

# ================================ n_grid_fun =================================

#' @keywords internal
#' @rdname rust-internal
n_grid_fn <- function(d) {
  # Sets a default values value of n_grid.
  #
  # Args:
  #   d          : A numeric scalar.  Dimension of the target density.
  #
  # Returns:
  #   A numeric scalar.  The value of n_grid.
  #
  return(ceiling(2501 ^ (1 / d)))
}

# =============================== init_ps_calc ================================

#' @keywords internal
#' @rdname rust-internal
init_psi_calc <- function(phi_to_psi, phi, lambda, gm, w, which_lam){
  # Estimates the mode and standard deviation of the Box-Cox transformed
  # target density.
  #
  # Args:
  #   phi_to_psi : A function.  The Box-Cox transformation function.
  #   phi        : A numeric matrix.  n_grid by d matrix of values of phi
  #   lambda     : A numeric vector.  The value of the Box-Cox transformation
  #                parameter lambda.
  #   gm         : A numeric vector.  The value of the Box-Cox scale
  #                parameter gm.
  #   w          : A numeric vector. The values of the target density at each
  #                value of the variables in phi.
  #   which_lam  : A numeric vector.  Indicates which variables to Box-Cox
  #                transform.
  #
  # Returns: A list containing
  #   init_psi : estimate of the mode of the Box-Cox transformed target
  #              density.
  #   sd_psi   : estimate of the standard deviation of the Box-Cox transformed
  #              target density.
  #
  psi <- phi_to_psi(phi)
  d <- ncol(phi)
  log_jac <- matrix(0, nrow(phi), d)
  for (j in which_lam) {
    log_jac[, j] <- (lambda[j] - 1) * log(phi[, j] / gm[j])
  }
  log_jac <- rowSums(log_jac)
  w_psi <- w * exp(-log_jac)
  temp_cov <- stats::cov.wt(psi, w_psi, method="ML")
  init_psi <- as.numeric(psi[which.max(w_psi), ])
  sd_psi <- sqrt(diag(temp_cov$cov))
  list(init_psi = init_psi, sd_psi = sd_psi)
}

# ============================= log_gpd_mdi_prior =============================

#' @keywords internal
#' @rdname rust-internal
log_gpd_mdi_prior <- function(pars) {
  # Generalized Pareto MDI prior log-density
  #
  # Calculates the log-density of a Maximal Data Information (MDI) prior
  # for generalized Pareto parameters.  The prior density is truncated,
  # (to xi >= -1), to produce a posterior density that is proper.
  #
  # Args:
  #   pars : A numeric vector containing the values of the generalized Pareto
  #          parameters sigma and xi.
  #
  # Returns:
  #   A numeric scalar. The value of the prior log-density.
  #
  if (pars[1] <= 0 | pars[2] < -1) {
    return(-Inf)
  }
  return(-log(pars[1]) - pars[2] - 1)
}

# ================================ gpd_loglik =================================

#' @keywords internal
#' @rdname rust-internal
gpd_loglik <- function(pars, gpd_data, m, xm, sum_gp) {
  # Generalized Pareto log-likelihood
  #
  # Calculates the log-likelihood for a random sample from a Generalized Pareto.
  #
  # Args:
  #   pars     : A numeric vector containing the values of the generalized
  #              Pareto parameters sigma and xi.
  #   gpd_data : A numeric vector containing positive sample values.
  #   m        : A numeric scalar.  The sample size: the length of gpd_data.
  #   xm       : A numeric scalar. The sample maximum.
  #   sum_gp   : The sum of the sample values.
  #
  # Returns:
  #   A numeric scalar. The value of the log-likelihood.
  #
  if (pars[1] <= 0 | pars[2] <= -pars[1] / xm)
    return(-Inf)
  sdat <- gpd_data / pars[1]
  zz <- 1 + pars[2] * sdat
  if (abs(pars[2]) > 1e-6) {
    val <- -m * log(pars[1]) - (1 + 1 / pars[2]) * sum(log(zz))
  } else {
    i <- 1:4
    g_fun <- function(x) {
      t1 <- x ^ i
      t2 <- (i * x - i - 1)
      sum((-1) ^ i * t1 * t2 * pars[2] ^ i / i / (i + 1))
    }
    val <- -m * log(pars[1]) - sum_gp / pars[1] - sum(sapply(sdat, g_fun))
  }
  return(val)
}

# ================================== gpd_mle ==================================

#' @keywords internal
#' @rdname rust-internal
gpd_mle <- function(gpd_data) {
  # Maximum likelihood estimation for the generalized Pareto distribution
  #
  # Performs maximum likelihood estimation for the generalized Pareto
  # distribution.  Uses the function \code{gpdmle} associated with
  # Grimshaw (1993), which returns MLEs of sigma and k = - \xi.
  #
  # Grimshaw, S. D. (1993) Computing Maximum Likelihood Estimates
  #   for the Generalized Pareto Distribution.  Technometrics, 35(2), 185-191.
  #   and Computing (1991) 1, 129-133. https://doi.org/10.1007/BF01889987.
  #
  # Args:
  #   gpd_data : A numeric vector containing positive values, assumed to be a
  #             random sample from a generalized Pareto distribution.
  #
  # Returns:
  #   A list with components
  #     mle  : A numeric vector.  MLEs of GP parameters sigma and xi.
  #     nllh : A numeric scalar.  The negated log-likelihood at the MLE.
  #
  temp <- list()
  # Use revdbayes function if revdbayes is available, otherwise use fallback
  if (requireNamespace("revdbayes", quietly = TRUE)) {
    # Call Grimshaw (1993) function, note: k is -xi, a is sigma
    pjn <- revdbayes::grimshaw_gp_mle(gpd_data)
    temp$mle <- c(pjn$a, -pjn$k)  # mle for (sigma,xi)
  } else {
    temp <- fallback_gp_mle(init = c(mean(gpd_data), 0), gpd_data = gpd_data,
                            m = length(gpd_data), xm = max(gpd_data),
                            sum_gp = sum(gpd_data))
  }
  sc <- rep(temp$mle[1], length(gpd_data))
  xi <- temp$mle[2]
  temp$nllh <- sum(log(sc)) + sum(log(1 + xi * gpd_data / sc) * (1 / xi + 1))
  return(temp)
}

# ============================== fallback_gp_mle ==============================

#' @keywords internal
#' @rdname rust-internal
fallback_gp_mle <- function(init, ...){
  x <- stats::optim(init, gpd_loglik, ..., control = list(fnscale = -1),
                    hessian = FALSE)
  temp <- list()
  temp$mle <- x$par
  temp$nllh <- -x$value
  return(temp)
}

# =============================== gpd_obs_info ================================

#' @keywords internal
#' @rdname rust-internal
gpd_obs_info <- function(gpd_pars, y) {
  # Observed information for the generalized Pareto distribution
  #
  # Calculates the observed information matrix for a random sample \code{y}
  # from the generalized Pareto distribution, i.e. the negated Hessian matrix of
  # the generalized Pareto log-likelihood, evaluated at \code{gpd_pars}.
  #
  # Args:
  #   gpd_pars : A numeric vector. Parameters sigma and xi of the
  #   generalized Pareto distribution.
  #   y       : A numeric vector. A sample of positive data values.
  #
  # Returns:
  #   A 2 by 2 numeric matrix.  The observed information matrix.
  #
  s <- gpd_pars[1]
  x <- gpd_pars[2]
  i <- matrix(NA, 2, 2)
  i[1,1] <- -sum((1 - (1 + x) * y * (2 * s + x * y) / (s + x * y) ^ 2) / s ^ 2)
  i[1,2] <- i[2,1] <- -sum(y * (1 - y / s) / (1 + x * y / s) ^ 2 / s ^ 2)
  i[2,2] <- sum(2 * log(1 + x * y / s) / x ^ 3 - 2 * y / (s + x * y) / x ^ 2 -
                  (1 + 1 / x) * y ^ 2 / (s + x * y) ^ 2)
  return(i)
}

# =================================== find_a ==================================

#' @keywords internal
#' @rdname rust-internal
find_a <-  function(neg_logf_rho, init_psi, d, r, lower, upper, algor,
                    method, control, shoof, ...) {
  #
  # Finds the value of a(r).
  #
  # Args:
  #   neg_logf_rho : A function. Negated target log-density function.
  #   init_psi     : A numeric scalar.  Initial value of psi.
  #   d            : A numeric scalar. Dimension of f.
  #   r            : A numeric scalar. Parameter of generalized
  #                  ratio-of-uniforms.
  #   lower        : A numeric vector.  Lower bounds on the arguments of logf.
  #   upper        : A numeric vector.  Upper bounds on the arguments of logf.
  #   algor        : A character scalar.  Algorithm ("optim" or "nlminb").
  #   method       : A character scalar.  Only relevant if algorithm = "optim".
  #   control      : A numeric list.  Control arguments to algor.
  #
  # Returns: a list containing
  #   the standard returns from optim or nlminb
  #   hessian: the estimated hessian of neg_logf_rho/(d*r+1) at its minimum.
  #
  big_val <- Inf
  big_finite_val <- 10 ^ 10
  #
  # Function to minimize to find a(r).
  # Use a weird argument name (._psi) to avoid a potential argument-matching
  # problem when nlminb() is used.
  a_obj <- function(._psi, ...) {
    # Avoid possible problems with nlminb calling function with NaNs.
    # Indicative of solution on boundary, which is OK in the current context.
    # See https://stat_ethz.ch/pipermail/r-help/2015-May/428488.html
    if (any(is.na(._psi))) return(big_val)
    neg_logf_rho(._psi, ...) / (d * r + 1)
  }
  # Function for L-BFGS-B and Brent, which don't like Inf or NA
  a_obj_no_inf <- function(._psi, ...) {
    check <- neg_logf_rho(._psi, ...) / (d * r + 1)
    if (is.infinite(check)) {
      check <- big_finite_val
    }
    check
  }
  #
  if (algor == "optim") {
    if (method %in% c("L-BFGS-B","Brent")) {
      temp <- stats::optim(par = init_psi, fn = a_obj_no_inf, ...,
                           control = control, hessian = FALSE, method = method,
                           lower = lower, upper = upper)
    } else {
      temp <- stats::optim(par = init_psi, fn = a_obj, ..., control = control,
                           hessian = FALSE, method = method)
      # Sometimes Nelder-Mead fails if the initial estimate is too good.
      # ... so avoid non-zero convergence indicator by using L-BFGS-B instead.
      if (temp$convergence > 0) {
        # Start a little away from the optimum, to avoid erroneous
        # convergence warnings, using init_psi as a benchmark
        # If init_psi = temp$par then multiply temp$par by 1 - shoof
        if (sum(abs(init_psi - temp$par)) > .Machine$double.eps) {
          new_start <- shoof * init_psi + (1 - shoof) * temp$par
        } else {
          new_start <- temp$par * (1 - shoof)
        }
        temp <- stats::optim(par = new_start, fn = a_obj_no_inf, ...,
                             control = control, hessian = FALSE,
                             method = "L-BFGS-B", lower = lower, upper = upper)
      }
      # In some cases optim with method = "L-BFGS-B" may reach its iteration
      # limit without the convergence criteria being satisfied.  Then try
      # nlminb as a further check, but don't use the control argument in
      # case of conflict between optim() and nlminb().
      if (temp$convergence > 0) {
        temp <- stats::nlminb(start = new_start, objective = a_obj, ...,
                              lower = lower, upper = upper)
      }
    }
    # Try to calculate Hessian, but avoid problems if an error is produced.
    # An error may occur if the MAP estimate is very close to a parameter
    # boundary.
    temp$hessian <- try(stats::optimHess(par = temp$par, fn = a_obj, ...),
                        silent = TRUE)
    return(temp)
  }
  # If we get to here we are using nlminb() ...
  temp <- stats::nlminb(start = init_psi, objective = a_obj, ...,
                        lower = lower, upper = upper, control = control)
  # Sometimes nlminb isn't sure that it has found the minimum when in fact
  # it has.  Try to check this, and avoid a non-zero convergence indicator
  # by using optim with method="L-BFGS-B", again starting from new_start,
  # but don't use the control argument in case of conflict between
  # optim() and nlminb().
  if (temp$convergence > 0) {
    if (sum(abs(init_psi - temp$par)) > .Machine$double.eps) {
      new_start <- shoof * init_psi + (1 - shoof) * temp$par
    } else {
      new_start <- temp$par * (1 - shoof)
    }
    temp <- stats::optim(par = new_start, fn = a_obj_no_inf, ...,
                         hessian = FALSE, method = "L-BFGS-B", lower = lower,
                         upper = upper)
  }
  # Try to calculate Hessian, but avoid problems if an error is produced.
  # An error may occur if the MAP estimate is very close to a parameter
  # boundary.
  temp$hessian <- try(stats::optimHess(par = temp$par, fn = a_obj, ...),
                      silent = TRUE)
  return(temp)
}

# ================================== find_bs ==================================

#' @keywords internal
#' @rdname rust-internal
find_bs <-  function(f_rho, d, r, lower, upper, f_mode, ep, vals, conv, algor,
                     method, control, shoof, ...) {
  # Finds the values of b-(r) and b+(r).
  #
  # Args:
  #   f_rho        : A function.  Target probability density function.
  #   d            : A numeric scalar. Dimension of f.
  #   r            : A numeric scalar. Parameter of generalized
  #                  ratio-of-uniforms.
  #   lower        : A numeric vector.  Lower bounds on the arguments of logf.
  #   upper        : A numeric vector.  Upper bounds on the arguments of logf.
  #   f_mode       : A numeric scalar.  The estimated mode of the target
  #                  log-density logf.
  #   ep           : A numeric scalar.  Controls initial estimates for
  #                  optimizations to find the b-bounding box parameters.
  #                  The default (ep=0) corresponds to starting at the mode of
  #                  logf small positive values of ep move the constrained
  #                  variable slightly away from the mode in the correct
  #                  direction.  If ep is negative its absolute value is used,
  #                  with no warning given.
  #   vals         : A numeric matrix.  Will contain the values of the
  #                  variables at which the ru box dimensions occur.
  #                  Row 1 already contains the values for a(r).
  #   conv         : A numeric scalar.  Will contain the covergence
  #                  indicators returned by the optimisation algorithms.
  #                  Row 1 already contains the values for a(r).
  #   algor        : A character scalar. Algorithm ("optim" or "nlminb").
  #   method       : A character scalar.  Only relevant if algorithm = "optim".
  #   control      : A numeric list. Control arguments to algor.
  #
  # Returns: a list containing
  #   l_box : A numeric vector.  Values of biminus(r), i = 1, ...d.
  #   u_box : A numeric vector.  Values of biplus(r), i = 1, ...d.
  #   vals  : as described above in Args.
  #   conv  : as described above in Args.
  #
  big_val <- Inf
  big_finite_val <- 10^10
  #
  # Functions to minimize to find biminus(r) and biplus(s), i = 1, ...,d.
  # Use a weird argument name (._rho) to avoid a potential argument-matching
  # problem when nlminb() is used.
  lower_box <- function(._rho, j, ...) {
    # Avoid possible problems with nlminb calling function with NaNs.
    # Indicative of solution on boundary, which is OK in the current context.
    # See https://stat_ethz.ch/pipermail/r-help/2015-May/428488.html
    if (any(is.na(._rho))) return(big_val)
    if (._rho[j] == 0L) return(0L)
    if (._rho[j] > 0L) return(big_val)
    if (f_rho(._rho, ...) == 0L) return(big_val)
    ._rho[j] * f_rho(._rho, ...) ^ (r / (d * r + 1))
  }
  upper_box <- function(._rho, j, ...) {
    if (any(is.na(._rho))) return(big_val)
    if (._rho[j] == 0) return(0)
    if (._rho[j] < 0) return(big_val)
    if (f_rho(._rho, ...) == 0) return(big_val)
    -._rho[j] * f_rho(._rho, ...) ^ (r / (d * r + 1))
  }
  # Functions for L-BFGS-B and Brent, which don't like Inf or NA
  lower_box_no_inf <- function(._rho, j, ...) {
    if (any(is.na(._rho))) return(big_finite_val)
    if (._rho[j] == 0) return(0)
    if (._rho[j] > 0) return(big_finite_val)
    if (f_rho(._rho, ...) == 0) return(big_finite_val)
    check <- ._rho[j] * f_rho(._rho, ...) ^ (r / (d * r + 1))
    if (is.infinite(check)) check <- big_finite_val
    check
  }
  upper_box_no_inf <- function(._rho, j, ...) {
    if (any(is.na(._rho))) return(big_finite_val)
    if (._rho[j] == 0) return(0)
    if (._rho[j] < 0) return(big_finite_val)
    if (f_rho(._rho, ...) == 0) return(big_finite_val)
    check <- -._rho[j] * f_rho(._rho, ...) ^ (r / (d * r + 1))
    if (is.infinite(check)) check <- big_finite_val
    check
  }
  l_box <- u_box <- NULL
  zeros <- rep(0,d)
  #
  # Find biminus(r) and biplus(s), i = 1, ...,d.
  #
  for (j in 1:d) {
    #
    # Find biminus(r) ----------
    #
    rho_init <- zeros
    rho_init[j] <- -ep
    t_upper <- upper - f_mode
    t_upper[j] <- 0
    if (algor == "nlminb") {
      temp <- stats::nlminb(start = rho_init, objective = lower_box, j = j,
                            ..., upper = t_upper, lower = lower - f_mode,
                            control = control)
      l_box[j] <- temp$objective
      # Sometimes nlminb isn't sure that it has found the minimum when in fact
      # it has.  Try to check this, and avoid a non-zero convergence indicator
      # by using optim with method="L-BFGS-B", starting from nlminb's solution.
      if (temp$convergence > 0) {
        if (sum(abs(rho_init - temp$par)) > .Machine$double.eps) {
          new_start <- shoof * rho_init + (1 - shoof) * temp$par
        } else {
          new_start <- temp$par * (1 - shoof)
        }
        temp <- stats::optim(par = new_start, fn = lower_box_no_inf, j = j, ...,
                             hessian = FALSE, method = "L-BFGS-B",
                             upper = t_upper, lower = lower - f_mode)
        l_box[j] <- temp$value
      }
    }
    if (algor == "optim") {
      if (method == "L-BFGS-B" | method == "Brent") {
        temp <- stats::optim(par = rho_init, fn = lower_box_no_inf, j = j, ...,
                             upper = t_upper, lower = lower - f_mode,
                             control = control, method = method,
                             hessian = FALSE)
        l_box[j] <- temp$value
      } else {
        temp <- stats::optim(par = rho_init, fn = lower_box, j = j, ...,
                             control = control, method = method,
                             hessian = FALSE)
        l_box[j] <- temp$value
        # Sometimes Nelder-Mead fails if the initial estimate is too good.
        # ... so avoid non-zero convergence indicator using L-BFGS-B instead.
        if (temp$convergence > 0) {
          if (sum(abs(rho_init - temp$par)) > .Machine$double.eps) {
            new_start <- shoof * rho_init + (1 - shoof) * temp$par
          } else {
            new_start <- temp$par * (1 - shoof)
          }
          temp <- stats::optim(par = new_start, fn = lower_box_no_inf, j = j,
                               ..., control = control, method = "L-BFGS-B",
                               hessian = FALSE,
                               upper = t_upper, lower = lower - f_mode)
          l_box[j] <- temp$value
        }
        # Check using nlminb() if optim's iteration limit is reached.
        if (temp$convergence == 1) {
          if (sum(abs(rho_init - temp$par)) > .Machine$double.eps) {
            new_start <- shoof * rho_init + (1 - shoof) * temp$par
          } else {
            new_start <- temp$par * (1 - shoof)
          }
          temp <- stats::nlminb(start = new_start, objective = lower_box,
                                j = j, ..., upper = t_upper,
                                lower = lower - f_mode)
          l_box[j] <- temp$objective
        }
      }
    }
    vals[j+1, ] <- temp$par
    conv[j+1] <- temp$convergence
    #
    # Find biplus(r) --------------
    #
    rho_init <- zeros
    rho_init[j] <- ep
    t_lower <- lower - f_mode
    t_lower[j] <- 0
    if (algor == "nlminb") {
      temp <- stats::nlminb(start = rho_init, objective = upper_box, j = j,
                            ..., lower = t_lower, upper = upper - f_mode,
                            control = control)
      u_box[j] <- -temp$objective
      # Sometimes nlminb isn't sure that it has found the minimum when in fact
      # it has.  Try to check this, and avoid a non-zero convergence indicator
      # by using optim with method="L-BFGS-B", starting from nlminb's solution.
      if (temp$convergence > 0) {
        if (sum(abs(rho_init - temp$par)) > .Machine$double.eps) {
          new_start <- shoof * rho_init + (1 - shoof) * temp$par
        } else {
          new_start <- temp$par * (1 - shoof)
        }
        temp <- stats::optim(par = new_start, fn = upper_box_no_inf, j = j,
                             ..., hessian = FALSE, method = "L-BFGS-B",
                             lower = t_lower, upper = upper - f_mode)
        u_box[j] <- -temp$value
      }
    }
    if (algor == "optim") {
      if (method == "L-BFGS-B" | method == "Brent") {
        temp <- stats::optim(par = rho_init, fn = upper_box_no_inf, j = j, ...,
                             lower = t_lower, upper = upper - f_mode,
                             control = control, method = method)
        u_box[j] <- -temp$value
      } else {
        temp <- stats::optim(par = rho_init, fn = upper_box, j = j, ...,
                             control = control, method = method)
        u_box[j] <- -temp$value
        # Sometimes Nelder-Mead fails if the initial estimate is too good.
        # ... so avoid non-zero convergence indicator using L-BFGS-B instead.
        if (temp$convergence > 0) {
          if (sum(abs(rho_init - temp$par)) > .Machine$double.eps) {
            new_start <- shoof * rho_init + (1 - shoof) * temp$par
          } else {
            new_start <- temp$par * (1 - shoof)
          }
          temp <- stats::optim(par = new_start, fn = upper_box_no_inf, j = j,
                               ..., control = control, method = "L-BFGS-B",
                               lower = t_lower, upper = upper - f_mode)
          u_box[j] <- -temp$value
        }
        # Check using nlminb() if optim's iteration limit is reached.
        if (temp$convergence == 1) {
          if (sum(abs(rho_init - temp$par)) > .Machine$double.eps) {
            new_start <- shoof * rho_init + (1 - shoof) * temp$par
          } else {
            new_start <- temp$par * (1 - shoof)
          }
          temp <- stats::nlminb(start = new_start, objective = upper_box,
                                j = j, ..., lower = t_lower,
                                upper = upper - f_mode)
          u_box[j] <- -temp$objective
        }
      }
    }
    vals[j+d+1, ] <- temp$par
    conv[j+d+1] <- temp$convergence
  }
  return(list(l_box = l_box, u_box = u_box, vals = vals, conv = conv))
}

# ================================ cpp_find_a =================================

#' @keywords internal
#' @rdname rust-internal
cpp_find_a <-  function(init_psi, lower, upper, algor, method, control,
                        a_obj_fun, ru_args, shoof) {
  #
  # Finds the value of a(r).
  #
  # Args:
  #   init_psi  : A numeric scalar.  Initial value of psi.
  #   lower     : A numeric vector.  Lower bounds on the arguments of logf.
  #   upper     : A numeric vector.  Upper bounds on the arguments of logf.
  #   algor     : A character scalar.  Algorithm ("optim" or "nlminb").
  #   method    : A character scalar.  Only relevant if algorithm = "optim".
  #   control   : A numeric list.  Control arguments to algor.
  #   a_obj_fun : The function to be minimized to find a(r).
  #   ru_args   : A numeric list, containing:
  #     d         : A numeric scalar. Dimension of f.
  #     r         : A numeric scalar. Generalized ratio-of-uniforms parameter.
  #     psi_mode  : A numeric vector.  Latterly this will contain the estimated
  #                 mode of the target density after transformation but prior
  #                 to mode relation. Equal to rep(0, d) at this stage.
  #     rot_mat   : A numeric matrix.  Latterly this will contain the rotation
  #                 matrix.  Equal to diag(d) at this stage.
  #     hscale    : A numeric scalar.  Scales the target log-density.
  #                 Equal to logf evaluated at init at this stage.
  #     which_lam : A vector of integers indication which components of lambda
  #                 are Box-Cox transformed. Only present if trans = "BC".
  #     lambda    : Box-Cox transformation parameters.
  #                 Only present if trans = "BC".
  #     gm        : Box-Cox scale parameters.  Only present if trans = "BC".
  #     con       : lambda * gm ^(lambda - 1).  Only present if trans = "BC".
  #     logf      : A pointer to the (original) target log-density function.
  #     pars      : A numeric list of additional arguments to logf.
  #
  # Returns: a list containing
  #   the standard returns from optim or nlminb
  #   hessian: the estimated hessian of -cpp_logf_rho/(d*r+1) at its minimum.
  #
  big_val <- 10 ^ 10
  #
  if (algor == "optim") {
    if (method == "L-BFGS-B" | method == "Brent") {
      add_args <- list(par = init_psi, fn = a_obj_fun, method = method,
                       control = control, lower = lower, upper = upper,
                       big_val = big_val)
      temp <- do.call(stats::optim, c(ru_args, add_args))
    } else {
      add_args <- list(par = init_psi, fn = a_obj_fun, method = method,
                       control = control, big_val = Inf)
      temp <- do.call(stats::optim, c(ru_args, add_args))
      # Sometimes Nelder-Mead fails if the initial estimate is too good.
      # ... so avoid non-zero convergence indicator using L-BFGS-B instead.
      if (temp$convergence > 0) {
        # Start a little away from the optimum, to avoid erroneous
        # convergence warnings, using init_psi as a benchmark
        # If init_psi = temp$par then multiply temp$par by  1 - shoof
        if (sum(abs(init_psi - temp$par)) > .Machine$double.eps) {
          new_start <- shoof * init_psi + (1 - shoof) * temp$par
        } else {
          new_start <- temp$par * (1 - shoof)
        }
        add_args <- list(par = new_start, fn = a_obj_fun, method = "L-BFGS-B",
                         control = control, big_val = big_val,
                         lower = lower, upper = upper)
        temp <- do.call(stats::optim, c(ru_args, add_args))
      }
      # In some cases optim with method = "L-BFGS-B" may reach its iteration
      # limit without the convergence criteria being satisfied.  Then try
      # nlminb as a further check, but don't use the control argument in
      # case of conflict between optim() and nlminb().
      if (temp$convergence > 0) {
        add_args <- list(start = new_start, objective = a_obj_fun,
                         lower = lower, upper = upper, big_val = Inf)
        temp <- do.call(stats::nlminb, c(ru_args, add_args))
      }
    }
  } else {
    add_args <- list(start = init_psi, objective = a_obj_fun,
                     control = control, lower = lower, upper = upper,
                     big_val = Inf)
    temp <- do.call(stats::nlminb, c(ru_args, add_args))
    # Sometimes nlminb isn't sure that it has found the minimum when in fact
    # it has.  Try to check this, and avoid a non-zero convergence indicator
    # by using optim with method="L-BFGS-B", again starting from new_start.
    if (temp$convergence > 0) {
      if (sum(abs(init_psi - temp$par)) > .Machine$double.eps) {
        new_start <- shoof * init_psi + (1 - shoof) * temp$par
      } else {
        new_start <- temp$par * (1 - shoof)
      }
      add_args <- list(par = new_start, fn = a_obj_fun, hessian = FALSE,
                       method = "L-BFGS-B", big_val = big_val,
                       lower = lower, upper = upper)
      temp <- do.call(stats::optim, c(ru_args, add_args))
    }
  }
  # Try to calculate Hessian, but avoid problems if an error is produced.
  # An error may occur if the MAP estimate is very close to a parameter
  # boundary.
  add_args <- list(par = temp$par, fn = a_obj_fun, big_val = Inf)
  temp$hessian <- try(do.call(stats::optimHess, c(ru_args, add_args)),
                      silent = TRUE)
  return(temp)
}

# ================================ cpp_find_bs ================================

#' @keywords internal
#' @rdname rust-internal
cpp_find_bs <-  function(lower, upper, ep, vals, conv, algor, method,
                         control, lower_box_fun, upper_box_fun, ru_args,
                         shoof) {
  # Finds the values of b-(r) and b+(r).
  #
  # Args:
  #   lower         : A numeric vector.  Lower bounds on the arguments of logf.
  #   upper         : A numeric vector.  Upper bounds on the arguments of logf.
  #   ep            : A numeric scalar.  Controls initial estimates for
  #                   optimizations to find the b-bounding box parameters.
  #                   The default (ep=0) corresponds to starting at the mode of
  #                   logf small positive values of ep move the constrained
  #                   variable slightly away from the mode in the correct
  #                   direction.  If ep is negative its absolute value is used,
  #                   with no warning given.
  #   vals          : A numeric matrix.  Will contain the values of the
  #                   variables at which the ru box dimensions occur.
  #                   Row 1 already contains the values for a(r).
  #   conv          : A numeric scalar.  Will contain the covergence
  #                   indicators returned by the optimisation algorithms.
  #                   Row 1 already contains the values for a(r).
  #   algor         : A character scalar. Algorithm ("optim" or "nlminb").
  #   method        : A character scalar.  Only relevant if algor = "optim".
  #   control       : A numeric list. Control arguments to algor.
  #   lower_box_fun : The function to be minimized to find b-(r).
  #   upper_box_fun : The function to be minimized to find b+(r).
  #   ru_args       : A numeric list, containing:
  #     d         : A numeric scalar. Dimension of f.
  #     r         : A numeric scalar. Generalized ratio-of-uniforms parameter.
  #     psi_mode  : A numeric vector.  The estimated mode of the target
  #                 density after transformation but prior to mode relocation.
  #     rot_mat   : A numeric matrix.  Rotation matrix (equal to the identity
  #                 matrix if rotate = FALSE).
  #     hscale    : A numeric scalar.  Scales the target log-density.
  #     which_lam : A vector of integers indication which components of lambda
  #                 are Box-Cox transformed. Only present if trans = "BC".
  #     lambda    : Box-Cox transformation parameters.
  #                 Only present if trans = "BC".
  #     gm        : Box-Cox scale parameters.  Only present if trans = "BC".
  #     con       : lambda * gm ^(lambda - 1).  Only present if trans = "BC".
  #     logf      : A pointer to the (original) target log-density function.
  #     pars      : A numeric list of additional arguments to logf.
  #
  # Returns: a list containing
  #   l_box : A numeric vector.  Values of biminus(r), i = 1, ...d.
  #   u_box : A numeric vector.  Values of biplus(r), i = 1, ...d.
  #   vals  : as described above in Args.
  #   conv  : as described above in Args.
  #
  big_val <- 10 ^ 10
  f_mode <- ru_args$psi_mode
  d <- ru_args$d
  #
  l_box <- u_box <- NULL
  zeros <- rep(0, d)
  #
  # Find biminus(r) and biplus(s), i = 1, ...,d.
  #
  for (j in 1:d) {
    #
    # Find biminus(r) ----------
    #
    rho_init <- zeros
    rho_init[j] <- -ep
    t_upper <- upper - f_mode
    t_upper[j] <- 0
    if (algor == "nlminb") {
      add_args <- list(start = rho_init, objective = lower_box_fun,
                       upper = t_upper, lower = lower - f_mode, j = j - 1,
                       control = control, big_val = Inf)
      temp <- do.call(stats::nlminb, c(ru_args, add_args))
      l_box[j] <- temp$objective
      # Sometimes nlminb isn't sure that it has found the minimum when in fact
      # it has.  Try to check this, and avoid a non-zero convergence indicator
      # by using optim with method="L-BFGS-B", starting from nlminb's solution.
      if (temp$convergence > 0) {
        if (sum(abs(rho_init - temp$par)) > .Machine$double.eps) {
          new_start <- shoof * rho_init + (1 - shoof) * temp$par
        } else {
          new_start <- temp$par * (1 - shoof)
        }
        add_args <- list(par = new_start, fn = lower_box_fun, j = j - 1,
                         method = "L-BFGS-B", big_val = big_val,
                         upper = t_upper, lower = lower - f_mode)
        temp <- do.call(stats::optim, c(ru_args, add_args))
        l_box[j] <- temp$value
      }
    }
    if (algor == "optim") {
      # L-BFGS-B and Brent don't like Inf or NA
      if (method == "L-BFGS-B" | method == "Brent") {
        add_args <- list(par = rho_init, fn = lower_box_fun, upper = t_upper,
                         lower = lower - f_mode, j = j - 1, control = control,
                         method = method, big_val = big_val)
        temp <- do.call(stats::optim, c(ru_args, add_args))
        l_box[j] <- temp$value
      } else {
        add_args <- list(par = rho_init, fn = lower_box_fun, j = j - 1,
                         control = control, method = method, big_val = Inf)
        temp <- do.call(stats::optim, c(ru_args, add_args))
        l_box[j] <- temp$value
        # Sometimes Nelder-Mead fails if the initial estimate is too good.
        # ... so avoid non-zero convergence indicator using L-BFGS-B instead.
        if (temp$convergence > 0) {
          if (sum(abs(rho_init - temp$par)) > .Machine$double.eps) {
            new_start <- shoof * rho_init + (1 - shoof) * temp$par
          } else {
            new_start <- temp$par * (1 - shoof)
          }
          add_args <- list(par = new_start, fn = lower_box_fun, j = j - 1,
                           control = control, method = "L-BFGS-B",
                           big_val = big_val, upper = t_upper,
                           lower = lower - f_mode)
          temp <- do.call(stats::optim, c(ru_args, add_args))
          l_box[j] <- temp$value
        }
        # Check using nlminb() if optim's iteration limit is reached.
        if (temp$convergence == 1) {
          if (sum(abs(rho_init - temp$par)) > .Machine$double.eps) {
            new_start <- shoof * rho_init + (1 - shoof) * temp$par
          } else {
            new_start <- temp$par * (1 - shoof)
          }
          add_args <- list(start = new_start, objective = lower_box_fun,
                           upper = t_upper, lower = lower - f_mode, j = j - 1,
                           big_val = Inf)
          temp <- do.call(stats::nlminb, c(ru_args, add_args))
          l_box[j] <- temp$objective
        }
      }
    }
    vals[j+1, ] <- temp$par
    conv[j+1] <- temp$convergence
    #
    # Find biplus(r) --------------
    #
    rho_init <- zeros
    rho_init[j] <- ep
    t_lower <- lower - f_mode
    t_lower[j] <- 0
    if (algor == "nlminb") {
      add_args <- list(start = rho_init, objective = upper_box_fun,
                       lower = t_lower, upper = upper - f_mode, j = j - 1,
                       control = control, big_val = Inf)
      temp <- do.call(stats::nlminb, c(ru_args, add_args))
      u_box[j] <- -temp$objective
      # Sometimes nlminb isn't sure that it has found the minimum when in fact
      # it has.  Try to check this, and avoid a non-zero convergence indicator
      # by using optim with method="L-BFGS-B", starting from nlminb's solution.
      if (temp$convergence > 0) {
        if (sum(abs(rho_init - temp$par)) > .Machine$double.eps) {
          new_start <- shoof * rho_init + (1 - shoof) * temp$par
        } else {
          new_start <- temp$par * (1 - shoof)
        }
        add_args <- list(par = new_start, fn = upper_box_fun, j = j - 1,
                         method = "L-BFGS-B", big_val = big_val,
                         lower = t_lower, upper = upper - f_mode)
        temp <- do.call(stats::optim, c(ru_args, add_args))
        u_box[j] <- -temp$value
      }
    }
    if (algor == "optim") {
      # L-BFGS-B and Brent don't like Inf or NA
      if (method == "L-BFGS-B" | method == "Brent") {
        add_args <- list(par = rho_init, fn = upper_box_fun,
                         lower = t_lower, upper = upper - f_mode, j = j - 1,
                         control = control, method = method, big_val = big_val)
        temp <- do.call(stats::optim, c(ru_args, add_args))
        u_box[j] <- -temp$value
      } else {
        add_args <- list(par = rho_init, fn = upper_box_fun, j = j - 1,
                         control = control, method = method, big_val = Inf)
        temp <- do.call(stats::optim, c(ru_args, add_args))
        u_box[j] <- -temp$value
        # Sometimes Nelder-Mead fails if the initial estimate is too good.
        # ... so avoid non-zero convergence indicator using L-BFGS-B instead.
        if (temp$convergence > 0) {
          if (sum(abs(rho_init - temp$par)) > .Machine$double.eps) {
            new_start <- shoof * rho_init + (1 - shoof) * temp$par
          } else {
            new_start <- temp$par * (1 - shoof)
          }
          add_args <- list(par = new_start, fn = upper_box_fun, j = j - 1,
                           control = control, method = "L-BFGS-B",
                           big_val = big_val, lower = t_lower,
                           upper = upper - f_mode)
          temp <- do.call(stats::optim, c(ru_args, add_args))
          u_box[j] <- -temp$value
        }
        # Check using nlminb() if optim's iteration limit is reached.
        if (temp$convergence == 1) {
          if (sum(abs(rho_init - temp$par)) > .Machine$double.eps) {
            new_start <- shoof * rho_init + (1 - shoof) * temp$par
          } else {
            new_start <- temp$par * (1 - shoof)
          }
          add_args <- list(start = new_start, objective = upper_box_fun,
                           lower = t_lower, upper = upper - f_mode, j = j - 1,
                           big_val = Inf)
          temp <- do.call(stats::nlminb, c(ru_args, add_args))
          u_box[j] <- -temp$objective
        }
      }
    }
    vals[j+d+1, ] <- temp$par
    conv[j+d+1] <- temp$convergence
  }
  return(list(l_box = l_box, u_box = u_box, vals = vals, conv = conv))
}

# ================================= wecdf =====================================

wecdf <- function (x, weights = rep(1, length(x))) {
  # Sort x and reorder the weights
  w <- weights[order(x)]
  x <- sort(x)
  # The sum of the weights for each unique value in xs
  ws <- tapply(w, x, sum)
  vals <- unique(x)
  rval <- stats::approxfun(vals, cumsum(ws) / sum(ws), method = "constant",
                           yleft = 0, yright = 1, f = 0, ties = "ordered")
  class(rval) <- c("ecdf", "stepfun", class(rval))
  assign("nobs", length(x), envir = environment(rval))
  attr(rval, "call") <- sys.call()
  rval
}
