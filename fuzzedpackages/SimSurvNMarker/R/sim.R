.time_eps <- 1e-16

#' Faster Pointwise Function than ns
#'
#' @param knots sorted numeric vector with boundary and interior knots.
#' @param intercept logical for whether to include an intercept.
#' @param do_log logical for whether to evaluate the spline at \code{log(x)}
#'               or \code{x}.
#'
#' @description
#' Creates a function which can evaluate a natural cubic spline like
#' \code{\link{ns}}.
#'
#' The result may differ between different BLAS and LAPACK
#' implementations as the QR decomposition is not unique. However, the column space
#' of the returned matrix will always be the same regardless of the BLAS and LAPACK
#' implementation.
#'
#' @examples
#' # compare with splines
#' library(splines)
#' library(SimSurvNMarker)
#' xs <- seq(1, 5, length.out = 10L)
#' bks <- c(1, 5)
#' iks <- 2:4
#'
#' # we get the same
#' if(require(Matrix)){
#'   r1 <- unclass(ns(xs, knots = iks, Boundary.knots = bks, intercept = TRUE))
#'   r2 <- get_ns_spline(knots = sort(c(iks, bks)), intercept = TRUE,
#'                       do_log = FALSE)(xs)
#'
#'   cat("Rank is correct:      ", rankMatrix(cbind(r1, r2)) == NCOL(r1), "\n")
#'
#'   r1 <- unclass(ns(log(xs), knots = log(iks), Boundary.knots = log(bks),
#'                    intercept = TRUE))
#'   r2 <- get_ns_spline(knots = log(sort(c(iks, bks))), intercept = TRUE,
#'                       do_log = TRUE)(xs)
#'   cat("Rank is correct (log):", rankMatrix(cbind(r1, r2)) == NCOL(r1), "\n")
#' }
#'
#' # the latter is faster
#' system.time(
#'   replicate(100,
#'             ns(xs, knots = iks, Boundary.knots = bks, intercept = TRUE)))
#' system.time(
#'   replicate(100,
#'             get_ns_spline(knots = sort(c(iks, bks)), intercept = TRUE,
#'                           do_log = FALSE)(xs)))
#' func <- get_ns_spline(knots = sort(c(iks, bks)), intercept = TRUE,
#'                       do_log = FALSE)
#' system.time(replicate(100, func(xs)))
#'
#' @export
get_ns_spline <- function(knots, intercept = TRUE, do_log = TRUE)
  with(new.env(), {
    if(length(knots) > 1L){
      ptr_obj <- get_ns_ptr(knots          = knots[-c(1L, length(knots))],
                            boundary_knots = knots[ c(1L, length(knots))],
                            intercept = intercept)
      out <- if(do_log)
        function(x)
          ns_cpp(x = log(x), ns_ptr = ptr_obj)
      else
        function(x)
          ns_cpp(x =     x , ns_ptr = ptr_obj)

      attributes(out) <- list(
        knots = knots, intercept = intercept, do_log = do_log)
      return(out)

    }

    function(x)
      matrix(0, nrow = length(x), ncol = 0L)
  })

.surv_func_inner <- function(ti, omega, b_func, gl_dat){
  lb <- .time_eps
  if(ti < lb)
    return(0)
  ub <- ti


  nodes <- (ub - lb) / 2 * gl_dat$node + (ub + lb) / 2
  f <- exp(drop(b_func(nodes) %*% omega))
  -(ub - lb) / 2 * drop(gl_dat$weight %*% f)
}

#' Evaluates the Survival Function without a Marker
#'
#' @description
#' Evaluates the survival function at given points where the hazard is
#' given by
#'
#' \deqn{h(t) = \exp(\vec\omega^\top\vec b(t) + \delta).}
#'
#' @param ti numeric vector with time points.
#' @param omega numeric vector with coefficients for the baseline hazard.
#' @param b_func basis function for the baseline hazard like \code{\link{poly}}.
#' @param gl_dat Gauss–Legendre quadrature data.
#'               See \code{\link{get_gl_rule}}.
#' @param delta offset on the log hazard scale. Use \code{NULL} if
#'              there is no effect.
#'
#' @examples
#' # Example of a hazard function
#' b_func <- function(x)
#'   cbind(1, sin(2 * pi * x), x)
#' omega <- c(-3, 3, .25)
#' haz_fun <- function(x)
#'   exp(drop(b_func(x) %*% omega))
#'
#' plot(haz_fun, xlim = c(0, 10))
#'
#' # plot the hazard
#' library(SimSurvNMarker)
#' gl_dat <- get_gl_rule(60L)
#' plot(function(x) eval_surv_base_fun(ti = x, omega = omega,
#'                                     b_func = b_func, gl_dat = gl_dat),
#'      xlim = c(1e-4, 10), ylim = c(0, 1), bty = "l", xlab = "time",
#'      ylab = "Survival", yaxs = "i")
#'
#' # using to few nodes gives a wrong result in this case!
#' gl_dat <- get_gl_rule(15L)
#' plot(function(x) eval_surv_base_fun(ti = x, omega = omega,
#'                                     b_func = b_func, gl_dat = gl_dat),
#'      xlim = c(1e-4, 10), ylim = c(0, 1), bty = "l", xlab = "time",
#'      ylab = "Survival", yaxs = "i")
#' @export
eval_surv_base_fun <- function(
  ti, omega, b_func, gl_dat = get_gl_rule(30L), delta = NULL){
  cum_haz_fac <- vapply(
    ti, .surv_func_inner, FUN.VALUE = numeric(1L),
    omega = omega, b_func = b_func, gl_dat = gl_dat)

  if(length(delta) > 0)
    exp(exp(delta) * cum_haz_fac)
  else
    exp(cum_haz_fac)
}

#' Fast Evaluation of Time-varying Marker Mean Term
#'
#' @description
#' Evaluates the marker mean given by
#'
#' \deqn{\vec\mu(s, \vec u) = \vec o + B^\top\vec g(s) + U^\top\vec m(s).}
#'
#' @param ti numeric vector with time points.
#' @param B coefficient matrix for time-varying fixed effects.
#'          Use \code{NULL} if there is no effect.
#' @param g_func basis function for \code{B} like \code{\link{poly}}.
#' @param U random effects matrix for time-varying random effects.
#'          Use \code{NULL} if there is no effects.
#' @param m_func basis function for \code{U} like \code{\link{poly}}.
#' @param offset numeric vector with non-time-varying fixed effects.
#'
#' @examples
#' # compare R version with this function
#' library(SimSurvNMarker)
#' set.seed(1)
#' n <- 100L
#' n_y <- 3L
#'
#' ti <- seq(0, 1, length.out = n)
#' offset <- runif(n_y)
#' B <- matrix(runif(5L * n_y), nr = 5L)
#' g_func <- function(x)
#'   cbind(1, x, x^2, x^3, x^4)
#' U <- matrix(runif(3L * n_y), nr = 3L)
#' m_func <- function(x)
#'   cbind(1, x, x^2)
#'
#' r_version <- function(ti, B, g_func, U, m_func, offset){
#'   func <- function(ti)
#'     drop(crossprod(B, drop(g_func(ti))) + crossprod(U, drop(m_func(ti))))
#'
#'   vapply(ti, func, numeric(n_y)) + offset
#' }
#'
#' # check that we get the same
#' stopifnot(isTRUE(all.equal(
#'   c(r_version  (ti[1], B, g_func, U, m_func, offset)),
#'     eval_marker(ti[1], B, g_func, U, m_func, offset))))
#' stopifnot(isTRUE(all.equal(
#'   r_version  (ti, B, g_func, U, m_func, offset),
#'   eval_marker(ti, B, g_func, U, m_func, offset))))
#'
#' # check the computation time
#' system.time(replicate(100, r_version  (ti, B, g_func, U, m_func, offset)))
#' system.time(replicate(100, eval_marker(ti, B, g_func, U, m_func, offset)))
#'
#' @export
eval_marker <- function(ti, B, g_func, U, m_func, offset){
  out <- matrix(0., max(NCOL(B), NCOL(U), 1L), length(ti))
  if(length(B) > 0)
    eval_marker_cpp(B = B, m = g_func(ti), Sout = out)
  if(length(U) > 0)
    eval_marker_cpp(B = U, m = m_func(ti), Sout = out)
  if(length(offset) > 0)
    out <- out + offset

  drop(out)
}

#' Samples from a Multivariate Normal Distribution
#'
#' @description
#' Simulates from a multivariate normal distribution and returns a
#' matrix with appropriate dimensions.
#'
#' @param Psi_chol Cholesky decomposition of the covariance matrix.
#' @param n_y number of markers.
#'
#' @examples
#' library(SimSurvNMarker)
#' set.seed(1)
#' n_y <- 2L
#' K <- 3L * n_y
#' Psi <- drop(rWishart(1, K, diag(K)))
#' Psi_chol <- chol(Psi)
#'
#' # example
#' dim(draw_U(Psi_chol, n_y))
#' samples <- replicate(100, draw_U(Psi_chol, n_y))
#' samples <- t(apply(samples, 3, c))
#'
#' colMeans(samples) # ~ zeroes
#' cov(samples) # ~ Psi
#'
#' @importFrom stats rnorm
#' @export
draw_U <- function(Psi_chol, n_y){
  K <- NCOL(Psi_chol)
  structure(rnorm(K) %*%  Psi_chol, dim = c(K %/% n_y, n_y))
}

#' Simulate a Number of Observed Marker for an Individual
#'
#' @description
#' Simulates from
#'
#' \deqn{\vec U_i \sim N^{(K)}(\vec 0, \Psi)}
#' \deqn{\vec Y_{ij} \mid \vec U_i = \vec u_i \sim N^{(r)}(\vec \mu(s_{ij}, \vec u_i), \Sigma)}
#'
#' with
#'
#' \deqn{\vec\mu(s, \vec u) = \vec o + \left(I \otimes \vec g(s)^\top\right)vec(B) + \left(I \otimes \vec m(s)^\top\right) \vec u.}
#'
#' The number of observations and the observations times, \eqn{s_{ij}}s, are
#' determined from the passed generating functions.
#'
#' @param sigma_chol Cholesky decomposition of the noise's covariance matrix.
#' @param r_n_marker function to generate the number of observed markers.
#'                   Takes an integer for the individual's id.
#' @param r_obs_time function to generate the observations times given the
#'                   number of observed markers. Takes an
#'                   integer for the number of markers and an integer
#'                   for the individual's id.
#' @param id integer with id passed to \code{r_n_marker} and
#'           \code{r_obs_time}.
#' @inheritParams eval_marker
#'
#' @seealso
#' \code{\link{draw_U}}, \code{\link{eval_marker}}
#'
#' @examples
#' #####
#' # example with polynomial basis functions
#' g_func <- function(x){
#'   x <- x - 1
#'   cbind(x^3, x^2, x)
#' }
#' m_func <- function(x){
#'   x <- x - 1
#'   cbind(x^2, x, 1)
#' }
#'
#' # parameters
#' gamma <- matrix(c(.25, .5, 0, -.4, 0, .3), 3, 2)
#' Psi <- structure(c(0.18, 0.05, -0.05, 0.1, -0.02, 0.06, 0.05, 0.34, -0.25,
#'                    -0.06, -0.03, 0.29, -0.05, -0.25, 0.24, 0.04, 0.04,
#'                    -0.12, 0.1, -0.06, 0.04, 0.34, 0, -0.04, -0.02, -0.03,
#'                    0.04, 0, 0.1, -0.08, 0.06, 0.29, -0.12, -0.04, -0.08,
#'                    0.51), .Dim = c(6L, 6L))
#' B <- structure(c(-0.57, 0.17, -0.48, 0.58, 1, 0.86), .Dim = 3:2)
#' sig <- diag(c(.6, .3)^2)
#'
#' # generator functions
#' r_n_marker <- function(id){
#'   cat(sprintf("r_n_marker: passed id is %d\n", id))
#'   # the number of markers is Poisson distributed
#'   rpois(1, 10) + 1L
#' }
#' r_obs_time <- function(id, n_markes){
#'   cat(sprintf("r_obs_time: passed id is %d\n", id))
#'   # the observations times are uniform distributed
#'   sort(runif(n_markes, 0, 2))
#' }
#'
#' # simulate marker
#' set.seed(1)
#' U <- draw_U(chol(Psi), NCOL(B))
#' sim_marker(B = B, U = U, sigma_chol = chol(sig), r_n_marker = r_n_marker,
#'            r_obs_time = r_obs_time, m_func = m_func, g_func = g_func,
#'            offset = NULL, id = 1L)
#'
#' @importFrom stats rnorm
#' @export
sim_marker <- function(B, U, sigma_chol, r_n_marker, r_obs_time, m_func,
                       g_func, offset, id = 1L){
  n_markes <- r_n_marker(id)
  obs_time <- r_obs_time(id, n_markes)
  n_y <- NCOL(sigma_chol)

  noise <- matrix(rnorm(n_markes * n_y), ncol = n_y)  %*% sigma_chol
  mu <- eval_marker(
    ti = obs_time, B = B, g_func = g_func, m_func = m_func, offset = offset,
    U = U)
  y_obs <- if(is.vector(mu))
      mu  + noise
  else
    t(mu) + noise

  list(obs_time = obs_time, y_obs = y_obs)
}

.surv_func_joint_inner <- function(ti, B, U, omega, alpha, b_func, m_func,
                                   g_func, gl_dat){
  lb <- .time_eps
  if(ti < lb)
    return(0)
  ub <- ti

  nodes <- (ub - lb) / 2 * gl_dat$node + (ub + lb) / 2
  f <- exp(
    drop(b_func(nodes) %*% omega +
           drop(alpha %*% eval_marker(
             ti = nodes, B = B, U = U, m_func = m_func, g_func = g_func,
             offset = NULL))))
  -(ub - lb) / 2 * drop(gl_dat$weight %*% f)
}

#' Evaluates the Conditional Survival Function Given the Random Effects
#'
#' @description
#' Evaluates the conditional survival function given the random effects,
#' \eqn{\vec U}. The conditional hazard function is
#'
#' \deqn{h(t \mid \vec u) = \exp(\vec\omega^\top\vec b(t) + \delta +
#'   \vec\alpha^\top\vec o +
#'   \vec 1^\top(diag(\vec \alpha) \otimes \vec g(t)^\top)vec(B) +
#'   \vec 1^\top(diag(\vec \alpha) \otimes \vec m(t)^\top)\vec u).}
#'
#' @inheritParams eval_surv_base_fun
#' @inheritParams eval_marker
#' @param alpha numeric vector with association parameters.
#'
#' @seealso
#' \code{\link{sim_marker}}, \code{\link{draw_U}},
#' \code{\link{eval_surv_base_fun}}
#'
#' @examples
#' #####
#' # example with polynomial basis functions
#' b_func <- function(x){
#'   x <- x - 1
#'   cbind(x^3, x^2, x)
#' }
#' g_func <- function(x){
#'   x <- x - 1
#'   cbind(x^3, x^2, x)
#' }
#' m_func <- function(x){
#'   x <- x - 1
#'   cbind(x^2, x, 1)
#' }
#'
#' # parameters
#' omega <- c(1.4, -1.2, -2.1)
#' Psi <- structure(c(0.18, 0.05, -0.05, 0.1, -0.02, 0.06, 0.05, 0.34, -0.25,
#'                    -0.06, -0.03, 0.29, -0.05, -0.25, 0.24, 0.04, 0.04,
#'                    -0.12, 0.1, -0.06, 0.04, 0.34, 0, -0.04, -0.02, -0.03,
#'                    0.04, 0, 0.1, -0.08, 0.06, 0.29, -0.12, -0.04, -0.08,
#'                    0.51), .Dim = c(6L, 6L))
#' B <- structure(c(-0.57, 0.17, -0.48, 0.58, 1, 0.86), .Dim = 3:2)
#' alpha <- c(.5, .9)
#'
#' # simulate and draw survival curve
#' gl_dat <- get_gl_rule(30L)
#' set.seed(1)
#' U <- draw_U(chol(Psi), NCOL(B))
#' tis <- seq(0, 2, length.out = 100)
#' Survs <- surv_func_joint(ti = tis, B = B, U = U, omega = omega,
#'                          delta = NULL, alpha = alpha, b_func = b_func,
#'                          m_func = m_func, gl_dat = gl_dat, g_func = g_func,
#'                          offset = NULL)
#' par_old <- par(mar = c(5, 5, 1, 1))
#' plot(tis, Survs, xlab = "Time", ylab = "Survival", type = "l",
#'      ylim = c(0, 1), bty = "l", xaxs = "i", yaxs = "i")
#' par(par_old)
#'
#' @export
surv_func_joint <- function(ti, B, U, omega, delta, alpha, b_func, m_func,
                            gl_dat = get_gl_rule(30L), g_func, offset){
  cum_haz_fac <- vapply(
    ti, .surv_func_joint_inner, FUN.VALUE = numeric(1L),
    B = B, U = U, omega = omega, alpha = alpha, b_func = b_func,
    m_func = m_func, g_func = g_func, gl_dat = gl_dat)

  if(length(delta) > 0)
    cum_haz_fac <- cum_haz_fac * exp(delta)
  if(length(offset) > 0)
    cum_haz_fac <- cum_haz_fac * exp(drop(offset %*% alpha))
  exp(cum_haz_fac)
}

list_of_lists_to_data_frame <- function(dat)
  as.data.frame(do.call(mapply, c(list(FUN = function(...){
    if(is.matrix(...elt(1L)))
      do.call(rbind, list(...))
    else
      do.call(c, list(...))
  }, SIMPLIFY = FALSE), dat)))

#' Simulate Individuals from a Joint Survival and Marker Model
#'
#' @description
#' Simulates individuals from the following model
#'
#' \deqn{\vec U_i \sim N^{(K)}(\vec 0, \Psi)}
#' \deqn{\vec Y_{ij} \mid \vec U_i = \vec u_i \sim N^{(r)}(\vec \mu_i(s_{ij}, \vec u_i), \Sigma)}
#' \deqn{h(t \mid \vec u) = \exp(\vec\omega^\top\vec b(t) + \vec\delta^\top\vec z_i +
#'   \vec 1^\top(diag(\vec \alpha) \otimes \vec x_i^\top)vec(\Gamma) +
#'   \vec 1^\top(diag(\vec \alpha) \otimes \vec g(t)^\top)vec(B) +
#'   \vec 1^\top(diag(\vec \alpha) \otimes \vec m(t)^\top)\vec u)}
#'
#' with
#'
#' \deqn{\vec\mu_i(s, \vec u) = (I \otimes \vec x_i^\top)vec(\Gamma) + \left(I \otimes \vec g(s)^\top\right)vec(B) + \left(I \otimes \vec m(s)^\top\right) \vec u}
#'
#' where \eqn{h(t \mid \vec u)} is the conditional hazard function.
#'
#' @param delta coefficients for fixed effects in the log hazard.
#' @param n_obs integer with the number of individuals to draw.
#' @param Psi the random effects' covariance matrix.
#' @param sigma the noise's covariance matrix.
#' @param gamma coefficient matrix for the non-time-varying fixed effects.
#'              Use \code{NULL} if there is no effect.
#' @param r_z generator for the covariates in the log hazard. Takes an
#'            integer for the individual's id.
#' @param r_left_trunc generator for the left-truncation time.
#'                     Takes an integer for the individual's id.
#' @param r_right_cens generator for the right-censoring time. Takes an
#'                     integer for the individual's id.
#' @param r_x generator for the covariates in for the markers. Takes an
#'            integer for the individual's id.
#' @param y_max maximum survival time before administrative censoring.
#' @param tol convergence tolerance passed to \code{\link{uniroot}}.
#' @param use_fixed_latent logical for whether to include the
#'                         \eqn{\vec 1^\top(diag(\vec \alpha) \otimes \vec x_i^\top)vec(\Gamma)}
#'                         term in the log hazard. Useful if derivatives of
#'                         the latent mean should be used.
#' @param m_func_surv basis function for \code{U} like \code{\link{poly}}
#'                    in the log hazard. Can be different from
#'                    \code{m_func}. Useful if derivatives of the latent
#'                    mean should be used.
#' @param g_func_surv basis function for \code{B} like \code{\link{poly}}
#'                    in the log hazard. Can be different from
#'                    \code{g_func}. Useful if derivatives of the latent
#'                    mean should be used.
#' @inheritParams surv_func_joint
#' @inheritParams sim_marker
#'
#' @seealso
#' See the examples on Github at
#' \url{https://github.com/boennecd/SimSurvNMarker/tree/master/inst/test-data}
#' or this vignette
#' \code{vignette("SimSurvNMarker", package = "SimSurvNMarker")}.
#'
#' \code{\link{sim_marker}} and \code{\link{surv_func_joint}}
#'
#' @examples
#' #####
#' # example with polynomial basis functions
#' b_func <- function(x){
#'   x <- x - 1
#'   cbind(x^3, x^2, x)
#' }
#' g_func <- function(x){
#'   x <- x - 1
#'   cbind(x^3, x^2, x)
#' }
#' m_func <- function(x){
#'   x <- x - 1
#'   cbind(x^2, x, 1)
#' }
#'
#' # parameters
#' delta <- c(-.5, -.5, .5)
#' gamma <- matrix(c(.25, .5, 0, -.4, 0, .3), 3, 2)
#' omega <- c(1.4, -1.2, -2.1)
#' Psi <- structure(c(0.18, 0.05, -0.05, 0.1, -0.02, 0.06, 0.05, 0.34, -0.25,
#'                    -0.06, -0.03, 0.29, -0.05, -0.25, 0.24, 0.04, 0.04,
#'                    -0.12, 0.1, -0.06, 0.04, 0.34, 0, -0.04, -0.02, -0.03,
#'                    0.04, 0, 0.1, -0.08, 0.06, 0.29, -0.12, -0.04, -0.08,
#'                    0.51), .Dim = c(6L, 6L))
#' B <- structure(c(-0.57, 0.17, -0.48, 0.58, 1, 0.86), .Dim = 3:2)
#' sig <- diag(c(.6, .3)^2)
#' alpha <- c(.5, .9)
#'
#' # generator functions
#' r_n_marker <- function(id)
#'   # the number of markers is Poisson distributed
#'   rpois(1, 10) + 1L
#' r_obs_time <- function(id, n_markes)
#'   # the observations times are uniform distributed
#'   sort(runif(n_markes, 0, 2))
#' r_z <- function(id)
#'   # return a design matrix for a dummy setup
#'   cbind(1, (id %% 3) == 1, (id %% 3) == 2)
#' r_x <- r_z # same covariates for the fixed effects
#' r_left_trunc <- function(id)
#'   # no left-truncation
#'   0
#' r_right_cens <- function(id)
#'   # right-censoring time is exponentially distributed
#'   rexp(1, rate = .5)
#'
#' # simulate
#' gl_dat <- get_gl_rule(30L)
#' y_max <- 2
#' n_obs <- 100L
#' set.seed(1)
#' dat <- sim_joint_data_set(
#'   n_obs = n_obs, B = B, Psi = Psi, omega = omega, delta = delta,
#'   alpha = alpha, sigma = sig, gamma = gamma, b_func = b_func,
#'   m_func = m_func, g_func = g_func, r_z = r_z, r_left_trunc = r_left_trunc,
#'   r_right_cens = r_right_cens, r_n_marker = r_n_marker, r_x = r_x,
#'   r_obs_time = r_obs_time, y_max = y_max)
#'
#' # checks
#' stopifnot(
#'   NROW(dat$survival_data) == n_obs,
#'   NROW(dat$marker_data) >= n_obs,
#'   all(dat$survival_data$y <= y_max))
#'
#' @importFrom stats uniroot runif
#' @export
sim_joint_data_set <- function(
  n_obs, B, Psi, omega, delta, alpha, sigma, gamma, b_func, m_func, g_func,
  r_z, r_left_trunc, r_right_cens,
  r_n_marker, r_x, r_obs_time, y_max, use_fixed_latent = TRUE,
  m_func_surv = m_func, g_func_surv = g_func,
  gl_dat = get_gl_rule(30L),
  tol = .Machine$double.eps^(1/4)){
  Psi_chol <- chol(Psi)
  y_min <- 0
  sigma_chol <- chol(sigma)
  max_it <- 1000000L

  # checks
  n_y <- NCOL(sigma)
  d_m <- NROW(B)
  d_b <- length(omega)
  d_z <- length(delta)
  d_x <- length(gamma) / n_y
  K <- length(m_func(1L)) * n_y

  stopifnot(
    length(B) == 0L || is.matrix(B),
    length(B) == 0L || is.numeric(B),
    length(B) == 0L || NCOL(B) == n_y,
    is.matrix(Psi), is.numeric(Psi), NCOL(Psi) == K,
    is.matrix(sigma), is.numeric(sigma),
    is.numeric(alpha), length(alpha) == n_y,
    is.numeric(b_func(1)), length(b_func(1)) == d_b,
    is.numeric(m_func(1)),
    is.numeric(m_func_surv(1)),
    length(m_func_surv(1)) == length(m_func(1)),
    is.numeric(g_func(1)),
    is.numeric(g_func_surv(1)),
    length(g_func_surv(1)) == length(g_func(1)),
    is.numeric(gl_dat$node), is.numeric(gl_dat$weight),
    length(gl_dat$node) == length(gl_dat$weight),
    is.numeric(r_z(1L)), length(r_z(1L)) == d_z,
    is.numeric(r_x(1L)), length(r_x(1L)) == d_x,
    length(gamma) == 0L || is.matrix(gamma),
    length(gamma) == 0L || is.numeric(gamma),
    length(gamma) == 0L || NCOL(gamma) == n_y,
    is.numeric(r_left_trunc(1L)),
    is.numeric(r_right_cens(1L)),
    is.integer(r_n_marker(1L)),
    is.numeric(r_obs_time(1L, 1L)),
    is.logical(use_fixed_latent), length(use_fixed_latent) == 1L)

  out <- lapply(1:n_obs, function(i){
    z_delta <- if(d_z == 0L){
      z <- numeric()
      NULL
    }
    else {
      z <- r_z(i)
      drop(z %*% delta)
    }

    if(d_x == 0L){
      x <- numeric()
      mu_offest <- NULL
    } else {
      x <- r_x(i)
      mu_offest <- matrix(0., NCOL(gamma), 1L)
      eval_marker_cpp(B = gamma, m = x, Sout = mu_offest)
      mu_offest <- drop(mu_offest)
    }

    mu_offest_surv <- if(use_fixed_latent)
      mu_offest else NULL

    #####
    # start drawing
    U <- draw_U(Psi_chol, n_y = n_y)
    unif <- runif(1)
    fun <- function(x){
      surv <- surv_func_joint(
        ti = x, B = B, U = U, omega = omega,
        delta = z_delta, alpha = alpha, b_func = b_func,
        m_func = m_func_surv,
        gl_dat = gl_dat, g_func = g_func_surv, offset = mu_offest_surv)

      surv - unif
    }

    # simulate left-truncation time and right-censoring time
    left_trunc <- r_left_trunc(i)
    for(it in 1:max_it){
      if(it == max_it)
        stop("failed to find right-censoring time")
      right_cens <- r_right_cens(i)
      if(right_cens > left_trunc)
        break
    }

    y_min_use <- max(y_min, left_trunc)
    for(it in 1:max_it){
      f_lower <- fun(y_min_use)
      if(f_lower > 0)
        break
      if(it == max_it)
        stop("failed to find event time greater than left-truncation time")

      # the survival time is before the left-censoring time
      unif <- runif(1)
    }

    y_max_use <- min(y_max, right_cens * 1.00001)
    f_upper <- fun(y_max_use)

    y <- if(f_upper < 0){
      root <- uniroot(fun, interval = c(y_min_use, y_max_use),
                      f.lower = f_lower, f.upper = f_upper,
                      tol = tol)
      root$root
    } else
      y_max_use

    # simulate the marker
    for(it in 1:max_it){
      if(it == max_it)
        stop("failed to find markers between left-truncation time and the observed time")

      markers <- sim_marker(
        B = B, U = U, sigma_chol = sigma_chol, r_n_marker = r_n_marker,
        r_obs_time = r_obs_time, m_func = m_func, g_func = g_func,
        offset = mu_offest, id = i)
      keep <- markers$obs_time <= min(y, right_cens) &
        markers$obs_time >= left_trunc

      if(sum(keep) > 0L)
        break
    }

    markers <- cbind(obs_time = markers$obs_time, markers$y_obs)
    markers <- markers[keep, , drop = FALSE]
    colnames(markers)[-1] <- paste0("Y", 1:(NCOL(markers) - 1L))

    list(z = z, left_trunc = left_trunc, y = min(right_cens, y),
         event = y < min(y_max, y_max_use), U = U, markers = markers,
         x = x)
  })

  # form data frames for estimation
  ids <- seq_along(out)
  survival_dat <- lapply(ids, function(i){
      dat_i <- out[[i]][c("z", "left_trunc", "y", "event")]
      if(d_z > 0){
        dat_i$z <- matrix(
          dat_i$z, nrow = 1L,
          dimnames = list(NULL, paste0("Z", seq_along(dat_i$z))))
        names(dat_i)[1L] <- ""
      } else
        dat_i$z <- NULL
      dat_i$id <- i
      dat_i
    })
  survival_dat <- list_of_lists_to_data_frame(survival_dat)

  marker_dat <- lapply(ids, function(i){
      markers <- out[[i]][["markers"]]
      n_y <- NROW(markers)
      if(d_x > 0){
        x <- structure(rep(out[[i]]$x, each = n_y),
                       dimnames = list(NULL, paste0("X", seq_len(d_x))),
                       dim = c(n_y, d_x))
        markers <- cbind(markers, x)
      }

      out <- list(markers)
      out$id <- rep(i, NROW(markers))
      out
    })
  marker_dat <- list_of_lists_to_data_frame(marker_dat)

  list(survival_data = survival_dat,
       marker_data   = marker_dat,
       complete_data = out,
       params        = list(
         gamma = gamma, B = B, Psi = Psi, omega = omega, delta = delta,
         alpha = alpha, sigma = sigma, b_attr = attributes(b_func),
         m_attr = attributes(m_func), g_attr = attributes(g_func),
         m_surv_attr = attributes(m_func_surv),
         g_surv_attr = attributes(g_func_surv)))
}
