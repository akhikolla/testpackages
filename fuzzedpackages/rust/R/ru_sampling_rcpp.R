# =========================== ru_rcpp ===========================

#' Generalized ratio-of-uniforms sampling using C++ via Rcpp
#'
#' Uses the generalized ratio-of-uniforms method to simulate from a
#' distribution with log-density \eqn{log f} (up to an additive constant).
#' \eqn{f} must be bounded, perhaps after a transformation of variable.
#' The file file `user_fns.cpp` that is sourced before running the examples
#' below is available at the rust Github page at
#' \url{https://github.com/paulnorthrop/rust/blob/master/src/user_fns.cpp}.
#'
#' @param logf An external pointer to a compiled C++ function returning the
#'   log of the target density \eqn{f}.
#'   This function should return \code{-Inf} when the density is zero.
#'   See the \code{vignette("rust-c-using-rcpp-vignette", package = "rust")},
#'   particularly the Section
#'   \strong{Providing a C++ function to \code{ru_rcpp}}, for details.
#' @param ... Further arguments to be passed to \code{logf} and related
#'   functions.
#' @param n A numeric scalar.  Number of simulated values required.
#' @param d A numeric scalar. Dimension of f.
#' @param init A numeric vector. Initial estimates of the mode of \code{logf}.
#'   If \code{trans = "BC"} or \code{trans = "user"} this is \emph{after}
#'   Box-Cox transformation or user-defined transformation, but \emph{before}
#'   any rotation of axes.
#' @param trans A character scalar. "none" for no transformation, "BC" for
#'   Box-Cox transformation, "user" for a user-defined transformation.
#'   If \code{trans = "user"} then the transformation should be specified
#'   using \code{phi_to_theta} and \code{log_j} and \code{user_args} may be
#'   used to pass arguments to \code{phi_to_theta} and \code{log_j}.
#' @param phi_to_theta An external pointer to a compiled C++ function returning
#'   (the inverse) of the transformation from theta to phi used to ensure
#'   positivity of phi prior to Box-Cox transformation.  The argument is
#'   phi and the returned value is theta.  If \code{phi_to_theta}
#'   is undefined at the input value then the function should return NA.
#' @param log_j An external pointer to a compiled C++ function returning the
#'  log of the Jacobian of the transformation from theta to phi, i.e. based on
#'  derivatives of phi with respect to theta. Takes theta as its argument.
#' @param user_args A list of numeric components. If \code{trans = ``user''}
#'   then \code{user_args} is a list providing arguments to the user-supplied
#'   functions \code{phi_to_theta} and \code{log_j}.
#' @param lambda Either
#' \itemize{
#'   \item {A numeric vector.  Box-Cox transformation parameters, or}
#'   \item {A list with components}
#'   \describe{
#'     \item{lambda}{A numeric vector.  Box-Cox parameters (required).}
#'     \item{gm}{A numeric vector.  Box-cox scaling parameters (optional).
#'       If supplied this overrides any \code{gm} supplied by the individual
#'       \code{gm} argument described below.}
#'     \item{init_psi}{A numeric vector.  Initial estimate of mode \emph{after}
#'       Box-Cox transformation (optional).}
#'     \item{sd_psi}{A numeric vector.  Estimates of the marginal standard
#'       deviations of the Box-Cox transformed variables (optional).}
#'     \item{phi_to_theta}{as above (optional).}
#'     \item{log_j}{as above (optional).}
#'     \item{user_args}{as above (optional).}
#'   }
#'   This list may be created using \code{\link{find_lambda_one_d_rcpp}}
#'   (for \code{d} = 1) or \code{\link{find_lambda_rcpp}} (for any \code{d}).
#' }
#' @param lambda_tol A numeric scalar.  Any values in lambda that are less
#'  than \code{lambda_tol} in magnitude are set to zero.
#' @param gm A numeric vector. Box-cox scaling parameters (optional). If
#'   \code{lambda$gm} is supplied in input list \code{lambda} then
#'   \code{lambda$gm} is used, not \code{gm}.
#' @param rotate A logical scalar. If TRUE (\code{d} > 1 only) use Choleski
#'   rotation.  If d = 1 and \code{rotate} = TRUE then rotate will be set to
#'   FALSE with a warning.
#' @param lower,upper Numeric vectors.  Lower/upper bounds on the arguments of
#'   the function \emph{after} any transformation from theta to phi implied by
#'   the inverse of \code{phi_to_theta}. If \code{rotate = FALSE} these
#'   are used in all of the optimizations used to construct the bounding box.
#'   If \code{rotate = TRUE} then they are use only in the first optimisation
#'   to maximise the target density.`
#'   If \code{trans = "BC"} components of \code{lower} that are negative are
#'   set to zero without warning and the bounds implied after the Box-Cox
#'   transformation are calculated inside \code{ru}.
#' @param r A numeric scalar.  Parameter of generalized ratio-of-uniforms.
#' @param ep A numeric scalar.  Controls initial estimates for optimizations
#'   to find the b-bounding box parameters.  The default (\code{ep}=0)
#'   corresponds to starting at the mode of \code{logf} small positive values
#'   of \code{ep} move the constrained variable slightly away from the mode in
#'   the correct direction.  If \code{ep} is negative its absolute value is
#'   used, with no warning given.
#' @param a_algor,b_algor Character scalars.  Either "nlminb" or "optim".
#'   Respective optimization algorithms used to find a(r) and (bi-(r), bi+(r)).
#' @param a_method,b_method Character scalars.  Respective methods used by
#'   \code{optim} to find a(r) and (bi-(r), bi+(r)).  Only used if \code{optim}
#'   is the chosen algorithm.  If \code{d} = 1 then \code{a_method} and
#'   \code{b_method} are set to \code{"Brent"} without warning.
#' @param a_control,b_control  Lists of control arguments to \code{optim} or
#'   \code{nlminb} to find a(r) and (bi-(r), bi+(r)) respectively.
#' @param var_names A character vector.  Names to give to the column(s) of
#'   the simulated values.
#' @param shoof A numeric scalar in [0, 1].  Sometimes a spurious
#'   non-zero convergence indicator is returned from
#'   \code{\link[stats]{optim}} or \code{\link[stats]{nlminb}}).
#'   In this event we try to check that a minimum has indeed been found using
#'   different algorithm.  \code{shoof} controls the starting value provided
#'   to this algorithm.
#'   If \code{shoof = 0} then we start from the current solution.
#'   If \code{shoof = 1} then we start from the initial estimate provided
#'   to the previous minimisation.  Otherwise, \code{shoof} interpolates
#'   between these two extremes, with a value close to zero giving a starting
#'   value that is close to the current solution.
#'   The exception to this is when the initial and current solutions are equal.
#'   Then we start from the current solution multiplied by \code{1 - shoof}.
#' @details If \code{trans = "none"} and \code{rotate = FALSE} then \code{ru}
#'   implements the (multivariate) generalized ratio of uniforms method
#'   described in Wakefield, Gelfand and Smith (1991) using a target
#'   density whose mode is relocated to the origin (`mode relocation') in the
#'   hope of increasing efficiency.
#'
#'   If \code{trans = "BC"} then marginal Box-Cox transformations of each of
#'   the \code{d} variables is performed, with parameters supplied in
#'   \code{lambda}.  The function \code{phi_to_theta} may be used, if
#'   necessary, to ensure positivity of the variables prior to Box-Cox
#'   transformation.
#'
#'   If \code{trans = "user"} then the function \code{phi_to_theta} enables
#'   the user to specify their own transformation.
#'
#'   In all cases the mode of the target function is relocated to the origin
#'   \emph{after} any user-supplied transformation and/or Box-Cox
#'   transformation.
#'
#'   If \code{d} is greater than one and \code{rotate = TRUE} then a rotation
#'   of the variable axes is performed \emph{after} mode relocation.  The
#'   rotation is based on the Choleski decomposition (see \link{chol}) of the
#'   estimated Hessian (computed using \code{\link[stats:optim]{optimHess}}
#'   of the negated
#'   log-density after any user-supplied transformation or Box-Cox
#'   transformation.  If any of the eigenvalues of the estimated Hessian are
#'   non-positive (which may indicate that the estimated mode of \code{logf}
#'   is close to a variable boundary) then \code{rotate} is set to \code{FALSE}
#'   with a warning.  A warning is also given if this happens when
#'   \code{d} = 1.
#'
#'   The default value of the tuning parameter \code{r} is 1/2, which is
#'   likely to be close to optimal in many cases, particularly if
#'   \code{trans = "BC"}.
#'
#' See \code{vignette("rust-b-using-rcpp-vignette", package = "rust")} and
#' \code{vignette("rust-a-vignette", package = "rust")} for full details.
#'
#' @return An object of class "ru" is a list containing the following
#'   components:
#'     \item{sim_vals}{An \code{n} by \code{d} matrix of simulated values.}
#'     \item{box}{A (2 * \code{d} + 1) by \code{d} + 2 matrix of
#'       ratio-of-uniforms bounding box information, with row names indicating
#'       the box parameter.  The columns contain
#'       \describe{
#'         \item{column 1}{values of box parameters.}
#'         \item{columns 2 to (2+\code{d}-1)}{values of variables at which
#'          these box parameters are obtained.}
#'         \item{column 2+\code{d}}{convergence indicators.}
#'       }
#'       Scaling of f within \code{ru} and relocation of the
#'       mode to the origin means that the first row of \code{box} will always
#'       be \code{c(1, rep(0, d))}.
#'     }
#'     \item{pa}{A numeric scalar.  An estimate of the probability of
#'       acceptance.}
#'     \item{d}{A numeric scalar.  The dimension of \code{logf}.}
#'     \item{logf}{A function. \code{logf} supplied by the user, but
#'       with f scaled by the maximum of the target density used in the
#'       ratio-of-uniforms method (i.e. \code{logf_rho}), to avoid numerical
#'       problems in contouring f in \code{\link{plot.ru}} when
#'       \code{d = 2}.}
#'     \item{logf_rho}{A function. The target function actually used in the
#'       ratio-of-uniforms algorithm.}
#'     \item{sim_vals_rho}{An \code{n} by \code{d} matrix of values simulated
#'       from the function used in the ratio-of-uniforms algorithm.}
#'     \item{logf_args}{A list of further arguments to \code{logf}.}
#'     \item{logf_rho_args}{A list of further arguments to \code{logf_rho}.
#'       Note: this component is returned by \code{ru_rcpp} but not
#'       by \code{ru}.}
#'     \item{f_mode}{The estimated mode of the target density f, after any
#'       Box-Cox transformation and/or user supplied transformation, but before
#'       mode relocation.}
#' @references Wakefield, J. C., Gelfand, A. E. and Smith, A. F. M. (1991)
#'  Efficient generation of random variates via the ratio-of-uniforms method.
#'  \emph{Statistics and Computing} (1991), \strong{1}, 129-133.
#'  \url{https://doi.org/10.1007/BF01889987}.
#' @references Eddelbuettel, D. and Francois, R. (2011). Rcpp: Seamless
#'  R and C++ Integration. \emph{Journal of Statistical Software},
#'  \strong{40}(8), 1-18.
#'  \url{https://www.jstatsoft.org/v40/i08/}.
#' @references Eddelbuettel, D. (2013). \emph{Seamless R and C++ Integration
#'  with Rcpp}, Springer, New York. ISBN 978-1-4614-6867-7.
#' @examples
#' n <- 1000
#'
#' # Normal density ===================
#'
#' # One-dimensional standard normal ----------------
#' ptr_N01 <- create_xptr("logdN01")
#' x <- ru_rcpp(logf = ptr_N01, d = 1, n = n, init = 0.1)
#'
#' # Two-dimensional standard normal ----------------
#' ptr_bvn <- create_xptr("logdnorm2")
#' rho <- 0
#' x <- ru_rcpp(logf = ptr_bvn, rho = rho, d = 2, n = n,
#'   init = c(0, 0))
#'
#' # Two-dimensional normal with positive association ===================
#' rho <- 0.9
#' # No rotation.
#' x <- ru_rcpp(logf = ptr_bvn, rho = rho, d = 2, n = n, init = c(0, 0),
#'              rotate = FALSE)
#'
#' # With rotation.
#' x <- ru_rcpp(logf = ptr_bvn, rho = rho, d = 2, n = n, init = c(0, 0))
#'
#' # Using general multivariate normal function.
#' ptr_mvn <- create_xptr("logdmvnorm")
#' covmat <- matrix(rho, 2, 2) + diag(1 - rho, 2)
#' x <- ru_rcpp(logf = ptr_mvn, sigma = covmat, d = 2, n = n, init = c(0, 0))
#'
#' # Three-dimensional normal with positive association ----------------
#' covmat <- matrix(rho, 3, 3) + diag(1 - rho, 3)
#'
#' # No rotation.
#' x <- ru_rcpp(logf = ptr_mvn, sigma = covmat, d = 3, n = n,
#'              init = c(0, 0, 0), rotate = FALSE)
#'
#' # With rotation.
#' x <- ru_rcpp(logf = ptr_mvn, sigma = covmat, d = 3, n = n,
#'              init = c(0, 0, 0))
#'
#' # Log-normal density ===================
#'
#' ptr_lnorm <- create_xptr("logdlnorm")
#' mu <- 0
#' sigma <- 1
#' # Sampling on original scale ----------------
#' x <- ru_rcpp(logf = ptr_lnorm, mu = mu, sigma = sigma, d = 1, n = n,
#'              lower = 0, init = exp(mu))
#'
#' # Box-Cox transform with lambda = 0 ----------------
#' lambda <- 0
#' x <- ru_rcpp(logf = ptr_lnorm, mu = mu, sigma = sigma, d = 1, n = n,
#'              lower = 0, init = exp(mu), trans = "BC", lambda = lambda)
#'
#' # Equivalently, we could use trans = "user" and supply the (inverse) Box-Cox
#' # transformation and the log-Jacobian by hand
#' ptr_phi_to_theta_lnorm <- create_phi_to_theta_xptr("exponential")
#' ptr_log_j_lnorm <- create_log_j_xptr("neglog")
#' x <- ru_rcpp(logf = ptr_lnorm, mu = mu, sigma = sigma, d = 1, n = n,
#'   init = 0.1, trans = "user", phi_to_theta = ptr_phi_to_theta_lnorm,
#'   log_j = ptr_log_j_lnorm)
#'
#' # Gamma (alpha, 1) density ===================
#'
#' # Note: the gamma density in unbounded when its shape parameter is < 1.
#' # Therefore, we can only use trans="none" if the shape parameter is >= 1.
#'
#' # Sampling on original scale ----------------
#'
#' ptr_gam <- create_xptr("logdgamma")
#' alpha <- 10
#' x <- ru_rcpp(logf = ptr_gam, alpha = alpha, d = 1, n = n,
#'   lower = 0, init = alpha)
#'
#' alpha <- 1
#' x <- ru_rcpp(logf = ptr_gam, alpha = alpha, d = 1, n = n,
#'   lower = 0, init = alpha)
#'
#' # Box-Cox transform with lambda = 1/3 works well for shape >= 1. -----------
#'
#' alpha <- 1
#' x <- ru_rcpp(logf = ptr_gam, alpha = alpha, d = 1, n = n,
#'   trans = "BC", lambda = 1/3, init = alpha)
#' summary(x)
#'
#' # Equivalently, we could use trans = "user" and supply the (inverse) Box-Cox
#' # transformation and the log-Jacobian by hand
#'
#' lambda <- 1/3
#' ptr_phi_to_theta_bc <- create_phi_to_theta_xptr("bc")
#' ptr_log_j_bc <- create_log_j_xptr("bc")
#' x <- ru_rcpp(logf = ptr_gam, alpha = alpha, d = 1, n = n,
#'   trans = "user", phi_to_theta = ptr_phi_to_theta_bc, log_j = ptr_log_j_bc,
#'   user_args = list(lambda = lambda), init = alpha)
#' summary(x)
#'
#' \donttest{
#' # Generalized Pareto posterior distribution ===================
#'
#' # Sample data from a GP(sigma, xi) distribution
#' gpd_data <- rgpd(m = 100, xi = -0.5, sigma = 1)
#' # Calculate summary statistics for use in the log-likelihood
#' ss <- gpd_sum_stats(gpd_data)
#' # Calculate an initial estimate
#' init <- c(mean(gpd_data), 0)
#'
#' n <- 1000
#' # Mode relocation only ----------------
#' ptr_gp <- create_xptr("loggp")
#' for_ru_rcpp <- c(list(logf = ptr_gp, init = init, d = 2, n = n,
#'                  lower = c(0, -Inf)), ss, rotate = FALSE)
#' x1 <- do.call(ru_rcpp, for_ru_rcpp)
#' plot(x1, xlab = "sigma", ylab = "xi")
#' # Parameter constraint line xi > -sigma/max(data)
#' # [This may not appear if the sample is far from the constraint.]
#' abline(a = 0, b = -1 / ss$xm)
#' summary(x1)
#'
#' # Rotation of axes plus mode relocation ----------------
#' for_ru_rcpp <- c(list(logf = ptr_gp, init = init, d = 2, n = n,
#'                  lower = c(0, -Inf)), ss)
#' x2 <- do.call(ru_rcpp, for_ru_rcpp)
#' plot(x2, xlab = "sigma", ylab = "xi")
#' abline(a = 0, b = -1 / ss$xm)
#' summary(x2)
#'
#' # Cauchy ========================
#'
#' ptr_c <- create_xptr("logcauchy")
#'
#' # The bounding box cannot be constructed if r < 1.  For r = 1 the
#' # bounding box parameters b1-(r) and b1+(r) are attained in the limits
#' # as x decreases/increases to infinity respectively.  This is fine in
#' # theory but using r > 1 avoids this problem and the largest probability
#' # of acceptance is obtained for r approximately equal to 1.26.
#'
#' res <- ru_rcpp(logf = ptr_c, log = TRUE, init = 0, r = 1.26, n = 1000)
#'
#' # Half-Cauchy ===================
#'
#' ptr_hc <- create_xptr("loghalfcauchy")
#'
#' # Like the Cauchy case the bounding box cannot be constructed if r < 1.
#' # We could use r > 1 but the mode is on the edge of the support of the
#' # density so as an alternative we use a log transformation.
#'
#' x <- ru_rcpp(logf = ptr_hc, init = 0, trans = "BC", lambda = 0, n = 1000)
#' x$pa
#' plot(x, ru_scale = TRUE)
#'
#' # Example 4 from Wakefield et al. (1991) ===================
#' # Bivariate normal x bivariate student-t
#'
#' ptr_normt <- create_xptr("lognormt")
#' rho <- 0.9
#' covmat <- matrix(c(1, rho, rho, 1), 2, 2)
#' y <- c(0, 0)
#'
#' # Case in the top right corner of Table 3
#' x <- ru_rcpp(logf = ptr_normt, mean = y, sigma1 = covmat, sigma2 = covmat,
#'   d = 2, n = 10000, init = y, rotate = FALSE)
#' x$pa
#'
#' # Rotation increases the probability of acceptance
#' x <- ru_rcpp(logf = ptr_normt, mean = y, sigma1 = covmat, sigma2 = covmat,
#'   d = 2, n = 10000, init = y, rotate = TRUE)
#' x$pa
#' }
#'
#' @seealso \code{\link{ru}} for a version of \code{\link{ru_rcpp}} that
#'   accepts R functions as arguments.
#' @seealso \code{\link{summary.ru}} for summaries of the simulated values
#'   and properties of the ratio-of-uniforms algorithm.
#' @seealso \code{\link{plot.ru}} for a diagnostic plot.
#' @seealso \code{\link{find_lambda_one_d_rcpp}} to produce (somewhat)
#'   automatically a list for the argument \code{lambda} of \code{ru} for the
#'   \code{d} = 1 case.
#' @seealso \code{\link{find_lambda_rcpp}} to produce (somewhat) automatically
#'   a list for the argument \code{lambda} of \code{ru} for any value of
#'   \code{d}.
#' @seealso \code{\link[stats]{optim}} for choices of the arguments
#'   \code{a_method}, \code{b_method}, \code{a_control} and \code{b_control}.
#' @seealso \code{\link[stats]{nlminb}} for choices of the arguments
#'   \code{a_control} and \code{b_control}.
#' @seealso \code{\link[stats:optim]{optimHess}} for Hessian estimation.
#' @seealso \code{\link[base]{chol}} for the Choleski decomposition.
#'
#' @export
ru_rcpp <- function(logf, ..., n = 1, d = 1, init = NULL,
               trans = c("none", "BC", "user"),  phi_to_theta = NULL,
               log_j = NULL, user_args = list(), lambda = rep(1L, d),
               lambda_tol = 1e-6, gm = NULL,
               rotate = ifelse(d == 1, FALSE, TRUE),
               lower = rep(-Inf, d),
               upper = rep(Inf, d), r = 1 / 2, ep = 0L,
               a_algor = if (d == 1) "nlminb" else "optim",
               b_algor = c("nlminb", "optim"),
               a_method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
                            "Brent"),
               b_method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
                            "Brent"),
               a_control = list(), b_control = list(), var_names = NULL,
               shoof = 0.2) {
  #
  Call <- match.call(expand.dots = TRUE)
  # Check that shoof is in [0, 1]
  if (shoof < 0 || shoof > 1) {
    stop("''shoof'' must be in [0, 1]")
  }
  # Check that logf is an external pointer.
  is_pointer <- (class(logf) == "externalptr")
  if (!is_pointer) {
    stop("logf must be an external pointer to a function")
  }
  # Extract list of parameters for logf
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
  # Check that the values of key arguments are suitable
  if (r < 0) {
    stop("r must be non-negative")
  }
  #
  a_algor <- match.arg(a_algor, c("optim", "nlminb"))
  a_method <- match.arg(a_method)
  b_algor <- match.arg(b_algor)
  b_method <- match.arg(b_method)
  if (any(upper <= lower)) {
    stop("upper must be greater than lower, componentwise.")
  }
  #
  trans <- match.arg(trans)
  # If Box-Cox scale parameter is not supplied (directly) set it to 1,
  # at least for the moment
  if (is.null(gm)) {
    gm <- rep(1, d)
  }
  # Set up Box-Cox transformation parameters (if necessary)
  if (trans == "BC") {
    lambda_type <- "numeric"
    # If lambda is a list then extract information from it
    if (is.list(lambda)) {
      lambda_type <- "list"
      if (is.null(lambda$lambda)) {
        stop("The list lambda must contain the object lambda$lambda")
      }
      if (!is.null(lambda$gm)) {
        gm <- lambda$gm
      }
      if (!is.null(lambda$init_psi)) {
        init <- lambda$init_psi
      }
      if (a_algor == "optim" & is.null(a_control$parscale)) {
        a_control <- c(a_control, list(parscale = lambda$sd_psi))
      }
      if (!is.null(lambda$phi_to_theta)) {
        phi_to_theta <- lambda$phi_to_theta
      }
      if (!is.null(lambda$log_j)) {
        log_j <- lambda$log_j
      }
      if (!is.null(lambda$user_args)) {
        user_args <- lambda$user_args
      }
      lambda <- lambda$lambda
    }
    # Set to zero values of lambda that are close to zero.
    lambda <- ifelse(abs(lambda) < lambda_tol, 0L, lambda)
    # Check that lambda is a vector of the correct length.
    if (!is.vector(lambda)) {
      stop("lambda must be a numeric vector")
    }
    if (!(length(lambda) %in% c(1, d))) {
      if (lambda_type == "numeric") {
        stop("lambda must be a numeric vector of length d")
      }
      if (lambda_type == "list") {
        stop("lambda$lambda must be a numeric vector of length d")
      }
    }
    if (length(lambda) == 1) {
      lambda <- rep(lambda, d)
    }
    # Adjust lower and upper for the value of lambda
    # Check that all components of upper are positive
    if (any(upper <= 0)) {
      stop("when trans = ``BC'' all elements of upper must be positive")
    }
    # If any components of lower or upper are negative then set them to zero.
    lower <- pmax(0, lower)
    lower <- ifelse(lambda == 0, gm * log(lower),
                    (lower^lambda - 1) / (lambda * gm ^ (lambda -1)))
    upper <- ifelse(lambda == 0, gm * log(upper),
                    (upper^lambda - 1) / (lambda * gm ^ (lambda -1)))
  }
  # Check that the optimization algorithm is appropriate given the bounds in
  # lower and upper.  If not then change it, with a warning.
  if (d == 1 & a_algor == "optim" & any(is.infinite(c(lower,upper)))) {
    a_algor = "nlminb"
    warning("For d = 1 finite lower and upper bounds must be supplied when
            using a_algor = `optim'.  a_algor has been changed to `nlminb'")
  }
  if (d == 1 & b_algor == "optim" & any(is.infinite(c(lower,upper)))) {
    b_algor = "nlminb"
    warning("For d = 1 finite lower and upper bounds must be supplied when
            using b_algor = `optim'.  b_algor has been changed to `nlminb'")
  }
  if (b_algor == "optim") {
    if (b_method == "BFGS" | b_method == "CG") {
      warning("Using optim with b_method==`BFGS' or `CG' can produce the error
              message `non-finite finite-difference value'.  If you really want
              to use BFGS or CG try setting ep to be positive but small, e.g.
              ep=0.001.", immediate. = TRUE, noBreaks. = TRUE)
    }
  }
  # If d = 1 then set a_method and b_method to "Brent", just in case optim is
  # being used.
  if (d == 1) {
    a_method <- "Brent"
    b_method <- "Brent"
  }
  # Determine which values in lambda are not equal to 1.
  if (d == 1) {
    which_lam <- 1L
  } else {
    which_lam <- which(lambda != 1L)
  }
  # If no initial estimates have been supplied then use a vector of ones.
  if (is.null(init)) {
    init <- rep(1, d)
    warning("No initial estimate of the mode given: a vector of ones has
            been used", noBreaks. = TRUE)
  }
  len_init <- length(init)
  if (len_init == 1 & d > 1) {
    init <- rep(init, length.out = d)
    warning("d > 1 but init has length 1: a d-vector of inits has been used")
  }
  if (len_init != d & len_init != 1) {
      stop("the length of init is incompatible with d")
  }
  # The option rotate = TRUE is not relevant when d = 1
  if (d == 1 & rotate) {
    rotate <- FALSE
    warning("rotation is not relevant when d=1: no rotation is used")
  }
  ep <- abs(ep)
  # Objects to store the values underlying the ru box and convergence
  # indicators.
  #   vals: will contain the values of the variables at which
  #         the ru box dimensions occur.
  #   conv: will contain the corresponding covergence indicators returned by
  #         the optimisation algorithms.
  vals <- matrix(NA, ncol = d, nrow = 2 * d + 1)
  colnames(vals) <- paste("vals", 1:d, sep="")
  conv <- rep(NA, 2 * d + 1)
  # Large value to return when parameters are out of bounds
  big_val <- Inf
  # f will be scaled by exp(hscale).  ... but no scaling yet_
  hscale <- 0
  # Mode of (transformed) target density: set to zero initially
  psi_mode <- rep(0, d)
  is_pointer <- (class(phi_to_theta) == "externalptr")
  if (trans == "none" & is_pointer) {
    warning("phi_to_theta() not used when trans = ``none'': identity fn used")
  }
  if (!is_pointer & !is.null(phi_to_theta)) {
    stop("phi_to_theta must be an external pointer to a function or NULL")
  }
  if (trans == "user" & is.null(phi_to_theta)) {
    stop("When trans = ``user'' phi_to_theta must be supplied")
  }
  is_pointer <- (class(log_j) == "externalptr")
  if (!is_pointer & !is.null(log_j)) {
    stop("log_j must be an external pointer to a function or NULL")
  }
  # variable rotation matrix: set to matrix of ones initially
  rot_mat <- diag(d)
  init_psi <- init
  #
  # Determine which case applies:
  #
  if (trans == "none") {
    logf_fun <- cpp_logf_rho
    a_obj_fun <- cpp_a_obj
    lower_box_fun <- cpp_lower_box
    upper_box_fun <- cpp_upper_box
    ru_fun <- ru_cpp
    logf_args <- list(psi_mode = rep(0, d), rot_mat = diag(d), hscale = 0,
                      logf = logf, pars = pars)
    ru_args <- list(d = d, r = r)
  } else if (trans == "BC" & is.null(phi_to_theta)) {
    logf_fun <- cpp_logf_rho_2
    a_obj_fun <- cpp_a_obj_2
    lower_box_fun <- cpp_lower_box_2
    upper_box_fun <- cpp_upper_box_2
    ru_fun <- ru_cpp_2
    con <- lambda * gm ^ (lambda - 1)
    tpars <- list(which_lam = which_lam - 1, lambda = lambda, gm = gm,
                  con = con)
    tfun <- create_trans_xptr("case_2")
    if (all(lambda != 0)) {
      ptpfun <- create_psi_to_phi_xptr("no_zero")
    } else {
      ptpfun <- create_psi_to_phi_xptr("has_zero")
    }
    phi_to_theta <- null_phi_to_theta_xptr("no_trans")
    log_j <- create_log_jac_xptr("log_none_jac")
    logf_args <- list(psi_mode = rep(0, d), rot_mat = diag(d), hscale = 0,
                      logf = logf, pars = pars, tpars = tpars, ptpfun = ptpfun,
                      phi_to_theta = phi_to_theta, log_j = log_j,
                      user_args = user_args)
    ru_args <- list(d = d, r = r, tfun = tfun)
  } else if (trans == "BC" & !is.null(phi_to_theta)) {
    logf_fun <- cpp_logf_rho_3
    a_obj_fun <- cpp_a_obj_2
    lower_box_fun <- cpp_lower_box_2
    upper_box_fun <- cpp_upper_box_2
    ru_fun <- ru_cpp_3
    con <- lambda * gm ^ (lambda - 1)
    tpars <- list(which_lam = which_lam - 1, lambda = lambda, gm = gm,
                  con = con)
    tfun <- create_trans_xptr("case_3")
    if (all(lambda != 0)) {
      ptpfun <- create_psi_to_phi_xptr("no_zero")
    } else {
      ptpfun <- create_psi_to_phi_xptr("has_zero")
    }
    if (is.null(log_j)) {
      log_j <- create_log_jac_xptr("case_3")
    }
    logf_args <- list(psi_mode = rep(0, d), rot_mat = diag(d), hscale = 0,
                      logf = logf, pars = pars, tpars = tpars, ptpfun = ptpfun,
                      phi_to_theta = phi_to_theta, log_j = log_j,
                      user_args = user_args)
    ru_args <- list(d = d, r = r, tfun = tfun)
  } else {
    logf_fun <- cpp_logf_rho_4
    a_obj_fun <- cpp_a_obj_2
    lower_box_fun <- cpp_lower_box_2
    upper_box_fun <- cpp_upper_box_2
    ru_fun <- ru_cpp_4
    tpars <- list()
    tfun <- create_trans_xptr("case_4")
    ptpfun <- create_psi_to_phi_xptr("no_trans")
    if (is.null(log_j)) {
      log_j <- create_log_jac_xptr("case_4")
    }
    logf_args <- list(psi_mode = rep(0, d), rot_mat = diag(d), hscale = 0,
                      logf = logf, pars = pars, tpars = tpars, ptpfun = ptpfun,
                      phi_to_theta = phi_to_theta, log_j = log_j,
                      user_args = user_args)
    ru_args <- list(d = d, r = r, tfun = tfun)
  }
  #
  # Scale logf ---------------------------------
  #
  logf_args$hscale <- do.call(logf_fun, c(list(rho = init_psi), logf_args))
  if (is.infinite(logf_args$hscale)) {
    stop("The target density is zero at initial parameter values")
  }
  #
  # Calculate a(r) ----------------------------------
  # Create list of arguments for find_a()
  ru_args <- c(ru_args, logf_args)
  for_find_a <- list(init_psi = init_psi, lower = lower, upper = upper,
                     algor = a_algor, method = a_method, control = a_control,
                     a_obj_fun = a_obj_fun, ru_args = ru_args,
                     shoof = shoof)
  temp <- do.call("cpp_find_a", for_find_a)
  #
  # Check that logf is finite at 0
  #
  check_finite <- do.call(logf_fun, c(list(rho = temp$par), logf_args))
  if (!is.finite(check_finite)) {
    stop(paste("The target log-density is not finite at its mode: mode = ",
               paste(temp$par, collapse = ","), ",
               function value = ", check_finite, ".", sep=""))
  }
  #
  # Scale logf to have a maximum at 0, i.e. a=1 ------------
  #
  ru_args$hscale <- check_finite + logf_args$hscale
  logf_args$hscale <- ru_args$hscale
  a_box <- 1
  f_mode <- temp$par
  vals[1, ] <- rep(0, d)
  conv[1] <- temp$convergence
  pos_def <- TRUE
  if (inherits(temp$hessian, "try-error")) {
    pos_def <- FALSE
  } else {
    hess_mat <- temp$hessian
    e_vals <- eigen(hess_mat, symmetric = TRUE, only.values = TRUE)$values
    if (any(e_vals < 1e-6)) {
      pos_def <- FALSE
    }
  }
  # We check the eigenvalues of the estimated Hessian, hess_mat, of -logf at
  # its minimum.  This has two purposes:
  #
  #   1. If rotate = TRUE the transformation uses the Cholesky decomposition
  #      of hess_mat so hess_mat must be symmetric and positive definite, with
  #      all eigenvalues positive.
  #   2. Even if rotate = FALSE it is worth checking the positive-definiteness
  #      of hess_mat.  Lack of positive-definiteness of hess_mat could result
  #      from the mode, f_mode, of logf could being at or near a variable
  #      boundary.  This may be fine if logf is finite at the boundary, i.e.
  #      we have found the maximum of logf, but it is possible that logf is
  #      unbounded.
  #
  # In both cases, i.e. regardless of rotate, a warning is given.
  if (!pos_def) {
    warning("The Hessian of the target log-density at its mode is not positive
            definite. This may not be a problem, but it may be that a mode
            at/near a parameter boundary has been found and/or that the target
            function is unbounded.", immediate. = TRUE, noBreaks. = TRUE)
    if (trans != "BC") {
      cat("  It might be worth using the option trans = ``BC''.", "\n")
    }
    if (rotate) {
      rotate <- FALSE
      warning("rotate has been changed to FALSE.", immediate. = TRUE)
    }
  }
  if (rotate) {
    rot_mat <- solve(t(chol(hess_mat)))
    # We standardize so that the determinant of rot_mat is 1, i.e. the
    # transformation rotates but doesn't scale.  To do this we divide
    # by det(rot_mat)^(1/d).  The determinant is the product of the
    # eigenvalues of rot_mat.  These eigenvalues are the eigenvalues of
    # hess_mat in e_vals, raised to the power -1/2.
    rot_mat <- rot_mat / exp(-mean(log(e_vals)) / 2)
  }
  # In the C++ function cpp_rho_to_psi() the vectors are column vectors
  # so we need to transpose rot_mat.
  ru_args$rot_mat <- t(rot_mat)
  ru_args$psi_mode <- f_mode
  logf_args$rot_mat <- t(rot_mat)
  logf_args$psi_mode <- f_mode
  #
  # If rotate = TRUE then don't use impose any (finite) bounds
  if (rotate) {
    lower <- rep(-Inf, d)
    upper <- rep(Inf, d)
  }
  #
  # Calculate biminus(r) and biplus(r), i = 1, ...d -----------
  # Create list of arguments for find_bs()
  for_find_bs <- list(lower = lower, upper = upper, ep = ep, vals = vals,
                      conv = conv,
                      algor = b_algor, method = b_method, control = b_control,
                      lower_box_fun = lower_box_fun,
                      upper_box_fun = upper_box_fun, ru_args = ru_args,
                      shoof = shoof)
  temp <- do.call("cpp_find_bs", for_find_bs)
  vals <- temp$vals
  conv <- temp$conv
  l_box <- temp$l_box
  u_box <- temp$u_box
  #
  # Perform ratio-of-uniforms rejection samping ---------------
  # Call C++ function to do this.
  #
  box_args <- list(n = n, a_box = a_box, l_box = l_box, u_box = u_box)
  ru_args <- c(box_args, ru_args)
  ru_args$tfun <- NULL
  res <- do.call(ru_fun, ru_args)
  res$pa <- n / res$ntry
  res$ntry <- NULL
  colnames(res$sim_vals) <- var_names
  colnames(res$sim_vals_rho) <- paste("rho[", 1:d, "]", sep="")
  box <- c(a_box, l_box, u_box)
  res$box <- cbind(box, vals, conv)
  bs <- paste(paste("b", 1:d, sep=""), rep(c("minus", "plus"), each=d), sep="")
  rownames(res$box) <- c("a", bs)
  if (any(conv != 0)) {
    warning("One or more convergence indicators are non-zero.",
            immediate.=TRUE)
    print(res$box)
  }
  res$d <- d
  # Add hscale to pars (and hence res$logf_args) and return cpp_logf_scaled
  # (which is cpp_logf - hscale) rather than cpp_logf.
  # This is to avoid over/under-flow in plot.ru() when d = 2.
  pars$hscale <- logf_args$hscale
  res$logf <- cpp_logf_scaled
  res$logf_args <- list(logf = logf, pars = pars)
  res$logf_rho <- logf_fun
  res$logf_rho_args <- logf_args
  res$f_mode <- f_mode
  res$call <- Call
  class(res) <- "ru"
  return(res)
}
