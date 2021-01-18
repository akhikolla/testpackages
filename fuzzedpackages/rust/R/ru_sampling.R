# =========================== ru ===========================

#' Generalized ratio-of-uniforms sampling
#'
#' Uses the generalized ratio-of-uniforms method to simulate from a
#' distribution with log-density \eqn{log f} (up to an additive constant).
#' The density \eqn{f} must be bounded, perhaps after a transformation of variable.
#'
#' @param logf A function returning the log of the target density \eqn{f}.
#'   This function should return \code{-Inf} when the density is zero.
#' @param ... Further arguments to be passed to \code{logf} and related
#'   functions.
#' @param n A numeric scalar.  Number of simulated values required.
#' @param d A numeric scalar. Dimension of f.
#' @param init A numeric vector. Initial estimates of the mode of \code{logf}.
#'   If \code{trans="BC"} or \code{trans = "user"} this is \emph{after} Box-Cox
#'   transformation or user-defined transformation, but \emph{before} any
#'   rotation of axes.
#' @param trans A character scalar. "none" for no transformation, "BC" for
#'   Box-Cox transformation, "user" for a user-defined transformation.
#'   If \code{trans = "user"} then the transformation should be specified
#'   using \code{phi_to_theta} and \code{log_j} and \code{user_args} may be
#'   used to pass arguments to \code{phi_to_theta} and \code{log_j}.
#' @param phi_to_theta A function returning (the inverse) of the transformation
#'   from theta to phi used to ensure positivity of phi prior to Box-Cox
#'   transformation.  The argument is phi and the returned value is theta.
#'   If \code{phi_to_theta} is undefined at the input value then the
#'   function should return NA.
#' @param log_j A function returning the log of the Jacobian of the
#'  transformation from theta to phi, i.e. based on derivatives of phi with
#'  respect to theta. Takes theta as its argument.
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
#'   }
#'   This list may be created using \link{find_lambda_one_d} (for \code{d} = 1)
#'   or \link{find_lambda} (for any \code{d}).
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
#' See \code{vignette("rust-a-vignette", package = "rust")} for full details.
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
#'       problems in contouring fin \code{\link{plot.ru}} when
#'       \code{d = 2}.}
#'     \item{logf_rho}{A function. The target function actually used in the
#'       ratio-of-uniforms algorithm.}
#'     \item{sim_vals_rho}{An \code{n} by \code{d} matrix of values simulated
#'       from the function used in the ratio-of-uniforms algorithm.}
#'     \item{logf_args}{A list of further arguments to \code{logf}.}
#'     \item{f_mode}{The estimated mode of the target density f, after any
#'       Box-Cox transformation and/or user supplied transformation, but before
#'       mode relocation.}
#' @references Wakefield, J. C., Gelfand, A. E. and Smith, A. F. M. (1991)
#'  Efficient generation of random variates via the ratio-of-uniforms method.
#'  \emph{Statistics and Computing} (1991), \strong{1}, 129-133.
#'  \url{https://doi.org/10.1007/BF01889987}.
#' @examples
#' # Normal density ===================
#'
#' # One-dimensional standard normal ----------------
#' x <- ru(logf = function(x) -x ^ 2 / 2, d = 1, n = 1000, init = 0.1)
#'
#' # Two-dimensional standard normal ----------------
#' x <- ru(logf = function(x) -(x[1]^2 + x[2]^2) / 2, d = 2, n = 1000,
#'         init = c(0, 0))
#'
#' # Two-dimensional normal with positive association ----------------
#' rho <- 0.9
#' covmat <- matrix(c(1, rho, rho, 1), 2, 2)
#' log_dmvnorm <- function(x, mean = rep(0, d), sigma = diag(d)) {
#'   x <- matrix(x, ncol = length(x))
#'   d <- ncol(x)
#'   - 0.5 * (x - mean) %*% solve(sigma) %*% t(x - mean)
#' }
#'
#' # No rotation.
#' x <- ru(logf = log_dmvnorm, sigma = covmat, d = 2, n = 1000, init = c(0, 0),
#'         rotate = FALSE)
#'
#' # With rotation.
#' x <- ru(logf = log_dmvnorm, sigma = covmat, d = 2, n = 1000, init = c(0, 0))
#'
#' # three-dimensional normal with positive association ----------------
#' covmat <- matrix(rho, 3, 3) + diag(1 - rho, 3)
#'
#' # No rotation.  Slow !
#' x <- ru(logf = log_dmvnorm, sigma = covmat, d = 3, n = 1000,
#'         init = c(0, 0, 0), rotate = FALSE)
#'
#' # With rotation.
#' x <- ru(logf = log_dmvnorm, sigma = covmat, d = 3, n = 1000,
#'         init = c(0, 0, 0))
#'
#' # Log-normal density ===================
#'
#' # Sampling on original scale ----------------
#' x <- ru(logf = dlnorm, log = TRUE, d = 1, n = 1000, lower = 0, init = 1)
#'
#' # Box-Cox transform with lambda = 0 ----------------
#' lambda <- 0
#' x <- ru(logf = dlnorm, log = TRUE, d = 1, n = 1000, lower = 0, init = 0.1,
#'         trans = "BC", lambda = lambda)
#'
#' # Equivalently, we could use trans = "user" and supply the (inverse) Box-Cox
#' # transformation and the log-Jacobian by hand
#' x <- ru(logf = dlnorm, log = TRUE, d = 1, n = 1000, init = 0.1,
#'         trans = "user", phi_to_theta = function(x) exp(x),
#'         log_j = function(x) -log(x))
#'
#' # Gamma(alpha, 1) density ===================
#'
#' # Note: the gamma density in unbounded when its shape parameter is < 1.
#' # Therefore, we can only use trans="none" if the shape parameter is >= 1.
#'
#' # Sampling on original scale ----------------
#'
#' alpha <- 10
#' x <- ru(logf = dgamma, shape = alpha, log = TRUE, d = 1, n = 1000,
#'         lower = 0, init = alpha)
#'
#' alpha <- 1
#' x <- ru(logf = dgamma, shape = alpha, log = TRUE, d = 1, n = 1000,
#'         lower = 0, init = alpha)
#'
#' # Box-Cox transform with lambda = 1/3 works well for shape >= 1. -----------
#'
#' alpha <- 1
#' x <- ru(logf = dgamma, shape = alpha, log = TRUE, d = 1, n = 1000,
#'         trans = "BC", lambda = 1/3, init = alpha)
#' summary(x)
#'
#' # Equivalently, we could use trans = "user" and supply the (inverse) Box-Cox
#' # transformation and the log-Jacobian by hand
#'
#' # Note: when phi_to_theta is undefined at x this function returns NA
#' phi_to_theta  <- function(x, lambda) {
#'   ifelse(x * lambda + 1 > 0, (x * lambda + 1) ^ (1 / lambda), NA)
#' }
#' log_j <- function(x, lambda) (lambda - 1) * log(x)
#' lambda <- 1/3
#' x <- ru(logf = dgamma, shape = alpha, log = TRUE, d = 1, n = 1000,
#'         trans = "user", phi_to_theta = phi_to_theta, log_j = log_j,
#'         user_args = list(lambda = lambda), init = alpha)
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
#' # Mode relocation only ----------------
#' n <- 1000
#' x1 <- ru(logf = gpd_logpost, ss = ss, d = 2, n = n, init = init,
#'          lower = c(0, -Inf), rotate = FALSE)
#' plot(x1, xlab = "sigma", ylab = "xi")
#' # Parameter constraint line xi > -sigma/max(data)
#' # [This may not appear if the sample is far from the constraint.]
#' abline(a = 0, b = -1 / ss$xm)
#' summary(x1)
#'
#' # Rotation of axes plus mode relocation ----------------
#' x2 <- ru(logf = gpd_logpost, ss = ss, d = 2, n = n, init = init,
#'          lower = c(0, -Inf))
#' plot(x2, xlab = "sigma", ylab = "xi")
#' abline(a = 0, b = -1 / ss$xm)
#' summary(x2)
#'
#' # Cauchy ========================
#'
#' # The bounding box cannot be constructed if r < 1.  For r = 1 the
#' # bounding box parameters b1-(r) and b1+(r) are attained in the limits
#' # as x decreases/increases to infinity respectively.  This is fine in
#' # theory but using r > 1 avoids this problem and the largest probability
#' # of acceptance is obtained for r approximately equal to 1.26.
#'
#' res <- ru(logf = dcauchy, log = TRUE, init = 0, r = 1.26, n = 1000)
#'
#' # Half-Cauchy ===================
#'
#' log_halfcauchy <- function(x) {
#'   return(ifelse(x < 0, -Inf, dcauchy(x, log = TRUE)))
#' }
#'
#' # Like the Cauchy case the bounding box cannot be constructed if r < 1.
#' # We could use r > 1 but the mode is on the edge of the support of the
#' # density so as an alternative we use a log transformation.
#'
#' x <- ru(logf = log_halfcauchy, init = 0, trans = "BC", lambda = 0, n = 1000)
#' x$pa
#' plot(x, ru_scale = TRUE)
#'
#' # Example 4 from Wakefield et al. (1991) ===================
#'
#' # Bivariate normal x bivariate student-t
#' log_norm_t <- function(x, mean = rep(0, d), sigma1 = diag(d), sigma2 = diag(d)) {
#'   x <- matrix(x, ncol = length(x))
#'   d <- ncol(x)
#'   log_h1 <- -0.5 * (x - mean) %*% solve(sigma1) %*% t(x - mean)
#'   log_h2 <- -2 * log(1 + 0.5 * x %*% solve(sigma2) %*% t(x))
#'   return(log_h1 + log_h2)
#' }
#'
#' rho <- 0.9
#' covmat <- matrix(c(1, rho, rho, 1), 2, 2)
#' y <- c(0, 0)
#'
#' # Case in the top right corner of Table 3
#' x <- ru(logf = log_norm_t, mean = y, sigma1 = covmat, sigma2 = covmat,
#'   d = 2, n = 10000, init = y, rotate = FALSE)
#' x$pa
#'
#' # Rotation increases the probability of acceptance
#' x <- ru(logf = log_norm_t, mean = y, sigma1 = covmat, sigma2 = covmat,
#'   d = 2, n = 10000, init = y, rotate = TRUE)
#' x$pa
#' }
#' @seealso \code{\link{ru_rcpp}} for a version of \code{\link{ru}} that uses
#'   the Rcpp package to improve efficiency.
#' @seealso \code{\link{summary.ru}} for summaries of the simulated values
#'   and properties of the ratio-of-uniforms algorithm.
#' @seealso \code{\link{plot.ru}} for a diagnostic plot.
#' @seealso \code{\link{find_lambda_one_d}} to produce (somewhat) automatically
#'   a list for the argument \code{lambda} of \code{ru} for the
#'   \code{d} = 1 case.
#' @seealso \code{\link{find_lambda}} to produce (somewhat) automatically
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
ru <- function(logf, ..., n = 1, d = 1, init = NULL,
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
  if (trans == "none" & is.function(phi_to_theta)) {
    warning("phi_to_theta() not used when trans = ``none'': identity fn used")
  }
  if (!is.function(phi_to_theta) & !is.null(phi_to_theta)) {
    stop("phi_to_theta must be a function or NULL")
  }
  if (trans == "user" & is.null(phi_to_theta)) {
    stop("When trans = ``user'' phi_to_theta must be supplied")
  }
  if (is.function(phi_to_theta) & is.null(log_j)) {
    log_j <- function(x) 0
    warning("No Jacobian for phi_to_theta(): constant Jacobian has been used")
  }
  if (is.null(phi_to_theta)) {
    phi_to_theta <- identity
    log_j <- function(x) 0
  }
  # variable rotation matrix: set to matrix of ones initially
  rot_mat <- diag(1, d)
  if (rotate) {
    rho_to_psi <- function(rho) psi_mode + rho %*% rot_mat
  } else {
    rho_to_psi <- function(rho) psi_mode + rho
  }
  if (trans == "BC") {
    const <- lambda * gm ^ (lambda - 1)
    if (any(lambda == 0L)) {
      psi_to_phi <- function(psi) {
        ifelse(lambda == 0, exp(psi / gm), (psi * const + 1) ^ (1 / lambda))
      }
    } else {
      psi_to_phi <- function(psi) (psi * const + 1) ^ (1 / lambda)
    }
  }
  init_psi <- init
  #
  if (trans == "none") {
    logf_rho <- function(._rho, ...) {
      theta <- rho_to_psi(._rho)
      val <- logf(theta, ...) - hscale
      structure(val, theta = theta)
    }
  }
  if (trans == "BC") {
    logf_rho <- function(._rho, ...) {
      psi <- rho_to_psi(._rho)
      # When lambda is not equal to one (psi * const + 1) must be positive
      test <- (psi * const + 1)[which_lam]
      if (any(test <= 0)) return(-Inf)
      phi <- psi_to_phi(psi)
      theta <- phi_to_theta(phi)
      log_bc_jac <- sum((lambda - 1)[which_lam] * log(phi[which_lam]))
      val <- logf(theta, ...) - log_bc_jac - log_j(theta) - hscale
      structure(val, theta = theta)
    }
  }
  if (trans == "user") {
    logf_rho <- function(._rho, ...) {
      phi <- rho_to_psi(._rho)
      theta <- do.call(phi_to_theta, c(list(phi), user_args))
      if (any(!is.finite(theta))) return(-Inf)
      logj <- do.call(log_j, c(list(theta), user_args))
      val <- logf(theta, ...) - logj - hscale
      structure(val, theta = theta)
    }
  }
  neg_logf_rho <- function(._rho, ...) -logf_rho(._rho, ...)
  f_rho <- function(._rho, ...) exp(logf_rho(._rho, ...))
  #
  # Scale logf ---------------------------------
  #
  hscale <- logf_rho(init_psi, ...)
  if (is.infinite(hscale)) {
      stop("The target density is zero at initial parameter values")
  }
  #
  # Calculate a(r) ----------------------------------
  # Create list of arguments for find_a()
  for_find_a <- list(neg_logf_rho = neg_logf_rho, init_psi = init_psi, d = d,
                     r = r, lower = lower, upper = upper, algor = a_algor,
                     method = a_method, control = a_control, shoof = shoof,
                     ...)
  temp <- do.call("find_a", for_find_a)
  #
  # Check that logf is finite at 0
  #
  check_finite <- logf_rho(temp$par, ...)
  if (!is.finite(check_finite)) {
    stop(paste("The target log-density is not finite at its mode: mode = ",
               paste(temp$par, collapse = ","), ",
               function value = ", check_finite, ".", sep=""))
  }
  #
  # Scale logf to have a maximum at 0, i.e. a=1 ------------
  #
  hscale <- check_finite + hscale
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
  psi_mode <- f_mode
  #
  # If rotate = TRUE then don't impose any (finite) bounds
  if (rotate) {
    lower <- rep(-Inf, d)
    upper <- rep(Inf, d)
  }
  #
  # Calculate biminus(r) and biplus(r), i = 1, ...d -----------
  # Create list of arguments for find_bs()
  for_find_bs <- list(f_rho = f_rho, d = d, r = r, lower = lower,
                      upper = upper, f_mode = f_mode, ep = ep, vals = vals,
                      conv = conv, algor = b_algor, method = b_method,
                      control = b_control, shoof = shoof, ...)
  temp <- do.call("find_bs", for_find_bs)
  vals <- temp$vals
  conv <- temp$conv
  l_box <- temp$l_box
  u_box <- temp$u_box
  #
  # Perform ratio-of-uniforms rejection samping ---------------
  #
  n_try <- n_acc <- 0L     # initialize number of tries and accepted values
  res <- list()            # list to return posterior sample and other stuff
  res$sim_vals <- matrix(NA, ncol = d, nrow = n)
  res$sim_vals_rho <- matrix(NA, ncol = d, nrow = n)
  colnames(res$sim_vals) <- var_names
  colnames(res$sim_vals_rho) <- paste("rho[", 1:d, "]", sep="")
  #
  d_box <- u_box - l_box
  d_r <- d * r + 1
  #----------------------------------# start of while loop  while (n_acc < n) {
  while (n_acc < n) {
    u <- runif(1, 0, a_box)
    vs <- d_box * runif(d) + l_box
    rho <- vs / u ^ r
    rhs <- logf_rho(rho, ...)
    n_try <- n_try + 1L
    if (d_r * log(u) < rhs) {
      n_acc <- n_acc + 1L
      res$sim_vals[n_acc, ] <- attr(rhs, "theta")
      res$sim_vals_rho[n_acc, ] <- rho
    }
  }
  #-----------------------------------------# end of while (n_acc < n_sim) loop
  box <- c(a_box, l_box, u_box)
  res$box <- cbind(box, vals, conv)
  bs <- paste(paste("b", 1:d, sep=""),rep(c("minus", "plus"), each=d), sep="")
  rownames(res$box) <- c("a", bs)
  res$pa <- n_acc / n_try
  if (any(conv != 0)) {
    warning("One or more convergence indicators are non-zero.",
            immediate.=TRUE)
    print(res$box)
  }
  res$d <- d
  # Return logf - hscale to avoid over/under-flow in plot.ru() when d = 2.
  logf_scaled <- function(theta, ...) {
    return(logf(theta, ...) - hscale)
  }
  res$logf <- logf_scaled
  res$logf_args <- list(...)
  res$logf_rho <- logf_rho
  res$f_mode <- f_mode
  res$call <- Call
  class(res) <- "ru"
  return(res)
}
