.make_initial <- function(initial, lower, upper, n_initial, is_missing) {
  if (is_missing) {
    n_initial_length <- max(length(lower), length(upper), length(n_initial))
    initial <- list()
    for (i in 1 : n_initial_length) {
      lower_i <- lower[1 + (i - 1) %% length(lower)]
      upper_i <- upper[1 + (i - 1) %% length(upper)]
      n_initial_i <- n_initial[1 + (i - 1) %% length(n_initial)]
      initial[[i]] <- (
        lower_i +
        (1 : n_initial_i) * (upper_i - lower_i) / (n_initial_i + 1)
      )
    }
  }
  if (!is.list(initial)) {
    initial <- list(initial)
  }
  initial
}

#' Adaptive Rejection Metropolis Sampling
#'
#' This function performs Adaptive Rejection Metropolis Sampling to sample from
#' a target distribution specified by its (potentially unnormalised) log
#' density. The function constructs a rejection distribution based on piecewise
#' linear functions that envelop the log density of the target.
#' \cr\cr
#' If the target is log-concave, the \code{metropolis} parameter can be set to
#' \code{FALSE}, and an accept-reject sampling scheme is used which yields
#' independent samples.
#' \cr\cr
#' Otherwise, if \code{metropolis} is \code{TRUE}, a Metropolis-Hastings step is
#' used to construct a Markov chain with a stationary distribution matching the
#' target. It is possible in this case for the rejection distribution to be a
#' poor proposal, so users should be careful to check the output matches the
#' desired distribution.
#' \cr\cr
#' All arguments other than \code{n_samples}, \code{include_n_evaluations} and
#' \code{arguments} can be either vectors or lists as appropriate. If they are
#' vectors, they will be recycled in the same manner as, e.g., rnorm. The
#' entries of \code{arguments} may be vectors/lists and will also be recycled
#' (see examples).
#'
#' @param n_samples Number of samples to return.
#' @param log_pdf Potentially unnormalised log density of target distribution.
#' Can also be a list of functions.
#' @param lower Lower bound of the support of the target distribution.
#' @param upper Upper bound of the support of the target distribution.
#' @param initial Initial points with which to build the rejection distribution.
#' @param n_initial Number of points used to form \code{initial}; ignored if
#' \code{initial} provided.
#' @param convex Convexity adjustment.
#' @param max_points Maximum number of points to allow in the rejection
#' distribution.
#' @param metropolis Whether to use a Metropolis-Hastings step after rejection
#' sampling. Not necessary if the target distribution is log concave.
#' @param previous The previous value of the Markov chain to be used if
#' \code{metropolis = TRUE}.
#' @param include_n_evaluations Whether to return an object specifying the
#' number of function evaluations used.
#' @param arguments List of additional arguments to be passed to log_pdf
#' @return Vector or matrix of samples if \code{include_n_evaluations} is
#' \code{FALSE}, otherwise a list.
#' @seealso \url{http://www1.maths.leeds.ac.uk/~wally.gilks/adaptive.rejection/web_page/Welcome.html}
#' @references Gilks, W. R., Best, N. G. and Tan, K. K. C. (1995) Adaptive
#' rejection Metropolis sampling. Applied Statistics, 44, 455-472.
#' @examples
#' # The normal distribution, which is log concave, so metropolis can be FALSE
#' result <- arms(
#'   1000, dnorm, -1000, 1000, metropolis = FALSE,
#'   arguments = list(log = TRUE), include_n_evaluations = TRUE
#' )
#' print(result$n_evaluations)
#' hist(result$samples, freq = FALSE, br = 20)
#' curve(dnorm(x), min(result$samples), max(result$samples), col = 'red', add = TRUE)
#'
#' # Mixture of normals: 0.4 N(-1, 1) + 0.6 N(4, 1). Not log concave.
#' dnormmixture <- function(x) {
#'   parts <- log(c(0.4, 0.6)) + dnorm(x, mean = c(-1, 4), log = TRUE)
#'   log(sum(exp(parts - max(parts)))) + max(parts)
#' }
#' samples <- arms(1000, dnormmixture, -1000, 1000)
#' hist(samples, freq = FALSE)
#'
#' # List of log pdfs, demonstrating recycling of log_pdf argument
#' samples <- arms(
#'   10,
#'   list(
#'     function(x) -x ^ 2 / 2,
#'     function(x) -(x - 10) ^ 2 / 2
#'   ),
#'   -1000,
#'   1000
#' )
#' print(samples)
#'
#' # Another way to achieve the above, this time with recycling in arguments
#' samples <- arms(
#'   10, dnorm, -1000, 1000,
#'   arguments = list(
#'     mean = c(0, 10), sd = 1, log = TRUE
#'   )
#' )
#' print(samples)
#' @export
arms <- function(
  n_samples,
  log_pdf,
  lower,
  upper,
  previous = (upper + lower) / 2,
  initial = lower + (1 : n_initial) * (upper - lower) / (n_initial + 1),
  n_initial = 10,
  convex = 0,
  max_points = max(2 * n_initial + 1, 100),
  metropolis = TRUE,
  include_n_evaluations = FALSE,
  arguments = list()
) {
  initial <- .make_initial(
    initial, lower, upper, n_initial, missing(initial)
  )
  if (missing(previous)) {
    previous_length <- max(length(lower), length(upper))
    previous <- rep(0, previous_length)
    for (i in 1 : previous_length) {
      lower_i <- lower[1 + (i - 1) %% length(lower)]
      upper_i <- upper[1 + (i - 1) %% length(upper)]
      previous[i] <- (upper_i + lower_i) / 2
    }
  }
  if (!is.list(log_pdf)) {
    log_pdf <- list(log_pdf)
  }
  .arms(
    n_samples,
    log_pdf,
    lower,
    upper,
    initial,
    convex,
    max_points,
    metropolis,
    previous,
    arguments,
    include_n_evaluations
  )
}


#' Gibbs sampling using ARMS
#'
#' This function uses ARMS (see also \code{\link{arms}}) to sample from
#' a multivariate target distribution specified by its (potentially
#' unnormalised) log density using Gibbs sampling. The function updates each
#' argument to the log pdf in turn using ARMS, returning a matrix of samples.
#' \cr\cr
#' The arguments to this function have the same meaning as for
#' \code{\link{arms}}, except here they are recycled along the dimension of
#' \code{previous}, rather than from sample to sample.
#'
#' @param n_samples Number of samples to return.
#' @param log_pdf Potentially unnormalised log density of target distribution.
#' @param previous Starting value for the Gibbs sampler. Vectors of this length
#' are passed to \code{log_pdf} as the first argument.
#' @param lower Lower bound of the support of the target distribution.
#' @param upper Upper bound of the support of the target distribution.
#' @param initial Initial points with which to build the rejection distribution.
#' @param n_initial Number of points used to form \code{initial}; ignored if
#' \code{initial} provided.
#' @param convex Convexity adjustment.
#' @param max_points Maximum number of points to allow in the rejection
#' distribution.
#' @param metropolis Whether to use a Metropolis-Hastings step after rejection
#' sampling. Not necessary if the target distribution is log concave.
#' @param include_n_evaluations Whether to return an object specifying the
#' number of function evaluations used.
#' @param show_progress If \code{TRUE}, a progress bar is shown.
#' @return Matrix of samples if \code{include_n_evaluations} is \code{FALSE},
#' otherwise a list.
#' @seealso \url{http://www1.maths.leeds.ac.uk/~wally.gilks/adaptive.rejection/web_page/Welcome.html}
#' @references Gilks, W. R., Best, N. G. and Tan, K. K. C. (1995) Adaptive
#' rejection Metropolis sampling. Applied Statistics, 44, 455-472.
#' @examples
#' # The classic 8schools example from RStan
#' # https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
#' schools_data <- list(
#'   J = 8,
#'   y = c(28,  8, -3,  7, -1,  1, 18, 12),
#'   sigma = c(15, 10, 16, 11,  9, 11, 10, 18)
#' )
#'
#' log_pdf <- function(params, p) {
#'   mu <- params[1]
#'   tau <- params[2]
#'   eta <- params[3 : (3 + schools_data$J - 1)]
#
#'   output <- (
#'     sum(dnorm(eta, 0, 1, log = TRUE)) +
#'     sum(dnorm(
#'       schools_data$y,
#'       mu + tau * eta,
#'       schools_data$sigma,
#'       log = TRUE
#'     ))
#'   )
#'   return(output)
#' }
#'
#' x_start <- c(0, 0, rep(0, schools_data$J))
#' names(x_start) <- c(
#'   'mu',
#'   'tau',
#'   paste0('eta', 1 : schools_data$J)
#' )
#'
#' samples <- arms_gibbs(
#'   250,
#'   x_start,
#'   log_pdf,
#'   c(-1000, 0, rep(-1000, schools_data$J)),
#'   1000,
#'   metropolis = FALSE
#' )
#' print(colMeans(samples[51 : 250, ]))
#' @export
arms_gibbs <- function(
  n_samples,
  previous,
  log_pdf,
  lower,
  upper,
  initial = NULL,
  n_initial = 10,
  convex = 0,
  max_points = 100,
  metropolis = TRUE,
  include_n_evaluations = FALSE,
  show_progress = FALSE
) {
  initial <- .make_initial(
    initial, lower, upper, n_initial, missing(initial)
  )

  .arms_gibbs(
    n_samples,
    previous,
    log_pdf,
    lower, upper,
    initial,
    convex,
    max_points,
    metropolis,
    include_n_evaluations,
    show_progress
  )
}
