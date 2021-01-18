#' Auxiliary Function for Setting DMBC Model Priors
#' 
#' @description{
#' \code{dmbc_prior()} is an auxiliary function as user interface for
#'   \code{dmbc()} fitting. Typically only used when calling the \code{dmbc()}
#'   function. It is used to set prior hyperparameters.
#' 
#' \code{prior_dmbc()} is an alias for \code{dmbc_prior()}.
#' 
#' \code{check_prior()} is an auxiliary function that verifies the
#'   correctness of the prior hyperparameters provided before a DMBC is fitted
#'   with \code{\link{dmbc}()}.
#' 
#' \code{update_prior()} is an auxiliary function to modify a set of prior
#'   choices using a new value of \emph{p} and \emph{G}. It is intended for
#'   internal use mainly in the \code{\link{dmbc_ic}()} function.
#' }
#'
#' @param eta A named list containing the hyperparameters for the prior
#'   distribution of the \eqn{\eta_1,\ldots,\eta_G} parameters. It should
#'   contain two numeric vectors, namely \code{a} and \code{b}.
#' @param sigma2 A named list containing the hyperparameters for the prior
#'   distributions of the \eqn{\sigma^2_1,\ldots,\sigma^2_G} parameters. It
#'   should contain two numeric scalars, namely \code{a} and \code{b}.
#' @param lambda A list containing the hyperparameters for the prior
#'   distribution of the \eqn{\lambda_1,\ldots,\lambda_G} parameters. It should
#'   contain a single numeric vector.
#' @param prior A named list of prior hyperparameters.
#' @return A list with the prior hyperparameters as components.
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#' @seealso \code{\link{dmbc}()}
#' @keywords model based clustering
#' @examples
#' \dontrun{
#' data(simdiss, package = "dmbc")
#' # Shorter run than default.
#' sim.fit <- dmbc(simdiss,
#'   control = dmbc_control(burnin = 1000, nsim = 2000, thin = 1, verbose = TRUE),
#'   prior = dmbc_prior(sigma2 = list(a = 1, b = 4)))
#' }
#'
#' @export
dmbc_prior <- function(eta = list(a = rep(1.5, .dmbcEnv$current_G), b = rep(.5, .dmbcEnv$current_G)),
                       sigma2 = list(a = 1e-1, b = 1e-1),
                       lambda = rep(1, .dmbcEnv$current_G)){
  prior <- list()
  for (arg in names(formals(sys.function())))
    prior[[arg]] <- get(arg)
  prior
}

#' @rdname dmbc_prior
#' @export
prior_dmbc <- dmbc_prior


#' @rdname dmbc_prior
#' @export
check_prior <- function(prior) {
  prior_ok <- TRUE

  # check prior list
  if (!is.list(prior)) {
    prior_ok <- FALSE
    return(prior_ok)
  }

  # check eta prior
  if (!is.list(prior[["eta"]])) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (length(prior[["eta"]][["a"]]) != .dmbcEnv$current_G) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (length(prior[["eta"]][["b"]]) != .dmbcEnv$current_G) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (any(prior[["eta"]][["a"]] < 0)) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (any(prior[["eta"]][["b"]] < 0)) {
    prior_ok <- FALSE
    return(prior_ok)
  }

  # check sigma2 prior
  if (!is.list(prior[["sigma2"]])) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (length(prior[["sigma2"]]) != 2) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (any(prior[["sigma2"]] < 0)) {
    prior_ok <- FALSE
    return(prior_ok)
  }

  # check lambda prior
  if (length(prior[["lambda"]]) != .dmbcEnv$current_G) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (any(prior[["lambda"]] < 0)) {
    prior_ok <- FALSE
    return(prior_ok)
  }

  return(prior_ok)
}

#' @rdname dmbc_prior
#' @param p A length-one numeric vector indicating the number of dimensions of the
#'   latent space.
#' @param G A length-one numeric vector indicating the number of cluster to
#'   partition the \emph{S} subjects.
#' @export
update_prior <- function(prior, p, G) {
  out <- dmbc_prior(eta = list(a = rep(prior[["eta"]][["a"]][1], G), b = rep(prior[["eta"]][["b"]][1], G)),
                    sigma2 = prior[["sigma2"]],
                    lambda = rep(prior[["lambda"]][1], G))

  return(out)
}
