#' Auxiliary Function for Controlling DMBC Model Fitting
#' 
#' @description{
#' \code{dmbc_control()} is an auxiliary function as user interface for
#'   \code{dmbc()} fitting. Typically only used when calling the \code{dmbc()}
#'   function. It is used to set parameters that affect the sampling but do
#'   not affect the posterior distribution.
#' 
#' \code{control_dmbc()} is an alias for \code{dmbc_control()}.
#' 
#' \code{check_control()} is an auxiliary function that verifies the
#'   correctness of the controls provided before a DMBC is fitted with
#'   \code{\link{dmbc}()}.
#' }
#'
#' @param nsim A length-one numeric vector for the number of draws to be taken
#'   from the posterior distribution.
#' @param burnin A length-one numeric vector for the number of initial MCMC
#'   iterations (usually to be discarded).
#' @param thin A length-one numeric vector for the number of iterations between
#'   consecutive draws.
#' @param nchains A length-one numeric vector for the number of parallel chains to run.
#' @param threads A length-one numeric vector for the number of chains to run.
#'   If greater than 1, package \pkg{\link{parallel}} is used to take advantage of any
#'   multiprocessing or distributed computing capabilities that may be available.
#' @param seed An integer scalar. If supplied, provides the random number seed.
#' @param parallel A length-one character vector indicating the type of parallel
#'   operation to be used (if any). Possible values are \code{multicore}
#'   (which worksonly on Unix/mcOS), \code{snow} and \code{no} (i.e. serial
#'   instead of parallel computing).
#' @param z.prop A length-one numeric vector providing the standard deviation of the
#'   proposal distribution for the jump in the individual latent space
#'   position.
#' @param alpha.prop A length-one numeric vector providing the standard deviation
#'   of the proposal distribution for the jump in the individual random effect value.
#' @param random.start A length-one logical vector. If \code{TRUE} the starting
#'   values are drawn randomly, otherwise a user-defined starting partition must
#'   be provided through the \code{partition} argument.
#' @param partition A length-one numeric vector providing the user-defined
#'   starting partition.
#' @param method A length-one character vector that specifies the distance
#'   measure to use in determining the initial partition. Allowed values are
#'   those from the \code{\link{dist}()} function.
#' @param procrustes A length-one logical vector. If \code{TRUE} the simulated
#'   MCMC chains are post-processed through a Procrustes transformation.
#' @param relabel A length-one logical vector. If \code{TRUE} the simulated
#'   MCMC chains are relabelled to address the label-switching problem.
#' @param store.burnin A logical scalar. If \code{TRUE}, the samples from the
#'   burnin are also stored and returned.
#' @param verbose A logical scalar. If \code{TRUE}, causes information to be
#'   printed out about the progress of the fitting.
#' @param control A list of control options.
#' @return A named list with the control options as components.
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#' @seealso \code{\link{dmbc}()}
#' @keywords model based clustering
#' @examples
#' \dontrun{
#' data(simdiss, package = "dmbc")
#' # Shorter run than default.
#' sim.fit <- dmbc(simdiss,
#'   control = dmbc_control(burnin = 1000, nsim = 2000, thin = 5, verbose = TRUE))
#' }
#' 
#' @export
dmbc_control <- function(nsim = 5000,
                         burnin = 10000,
                         thin = 1,
                         nchains = 1,
                         threads = 1,
                         seed = NULL,
                         parallel = "no",
                         z.prop = 1.5,
                         alpha.prop = 0.75,
                         random.start = TRUE,
                         partition = NULL,
                         method = "manhattan",
                         procrustes = TRUE,
                         relabel = TRUE,
                         store.burnin = TRUE,
                         verbose = FALSE){
  control <- list()
  for (arg in names(formals(sys.function())))
    control[[arg]] <- get(arg)
  control
}

#' @rdname dmbc_control
#' @export
control_dmbc <- dmbc_control

#' @rdname dmbc_control
#' @export
check_control <- function(control) {
  all_method <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
  control_ok <- TRUE

  if (!is.list(control)) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (control[["nsim"]] < 1) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (control[["burnin"]] < 1) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (control[["thin"]] < 1) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (control[["nchains"]] < 1) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (control[["threads"]] < 1) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (!is.null(control[["seed"]])) {
    if (control[["seed"]] < 1) {
      control_ok <- FALSE
      return(control_ok)
    }
  }
  if (!(control[["parallel"]] %in% c("no", "snow", "multicore"))) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (control[["z.prop"]] < 0) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (control[["alpha.prop"]] < 0) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (!is.logical(control[["random.start"]])) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (!control[["random.start"]]) {
    if (is.null(control[["partition"]])) {
      if (!(control[["method"]] %in% all_method)) {
        control_ok <- FALSE
        return(control_ok)
      }
    } else {
      if (!is.numeric(control[["partition"]])) {
        control_ok <- FALSE
        return(control_ok)
      }
    }
  }
  if (!is.logical(control[["store.burnin"]])) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (!is.logical(control[["verbose"]])) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (!is.logical(control[["procrustes"]])) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (!is.logical(control[["relabel"]])) {
    control_ok <- FALSE
    return(control_ok)
  }

  return(control_ok)
}
