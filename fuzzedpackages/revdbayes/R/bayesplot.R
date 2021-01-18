# ========================= pp_check.evpost ===========================

#' Posterior predictive checks for an evpost object
#'
#' \code{pp_check} method for class "evpost".  This provides an interface
#' to the functions that perform posterior predictive checks in the
#' \strong{bayesplot} package.  See \link[bayesplot]{PPC-overview} for
#' details of these functions.
#'
#' @aliases pp_check
#'
#' @param object An object of class "evpost", a result of a call to
#'   \code{\link{rpost}} or \code{\link{rpost_rcpp}}.
#'   Currently \code{object$model = "gev"},
#'   \code{"gp"}, \code{"bingp"} and \code{"pp"} are supported.
#' @param ... Additional arguments passed on to bayesplot functions.
#' @param type A character vector.  The type of bayesplot plot required:
#' \itemize{
#'   \item "stat" for predictive test statistics
#'     (see \link[bayesplot]{PPC-test-statistics}),
#'   \item "overlaid" for comparison of observed data to predictive simulated
#'     datasets using overlaid density function or distribution functions
#'     (see \link[bayesplot]{PPC-distributions}),
#'   \item "multiple" for comparison of observed data to predictive simulated
#'     datasets using multiple summary plots
#'     (see \link[bayesplot]{PPC-distributions}),
#'   \item "intervals" for comparison of observed data to predictive simulated
#'     datasets using sample medians and a predictive interval,
#'     (see \link[bayesplot]{PPC-intervals}),
#'   \item "user" for direct access to the default bayesplot function
#'     \code{\link[bayesplot]{pp_check}}.  This requires the argument
#'     \code{fun} to be supplied
#'     (see \link[bayesplot]{pp_check}).
#' }
#' @param subtype A character scalar.  Specifies the form of the plot(s)
#'   produced.  Could be one of
#'   \code{"dens", "hist", "boxplot", "ribbon"} or \code{"intervals"}.
#'   If \code{subtype} is not supplied then the defaults are:
#'   \code{"ecdf"} if \code{type = overlaid},
#'   \code{"dens"} if \code{type = multiple},
#'   \code{"intervals"} if \code{type = intervals}.
#'   \code{subtype} is not relevant if \code{type = "stat"}.
#' @param stat See \link[bayesplot]{PPC-test-statistics}.
#' @param nrep If \code{type = "multiple"} the maximum number of
#'   summary plots of the predictive simulated datasets to include.
#'   If \code{nrep} is greater than \code{nrow(object$data_rep)} then
#'   \code{nrep} is set equal to \code{nrow(object$data_rep)}.
#' @param fun The plotting function to call.
#'   Only relevant if \code{type = "user"}.
#'   Can be any of the functions detailed at \link[bayesplot]{PPC-overview}.
#'   The "ppc_" prefix can optionally be dropped if fun is specified
#'   as a string.
#' @details For details of these functions see \link[bayesplot]{PPC-overview}.
#'   See also the vignette
#'   \href{https://CRAN.R-project.org/package=revdbayes}{Posterior Predictive Extreme Value Inference}
#'   and the \strong{bayesplot} vignette
#'   \href{https://CRAN.R-project.org/package=bayesplot}{Graphical posterior predictive checks}.
#'
#'   The general idea is to compare the observed data \code{object$data}
#'   with a matrix \code{object$data_rep} in which each row is a
#'   replication of the observed data simulated from the posterior predictive
#'   distribution.  For greater detail see Chapter 6 of Gelman et al. (2013).
#'
#'   The format of \code{object$data} depends on the model:
#'   \itemize{
#'     \item{\code{model = "gev"}.} A vector of block maxima.
#'     \item{\code{model = "gp"}.} Data that lie above the threshold,
#'       i.e. threshold exceedances.
#'     \item{\code{model = "bingp"} or \code{"pp"}} The input data are
#'       returned but any value lying below the threshold is set to
#'       \code{object$thresh}.
#'   }
#'   In all cases any missing values have been removed from the data.
#'
#'   If \code{model = "bingp"} or \code{"pp"} the rate of threshold exceedance
#'   is part of the inference.  Therefore, the number of values in
#'   \code{object$data_rep} that lie above the threshold varies between
#'   predictive replications, with values below the threshold being
#'   left-censored at the threshold.  This limits a little the posterior
#'   predictive checks that it is useful to perform.  In the examples below
#'   we have compared \code{object$data} and \code{object$data_rep} using
#'   only their sample maxima.
#' @return A ggplot object that can be further customized using the
#'   \strong{ggplot2} package.
#' @seealso \code{\link{rpost}} and \code{\link{rpost_rcpp}} for sampling
#'   from an extreme value posterior distribution.
#' @seealso \strong{bayesplot} functions \link[bayesplot]{PPC-overview},
#'   \link[bayesplot]{PPC-distributions},
#'   \link[bayesplot]{PPC-test-statistics},
#'   \link[bayesplot]{PPC-intervals},
#'   \link[bayesplot]{pp_check}.
#' @references Jonah Gabry (2016). bayesplot: Plotting for Bayesian
#' Models. R package version 1.1.0.
#' \url{https://CRAN.R-project.org/package=bayesplot}
#' @references Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B.,
#' Vehtari, A., and Rubin, D. B. (2013). \emph{Bayesian Data Analysis}.
#' Chapman & Hall/CRC Press, London, third edition. (Chapter 6)
#' @examples
#' \donttest{
#' # GEV model
#' data(portpirie)
#' mat <- diag(c(10000, 10000, 100))
#' pn <- set_prior(prior = "norm", model = "gev", mean = c(0,0,0), cov = mat)
#' gevp  <- rpost(1000, model = "gev", prior = pn, data = portpirie,
#'                nrep = 50)
#'
#' # Posterior predictive test statistics
#' pp_check(gevp)
#' pp_check(gevp, stat = "min")
#' pp_check(gevp, stat = c("min", "max"))
#' iqr <- function(y) diff(quantile(y, c(0.25, 0.75)))
#' pp_check(gevp, stat = "iqr")
#'
#' # Overlaid density and distributions functions
#' pp_check(gevp, type = "overlaid")
#' pp_check(gevp, type = "overlaid", subtype = "dens")
#'
#' # Multiple plots
#' pp_check(gevp, type = "multiple")
#' pp_check(gevp, type = "multiple", subtype = "hist")
#' pp_check(gevp, type = "multiple", subtype = "boxplot")
#'
#' # Intervals
#' pp_check(gevp, type = "intervals")
#' pp_check(gevp, type = "intervals", subtype = "ribbon")
#'
#' # User-supplied bayesplot function
#' # Equivalent to p_check(gevp, type = "overlaid")
#' pp_check(gevp, type = "user", fun = "dens_overlay")
#'
#' # GP model
#' data(gom)
#' u <- quantile(gom, probs = 0.65)
#' fp <- set_prior(prior = "flat", model = "gp", min_xi = -1)
#' gpg <- rpost(n = 1000, model = "gp", prior = fp, thresh = u,
#'              data = gom, nrep = 50)
#' pp_check(gpg)
#' pp_check(gpg, type = "overlaid")
#'
#' # bin-GP model
#' bp <- set_bin_prior(prior = "jeffreys")
#' bgpg <- rpost(n = 1000, model = "bingp", prior = fp, thresh = u,
#'               data = gom, bin_prior = bp, nrep = 50)
#' pp_check(bgpg, stat = "max")
#'
#' # PP model
#' data(rainfall)
#' rthresh <- 40
#' pf <- set_prior(prior = "flat", model = "gev", min_xi = -1)
#' ppr <- rpost(n = 1000, model = "pp", prior = pf, data = rainfall,
#'              thresh = rthresh, noy = 54, nrep = 50)
#' pp_check(ppr, stat = "max")
#' }
#' @export pp_check
#' @export
pp_check.evpost <- function(object, ...,
                            type = c("stat", "overlaid", "multiple",
                                     "intervals", "user"),
                            subtype = NULL, stat = "median", nrep = 8,
                            fun = NULL) {
  if (!inherits(object, "evpost")) {
    stop("use only with \"evpost\" objects")
  }
  if (is.null(object[["data_rep"]])) {
    stop("data_rep is NULL: call rpost() or rpost_rcpp() again supplying nrep")
  }
  type <- match.arg(type)
  if (is.null(subtype)) {
    subtype <- switch(type, overlaid = "ecdf", multiple = "dens",
                       intervals = "intervals", stat = NULL)
  }
  y <- object[["data"]]
  yrep <- object[["data_rep"]]
  if (type == "stat") {
    if (length(stat) == 1) {
      return(bayesplot::ppc_stat(y, yrep, stat = stat, ...))
    } else if (length(stat) == 2) {
      return(bayesplot::ppc_stat_2d(y, yrep, stat = stat, ...))
    } else {
      warning("stat cannot have length > 2, only first two elements used.")
      return(bayesplot::ppc_stat_2d(y, yrep, stat = stat[1:2], ...))
    }
  }
  if (type == "overlaid") {
    return(switch(subtype,
      dens = bayesplot::ppc_dens_overlay(y, yrep, ...),
      ecdf = bayesplot::ppc_ecdf_overlay(y, yrep, ...)
    ))
  }
  if (type == "multiple") {
    yrep <- yrep[1:min(nrep, nrow(yrep)), , drop = FALSE]
    return(switch(subtype,
      dens = bayesplot::ppc_dens(y, yrep, ...),
      hist = bayesplot::ppc_hist(y, yrep, ...),
      boxplot = bayesplot::ppc_boxplot(y, yrep, ...)
    ))
  }
  if (type == "intervals") {
    return(switch(subtype,
                  intervals = bayesplot::ppc_intervals(y, yrep, ...),
                  ribbon = bayesplot::ppc_ribbon(y, yrep, ...)
    ))
  }
  if (type == "user") {
    return(bayesplot::pp_check(y, yrep, fun = fun, ...))
  }
}
