# =========================== plot.evpost ===========================

#' Plot diagnostics for an evpost object
#'
#' \code{plot} method for class "evpost".  For \code{d = 1} a histogram of the
#'   simulated values is plotted with a the density function superimposed.
#'   The density is normalized crudely using the trapezium rule.  For
#'   \code{d = 2} a scatter plot of the simulated values is produced with
#'   density contours superimposed.  For \code{d > 2} pairwise plots of the
#'   simulated values are produced.
#'   An interface is also provided to the functions in the \strong{bayesplot}
#'   package that produce plots of Markov chain Monte Carlo (MCMC)
#'   simulations.  See \link[bayesplot]{MCMC-overview} for details of these
#'   functions.
#'
#' @param x An object of class "evpost", a result of a call to
#'   \code{\link{rpost}} or \code{\link{rpost_rcpp}}.
#' @param y Not used.
#' @param ... Additional arguments passed on to \code{hist}, \code{lines},
#'   \code{contour}, \code{points} or functions from the \strong{bayesplot}
#'   package.
#' @param n A numeric scalar.  Only relevant if \code{x$d = 1} or
#'   \code{x$d = 2}. The meaning depends on the value of x$d.
#' \itemize{
#'   \item {For d = 1 : n + 1 is the number of abscissae in the trapezium
#'      method used to normalize the density.}
#'   \item {For d = 2 : an n by n regular grid is used to contour the density.}
#' }
#' @param prob Numeric vector. Only relevant for d = 2.  The contour lines are
#'   drawn such that the respective probabilities that the variable lies
#'   within the contour are approximately prob.
#' @param ru_scale A logical scalar.  Should we plot data and density on the
#'   scale used in the ratio-of-uniforms algorithm (TRUE) or on the original
#'   scale (FALSE)?
#' @param rows A numeric scalar.  When \code{d} > 2 this sets the number of
#'   rows of plots.  If the user doesn't provide this then it is set
#'   internally.
#' @param xlabs,ylabs Numeric vectors.  When \code{d} > 2 these set the labels
#'   on the x and y axes respectively.  If the user doesn't provide these then
#'   the column names of the simulated data matrix to be plotted are used.
#' @param points_par A list of arguments to pass to
#'   \code{\link[graphics]{points}} to control the appearance of points
#'   depicting the simulated values. Only relevant when \code{d = 2}.
#' @param pu_only Only produce a plot relating to the posterior distribution
#'   for the threshold exceedance probability \eqn{p}. Only relevant when
#'   \code{model == "bingp"} was used in the call to \code{rpost} or
#'   \code{rpost_rcpp}.
#' @param add_pu Before producing the plots add the threshold exceedance
#'   probability \eqn{p} to the parameters of the extreme value model. Only
#'   relevant when \code{model == "bingp"} was used in the call to
#'   \code{rpost} or \code{rpost_rcpp}.
#' @param use_bayesplot A logical scalar. If \code{TRUE} the bayesplot
#'   function indicated by \code{fun_name} is called.  In principle \emph{any}
#'   bayesplot function (that starts with \code{mcmc_}) can be called but
#'   this may not always be successful because, for example, some of the
#'   bayesplot functions work only with multiple MCMC simulations.
#' @param fun_name A character scalar.  The name of the bayesplot function,
#'   with the initial \code{mcmc_} part removed.  See
#'   \link[bayesplot]{MCMC-overview} and links therein for the names of these
#'   functions. Some examples are given below.
#' @details For details of the \strong{bayesplot} functions available when
#'   \code{use_bayesplot = TRUE} see \link[bayesplot]{MCMC-overview} and
#'   the \strong{bayesplot} vignette
#'   \href{https://CRAN.R-project.org/package=bayesplot}{Plotting MCMC draws}.
#' @return Nothing is returned unless \code{use_bayesplot = TRUE} when a
#'   ggplot object, which can be further customized using the
#'   \strong{ggplot2} package, is returned.
#' @seealso \code{\link{summary.evpost}} for summaries of the simulated values
#'   and properties of the ratio-of-uniforms algorithm.
#' @seealso \code{\link[bayesplot]{MCMC-overview}},
#'   \code{\link[bayesplot]{MCMC-intervals}},
#'   \code{\link[bayesplot]{MCMC-distributions}}.
#' @references Jonah Gabry (2016). bayesplot: Plotting for Bayesian
#' Models. R package version 1.1.0.
#' \url{https://CRAN.R-project.org/package=bayesplot}
#' @examples
#' ## GP posterior
#' data(gom)
#' u <- stats::quantile(gom, probs = 0.65)
#' fp <- set_prior(prior = "flat", model = "gp", min_xi = -1)
#' gpg <- rpost(n = 1000, model = "gp", prior = fp, thresh = u, data = gom)
#' plot(gpg)
#'
#' \donttest{
#' # Using the bayesplot package
#' plot(gpg, use_bayesplot = TRUE)
#' plot(gpg, use_bayesplot = TRUE, pars = "xi", prob = 0.95)
#' plot(gpg, use_bayesplot = TRUE, fun_name = "intervals", pars = "xi")
#' plot(gpg, use_bayesplot = TRUE, fun_name = "hist")
#' plot(gpg, use_bayesplot = TRUE, fun_name = "dens")
#' plot(gpg, use_bayesplot = TRUE, fun_name = "scatter")
#' }
#'
#' ## bin-GP posterior
#' data(gom)
#' u <- quantile(gom, probs = 0.65)
#' fp <- set_prior(prior = "flat", model = "gp", min_xi = -1)
#' bp <- set_bin_prior(prior = "jeffreys")
#' npy_gom <- length(gom)/105
#' bgpg <- rpost(n = 1000, model = "bingp", prior = fp, thresh = u,
#'               data = gom, bin_prior = bp, npy = npy_gom)
#' plot(bgpg)
#' plot(bgpg, pu_only = TRUE)
#' plot(bgpg, add_pu = TRUE)
#'
#' \donttest{
#' # Using the bayesplot package
#' dimnames(bgpg$bin_sim_vals)
#' plot(bgpg, use_bayesplot = TRUE)
#' plot(bgpg, use_bayesplot = TRUE, fun_name = "hist")
#' plot(bgpg, use_bayesplot = TRUE, pars = "p[u]")
#' }
#' @export
plot.evpost <- function(x, y, ..., n = ifelse(x$d == 1, 1001, 101),
                        prob = c(0.5, 0.1, 0.25, 0.75, 0.95, 0.99),
                        ru_scale = FALSE, rows = NULL, xlabs = NULL,
                        ylabs = NULL, points_par = list(col = 8),
                        pu_only = FALSE, add_pu = FALSE, use_bayesplot = FALSE,
                        fun_name = c("areas", "intervals", "dens", "hist",
                                     "scatter")) {
  fun_name <- match.arg(fun_name)
  if (!inherits(x, "evpost")) {
    stop("use only with \"evpost\" objects")
  }
  if (n < 1) {
    stop("n must be no smaller than 1")
  }
  #
  if (use_bayesplot) {
    fun_name <- paste("mcmc_", fun_name, sep = "")
    bfun <- utils::getFromNamespace(fun_name, "bayesplot")
    x <- create_sim_vals(x)
    if (fun_name %in% c("mcmc_areas", "mcmc_intervals")) {
      return(bfun(x, prob = prob[1], ...))
    } else {
      return(bfun(x, ...))
    }
  }
  #
  prob <- sort(prob)
  if (ru_scale) {
    plot_data <- x$sim_vals_rho
    plot_density <- x$logf_rho
    if (is.null(x$logf_rho_args)) {
      density_args <- x$logf_args
    } else {
      density_args <- x$logf_rho_args
    }
  } else {
    plot_data <- x$sim_vals
    plot_density <- x$logf
    density_args <- x$logf_args
  }
  #
  if (pu_only & add_pu) {
    stop("At most one of 'pu_only' and 'add_pu' can be TRUE")
  }
  if (pu_only & is.null(x$bin_sim_vals)) {
    warning("pu_only = TRUE is not relevant and has been ignored",
            immediate. = TRUE)
    pu_only <- FALSE
  }
  if (add_pu & is.null(x$bin_sim_vals)) {
    warning("add_pu = TRUE is not relevant and has been ignored",
            immediate. = TRUE)
    add_pu <- FALSE
  }
  if (pu_only) {
    plot_data <- x$bin_sim_vals
    plot_density <- x$bin_logf
    density_args <- x$bin_logf_args
    x$d <- 1
  }
  if (add_pu) {
    plot_data <- cbind(x$bin_sim_vals, plot_data)
    x$d <- x$d + 1
  }
  if (x$d == 1) {
    temp <- graphics::hist(plot_data, plot = FALSE)
    a <- temp$breaks[1]
    b <- temp$breaks[length(temp$breaks)]
    h <- (b-a)/n
    xx <- seq(a, b, by = h)
    density_fun <- function(z) {
      density_list <- c(list(z), density_args)
      exp(do.call(plot_density, density_list))
    }
    yy <- sapply(xx, density_fun)
    # Remove any infinite, missing, or undefined values
    cond <- is.finite(yy)
    yy <- yy[cond]
    xx <- xx[cond]
    n <- length(yy) - 1
    #
    if (pu_only) {
      area <- 1
    } else {
      area <- h * (yy[1] / 2 + sum(yy[2:n]) + yy[n + 1] / 2)
    }
    yy <- yy / area
    max_y <- max(temp$density, yy)
    temp <- list(...)
    my_hist <- function(x, ..., type, lty, lwd, pch, lend, ljoin, lmitre) {
      temp_hist <- list(...)
      if (is.null(temp_hist$main)) {
        graphics::hist(x, main = "", ...)
      } else {
        graphics::hist(x, ...)
      }
    }
    if (is.null(temp$xlab)) {
      my_hist(plot_data, prob = TRUE, ylim = c(0, max_y), xlab = "", ...)
      if (!is.null(colnames(plot_data))) {
        graphics::title(xlab = parse(text = colnames(plot_data)[1]))
      }
    } else {
      my_hist(plot_data, prob = TRUE, ylim = c(0, max_y), ...)
    }
    my_lines <- function(x, y, ..., breaks, freq, probability, include.lowest,
                         right, density, angle, border, plot, labels, nclass) {
      graphics::lines(x, y, ...)
    }
    my_lines(xx, yy, ...)
  }
  if (x$d == 2) {
    r <- apply(plot_data, 2, range)
    xx <- seq(r[1,1], r[2,1], len = n)
    yy <- seq(r[1,2], r[2,2], len = n)
    zz <- matrix(NA, ncol = length(xx), nrow = length(yy))
    for (i in seq_along(xx)) {
      for (j in seq_along(yy)) {
        for_logf <- c(list(c(xx[i], yy[j])), density_args)
        zz[i, j] <- exp(do.call(plot_density, for_logf))
      }
    }
    zz[zz == -Inf] <- NA
    dx <- diff(xx[1:2]); dy <- diff(yy[1:2]); sz <- sort(zz)
    c1 <- cumsum(sz) * dx * dy; c1 <- c1/max(c1)
    # Suppress a warning that comes from stats:::regularize.values
    con_levs <- suppressWarnings(sapply(prob, function(x)
      stats::approx(c1, sz, xout = 1 - x)$y))
    #
    graphics::contour(xx, yy, zz, levels = con_levs, add = FALSE, ann = FALSE,
      labels = prob * 100, ...)
    do.call(graphics::points, c(list(x = plot_data), points_par))
    graphics::contour(xx, yy, zz, levels = con_levs, add = TRUE, ann = TRUE,
      labels = prob * 100, ...)
    temp <- list(...)
    if (is.null(temp$xlab)) {
      if (!is.null(colnames(plot_data))) {
        graphics::title(xlab = parse(text = colnames(plot_data)[1]))
      }
    }
    if (is.null(temp$ylab)) {
      if (!is.null(colnames(plot_data))) {
        graphics::title(ylab = parse(text = colnames(plot_data)[2]))
      }
    }
  }
  if (x$d > 2) {
    if (is.null(rows)) {
      rows <- x$d -2
    }
    cols <- ceiling(choose(x$d, 2) / rows)
    if (is.null(xlabs)) {
      if (!is.null(colnames(plot_data))) {
        xlabs <- colnames(plot_data)
      } else {
        xlabs <- rep(NA, x$d)
      }
    }
    if (is.null(ylabs)) {
      if (!is.null(colnames(plot_data))) {
        ylabs <- colnames(plot_data)
      } else {
        ylabs <- rep(NA, x$d)
      }
    }
    oldpar <- graphics::par(mfrow = c(rows, cols))
    on.exit(graphics::par(oldpar))
    pairwise_plots <- function(x) {
      for (i in 1:(ncol(x)-1)) {
        for (j in (i+1):ncol(x)) {
          graphics::plot(x[, i], x[, j], xlab = "", ylab = "", ...)
          graphics::title(xlab = parse(text = xlabs[i]), ylab =
                            parse(text = ylabs[j]))
        }
      }
    }
    pairwise_plots(plot_data)
  }
}

# =========================== summary.evpost ===========================

#' Summarizing an evpost object
#'
#' \code{summary} method for class "evpost"
#'
#' @param object An object of class "evpost", a result of a call to
#'   \code{\link{rpost}} or \code{\link{rpost_rcpp}}.
#' @param add_pu Includes in the summary of the simulated values the threshold
#'   exceedance probability \eqn{p}. Only relevant when \code{model == "bingp"}
#'   was used in the call to \code{rpost} or \code{rpost_rcpp}.
#' @param ... Additional arguments passed on to \code{print}.
#' @return Prints
#' \itemize{
#'   \item {information about the ratio-of-uniforms bounding box, i.e.
#'     \code{object$box}}
#'   \item {an estimate of the probability of acceptance, i.e.
#'     \code{object$pa}}
#'   \item {a summary of the simulated values, via
#'     \code{summary(object$sim_vals)}}
#' }
#' @examples
#' # GP posterior
#' data(gom)
#' u <- stats::quantile(gom, probs = 0.65)
#' fp <- set_prior(prior = "flat", model = "gp", min_xi = -1)
#' gpg <- rpost_rcpp(n = 1000, model = "gp", prior = fp, thresh = u,
#'                   data = gom)
#' summary(gpg)
#' @seealso \code{\link[rust]{ru}} or \code{\link[rust]{ru_rcpp}} for
#'   descriptions of \code{object$sim_vals} and \code{object$box}.
#' @seealso \code{\link{plot.evpost}} for a diagnostic plot.
#' @export
summary.evpost <- function(object, add_pu = FALSE, ...) {
  if (!inherits(object, "evpost")) {
    stop("use only with \"evpost\" objects")
  }
  if (!add_pu) {
    sim_vals <- object$sim_vals
  } else {
    sim_vals <- cbind(object$bin_sim_vals, object$sim_vals)
  }
  object <- object[c("box", "pa")]
  object$sim_vals <- sim_vals
  class(object) <- "summary.evpost"
  return(object)
}

# ========================== print.summary.evpost =============================

#' Print method for objects of class "summary.evpost"
#'
#' \code{print} method for an object \code{object} of class "summary.evpost".
#'
#' @param x An object of class "summary.evpost", a result of a call to
#'   \code{\link{summary.evpost}}.
#' @param ... Additional arguments passed on to \code{print}.
#' @return Prints
#' \itemize{
#'   \item {information about the ratio-of-uniforms bounding box, i.e.
#'     \code{object$box}}
#'   \item {an estimate of the probability of acceptance, i.e.
#'     \code{object$pa}}
#'   \item {a summary of the simulated values, via
#'     \code{summary(object$sim_vals)}}
#' }
#' @examples
#' # GP posterior
#' data(gom)
#' u <- stats::quantile(gom, probs = 0.65)
#' fp <- set_prior(prior = "flat", model = "gp", min_xi = -1)
#' gpg <- rpost_rcpp(n = 1000, model = "gp", prior = fp, thresh = u,
#'                   data = gom)
#' summary(gpg)
#' @seealso \code{\link[rust]{ru}} or \code{\link[rust]{ru_rcpp}} for
#'   descriptions of \code{object$sim_vals} and \code{$box}.
#' @seealso \code{\link{plot.evpost}} for a diagnostic plot.
#' @export
print.summary.evpost <- function(x, ...) {
  if (!inherits(x, "summary.evpost")) {
    stop("use only with \"summary.evpost\" objects")
  }
  cat("ru bounding box: ", "\n")
  print(x$box, ...)
  cat("\n")
  cat("estimated probability of acceptance: ", "\n")
  print(x$pa, ...)
  cat("\n")
  cat("sample summary", "\n")
  print(summary(x$sim_vals, ...), ...)
  invisible(x)
}

# ========================= create_sim_vals ===========================

create_sim_vals <- function(object) {
  x <- object$sim_vals
  if (object$model == "bingp") {
    x <- cbind(x, object$bin_sim_vals)
  }
  return(x)
}
