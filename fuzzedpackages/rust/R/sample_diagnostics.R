# =========================== plot.ru ===========================

#' Plot diagnostics for an ru object
#'
#' \code{plot} method for class "ru".  For \code{d = 1} a histogram of the
#'   simulated values is plotted with a the density function superimposed.
#'   The density is normalized crudely using the trapezium rule.  For
#'   \code{d = 2} a scatter plot of the simulated values is produced with
#'   density contours superimposed.  For \code{d > 2} pairwise plots of the
#'   simulated values are produced.
#'
#' @param x an object of class "ru", a result of a call to \code{ru}.
#' @param y Not used.
#' @param ... Additional arguments passed on to \code{hist}, \code{lines},
#'   \code{contour} or \code{points}.
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
#' @examples
#' # Log-normal density ----------------
#' x <- ru(logf = dlnorm, log = TRUE, d = 1, n = 1000, lower = 0, init = 1)
#'
#' \donttest{
#' plot(x)
#' }
#'
#' # Improve appearance using arguments to plot() and hist()
#' \donttest{
#' plot(x, breaks = seq(0, ceiling(max(x$sim_vals)), by = 0.25),
#'   xlim = c(0, 10))
#' }
#'
#' # Two-dimensional normal with positive association ----------------
#' rho <- 0.9
#' covmat <- matrix(c(1, rho, rho, 1), 2, 2)
#' log_dmvnorm <- function(x, mean = rep(0, d), sigma = diag(d)) {
#'   x <- matrix(x, ncol = length(x))
#'   d <- ncol(x)
#'   - 0.5 * (x - mean) %*% solve(sigma) %*% t(x - mean)
#' }
#' x <- ru(logf = log_dmvnorm, sigma = covmat, d = 2, n = 1000, init = c(0, 0))
#'
#' \donttest{
#' plot(x)
#' }
#' @seealso \code{\link{summary.ru}} for summaries of the simulated values
#'   and properties of the ratio-of-uniforms algorithm.
#' @export
plot.ru <- function(x, y, ..., n = ifelse(x$d == 1, 1001, 101),
                    prob = c(0.1, 0.25, 0.5, 0.75, 0.95, 0.99),
                    ru_scale = FALSE, rows = NULL, xlabs = NULL,
                    ylabs = NULL, points_par = list(col = 8)) {
  if (!inherits(x, "ru")) {
    stop("use only with \"ru\" objects")
  }
  if (n < 1) {
    stop("n must be no smaller than 1")
  }
  if (ru_scale) {
    plot_data <- x$sim_vals_rho
    plot_density <- x$logf_rho
  } else {
    plot_data <- x$sim_vals
    pmedian <- apply(x$sim_vals, 2, stats::median)
    hshift <- do.call(x$logf, c(list(pmedian), x$logf_args))
    xlogf <- x$logf
    plot_density <- function(x, ...) {
      xlogf(x, ...) - hshift
    }
  }
  if (!is.null(x$logf_rho_args) & ru_scale) {
    density_args <- x$logf_rho_args
  } else {
    density_args <- x$logf_args
  }
  if (x$d == 1) {
    temp <- graphics::hist(plot_data, plot = FALSE)
    a <- temp$breaks[1]
    b <- temp$breaks[length(temp$breaks)]
    h <- (b - a) / n
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
    area <- h * (yy[1] / 2 + sum(yy[2:n]) + yy[n + 1] / 2)
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
    xx <- seq(r[1, 1], r[2, 1], len = n)
    yy <- seq(r[1, 2], r[2, 2], len = n)
    zz <- matrix(NA, ncol = length(xx), nrow = length(yy))
    for (i in 1:length(xx)) {
      for (j in 1:length(yy)) {
        for_logf <- c(list(c(xx[i], yy[j])), density_args)
        zz[i, j] <- exp(do.call(plot_density, for_logf))
      }
    }
    zz[zz == -Inf] <- NA
    dx <- diff(xx[1:2])
    dy <- diff(yy[1:2])
    sz <- sort(zz)
    c1 <- cumsum(sz) * dx * dy
    c1 <- c1 / max(c1)
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
      rows <- x$d - 2
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
      for (i in 1:(ncol(x) - 1)) {
        for (j in (i + 1):ncol(x)) {
          graphics::plot(x[, i], x[, j], xlab = "", ylab = "", ...)
          graphics::title(xlab = parse(text = xlabs[i]),
                          ylab = parse(text = ylabs[j]))
        }
      }
    }
    pairwise_plots(plot_data)
  }
}

# =========================== summary.ru ===========================

#' Summarizing ratio-of-uniforms samples
#'
#' \code{summary} method for class "ru"
#'
#' @param object an object of class "ru", a result of a call to \code{ru}.
#' @param ... Additional arguments passed on to \code{summary}.
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
#' # one-dimensional standard normal ----------------
#' x <- ru(logf = function(x) -x ^ 2 / 2, d = 1, n = 1000, init = 0)
#' summary(x)
#'
#' # two-dimensional normal with positive association ----------------
#' rho <- 0.9
#' covmat <- matrix(c(1, rho, rho, 1), 2, 2)
#' log_dmvnorm <- function(x, mean = rep(0, d), sigma = diag(d)) {
#'   x <- matrix(x, ncol = length(x))
#'   d <- ncol(x)
#'   - 0.5 * (x - mean) %*% solve(sigma) %*% t(x - mean)
#' }
#' x <- ru(logf = log_dmvnorm, sigma = covmat, d = 2, n = 1000, init = c(0, 0))
#' summary(x)
#' @seealso \code{\link{ru}} for descriptions of \code{object$sim_vals} and
#'   \code{object$box}.
#' @seealso \code{\link{plot.ru}} for a diagnostic plot.
#' @export
summary.ru <- function(object, ...) {
  if (!inherits(object, "ru")) {
    stop("use only with \"ru\" objects")
  }
  sim_summary <- summary(object$sim_vals, ...)
  object <- object[c("box", "pa")]
  object$sim_summary <- sim_summary
  class(object) <- "summary.ru"
  return(object)
}

# =========================== print.summary.ru ===========================

#' Print method for objects of class "summary.ru"
#'
#' \code{print} method for an object \code{object} of class "summary.ru".
#'
#' @param x an object of class "summary.ru", a result of a call to
#'   \code{\link{summary.ru}}.
#' @param ... Additional optional arguments to be passed to
#'   \code{\link{print}}.
#'
#' @return Prints
#' \itemize{
#'   \item {a summary of the simulated values, via
#'     \code{summary(object$sim_vals)}}
#'   \item {an estimate of the probability of acceptance, i.e.
#'     \code{object$pa}}
#'   \item {information about the ratio-of-uniforms bounding box, i.e.
#'     \code{object$box}}
#' }
#' @examples
#' # one-dimensional standard normal ----------------
#' x <- ru(logf = function(x) -x ^ 2 / 2, d = 1, n = 1000, init = 0)
#' summary(x)
#'
#' # two-dimensional normal with positive association ----------------
#' rho <- 0.9
#' covmat <- matrix(c(1, rho, rho, 1), 2, 2)
#' log_dmvnorm <- function(x, mean = rep(0, d), sigma = diag(d)) {
#'   x <- matrix(x, ncol = length(x))
#'   d <- ncol(x)
#'   - 0.5 * (x - mean) %*% solve(sigma) %*% t(x - mean)
#' }
#' x <- ru(logf = log_dmvnorm, sigma = covmat, d = 2, n = 1000, init = c(0, 0))
#' summary(x)
#' @seealso \code{\link{summary.ru}} for summaries of the simulated values
#'   and properties of the ratio-of-uniforms algorithm.
#' @seealso \code{\link{plot.ru}} for a diagnostic plot.
#' @seealso \code{\link{ru}} for descriptions of \code{object$sim_vals} and
#'   \code{object$box}.
#' @export
print.summary.ru <- function(x, ...) {
  if (!inherits(x, "summary.ru")) {
    stop("use only with \"summary.ru\" objects")
  }
  cat("ru bounding box: ", "\n")
  print(x$box, ...)
  cat("\n")
  cat("estimated probability of acceptance: ", "\n")
  print(x$pa, ...)
  cat("\n")
  cat("sample summary", "\n")
  print(x$sim_summary, ...)
  invisible(x)
}

# ============================= print.ru() ================================== #

#' Print method for an \code{"ru"} object
#'
#' \code{print} method for class \code{"ru"}.
#'
#' @param x an object of class \code{"ru"}, a result of a call to
#'   \code{\link{ru}} or \code{\link{ru_rcpp}}.
#' @param ... Additional arguments.  None are used in this function.
#' @details Simply prints the call to \code{ru} or \code{ru_rcpp}.
#' @return The argument \code{x}, invisibly, as for all
#'   \code{\link[base]{print}} methods.
#' @seealso \code{\link{summary.ru}} for summaries of the simulated values
#'   and properties of the ratio-of-uniforms algorithm.
#' @seealso \code{\link{plot.ru}} for a diagnostic plot.
#' @export
print.ru <- function(x, ...) {
  if (!inherits(x, "ru")) {
    stop("use only with \"ru\" objects")
  }
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  return(invisible(x))
}
