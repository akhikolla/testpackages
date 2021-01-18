# =========================== plot.evpred ===========================

#' Plot diagnostics for an evpred object
#'
#' \code{plot} method for class "evpred".  Plots summarising the predictive
#' distribution of the largest value to be observed in N years are produced.
#' The plot produced depends on \code{x$type}.
#' If \code{x$type = "d", "p"} or \code{"q"} then
#' \code{\link[graphics]{matplot}} is used to produce a line plot of the
#' predictive density, distribution or quantile function, respectively,
#' with a line for each value of N in \code{x$n_years}.
#' If \code{x$type = "r"} then estimates of the predictive density
#' (from \code{\link[stats]{density}}) are plotted with a line for each N.
#' If \code{x$type = "i"} then lines representing estimated predictive
#' intervals are plotted, with the level of the interval indicated next to
#' the line.
#'
#' @param x An object of class "evpost", a result of a call to
#'   \code{\link{rpost}}.
#' @param ... Additional arguments passed on to \code{\link{matplot}}.
#' @param leg_pos A character scalar. Keyword for the position of legend.
#'   See \code{\link{legend}}.
#' @param leg_text A character or expression vector.  Text for legend.
#'   See \code{\link{legend}}.
#' @param which_int  A character scalar.  If \code{x$type = "i"} which
#'   intervals should be plotted?  \code{"long"} for equi-tailed intervals,
#'   \code{"short"} for the shortest possible intervals, \code{"both"} for
#'   both.
#' @seealso \code{\link{predict.evpost}} for the S3 \code{predict} method
#'  for objects of class \code{evpost}.
#' @examples
#' data(portpirie)
#' mat <- diag(c(10000, 10000, 100))
#' pn <- set_prior(prior = "norm", model = "gev", mean = c(0,0,0), cov = mat)
#' gevp  <- rpost(n = 1000, model = "gev", prior = pn, data = portpirie)
#'
#' # Predictive density function
#' d_gevp <- predict(gevp, type = "d", n_years = c(100, 1000))
#' plot(d_gevp)
#'
#' # Predictive distribution function
#' p_gevp <- predict(gevp, type = "p", n_years = c(100, 1000))
#' plot(p_gevp)
#'
#' # Predictive quantiles
#' q_gevp <- predict(gevp, type = "q", n_years = c(100, 1000))
#' plot(q_gevp)
#'
#' # Predictive intervals
#' i_gevp <- predict(gevp, type = "i", n_years = c(100, 1000), hpd = TRUE)
#' plot(i_gevp, which_int = "both")
#'
#' # Sample from predictive distribution
#' r_gevp <- predict(gevp, type = "r", n_years = c(100, 1000))
#' plot(r_gevp)
#' plot(r_gevp, xlim = c(4, 10))
#'
#' @export
plot.evpred <- function(x, ..., leg_pos = NULL, leg_text = NULL,
                        which_int = c("long", "short", "both")) {
  object <- x
  if (!inherits(object, "evpred")) {
    stop("use only with \"evpred\" objects, produced by predict.evpost")
  }
  which_int <- match.arg(which_int)
  my_title <- function(..., asp) {
    graphics::title(...)
  }
  my_segments <- function(..., lty = 2) {
    graphics::segments(..., lty = lty)
  }
  my_legend <- function(x, legend, ..., xlab, ylab, main, sub, line, outer,
                        asp, cex.lab, cex.axis, cex.main, cex.sub, col.axis,
                        col.lab, col.main, col.sub, font, font.axis, font.lab,
                        font.main, font.sub, lab, las, lend, tck, xaxp, xaxs,
                        xaxt, xlog, xpd, yaxp, yaxs, yaxt, xlim, ylim) {
    graphics::legend(x = x, legend = legend, ...)
  }
  if (is.null(leg_pos)) {
    if (object$type %in% c("d", "r")) {
      leg_pos <- "topright"
    } else {
      leg_pos <- "bottomright"
    }
  }
  temp <- list(...)
  type <- object$type
  n_years <- object$n_years
  level <- object$level
  n_y <- length(n_years)
  n_l <- length(level)
  x <- object$x
  y <- object$y
  if (object$type == "i") {
    if ((which_int == "short" | which_int == "both") &
        all(is.na(object$short[, 1:2]))) {
      which_int <- "long"
      warning("All hpd intervals missing, plot has equi-tailed intervals only")
    }
    if (which_int == "long") {
      x <- object$long[, 1:2, drop = FALSE]
    } else if (which_int == "short") {
      x <- object$short[, 1:2, drop = FALSE]
    } else {
      x <- object$long[, 1:2, drop = FALSE]
      x2 <- object$short[, 1:2, drop = FALSE]
    }
    y <- object$long[, 3]
    y <- factor(y)
    levels(y) <- 1:(n_y * n_l)
    y <- as.numeric(y)
    epy <- 0.5 / n_l
    epy2 <- epy / 4
    shoof <- ((n_l - 1) : 0) * epy
    y <- y + shoof
    if (which_int == "both") {
      y2 <- y + epy2
      y_all <- cbind(y, y2)
      x_all <- cbind(x, x2)
    } else {
      y_all <- y
      x_all <- x
    }
    y_lab <- y + epy2 / 2
    graphics::matplot(range(x_all), range(y_all), type = "n", axes = FALSE,
                      ann = FALSE, ...)
    graphics::segments(x[, 1], y, x[, 2], y, ...)
    if (which_int == "both") {
      my_segments(x2[, 1], y2, x2[, 2], y2, ...)
    }
    u <- graphics::par("usr")
    epx <- (u[2] - u[1]) / 50
    if (is.null(temp$cex)) {
      graphics::text(x[, 2] + epx, y_lab, labels = as.character(level),
                     cex = 0.7, xpd = TRUE, ...)
    } else {
      graphics::text(x[, 2] + epx, y_lab, labels = as.character(level),
                     xpd = TRUE, ...)
    }
    graphics::axis(1, ...)
    if (which_int == "both") {
      my_at <- 1:n_y + epy * (n_l - 1) / 2 + epy2 / 2
    } else {
      my_at <- 1:n_y + epy * (n_l - 1) / 2
    }
    graphics::axis(2, at = my_at, labels = n_years, ...)
    graphics::box(bty ="l", ...)
    if (is.null(temp$xlab)) {
      graphics::title(xlab = "quantile")
    }
    if (is.null(temp$ylab)) {
      graphics::title(ylab = "time horizon N, in years")
    }
    my_title(...)
  }
  if (object$type == "r") {
    x_d <- y_d <- matrix(NA, ncol = n_y, nrow = 512)
    for (i in 1:n_y) {
      temp <- stats::density(y[, i])
      y_d[, i] <- temp$y
      x_d[, i] <- temp$x
    }
    x <- x_d
    y <- y_d
  }
  if (object$type %in% c("p", "d", "q", "r")) {
    if (is.null(temp$col)) {
      graphics::matplot(x, y, type = "l", ann = FALSE, col = 1, ...)
    } else {
      graphics::matplot(x, y, type = "l", ann = FALSE, ...)
    }
    if (is.null(temp$xlab)) {
      graphics::title(xlab = switch(type, p = "quantile", d = "quantile",
                      r = "quantile", q = "probability", "i" = "quantile"))
    }
    if (is.null(temp$ylab)) {
      graphics::title(ylab = switch(type, p = "probability", d = "density",
                      r = "density", q = "quantile", "i" = "quantile"))
    }
    my_title(...)
    if (is.null(leg_text)) {
      leg_text <- paste("N =", n_years)
    }
    if (is.null(temp$lty)) {
      my_legend(x = leg_pos, legend = leg_text, lty = 1:n_y, ...)
    } else {
      my_legend(x = leg_pos, legend = leg_text, ...)
    }
  }
}
