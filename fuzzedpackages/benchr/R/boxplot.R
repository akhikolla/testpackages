#' @title Boxplot method for the `benchmark` timings
#'
#' @description
#' Displays measurement results as box plots, with R expressions on X axis and
#' execution time on Y axis.
#'
#' @param x An object of class `benchmark`.
#' @param units Character. The units to be used in printing the timings.
#' The available units are nanoseconds (`"ns"`), microseconds
#' (`"us"`), milliseconds (`"ms"`), seconds (`"s"`).
#' @param log Logical. Should times be plotted on log scale?
#' @param xlab Character. X axis label.
#' @param ylab Character. Y axis label.
#' @param horizontal Logical. If set to `TRUE`, X and Y axes will be switched. Defaults to `FALSE`.
#' @param violin Logical. Use [ggplot2::geom_violin()] instead [ggplot2::geom_boxplot()] for the `ggplot` plots.
#' @param ... Arguments passed on to [boxplot.default()].
#'
#' @details
#' If `ggplot2` package is available, it will be used. In order to switch to default
#' `boxplot` from the `graphics` package set option `benchr.use_ggplot`
#' to `FALSE`: `options(benchr.use_ggplot = FALSE)`.
#'
#' @include utils.R
#' @importFrom graphics boxplot
#' @method boxplot benchmark
#' @export
#'
#' @author Artem Klevtsov \email{a.a.klevtsov@gmail.com}
#'
#' @seealso
#' [plot.benchmark()]
#'
#' @examples
#' timings <- benchmark(
#'   rchisq(100, 0), rchisq(100, 1), rchisq(100, 2),
#'   rchisq(100, 3), rchisq(100, 5),
#'   times = 1000L
#' )
#' boxplot(timings)
boxplot.benchmark <- function(x, units = "auto", log = TRUE,
                              xlab, ylab, horizontal = FALSE, violin = FALSE, ...) {
  x <- convert_units(x, units)
  if (missing(xlab)) {
    xlab <- "Expression"
  }
  if (missing(ylab)) {
    if (log) {
      ylab <- sprintf("Log Execution Time [%s]", attr(x, "units"))
    } else {
      ylab <- sprintf("Execution Time [%s]", attr(x, "units"))
    }
  }
  if (is_installed("ggplot2") && getOption("benchr.use_ggplot", TRUE)) {
    return(boxplot_ggplot(x, log, xlab, ylab, horizontal, violin))
  } else {
    boxplot_default(x, log, xlab, ylab, horizontal, ...)
  }
  invisible()
}

# ggplot2 variant of the boxplot
boxplot_ggplot <- function(x, log = TRUE, xlab, ylab, horizontal = FALSE, violin = FALSE) {
  expr <- as.symbol("expr")
  time <- as.symbol("time")
  aes <- ggplot2::aes(x = expr, y = time)
  p <- ggplot2::ggplot(x, aes) +
    ggplot2::labs(x = xlab, y = ylab)
  if (violin) {
    p <- p + ggplot2::geom_violin()
  } else {
    p <- p + ggplot2::geom_boxplot()
  }
  if (log) {
    p <- p + ggplot2::scale_y_log10()
  } else {
    p <- p + ggplot2::expand_limits(y = 0)
  }
  if (horizontal) {
    p <- p + ggplot2::coord_flip()
  }
  p
}

# graphics variant of the boxplot
#' @importFrom graphics boxplot.default
boxplot_default <- function(x, log = TRUE, xlab, ylab, horizontal = FALSE, ...) {
  if (horizontal) {
    s <- ylab
    ylab <- xlab
    xlab <- s
    ll <- if (log) "x" else ""
  } else {
    ll <- if (log) "y" else ""
  }
  boxplot.default(split(x$time, x$expr),
    xlab = xlab, ylab = ylab,
    log = ll, horizontal = horizontal, ...
  )
}
