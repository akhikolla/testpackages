#' Plot cbq object
#'
#' General plot function for \code{cbq} objects, which dispatches the chosen
#' type of plotting to the corresponding function.
#'
#' @param x A \code{cbq} object to be plotted.
#' @param type Character string giving the type of plotting. The options are
#'   \code{"trace"} for trace plots, \code{"coef"} for coefficient plots. The default is the traceplot.
#' @param ... Additional arguments to be passed to subsequent plot functions.
#'
#' @export
#' @return None.
#'
#'
plot.cbq <- function(x, type = "trace", ...) {
  printFunName <- paste0("plot_", type, ".cbq")
  do.call(printFunName, args = c(list(object = x), list(...)))
}



#' Make traceplots for cbq
#'
#' Plot traceplots from a \code{cbq} object.
#'
#' @param object A \code{cbq} object.
#' @param ... Additional parameters to be passed to the traceplot function.
#'
#' @export
#' @return None.
#'
#'
plot_trace.cbq <- function(object, ...) {
  names(object$stanfit)[1:object$npars] <- object$xnames
  rstan::traceplot(object$stanfit, ...)
}

#' Make coefficient plots for cbq
#'
#' Plot traceplots from a \code{cbq} object.
#'
#' @param object A \code{cbq} object.
#' @param ... Additional parameters to be passed to the plot function.
#'
#' @export
#' @return None.
#'
#'
plot_coef.cbq <- function(object, ...) {
  names(object$stanfit)[1:object$npars] <- object$xnames
  rstan::plot(object$stanfit, ...)
}
