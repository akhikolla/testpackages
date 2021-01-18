#' Plot qmix object
#'
#' General plot function for \code{qmix} objects, which dispatches the chosen
#' type of plotting to the corresponding function.
#'
#' @param x A \code{qmix} object to be plotted.
#' @param type Character string giving the type of plotting. The options are
#'   \code{"trace"} for trace plots, \code{"coef"} for coefficient plots. The default is "coef".
#' @param ... Additional arguments to be passed to subsequent plot functions (check the See Also section).
#'
#'
#' @return None.
#'
#' @seealso \code{\link{plot_trace.qmix}} and \code{\link{plot_coef.qmix}}.
#' @export
#'
#'
#'
plot.qmix <- function(x, type = "coef", ...) {
  printFunName <- paste0("plot_", type, ".qmix")
  do.call(printFunName, args = c(list(object = x), list(...)))
}



#' Trace plots for qmix
#'
#' \code{plot_trace.qmix} is used to produce trace plots from a \code{qmix} object from the main function \code{\link{qmix}}.
#'
#' @param object A \code{qmix} object from running the main function \code{\link{qmix}}.
#' @param ... Additional parameters to be passed to \code{\link[rstan]{traceplot}}.
#'
#' @return None.
#' @export
#'
#'
plot_trace.qmix <- function(object, ...) {
  if (object$binarylogic == TRUE & object$design == "fixed") {
    names(object$stanfit)[-length(names(object$stanfit))] <-
      c(paste0(
        rep(paste0("C", 1:object$nmix, ": "), object$nmix),
        rep(object$xnames, each = object$nmix)
      ),
      paste0("C", 1:object$nmix, ": proportion"))
  } else if (object$binarylogic == FALSE &
             object$design == "fixed") {
    names(object$stanfit)[-length(names(object$stanfit))] <-
      c(
        paste0(
          rep(paste0("C", 1:object$nmix, ": "), object$nmix),
          rep(object$xnames, each = object$nmix)
        ),
        paste0("C", 1:object$nmix, ": proportion"),
        paste0("C", 1:object$nmix, ": sigma")
      )
  } else if (object$binarylogic == TRUE &
             object$design == "random") {
    names(object$stanfit)[-length(names(object$stanfit))] <-
      c(
        paste0(
          rep(paste0("C", 1:object$nmix, ": "), object$nmix),
          rep(object$xnames, each = object$nmix)
        ),
        paste0("C", 1:object$nmix, ": proportion"),
        paste0("C", 1:object$nmix, ": quantile")
      )
  } else {
    names(object$stanfit)[-length(names(object$stanfit))] <-
      c(
        paste0(
          rep(paste0("C", 1:object$nmix, ": "), object$nmix),
          rep(object$xnames, each = object$nmix)
        ),
        paste0("C", 1:object$nmix, ": proportion"),
        paste0("C", 1:object$nmix, ": sigma"),
        paste0("C", 1:object$nmix, ": quantile")
      )
  }
  rstan::traceplot(object$stanfit, ...)
}

#' Make coefficient plots for a \code{qmix} object
#'
#' \code{plot_coef.qmix} is used to produce coefficient plots from a \code{qmix} object.
#'
#' @param object A \code{qmix} object from running the main function \code{\link{qmix}}.
#' @param ... Additional parameters to be passed to \code{\link[rstan]{stan_plot}}.
#'
#' @return None.
#' @export
#'
#'
plot_coef.qmix <- function(object, ...) {
  if (object$binarylogic == TRUE & object$design == "fixed") {
    names(object$stanfit)[-length(names(object$stanfit))] <-
      c(paste0(
        rep(paste0("C", 1:object$nmix, ": "), object$nmix),
        rep(object$xnames, each = object$nmix)
      ),
      paste0("C", 1:object$nmix, ": proportion"))
  } else if (object$binarylogic == FALSE &
             object$design == "fixed") {
    names(object$stanfit)[-length(names(object$stanfit))] <-
      c(
        paste0(
          rep(paste0("C", 1:object$nmix, ": "), object$nmix),
          rep(object$xnames, each = object$nmix)
        ),
        paste0("C", 1:object$nmix, ": proportion"),
        paste0("C", 1:object$nmix, ": sigma")
      )
  } else if (object$binarylogic == TRUE &
             object$design == "random") {
    names(object$stanfit)[-length(names(object$stanfit))] <-
      c(
        paste0(
          rep(paste0("C", 1:object$nmix, ": "), object$nmix),
          rep(object$xnames, each = object$nmix)
        ),
        paste0("C", 1:object$nmix, ": proportion"),
        paste0("C", 1:object$nmix, ": quantile")
      )
  } else {
    names(object$stanfit)[-length(names(object$stanfit))] <-
      c(
        paste0(
          rep(paste0("C", 1:object$nmix, ": "), object$nmix),
          rep(object$xnames, each = object$nmix)
        ),
        paste0("C", 1:object$nmix, ": proportion"),
        paste0("C", 1:object$nmix, ": sigma"),
        paste0("C", 1:object$nmix, ": quantile")
      )
  }
  rstan::plot(object$stanfit, ...)
}
