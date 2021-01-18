#' @importFrom Rcpp sourceCpp
#' @useDynLib psqn, .registration = TRUE
NULL

#' @details
#' The main method in the psqn package is the \code{\link{psqn}} function.
#' Notice that it is also possible to use the package from C++. This may
#' yield a large reduction in the computation time. See the vignette for
#' details e.g. by calling \code{vignette("psqn", package = "psqn")}.
#' A brief introduction is provided in the "quick-intro" vignette
#' (see \code{vignette("quick-intro", package = "psqn")}).
#'
#' This package is fairly new. Thus, results may change and
#' contributions and feedback is much appreciated.
#' @keywords internal
#' @aliases psqn-package
"_PACKAGE"
