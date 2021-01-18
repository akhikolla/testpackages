#' rust: Ratio-of-Uniforms Simulation with Transformation
#'
#' Uses the multivariate generalized ratio-of-uniforms method to simulate from a
#' distribution with log-density \code{logf} (up to an additive constant).
#' \code{logf} must be bounded, perhaps after a transformation of variable.
#'
#' @details
#'
#' The main functions in the rust package are \code{\link{ru}} and
#'   \code{\link{ru_rcpp}}, which implement the generalized ratio-of-uniforms
#'   algorithm.  The latter uses the Rcpp package to improve efficiency.
#'   Also provided are two functions,
#'   \code{find_lambda} and \code{find_lambda_one_d}, that may be used to
#'   set a suitable value for the parameter \code{lambda} if Box-Cox
#'   transformation is used prior to simulation.
#'   If \code{\link{ru_rcpp}} is used the equivalent functions are
#'   \code{find_lambda_rcpp} and \code{find_lambda_one_d_rcpp}
#'   Basic \code{plot} and \code{summary} methods are also provided.
#'
#' See \code{vignette("rust-a-vignette", package = "rust")} for an overview of
#' the package.
#'
#' @references Wakefield, J. C., Gelfand, A. E. and Smith, A. F. M. Efficient
#'  generation of random variates via the ratio-of-uniforms method. Statistics
#'  and Computing (1991) 1, 129-133.
#'  \url{https://doi.org/10.1007/BF01889987}.
#' @references Box, G. and Cox, D. R. (1964) An Analysis of Transformations.
#'  Journal of the Royal Statistical Society. Series B (Methodological), 26(2),
#'  211-252.
#' @references Eddelbuettel, D. and Francois, R. (2011). Rcpp: Seamless R and
#'   C++ Integration. Journal of Statistical Software, 40(8), 1-18.
#'   \url{https://www.jstatsoft.org/v40/i08/}.
#' @references Eddelbuettel, D. (2013) Seamless R and C++ Integration with
#'   Rcpp. Springer, New York. ISBN 978-1-4614-6867-7.
#'
#' @seealso \code{\link{ru}} and \code{\link{ru_rcpp}} to perform
#'   ratio-of-uniforms sampling.
#' @seealso \code{\link{summary.ru}} for summaries of the simulated values
#'   and properties of the ratio-of-uniforms algorithm.
#' @seealso \code{\link{plot.ru}} for a diagnostic plot.
#' @seealso \code{\link{find_lambda_one_d}} and
#'   \code{\link{find_lambda_one_d_rcpp}} to produce (somewhat) automatically
#'   a list for the argument \code{lambda} of \code{ru} for the
#'   \code{d} = 1 case.
#' @seealso \code{\link{find_lambda}} and \code{\link{find_lambda_rcpp}}
#'   to produce (somewhat) automatically
#'   a list for the argument \code{lambda} of \code{ru} for any value of
#'   \code{d}.
#'
#' @docType package
#' @name rust
#' @importFrom stats runif
#' @useDynLib rust, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
