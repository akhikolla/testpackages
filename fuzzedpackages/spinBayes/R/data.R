#' simulated data for demonstrating the features of BVCfit
#'
#' Simulated gene expression data for demonstrating the features of BVCfit.
#'
#' @docType data
#' @keywords datasets
#' @name data
#' @aliases gExp.new gExp.L spbayes X Y Z E clin X2 Y2 Z2 E2 clin2 X.new Y.new Z.new E.new clin.new
#' @usage data("gExp")
#' data("gExp.new")
#' data("gExp.L")
#' @format gExp consists of five components: X, Y, Z, E and clin. gExp.new contains the data of new observations
#' (X.new, Y.new, Z.new, E.new and clin.new) which can be used for evaluating the prediction performance.
#'
#' gExp.L contains larger datasets: X2, Y2, Z2, E2 and clin2
#'
#' @details the same true model is used for generating Y, Y.new and Y2
#' \deqn{Y = \beta_{0}(Z)+\beta_{1}(Z)X_{1} + \beta_{2}(Z)X_{2} + 1.5X_{3} - X_{5} + 1.3E - 1.2 EX_{2}+1.3 EX_{4} - clin_{1} + 1.5 clin_{2} + \epsilon}
#' where \eqn{\epsilon\sim N(0,1)}, \eqn{\beta_{0}=2\sin(2\pi*Z)}, \eqn{\beta_{1}=2\exp(2Z-1) } and \eqn{\beta_{2}=-6Z(1-Z) }
#'
#' @examples
#' data(gExp)
#' dim(X)
#'
#' data(gExp.L)
#' dim(X)
#'
#' @seealso \code{\link{BVCfit}}
NULL
