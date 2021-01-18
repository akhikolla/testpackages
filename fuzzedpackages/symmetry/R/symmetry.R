#' @useDynLib symmetry
#' @importFrom stats coef fitted model.matrix residuals rnorm rlogis
#' @importFrom utils tail
#' @importFrom Rdpack reprompt
#' @import parallel
#' @import Rcpp
#'
NULL

#' symmetry: A package which implements tests for symmetry of random samples,
#' linear models and generalized autoregressive conditional heteroskedasticity
#' (GARCH) models
#'
#' The package contains a large number of tests for symmetry (and their
#' bootstrap variants), which can be used to test the symmetry of random samples
#' or of model residuals. Currently, the supported models are linear models and
#' generalized autoregressive conditional heteroskedasticity (GARCH) models
#' (fitted with the 'fGarch' package). The tests are implemented using the
#' 'Rcpp' package which ensures great performance.
#'
#' To see the available tests, see \link{TestStatistics}
#'
#' For documentation on how to perform the tests, see \link{symmetry_test}
#'
#' @docType package
#' @name symmetry
NULL
