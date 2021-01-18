#' R6 class to obtain preliminary result for semi-confirmatory structural equation modeling
#'
#' @docType class
#' @useDynLib lslx
#' @import stats
#' @importFrom Rcpp sourceCpp
#' @importFrom R6 R6Class
#' @keywords NULL
#' @format NULL
#' @usage NULL
#' @return Object of \code{prelslx} R6 class.
#' @export
prelslx <-
  R6::R6Class(
    classname = "prelslx",
    private = list(
      model = "lslxModel",
      data = "lslxData",
      fitting = "lslxFitting"
    )
  )
