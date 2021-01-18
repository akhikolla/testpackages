#' @import methods
NULL

#' Class "jmcmMod" of Fitted Joint Mean-Covariance Models.
#'
#' @slot call the matched call
#' @slot opt the optimization result returned by optimizeJmcm
#' @slot args arguments m, Y, X, Z, W, time
#' @slot triple an integer vector of length three containing the degrees of the
#' three polynomial functions for the mean structure, the log innovation
#' -variances and the autoregressive or moving average coefficients when
#' 'mcd' or 'acd' is specified for cov.method. It refers to the mean structure,
#' variances and angles when 'hpc' is specified for cov.method.
#' @slot devcomp the deviance components list
#'
#' @exportClass jmcmMod
setClass("jmcmMod",
  representation(
    call = "call",
    opt = "list",
    args = "list",
    triple = "numeric",
    devcomp = "list"
    ))