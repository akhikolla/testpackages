#' @import methods
NULL

#' Class "geerMod" of Fitted GEE-MCD/WGEE-MCD  Models.
#'
#' @slot call the matched call
#' @slot opt the optimization result returned by optimizeGeer
#' @slot args arguments m, Y, X, Z, W, H, time
#' @slot triple an integer vector of length three containing the degrees of the
#' three polynomial functions for the mean structure, the log innovation
#' -variances and the autoregressive coefficients.
#' @slot rho a parameter used in the 'working' covariance structure.
#' @slot devcomp the deviance components list.
#'
#' @exportClass geerMod
setClass("geerMod",
  representation(
    call = "call",
    opt = "list",
    args = "list",
    triple = "numeric",
    rho = "numeric",
    devcomp = "list"
    ))

## #' Class gee_jmcm
## #'
## #' Class \code{gee_jmcm} defines a Modified Cholesky decomposition (MCD) based
## #' joint mean covariance model within the framework of (weighted) generalised
## #' estimating equations.
## #'
## #' @name gee_jmcm-class
## #'
## #' @exportClass gee_jmcm
## setClass("gee_jmcm", representation(pointer="externalptr"))

## gee_jmcm_method <- function(name) {
##   paste("gee_jmcm", name, sep = "__")
## }

## #' Extract parts of gee_jmcm.
## #'
## #' @param x an object of gee_jmcm class
## #' @param name member function of gee_jmcm class
## setMethod("$", "gee_jmcm", function(x, name) {
##   geefun <- gee_jmcm_method(name)
##   function(...) .Call(geefun, x@pointer, ...)
## })

## setMethod("initialize", "gee_jmcm", function(.Object, ...) {
##   .Object@pointer <- .Call(gee_jmcm_method("new"), ...)
##   .Object
## })

## #' Class ipw
## #'
## #' Class \code{ipw} defines the inverse probability weighting.
## #'
## #' @name ipw-class
## #'
## #' @exportClass ipw
## setClass("ipw", representation(pointer="externalptr"))

## ipw_method <- function(name) {
##   paste("ipw", name, sep = "__")
## }

## #' Extract parts of ipw.
## #'
## #' @param x an object of ipw class
## #' @param name member function of ipw class
## setMethod("$", "ipw", function(x, name) {
##   ipwfun <- ipw_method(name)
##   function(...) .Call(ipwfun, x@pointer, ...)
## })

## setMethod("initialize", "ipw", function(.Object, ...) {
##   .Object@pointer <- .Call(ipw_method("new"), ...)
##   .Object
## })
