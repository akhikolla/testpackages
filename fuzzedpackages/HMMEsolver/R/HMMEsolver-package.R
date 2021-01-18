#' HMMEsolver Package
#'
#' Consider the linear mixed model with normal random effects,
#' \deqn{Y = X\beta + Zv + \epsilon}
#' where \eqn{\beta} and \eqn{v} are vectors of fixed and random effects.
#' One of most popular methods to solve the Henderson's Mixed Model Equation
#' related to the problem is EM-type algorithm. Its drawback, however, comes from
#' repetitive matrix inversion during recursive estimation steps. Kim (2017) proposed
#' a novel method of avoiding such difficulty, letting the estimation more fast, stable, and
#' scalable.
#'
#' @docType package
#' @name HMMEsolver-package
#' @import Rdpack
#' @importFrom Rcpp evalCpp
#' @useDynLib HMMEsolver, .registration=TRUE
NULL

# NOTES
# tools::package_native_routine_registration_skeleton(".")
