#' rIsing: High-Dimensional Ising Model Selection.
#'
#' Fits an Ising model to a binary dataset using L1-regularized logistic regression and BIC.
#' Also includes a fast lasso logistic regression function for high-dimensional problems. Uses the
#' 'libLBFGS' optimization library by Naoki Okazaki.
#'
#' @section rIsing functions:
#' \itemize{
#'   \item \code{logreg}: L1-regularized logistic regression using OWL-QN L-BFGS-B optimization.
#'   \item \code{Ising}: Ising Model selection using L1-regularized logistic regression and extended BIC.
#' }
#'
#' @import Rcpp
#' @import data.table
#' @importFrom stats complete.cases
#' @useDynLib rIsing
#' @docType package
#' @name rIsing
NULL
