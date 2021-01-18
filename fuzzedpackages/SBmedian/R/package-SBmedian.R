#' Scalable Bayes with Median of Subset Posteriors
#' 
#' Median-of-means is a generic yet powerful framework for scalable and robust estimation. 
#' A framework for Bayesian analysis is called M-posterior, which estimates a median of subset posterior measures. 
#' 
#' @docType package
#' @name SBmedian
#' @aliases SBmedian-package
#' @import Rdpack
#' @importFrom expm expm logm
#' @importFrom stats rnorm rWishart density
#' @importFrom utils packageVersion
#' @importFrom Rcpp evalCpp
#' @useDynLib SBmedian
NULL