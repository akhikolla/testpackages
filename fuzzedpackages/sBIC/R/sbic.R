#' sBIC package documentation.
#'
#' Computes the sBIC for various model collections including:
#' \itemize{
#'  \item{Binomial mixtures}
#'  \item{Factor analyses}
#'  \item{Gaussian mixtures}
#'  \item{Latent forests}
#'  \item{Latent class analyses}
#'  \item{Reduced rank regressions}
#' }
#' The primary functionality of this package can be accessed through the
#' sBIC function.
#'
#' @importFrom Rcpp evalCpp
#' @importFrom R.oo Object setConstructorS3 extend
#' @importFrom R.methodsS3 setMethodS3 throw
#' @importFrom mclust Mclust mclustBIC meV meVVV
#' @importFrom stats as.formula cor cov factanal rexp runif
#' @useDynLib sBIC
#'
#' @name sBIC-package
NULL
