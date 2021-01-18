# Various imports
#' @importFrom Rcpp evalCpp
#' @importFrom MASS mvrnorm
#' @importFrom stats cov2cor
#' @importFrom stats var
#' @importFrom stats cov
#' @importClassesFrom  Matrix dgCMatrix
#' @useDynLib gif, .registration = TRUE
NULL

#' @title Synthetic multivariate Gaussian data
#' @name ar1
#' @docType data
#' @description
#' A synthetic dataset includes 200 samples under multivariate Gaussian distribution with 100 variables.
#'
#' @details
#' This synthetic dataset contains 200 samples, each of them is a vector following
#' multivariate Gaussian distribution with 100 variables.
#' The inverse covariance matrix of the distribution is as follows,
#' \itemize{
#'   \code{Omega[i, i] = 1}.
#'
#'   \code{Omega[i, i + 1] = Omega[i, i - 1] = 0.5}.
#'
#'   Otherwise: \code{Omega[i, j] = 0}.
#' }
#'
#' The corresponding graph structure is the so-called AR(1).
#'
#' @format
#' \describe{
#'   \item{ar1$x}{A numeric matrix with 200 rows and 100 variables where each row represents a sample.}
#'   \item{ar1$Omega}{The corresponding inverse covariance matrix of the Gaussian graphical model.}
#' }
NULL
