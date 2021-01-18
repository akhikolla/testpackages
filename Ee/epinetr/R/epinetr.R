#' epinetr
#'
#' epinetr is a package intended to aid in the investigation of the
#' contribution of epistatic networks to complex traits,
#'
#' This package provides a range of functions for running
#' forward-time simulation using epistatic networks, including
#' visualisation tools.
#'
#' For a complete list of functions, use library(help = "epinetr").
#'
#' @docType package
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#' @importFrom graphics title
#' @importFrom Rcpp sourceCpp
#' @importFrom stats cov rnorm runif sd var
#' @importFrom methods is
#' @useDynLib epinetr, .registration = TRUE
#' @name epinetr
NULL
