#' CLERE methodology for simultaneous variables clustering and regression
#' 
#' The methodology consists in creating clusters of variables involved in a
#' high dimensional linear regression model so as to reduce the dimensionality.
#' A model-based approach is proposed and fitted using a Stochastic EM-Gibbs
#' algorithm (SEM-Gibbs).
#' 
#' @name clere-package
#' @aliases clere clere-package
#' @docType package
#' 
#' @seealso Overview : \code{\link{clere-package}} \cr 
#' Classes : \code{\linkS4class{Clere}}, \code{\linkS4class{Pacs}} \cr 
#' Methods : \code{\link{plot}}, \code{\link{clusters}}, \code{\link{predict}}, \code{\link{summary}} \cr 
#' Functions : \code{\link{fitClere}}, \code{\link{fitPacs}} 
#' Datasets : \code{\link{numExpRealData}}, \code{\link{numExpSimData}}, \code{\link{algoComp}}
#' 
#' @examples
#'  
#'  # Simple example using simulated data
#'  # to see how to you the main function clere
#'  library(clere)
#'  x  <- matrix(rnorm(50 * 100), nrow = 50, ncol = 100)
#'  y  <- rnorm(50)
#'  model <- fitClere(y = y, x = x, g = 2, plotit = FALSE)
#'  plot(model) 
#'  clus <- clusters(model, threshold = NULL)
#'  predict(model, newx = x+1)
#'  summary(model)
#' 
NULL

## usethis namespace: start
#' @useDynLib clere, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
## usethis namespace: start
#' @import RcppEigen
## usethis namespace: end
NULL
