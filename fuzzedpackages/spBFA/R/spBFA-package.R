#' spBFA
#'
#' \code{spBFA} is a package for Bayesian spatial factor analysis. A corresponding manuscript is forthcoming.
#' 
#' @author Samuel I. Berchuck \email{sib2@duke.edu}
#'
#' @name spBFA
#' @docType package
#' @import Rcpp
#' @importFrom stats quantile rnorm runif median model.frame model.matrix
#' @importFrom graphics abline axis layout par plot points polygon title segments symbols rect text lines
#' @importFrom grDevices col2rgb colorRampPalette
#' @importFrom utils tail
#' @importFrom stats dnorm lm sd var pnorm
#' @importFrom msm rtnorm
#' @importFrom mvtnorm pmvnorm
#' @importFrom pgdraw pgdraw
#' @useDynLib spBFA
NULL
