#' womblR
#'
#' This package implements a spatiotemporal boundary detection
#' with a dissimilarity metric for areal data with inference in a Bayesian setting
#' using Markov chain Monte Carlo (MCMC). The response variable can be modeled as
#' Gaussian (no nugget), probit or Tobit link and spatial correlation is introduced
#' at each time point through a conditional autoregressive (CAR) prior. Temporal
#' correlation is introduced through a hierarchical structure and can be specified as
#' exponential or first-order autoregressive. Full details of the the package can be found
#' in the accompanying vignette. Furthermore, the details of the package can be found in
#' "Diagnosing Glaucoma Progression with Visual Field Data Using a Spatiotemporal Boundary
#' Detection Method", by Berchuck et al (2018), <arXiv:1805.11636>. The paper is in press
#' at the Journal of the American Statistical Association.
#'
#' @author Samuel I. Berchuck \email{sib2@duke.edu}
#'
#' @name womblR
#' @docType package
#' @import Rcpp
#' @importFrom graphics abline axis layout par plot points polygon title segments symbols rect text
#' @importFrom grDevices  col2rgb colorRampPalette
#' @importFrom utils tail
#' @importFrom stats lm sd var
#' @importFrom msm rtnorm
#' @importFrom mvtnorm pmvnorm
#' @useDynLib womblR
NULL
