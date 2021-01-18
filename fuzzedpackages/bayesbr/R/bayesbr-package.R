#' @title bayesbr: A package for Bayesian Beta Regression
#'
#' @description The package fits or the beta regression model under the view of Bayesian statistics using the No-U-Turn-Sampler (NUTS) method for computational calculations. In addition to showing the coefficients, the package also has functions for displaying residuals, checking the model's convergence, checking the quality of the model and other utilities that may be useful.
#'
#' @docType package
#' @name bayesbr-package
#' @aliases bayesbr-package
#' @useDynLib bayesbr, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom tidyr drop_na
#' @import dplyr
#' @import stringr
#' @import ggplot2
#' @import fdrtool
#' @import Formula
#' @importFrom stats as.formula cor dbeta median pbeta qlogis qnorm quantile sd model.frame model.matrix model.response
#' @import loo
#' @import coda
#' @importFrom RcppParallel RcppParallelLibs
#'@references
#'\href{https://arxiv.org/abs/1111.4246}{arXiv:1111.4246} Hoffman, M. D., and Gelman, A. (2014). The No-U-Turn sampler: adaptively setting path lengths in Hamiltonian Monte Carlo. \emph{Journal of Machine Learning Research}, \bold{15}(1), 1593-1623.
#'@references
#'  \doi{10.1080/0266476042000214501} Ferrari, S.L.P., and Cribari-Neto, F. (2004).
#'Beta Regression for Modeling Rates and Proportions. \emph{Journal of Applied Statistics}, \bold{31}(7), 799--815.
#'
NULL
