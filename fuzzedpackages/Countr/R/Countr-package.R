#' @docType package
#' @name Countr-package
#' @aliases  Countr-package Countr
#'
#' @title
#' Flexible Univariate Count Models Based on Renewal Processes
#'
#' @description
#' Flexible univariate count models based on renewal
#' processes. The models may include covariates and can be specified
#' with familiar formula syntax as in glm() and 'flexsurv'.
#'
#' @useDynLib Countr
#'
#' @import Matrix Rcpp Formula flexsurv dplyr
#' @importFrom stats nobs AIC  coef  confint  confint.default
#' @importFrom stats formula as.formula getCall  glm.fit  logLik  model.frame
#' @importFrom stats model.matrix  model.response  model.weights
#' @importFrom stats na.pass  pnorm  poisson  printCoefmat  quantile
#' @importFrom stats dnbinom dpois lm pchisq
#' @importFrom stats residuals update.formula  vcov
#' @importFrom stats density predict qqline qqnorm terms.formula
#' @importFrom graphics par plot
#' @importFrom RColorBrewer brewer.pal
#' @importFrom lattice barchart
#' @importFrom car Boot
#' @importFrom Rdpack reprompt
#' @importFrom utils capture.output
#' @importFrom standardize standardize
#' @importFrom lmtest lrtest
#' @importFrom xtable xtable
#'
#' @details
#'
#' The methodology is described in the forthcoming paper
#' \insertCiteOnly{CountrJssArticle}{Countr}
#' in the Journal of Statistical Software (included in the package as vignette
#' \code{vignette('Countr_guide_paper', package = "Countr")}).
#'
#' The main function is \code{\link{renewalCount}}, see its documentation for
#' examples.
#'
#' Goodness of fit chi-square (likelihood ratio and Pearson) tests for glm and
#' count renewal models are implemented in \code{\link{chiSq_gof}} and
#' \code{\link{chiSq_pearson}}.
#'
#' @references
#'
#' \insertRef{CountrJssArticle}{Countr}
#' 
#' \insertRef{baker2017event}{Countr}
#'
#' \insertRef{boshnakov2017bivariate}{Countr}
#'
#' \insertRef{cameron2013regression}{Countr}
#'
#' \insertRef{TarakEtAl2018jss}{Countr}
#'
#' \insertRef{mcshane2008count}{Countr}
#'
#' \insertRef{winkelmann1995duration}{Countr}
#'
NULL

## To stop the respective NOTE from R CMD check
utils::globalVariables(
    c("Counts", "Actual" # see compareToGLM()
      )
)
