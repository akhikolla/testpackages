#' A Package to Fit the Cox Regression from Right Truncated Data
#'
#' The method assumes that truncation is independent of covariates,
#' and of lifetime, and that there is no censoring.
#' The method uses Inverse-Probability-Weighting estimating equations with stabilized weights,
#'  IPW-S and IPW-SA, as described in Vakulenko-Lagun et al. (2018).
#'  Currently the code allows only time-independent covariates.
#'
#' The \pkg{coxrt} package provides two functions:
#' \code{\link{coxph.RT}} (IPW-S) that assumes positivity
#'  and \code{\link{coxph.RT.a0}} (IPW-SA) that allows
#'  for adjustment of estimation using plugged-in \code{a0}.
#'  The illustrative examples in these functions include analysis of AIDS
#'  latency data with age as a covariate, where the AIDS cases were retrospectively
#'  ascertained at June 30, 1986, and only those who developed AIDS by that time were
#'  included in the analysis (Kalbfeisch and Lawless, 1989).

#'
#'
#' @references Vakulenko-Lagun, B., Mandel, M., Betensky, R.A. Inverse probability weighting methods for Cox regression with right-truncated data. 2019, submitted to \emph{Biometrics}
#' @references Kalbfeisch, J.D. and Lawless, J.F. Inference based on retrospective ascertainment: an analysis
#' of the data on transfusion-related AIDS. Journal of the American Statistical Association,
#' 84 (406):360-372, 1989.

#' @import survival
#' @import BB
#' @import inline
#' @import gss
#' @import ggplot2
#' @importFrom stats pnorm qnorm quantile sd var
#' @importFrom Rcpp evalCpp
#' @useDynLib coxrt, .registration = TRUE
#'
#' @docType package
#' @name coxrt
NULL


