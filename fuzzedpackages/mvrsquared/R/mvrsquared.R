#' mvrsquared
#'
#' Compute the Coefficient of Determination for Vector or Matrix Outcomes
#'
#' @details
#' Welcome to the \code{mvrsquared} package! This package does one thing: calculate
#' the coefficient of determination or "R-squared". However, this implementation
#' is different from what you may be familiar with. In addition to the standard
#' R-squared used frequently in linear regression, `mvrsquared` calculates
#' R-squared for multivariate outcomes. (This is why there is an 'mv' in
#' \code{mvrsquared}).
#'
#' \code{mvrsquared} implements R-squared based on a derivation in this paper
#' (\url{https://arxiv.org/abs/1911.11061}). It's the same definition of R-squared
#' you're probably familiar with, i.e. 1 - SSE/SST but generalized to n-dimensions.
#'
#' In the standard case, your outcome and prediction are vectors. In other words,
#' each observation is a single number. This is fine if you are predicting a
#' single variable. But what if you are predicting multiple variables at once?
#' In that case, your outcome and prediction are matrices. This situation occurs
#' frequently in topic modeling or simultaneous equation modeling.
#'
#' @name mvrsquared
#' @docType package
NULL

#' @import Matrix
#' @import Rcpp
#' @importFrom methods as
#' @importFrom methods canCoerce
#' @useDynLib "mvrsquared", .registration=TRUE
NULL
