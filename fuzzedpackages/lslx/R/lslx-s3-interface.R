#' S3 interface for semi-confirmatory SEM via PL
#' 
#' \code{plsem()} is an \code{S3} interface for obaining a fitted \code{lslx} object. 
#' 
#' @return A fitted \code{lslx} object
#' @param model A \code{character} with length one to represent the model specification. 
#' @param data A \code{data.frame} of raw data. 
#' It must contains variables specified in \code{model} (and possibly the variables specified by \code{group_variable} and \code{weight_variable}).
#' @param penalty_method A \code{character} to specify the penalty method.
#'   The current version supports \code{"none"}, \code{"lasso"}, \code{"ridge"}, \code{"elastic"}, and \code{"mcp"}.
#' @param lambda_grid A non-negative \code{numeric} to specify penalty levels for both \code{"lasso"} and \code{"mcp"}.
#'   If it is set as \code{"default"}, its value will be generated automatically based on the variable scales.
#' @param delta_grid A non-negative \code{numeric} to specify the convexity level for \code{"mcp"}.
#'   If it is set as \code{"default"}, its value will be generated automatically based on the variable scales.
#' @param numeric_variable A \code{character} to specify which response variables should be transfromed into \code{numeric}.
#' @param ordered_variable A \code{character} to specify which response variables should be transfromed into \code{ordered}.
#' @param weight_variable A \code{character} with length one to specify what variable is used for sampling weight.
#' @param auxiliary_variable A \code{character} to specify what variable(s) is used as auxiliary variable(s) for estimating saturated moments when missing data presents and two-step method is implemented. 
#'   Auxiliary variable(s) must be numeric. If any categorical auxiliary is considered, please transform it into dummy variables before initialization.
#' @param group_variable A \code{character} with length one to specify what variable is used for labeling group.
#' @param reference_group A \code{character} with length one to specify which group is set as reference. 
#' @param sample_cov A numeric \code{matrix} (single group case) or a \code{list} of numeric \code{matrix} (multi-group case) to represent sample covariance matrixs. It must have row and column names that match the variable names specified in \code{model}.
#' @param sample_mean A \code{numeric} (single group case) or a \code{list} of \code{numeric} (multi-group case) to represent sample mean vectors.
#' @param sample_size A \code{numeric} (single group case) with length one or a \code{list} of \code{numeric} (multi-group case) to represent the sample sizes.
#' @param sample_moment_acov A numeric \code{matrix} (single group case) or a \code{list} of numeric \code{matrix} (multi-group case) to represent asymptotic covariance for moments. 
#' @param verbose A \code{logical} to specify whether messages made by \code{lslx} should be printed.
#' @param ... Other arguments. For details, please see the documentation of \code{lslx}.
#' @examples 
#' ## EXAMPLE: Semi-Confirmatory Factor Analysis with lavaan Style ##
#' # specify a factor analysis model with lavaan style
#' model_fa <- "visual  =~ x1 + x2 + x3
#'              textual =~ x4 + x5 + x6
#'              speed   =~ x7 + x8 + x9
#'              pen() * visual  =~ x4 + x5 + x6 + x7 + x8 + x9
#'              pen() * textual =~ x1 + x2 + x3 + x7 + x8 + x9
#'              pen() * speed   =~ x1 + x2 + x3 + x4 + x5 + x6
#'              visual  ~~ 1 * visual
#'              textual ~~ 1 * textual
#'              speed   ~~ 1 * speed"
#'              
#' # fit with mcp under specified penalty levels and convexity levels
#' lslx_fa <- plsem(model = model_fa, 
#'                  data = lavaan::HolzingerSwineford1939,
#'                  penalty_method = "mcp", 
#'                  lambda_grid = seq(.02, .60, .02), 
#'                  delta_grid = c(1.5, 3.0, Inf))
#' 
#' # summarize fitting result under the penalty level selected by 'bic'
#' summary(lslx_fa, selector = "bic")
#' 
#' @export
## \code{plsem()} is a wrapper for \code{lslx$()$fit()}. ##
plsem <- function(model, 
                  data,
                  penalty_method = "mcp",
                  lambda_grid = "default",
                  delta_grid = "default",
                  numeric_variable,
                  ordered_variable,
                  weight_variable,
                  auxiliary_variable,
                  group_variable,
                  reference_group,
                  sample_cov,
                  sample_mean,
                  sample_size,
                  sample_moment_acov,
                  verbose = TRUE,
                  ...) {
  r6_lslx <- lslx$new(model = model,
                      data = data,
                      numeric_variable = numeric_variable,
                      ordered_variable = ordered_variable,
                      weight_variable = weight_variable,
                      auxiliary_variable = auxiliary_variable,
                      group_variable = group_variable,
                      reference_group = reference_group,
                      sample_cov = sample_cov,
                      sample_mean = sample_mean,
                      sample_size = sample_size,
                      sample_moment_acov = sample_moment_acov,
                      verbose = verbose)
  r6_lslx$fit(penalty_method = penalty_method,
              lambda_grid = lambda_grid,
              delta_grid = delta_grid,
              verbose = verbose,
              ...)
  return(invisible(r6_lslx))
}


#' @export
## \code{print.lslx()} is an S3 method to summarize \code{lslx} fitting results.##
print.lslx <- function(x, ...) {
  x$print()
}



#' S3 method to summarize \code{lslx} fitting results
#' 
#' \code{summary.lslx()} is an \code{S3} interface for summarizing \code{lslx} fitting results. 
#' 
#' @param object A fitted \code{lslx} object. 
#' @param selector A \code{character} to specify a selector for determining an optimal penalty level.
#' @param lambda A \code{numeric} to specify a chosen optimal lambda value.
#' @param delta A \code{numeric} to specify a chosen optimal lambda value.
#' @param ... Other arguments. For details, please see the \code{$summarize()} method in \code{lslx}.
#' @export
## \code{summary.lslx()} is an S3 method to summarize \code{lslx} fitting results.##
summary.lslx <- function(object,
                         selector,
                         lambda,
                         delta,
                         ...) {
  object$summarize(selector = selector,
                   lambda = lambda,
                   delta = delta)
}

#' S3 method to extract parameter estimate from \code{lslx}
#' 
#' \code{coef.lslx()} is an \code{S3} interface for extracting parameter estimate from a \code{lslx} object. 
#'
#' @param object A fitted \code{lslx} object. 
#' @param selector A \code{character} to specify a selector for determining an optimal penalty level.
#' @param lambda A \code{numeric} to specify a chosen optimal lambda value.
#' @param delta A \code{numeric} to specify a chosen optimal lambda value.
#' @param ... Other arguments. For details, please see the \code{$extracted_coefficient()} method in \code{lslx}.
#' @export
## \code{coef.lslx()} is an S3 method to extract model parameters from \code{lslx}##
coef.lslx <- function(object,
                      selector,
                      lambda,
                      delta,
                      ...) {
  object$extract_coefficient(selector = selector,
                             lambda = lambda,
                             delta = delta)
}

#' S3 method to extract covariance matrix of estimates from \code{lslx}
#' 
#' \code{vcov.lslx()} is an \code{S3} interface for extracting covariance matrix of parameter estimate from a \code{lslx} object.
#'
#' @param object A fitted \code{lslx} object. 
#' @param selector A \code{character} to specify a selector for determining an optimal penalty level.
#'   Its value can be any one in \code{"aic"}, \code{"aic3"}, \code{"caic"}, \code{"bic"}, \code{"abic"}, \code{"hbic"},
#'   or their robust counterparts \code{"raic"}, \code{"raic3"}, \code{"rcaic"}, \code{"rbic"}, \code{"rabic"}, \code{"rhbic"} if raw data is available.
#' @param lambda A \code{numeric} to specific a chosen optimal penalty level. 
#'   If the specified \code{lambda} is not in \code{lambda_grid}, a nearest legitimate value will be used.
#' @param delta A \code{numeric} to specific a chosen optimal convexity level.
#'   If the specified \code{delta} is not in \code{delta_grid}, a nearest legitimate value will be used.
#' @param ... Other arguments. For details, please see the \code{$extracted_coefficient_acov()} method in \code{lslx}.
#' @export
## \code{vcov.lslx()} is an S3 method to extract asymptotic covariance matrix of parameter estimates from \code{lslx}.##
vcov.lslx <- function(object,
                      selector,
                      lambda,
                      delta,
                      ...) {
  object$extract_coefficient_acov(selector = selector,
                                  lambda = lambda,
                                  delta = delta,
                                  ...)
}

#' S3 method to extract model-implied moments from \code{lslx}
#' 
#' \code{fitted.lslx()} is an \code{S3} interface for extracting model-implied moments from a \code{lslx} object.
#'
#' @param object A fitted \code{lslx} object. 
#' @param selector A \code{character} to specify a selector for determining an optimal penalty level.
#'   Its value can be any one in \code{"aic"}, \code{"aic3"}, \code{"caic"}, \code{"bic"}, \code{"abic"}, \code{"hbic"},
#'   or their robust counterparts \code{"raic"}, \code{"raic3"}, \code{"rcaic"}, \code{"rbic"}, \code{"rabic"}, \code{"rhbic"} if raw data is available.
#' @param lambda A \code{numeric} to specific a chosen optimal penalty level. 
#'   If the specified \code{lambda} is not in \code{lambda_grid}, a nearest legitimate value will be used.
#' @param delta A \code{numeric} to specific a chosen optimal convexity level.
#'   If the specified \code{delta} is not in \code{delta_grid}, a nearest legitimate value will be used.
#' @param ... Other arguments. For details, please see the \code{$extracted_implied_mean()} and the \code{$extracted_implied_cov()} methods in \code{lslx}.
#' @export
## \code{fitted.lslx()} is an S3 method to extract model-implied moments from \code{lslx}.##
fitted.lslx <- function(object,
                        selector,
                        lambda,
                        delta,
                        ...) {
  implied_mean <- 
    object$extract_implied_mean(selector = selector,
                                lambda = lambda,
                                delta = delta,
                                ...)
  implied_cov <- 
    object$extract_implied_cov(selector = selector,
                               lambda = lambda,
                               delta = delta,
                               ...)
  return(list(mean = implied_mean,
              cov = implied_cov))
}

#' S3 method to extract residual moments from \code{lslx}
#' 
#' \code{residuals.lslx()} is an \code{S3} interface for extracting residuals from a \code{lslx} object.
#'
#' @param object A fitted \code{lslx} object. 
#' @param selector A \code{character} to specify a selector for determining an optimal penalty level.
#'   Its value can be any one in \code{"aic"}, \code{"aic3"}, \code{"caic"}, \code{"bic"}, \code{"abic"}, \code{"hbic"},
#'   or their robust counterparts \code{"raic"}, \code{"raic3"}, \code{"rcaic"}, \code{"rbic"}, \code{"rabic"}, \code{"rhbic"} if raw data is available.
#' @param lambda A \code{numeric} to specific a chosen optimal penalty level. 
#'   If the specified \code{lambda} is not in \code{lambda_grid}, a nearest legitimate value will be used.
#' @param delta A \code{numeric} to specific a chosen optimal convexity level.
#'   If the specified \code{delta} is not in \code{delta_grid}, a nearest legitimate value will be used.
#' @param ... Other arguments. For details, please see the \code{$extracted_residual_mean()} and the \code{$extracted_residual_cov()} methods in \code{lslx}.
#' @export
## \code{residuals.lslx()} is an S3 method to extract \residual moments from \code{lslx}.##
residuals.lslx <- function(object,
                           selector,
                           lambda,
                           delta,
                           ...) {
  residual_mean <- 
    object$extract_residual_mean(selector = selector,
                                 lambda = lambda,
                                 delta = delta,
                                 ...)
  residual_cov <- 
    object$extract_residual_cov(selector = selector,
                                lambda = lambda,
                                delta = delta,
                                ...)
  return(list(mean = residual_mean,
              cov = residual_cov))
}

