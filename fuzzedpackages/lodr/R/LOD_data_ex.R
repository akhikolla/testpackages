#' Simulated data with covariates subject to limits of detection
#'
#' A simulated dataset containing a generic outcome varible and three covariates,
#' two of which are subject to a lower limit of detection of 0, with a sample size
#' of 100.  See Details for information on how these data were generated.  
#' 
#' Each of the covariates were generated independently from 100 independent draws from 
#' the standard normal distributon.  The outcome variable was generated from
#' a linear model with these three covariates, along with an intercept of 1, 
#' a residual variance of 1, and regression coefficients of 1 for each covariates.
#' Then for two of the covariates, to reflect a lower limit of detection of 0, values below 
#' this limit were set to 0.  This results in a 50 percent probability of being below the limit of detection for each of the 
#' two corresponding covariates.
#'
#' @format A data frame with 100 rows and 4 variables:
#' \describe{
#'   \item{y}{Outcome}
#'   \item{x1}{First covariate , no limits of detection}
#'   \item{x2}{Second covariate, lower limit of detection of 0}
#'   \item{x3}{Third covariate, lower limit of detection of 0}
#' }
"lod_data_ex"