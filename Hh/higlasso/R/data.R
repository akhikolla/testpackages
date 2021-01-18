#' Synthetic Example Data For Higlasso
#'
#' This synthetic data is adapted from the linear interaction simulations from
#' the higlasso paper. The data generating model is:
    #' \deqn{Y = V_1 + V_2 + V_3 + V_7 + V_1 V_2 + V_2 V_7}
#' @format A data.frame with 300 observations on 8 variables:
#' \describe{
#'   \item{Y}{Continuous response.}
#'   \item{V1-V7}{Covariates.}
#' }
"higlasso.df"
