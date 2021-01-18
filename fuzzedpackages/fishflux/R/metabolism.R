#' A function to estimate f0 and alpha
#'
#' All model parameters below were estimated by Barneche & Allen 2018 Ecology
#' Letters doi: 10.1111/ele.12947. These parameters are for the best model
#' (Model 2 in the paper online supplementary material) of fish resting
#' metabolic rates reported in the paper, which also includes trophic level as
#' a covariate.
#'
#' @param family family fish
#' @param temp  Temperature in degrees Celsius
#' @param troph_m Trophic level mean (from 1 to 5)
#' @param troph_sd Trophic level sd (optional)
#' 
#' @details     All model parameters below were estimated by 
#' Barneche & Allen 2018 Ecology Letters doi: 10.1111/ele.12947. 
#' These parameters are for the best model 
#' (Model 2 in the paper online supplementary material) 
#' of fish resting metabolic rates reported in the paper, 
#' which also includes trophic level as a covariate.
#' 
#' @keywords    Fish metabolism
#' @importFrom rstan sampling summary
#' 
#' @examples
#' library(fishflux)
#' metabolism(family = "Pomacentridae", temp = 27, troph_m = 2)
#' 
#' @export
metabolism <- function (family, temp, troph_m, troph_sd = 0.0000000001) {

  ## get b0 and a from database
  metpar <- fishflux::metabolic_parameters
  metpars <- metpar[metpar$family == family, ]

  if (nrow(metpars) > 0) {
    message("values for b0 and alpha on family level")
  }

  if (nrow(metpars) == 0) {
    metpars <- metpar[metpar$family == "all", ]
    message("average values for b0 and alpha used")
  }

  metpars$troph_m <- troph_m
  metpars$troph_sd <- troph_sd
  metpars$temp <- temp
  metpars <- as.list(metpars[, -1])

  ## predict B0
  stanfit <-  sampling(stanmodels$getB0, data = metpars, iter = 1000,
                       algorithm = "Fixed_param", chains = 1)
  result <- as.data.frame(summary(stanfit)$summary)
  result <- result[1, c("mean", "sd")]
  colnames(result) <- c("B0_m", "B0_sd")
  result <- cbind(result, metpars)
  result <- result[, ! names(result) %in% c("troph_m", "troph_sd", "temp")]
  rownames(result) <- NULL
  colnames(result) <- c("f0_m", "f0_sd", "alpha_m", "alpha_sd", "b0_m", "b0_sd")
  result
}
