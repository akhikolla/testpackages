#' @export
residuals.evmOpt <- function(object,...){
    object$residuals
}

#' @export
residuals.evmSim <- function(object,...){
    resid(object$map)
}

#' @export
residuals.evmBoot <- function(object,...){
    resid(object$map)
}
