# Confidence Intervals on jmodelMult Parameters

confint.jmodelMult <- function(object, parm, level = 0.95, ...) {
  if (!inherits(object, "jmodelMult"))
    stop("Can only be used for 'jmodelMult' objects. \n")
  
  a <- (1 - level) / 2
  mult <- qnorm(1 - a)
  
  gammas <- object$coefficients$gamma
  Vcov <- object$Vcov
  ncb <- length(gammas)
  seG <- suppressWarnings(sqrt(diag(Vcov)[1 : ncb]))
  infoG <- cbind("Estimate" = gammas, "Lower" = gammas - mult * seG, "Upper" = gammas + mult * seG)
  
  ncz <- length(object$coefficients$phi)
  phis <- c(object$coefficients$phi, object$coefficients$alpha)
  seP <- suppressWarnings(sqrt(diag(Vcov)[(ncb + 1) : (ncb + ncz + 1)]))
  infoP <- cbind("Estimate" = phis, "Lower" = phis - mult * seP, "Upper" = phis + mult * seP)
  
  result <- list("infoLong" = infoG, "infoSurv" = infoP)
  result$level <- level
  class(result) <- "confint.jmodelMult"
  result
}