# Confidence Intervals on jmodelTM Parameters

confint.jmodelTM <- function(object, parm, level = 0.95, ...) {
  if (!inherits(object, "jmodelTM"))
    stop("Can only be used for 'jmodelTM' objects. \n")

  a <- (1 - level) / 2
  mult <- qnorm(1 - a)
  
  betas <- object$coefficients$beta
  Vcov <- object$Vcov
  ncx <- length(betas)
  seB <- suppressWarnings(sqrt(diag(Vcov)[1 : ncx]))
  infoB <- cbind("Estimate" = betas, "Lower" = betas - mult * seB, "Upper" = betas + mult * seB)
  
  ncw <- length(object$coefficients$phi)
  phis <- c(object$coefficients$phi, object$coefficients$alpha)
  seP <- suppressWarnings(sqrt(diag(Vcov)[(ncx + 1) : (ncx + ncw + 1)]))
  infoP <- cbind("Estimate" = phis, "Lower" = phis - mult * seP, "Upper" = phis + mult * seP)
  
  result <- list("infoLong" = infoB, "infoSurv" = infoP)
  result$level <- level
  class(result) <- "confint.jmodelTM"
  result
}