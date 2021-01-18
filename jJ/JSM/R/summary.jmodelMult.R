
#========== Function to compute summaries of objects in the joint model class with NMRE ==========#

summary.jmodelMult <- function (object, ...) 
{ 
  fit = object
  gammas <- fit$coefficients$gamma
  Vcov <- fit$Vcov
  ncb <- length(gammas)
  seG <- suppressWarnings(sqrt(diag(Vcov)[1 : ncb]))
  infoG <- cbind("Estimate" = gammas, "StdErr" = seG, "z.value" = gammas / seG, 
                 "p.value" = 2 * pnorm(abs(gammas / seG), lower.tail = FALSE))
  
  ncz <- length(fit$coefficients$phi)
  phis <- c(fit$coefficients$phi, fit$coefficients$alpha)
  seP <- suppressWarnings(sqrt(diag(Vcov)[(ncb + 1) : (ncb + ncz + 1)]))
  infoP <- cbind("Estimate" = phis, "StdErr" = seP, "z.value" = phis / seP,
                 "p.value" = 2 * pnorm(abs(phis / seP), lower.tail = FALSE))
  
  npar <- ncb + ncz + 3
  
  logLik <- fit$logLik
  AIC <- - 2 * logLik + 2 * npar
  BIC <- - 2 * logLik + log(fit$n) * npar
  
  result <- list("infoLong" = infoG, "infoSurv" = infoP, sigma.b = fit$coefficients$Bsigma,
                 sigma.e = fit$coefficients$Ysigma, logLik = logLik, AIC = AIC, BIC = BIC)
  result$N <- fit$N
  result$n <- fit$n
  result$d <- fit$d
  result$rho <- fit$rho
  result$control <- fit$control
  result$convergence <- fit$convergence
  result$numIter <- fit$numIter
  result$call <- fit$call
  class(result) <- "summary.jmodelMult"
  result
}
