
#========== Function to compute summaries of objects in the joint model class ==========#

summary.jmodelTM <- function (object, ...) 
{
  fit = object
  betas <- fit$coefficients$beta
  Vcov <- fit$Vcov
  ncx <- length(betas)
  seB <- suppressWarnings(sqrt(diag(Vcov)[1 : ncx]))
  infoB <- cbind("Estimate" = betas, "StdErr" = seB, "z.value" = betas / seB, 
                  "p.value" = 2 * pnorm(abs(betas / seB), lower.tail = FALSE))
  
  ncw <- length(fit$coefficients$phi)
  phis <- c(fit$coefficients$phi, fit$coefficients$alpha)
  seP <- suppressWarnings(sqrt(diag(Vcov)[(ncx + 1) : (ncx + ncw + 1)]))
  infoP <- cbind("Estimate" = phis, "StdErr" = seP, "z.value" = phis / seP,
                  "p.value" = 2 * pnorm(abs(phis / seP), lower.tail = FALSE))
  
  Bsigma <- fit$coefficients$Bsigma
  ncz <- if(length(Bsigma) == 1) 1 else nrow(Bsigma)
  npar <- ncx + ncw + ncz * (ncz + 1) / 2 + 2
  
  logLik <- fit$logLik
  AIC <- - 2 * logLik + 2 * npar
  BIC <- - 2 * logLik + log(fit$n) * npar
  
  result <- list("infoLong" = infoB, "infoSurv" = infoP, Bsigma = Bsigma,
                 sigma.e = fit$coefficients$Ysigma, logLik = logLik, AIC = AIC, BIC = BIC)
  result$N <- fit$N
  result$n <- fit$n
  result$d <- fit$d
  result$rho <- fit$rho
  result$control <- fit$control
  result$convergence <- fit$convergence
  result$numIter <- fit$numIter
  result$call <- fit$call
  class(result) <- "summary.jmodelTM"
  result
}
