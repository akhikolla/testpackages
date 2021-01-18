
#========== Function to print detailed summary of joint model with NMRE ==========#

print.summary.jmodelMult <- function (x, digits = max(4, getOption("digits") - 4), printKnots = FALSE, ...) 
{
  result = x
  cat("\nCall:\n")
  cat(paste(deparse(result$call), sep = "\n", collapse = "\n"), "\n")
  cat("\nData Descriptives:\n")
  percent <- round(100 * sum(result$d) / result$n, 1)
  cat("Longitudinal Process\t\tSurvival Process")
  cat("\nNumber of Observations: ", result$N, "\tNumber of Events: ", 
      sum(result$d), " (", percent, "%)", sep = "")
  cat("\nNumber of Groups:", result$n, "\n")
  model.sum <- data.frame(AIC = result$AIC, BIC = result$BIC, logLik = result$logLik, row.names = "")
  print(model.sum)
  
  cat("\nCoefficients:")
  cat("\nLongitudinal Process: Nonparametric multiplicative random effects model\n")
  printCoefmat(result$infoLong, P.values = TRUE, has.Pvalue = TRUE)
  cat("\nSurvival Process: ")
  if (result$rho == 0)
    cat("Proportional hazards model with unspecified baseline hazard function\n")
  else if (result$rho == 1)
    cat("Proportional odds model with unspecified baseline hazard function\n")
  else
    cat("Transformation model with rho=", result$rho, "and unspecified baseline hazard function\n")
  printCoefmat(result$infoSurv, P.values = TRUE, has.Pvalue = TRUE)
  
  cat("\nVariance Components:\n")
  SD <- c(result$sigma.b, result$sigma.e)
  Mat <- data.frame("StdDev" = SD, row.names = c("Random", "Residual"), 
                    check.rows = FALSE, check.names = FALSE)
  print(if(!is.numeric(Mat)) Mat else round(Mat, digits))
  
  cat("\nIntegration: (Adaptive Gauss-Hermite Quadrature)")
  cat("\nquadrature points:", result$control$nknot, "\n")
  
  cat("\nStdErr Estimation:")
  if (result$control$SE.method == "PLFD")
    cat("\nmethod: profile likelihood with forward difference\n")
  else if (result$control$SE.method == "PFDS")
    cat("\nmethod: profile Fisher score with forward difference\n")
  else if (result$control$SE.method == "PRES")
    cat("\nmethod: profile Fisher score with Richardson extrapolation\n")
  else
    cat("\nmethod: NA\n")
  
  cat("\nOptimization:")
  cat("\nconvergence:", result$convergence)
  cat("\niterations:", result$numIter, "\n")
  
  invisible(result)
}
