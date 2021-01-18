
#========== Function to print detailed summary of joint model ==========#

print.summary.jmodelTM <- function (x, digits = max(4, getOption("digits") - 4), printKnots = FALSE, ...) 
{
  result <- x
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
  cat("\nLongitudinal Process: Linear mixed-effects model\n")
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
  Bsigma <- result$Bsigma
  ncz <- if(length(Bsigma) == 1) 1 else nrow(Bsigma)
  SD <- if (ncz == 1) sqrt(Bsigma) else sqrt(diag(Bsigma))
  SD <- c(SD, result$sigma.e)
  if (ncz > 1) {
    corr <- cov2cor(Bsigma)
    corr[upper.tri(corr, TRUE)] <- 0
    corr <- rbind(corr, rep(0, ncz))
    tempMat <- round(cbind(SD, corr[ , - ncz]), digits)
    tempMat <- apply(tempMat, 2, sprintf, fmt = "% .4f")
    tempMat[tempMat == tempMat[1, 2]] <- ""
    tempMat[1, -1] <- abbreviate(colnames(tempMat)[- 1], 6)
    colnames(tempMat) <- c(colnames(tempMat)[1], rep("", ncz - 1))
    Mat <- data.frame(tempMat, check.rows = FALSE, check.names = FALSE)
    colnames(Mat) <- c("StdDev", "Corr", if (ncz > 2) rep(" ", ncz - 2) else NULL)
    rownames(Mat) <- c(dimnames(Bsigma)[[1]], "Residual")
  }else {
    Mat <- data.frame("StdDev" = SD, row.names = c("Random", "Residual"), 
                      check.rows = FALSE, check.names = FALSE)
  }
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
