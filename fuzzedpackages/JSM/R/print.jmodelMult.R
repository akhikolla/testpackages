
#========== Function to print basic summary of joint model with NMRE ==========#

print.jmodelMult <- function (x, digits = max(4, getOption("digits") - 4), ...) 
{
  result = x
  cat("Joint modeling fitting of longitudinal and survival data by ML\n")
  cat("Nonparametric multiplicative random effects are used\n")
  
  cat("\nCoefficients:\n")
  cat("\tLongitudinal Process\n")
  print(round(result$coefficients$gamma, digits))
  cat("\tSurvival Process\n")
  phis <- c(result$coefficients$phi, result$coefficients$alpha)
  print(round(phis, digits))
  
  cat("\nVariance Components:\n")
  SD <- c(result$coefficients$Bsigma, result$coefficients$Ysigma)
  Mat <- data.frame("StdDev" = SD, row.names = c("Random", "Residual"), 
                    check.rows = FALSE, check.names = FALSE)
  print(if(!is.numeric(Mat)) Mat else round(Mat, digits))
  
  cat("\nLog-likelihood:", result$logLik, "\n")
  cat("Number of Observations:", result$N, "\n")
  cat("Number of Subjects:", result$n, "\n")
  invisible(result)
}
