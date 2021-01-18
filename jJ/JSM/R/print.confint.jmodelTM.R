# Print Detailed Summary of Confidence Interval of jmodelTM

print.confint.jmodelTM <- function (x, digits = max(4, getOption("digits") - 4), printKnots = FALSE, ...) {
  result <- x
  cat(paste("Approximate ", result$level * 100, "%", " confidence intervals", sep = ""), "\n")
  
  cat("\nLongitudinal Process:\n")
  printCoefmat(result$infoLong)
  cat("\nSurvival Process:\n")
  printCoefmat(result$infoSurv)
  
  invisible(result)
}