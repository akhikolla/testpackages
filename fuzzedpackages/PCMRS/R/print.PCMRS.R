print.PCMRS <- function(x, ...){
  
  cat("Output of PCMRS estimation:","\n")
  
  cat("---","\n")
  
  cat("Estimates of item parameters delta","\n")
  print(x$delta, ...)
  
  cat("\n")
  
  cat("Estimates of covariance matrix Sigma","\n")
  print(x$Sigma, ...)
  
  cat("\n")
  cat("------------------","\n")
  cat("\n")
  
  cat("Output of simple PCM estimation:","\n")
  
  cat("---","\n")
  
  cat("Estimates of item parameters delta","\n")
  print(x$delta.PCM, ...)
  
  cat("\n")
  
  cat("Estimates of variance of theta","\n")
  print(x$sigma.PCM, ...)
  
  invisible(x)
}