
summary.altf4 <- function(object, ...)
  {
   x <- object
   cat("Mean coefficients: ")
   cat("\n")
   c <- colMeans(x$coeff.[[1]],na.rm=TRUE)
   c <- round(c,digits=4)
   print(c,quote=FALSE)
   cat("\n")
   cat("Forecast quality measures: ")
   cat("\n")
   print(x$summary)
   cat("\n")
  }
