
print.altf2 <- function(x, ...)
  {
   cat("Forecast quality measures: ")
   cat("\n")
   print(x$summary,quote=FALSE)
   cat("\n")
  }
