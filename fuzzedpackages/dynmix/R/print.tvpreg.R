

print.tvpreg <- function(x, ...)
  {
   cat("Models: ")
   cat("y ~ ")
   cat(paste(colnames(x$coef),collapse=" + "))
   cat("\n")
   cat("\n")
   cat("observations: ")
   cat(length(x$y.hat))
   cat("\n")
   cat("\n")
   cat("Regression coefficients from the last period: ")
   cat("\n")
   cat("\n")
   p <- x$coef[nrow(x$coef),,drop=FALSE]
   rownames(p) <- ""
   print(format(round(p,digits=4),scientific=FALSE),quote=FALSE)
   cat("\n")
  }

