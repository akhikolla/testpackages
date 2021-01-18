

print.mixest <- function(x, ...)
  {
   cat("Models: ")
   cat("y ~ ")
   cat(paste(colnames(x$components),collapse=" + "))
   cat("\n")
   cat("\n")
   cat("observations: ")
   cat(length(x$y.hat))
   cat("\n")
   cat("models: ")
   cat(nrow(x$components))
   cat("\n")
   cat("\n")
   cat("parameters: ")
   cat("\n")
   cat("====================")
   p <- x$parameters
   if (p[2]==0) { p[2] <- "mixture" }
   if (p[2]==1) { p[2] <- "averaging" }
   if (p[2]==2) { p[2] <- "selection" }
   if (p[2]==3) { p[2] <- "median probability" }
   p <- as.matrix(p)
   colnames(p) <- ""
   print(p,quote=FALSE)
   cat("\n")
  }

