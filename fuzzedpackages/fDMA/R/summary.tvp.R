
summary.tvp <- function(object, ...)
  {
   x <- object
   cat("Mean coefficients: ")
   c <- colMeans(x$thetas,na.rm=TRUE)
   cat("\n")
   print(round(c,digits=4),quote=FALSE)
   cat("\n")
   e <- (mean((as.vector(x$y)-as.vector(x$y.hat))^2,na.rm=TRUE))^0.5
   e <- round(e,digits=4)
   cat("RMSE: ",e)
   cat("\n")
   e <- mean(abs(as.vector(x$y)-as.vector(x$y.hat)),na.rm=TRUE)
   e <- round(e,digits=4)
   cat("MAE:  ",e)
   cat("\n")
  }
