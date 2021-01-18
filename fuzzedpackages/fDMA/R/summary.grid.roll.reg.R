
summary.grid.roll.reg <- function(object, ...)
  {
   x <- object
   
   print(round(x$fq,digits=4),quote=FALSE)
   cat("\n")
   cat("The model minimising RMSE:")
   temp <- which(x$fq[,1] == min(x$fq[,1]), arr.ind = TRUE)[1]
   names(temp) <- c("")
   print(temp,quote=FALSE)
   cat("\n")
   cat("The model minimising MAE:")
   temp <- which(x$fq[,2] == min(x$fq[,2]), arr.ind = TRUE)[1]
   names(temp) <- c("")
   print(temp,quote=FALSE)  
   cat("\n") 
  }
