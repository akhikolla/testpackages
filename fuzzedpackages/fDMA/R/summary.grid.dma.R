
summary.grid.dma <- function(object, ...)
  {
   x <- object
   
   cat("RMSE: ")
   cat("\n")
   print(round(x$RMSE,digits=4),quote=FALSE)
   cat("\n")
   cat("Indices of the model minimising RMSE:")
   temp <- which(x$RMSE == min(x$RMSE), arr.ind = TRUE)[1,]
   names(temp) <- c("","")
   print(temp,quote=FALSE)
   cat("\n")
   cat("MAE: ")
   cat("\n")
   print(round(x$MAE,digits=4))
   cat("\n")
   cat("Indices of the model minimising MAE:")
   temp <- which(x$MAE == min(x$MAE), arr.ind = TRUE)[1,]
   names(temp) <- c("","")
   print(temp,quote=FALSE)   
   cat("\n")
   cat("* alphas by columns, lambdas by rows")
   cat("\n")
  }
