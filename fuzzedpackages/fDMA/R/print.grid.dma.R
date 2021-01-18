
print.grid.dma <- function(x, ...)
  {
   cat("RMSE: ")
   cat("\n")
   print(round(x$RMSE,digits=4),quote=FALSE)
   cat("\n")
   cat("\n")
   cat("MAE: ")
   cat("\n")
   print(round(x$MAE,digits=4),quote=FALSE)
   cat("\n")
   cat("* alphas by columns, lambdas by rows")
   cat("\n")
  }
