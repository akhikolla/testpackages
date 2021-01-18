summary.qtlRes <- function(object, ...){
  
  cat("QTL Summary\n")
  cat("---------------\n")
  cat("Type of test        :",object$method,"\n")
  cat("Tested phenotypes   :",length(object$qtl),"\n")
  invisible(object)
} 
