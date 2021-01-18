#'@title Printing DiscreteFDR results
#'
#'@description
#'Prints the results of discrete FDR analysis, stored in a \code{DiscreteFDR} 
#'S3 class object.
#'
#'@param x          an object of class "\code{DiscreteFDR}".
#'@param ...        further arguments to be passed to or from other methods.
#'                  They are ignored in this function.
#'
#'@return
#'The respective input object is invisibly returned via \code{invisible(x)}. 
#'
#'@template example
#'@examples
#'
#'DBH.su.crit <- DBH(raw.pvalues, pCDFlist, direction = "su", ret.crit.consts = TRUE)
#'print(DBH.su.crit)
#'
#'@importFrom stats p.adjust
#'@method print DiscreteFDR
#'@export
## S3 method for class 'DiscreteFDR'
print.DiscreteFDR <- function(x, ...){
  m <- length(x$Data$raw.pvalues)
  k <- x$Num.rejected
  BH <- p.adjust(x$Data$raw.pvalues, "BH")
  
  # print title (i.e. algorithm)
  cat("\n")
  cat("\t", x$Method, "\n")
  # print dataset name(s)
  cat("\n")
  cat("data: ", x$Data$data.name, "\n")
  # print short results overview
  cat("number of tests =", m, "\n")
  cat("number of rejections =", k, "at global FDR level", x$Signif.level, "\n")
  cat("(Original BH rejections = ", sum(BH <= x$Signif.level), ")\n", sep = "")
  if(k) cat("largest rejected p value: ", max(x$Rejected), "\n")
  
  cat("\n")
  invisible(x)
}