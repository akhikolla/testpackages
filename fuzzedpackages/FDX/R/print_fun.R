#'@title Printing FDX results
#'
#'@description
#'Prints the results of discrete FDX analysis, stored in a \code{FDX}
#'S3 class object.
#'
#'@return
#'The respective input object is invisibly returned via \code{invisible(x)}. 
#'
#'@param x          an object of class "\code{FDX}".
#'@param ...        further arguments to be passed to or from other methods.
#'                  They are ignored in this function.
#'
#'@template example
#'@examples
#'
#'DPB.crit <- DPB(raw.pvalues, pCDFlist, critical.values = TRUE)
#'print(DPB.crit)
#'
#'@method print FDX
#'@importFrom stats p.adjust
#'@export
## S3 method for class 'FDX'
print.FDX <- function(x, ...){
  m <- length(x$Data$raw.pvalues)
  k <- x$Num.rejected
  if(grepl("Lehmann", x$Method)){
    n <- continuous.LR(x$Data$raw.pvalues, x$FDP.threshold, x$Exceedance.probability, TRUE, FALSE)$Num.rejected
    orig <- "Lehmann-Romano"
  }
  else{
    n <- continuous.GR(x$Data$raw.pvalues, x$FDP.threshold, x$Exceedance.probability, TRUE, FALSE)$Num.rejected
    orig <- "Guo-Romano"
  }
  
  # print title (i.e. algorithm)
  cat("\n")
  cat("\t", x$Method, "\n")
  
  # print dataset name(s)
  cat("\n")
  cat("Data: ", x$Data$data.name, "\n")
  
  # print short results overview
  cat("Number of tests =", m, "\n")
  cat("Number of rejections = ", k, " when controlling FDP at level ", x$FDP.threshold, " with probability ",
      x$Exceedance.probability, ",\n", paste(rep(" ", 24 + nchar(as.character(k))), collapse = ""),
      "i.e. P(FDP > ", x$FDP.threshold, ") <= ", x$Exceedance.probability, "\n", sep = "")
  
  if(!grepl("Continuous", x$Method))
    cat("Original", orig, "rejections =", n, "\n")
  
  cat("Original Benjamini-Hochberg rejections =", sum(p.adjust(x$Data$raw.pvalues, "BH") <= x$FDP.threshold),
      "at level", x$FDP.threshold, "\n")
  
  if(k){
    if(!grepl("Weighted", x$Method))
      cat("Largest rejected p value: ", max(x$Rejected), "\n")
    else
      cat("Largest rejected weighted p value: ", max(x$Weighted[x$Indices]), "\n")
  }
  
  cat("\n")
  invisible(x)
}