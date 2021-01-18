#'@name summary.FDX
#'
#'@title Summarizing Discrete FDX Results
#'
#'@description
#'\code{summary} method for class "\code{FDX}"
#'
#'@param object       an object of class "\code{FDX}".
#'@param x            an object of class "\code{summary.FDX}".
#'@param max          numeric or \code{NULL}, specifying the maximal number of
#'                    \emph{rows} of the p-value table to be printed. By
#'                    default, when \code{NULL}, \code{getOption("max.print")}
#'                    is used.
#'@param ...          further arguments passed to or from other methods.
#'
#'@details
#'\code{summary.FDX} objects include all data of an \code{FDX}
#'object, but also include an additional table which includes the raw p-values,
#'their indices, the respective critical values (if present), the adjusted
#'p-values (if present) and a logical column to indicate rejection. The table
#'is sorted in ascending order by the raw p-values.
#'
#'\code{print.summary.FDX} simply prints the same output as
#'\code{print.FDX}, but also prints the p-value table.
#'
#'@return
#'\code{summary.FDX} computes and returns a list that includes all the
#'data of an input \code{FDX}, plus
#'\item{Table}{a \code{data.frame}, sorted by the raw p-values, that contains
#'             the indices, that raw p-values themselves, their respective
#'             critical values (if present), their adjusted p-values (if
#'             present) and a logical column to indicate rejection.}
#'             
#'\code{print.summary.FDX} returns that object invisibly.
#'
#'@template example
#'@examples
#'
#'DGR.crit <- DGR(raw.pvalues, pCDFlist, critical.values = TRUE)
#'DGR.crit.summary <- summary(DGR.crit)
#'print(DGR.crit.summary)
#'
#'@method summary FDX
#'@export
## S3 method for class 'FDX'
summary.FDX <- function(object, ...){
  # include all data of x (FDX)
  out <- object
  
  # number of tests
  m <- length(object$Data$raw.pvalues)
  # determine order of raw p-values
  o <- order(object$Data$raw.pvalues)
  # ordered indices
  i <- (1:m)[o]
  # sort raw p-values
  y <- object$Data$raw.pvalues[o]
  # determine for each p-value if its corresponding null hypothesis is rejected
  r <- ((1:m) %in% object$Indices)[o]
  
  # create summary table
  out <- c(out, list(Table = data.frame('P.value' = y, 'Index' = i)))
  if(grepl("Weighted", out$Method)) out$Table <- data.frame(out$Table, 'Weights' = out$Data$weights[o], 'Weighted' = out$Weighted[o])
  if(exists('Critical.values', where = object)){
    if(grepl("Weighted", out$Method))
      out$Table <- data.frame(out$Table, 'Critical.value' = object$Critical.values[order(order(out$Weighted))][o])
    else
      out$Table <- data.frame(out$Table, 'Critical.value' = object$Critical.values)
  }
  if(exists('Adjusted', where = object)){
    out$Table <- data.frame(out$Table, 'Adjusted' = object$Adjusted[o])
  }
  out$Table <- data.frame(out$Table, 'Rejected' = r)
  
  # return output object
  class(out) <- "summary.FDX" # basically a 'FDX' object, but with a summary table (just like 'lm' and 'summary.lm' classes)
  return(out)
}

#'@rdname summary.FDX
#'@method print summary.FDX
#'@export
## S3 method for class 'summary.FDX'
print.summary.FDX <- function(x, max = NULL, ...){
  # determine number of tests
  m <- length(x$Data$raw.pvalues)
  
  # print 'FDX' part of the object
  print.FDX(x)
  
  # rows to print: number of rejections + 5 (if not requested otherwise)
  max <- if(!is.null(max)) ncol(x$Table) * max else getOption("max.print")
  
  # print additional summary table
  print(x$Table, max = max, ...)
  
  cat("\n")
  invisible(x)
}