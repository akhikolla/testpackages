#'@title Summarizing Discrete FDR Results
#'
#'@description
#'\code{summary} method for class "\code{DiscreteFDR}"
#'
#'@param object       an object of class "\code{DiscreteFDR}".
#'@param x            an object of class "\code{summary.DiscreteFDR}".
#'@param max          numeric or \code{NULL}, specifying the maximal number of
#'                    \emph{rows} of the p-value table to be printed. By default,
#'                    when \code{NULL}, \code{getOption("max.print")} is used.
#'@param ...          further arguments passed to or from other methods.
#'
#'@details
#'\code{summary.DiscreteFDR} objects include all data of an \code{DiscreteFDR}
#'object, but also include an additional table which includes the raw p-values,
#'their indices, the respective critical values (if present), the adjusted
#'p-values (if present) and a logical column to indicate rejection. The table
#'is sorted in ascending order by the raw p-values.
#'
#'\code{print.summary.DiscreteFDR} simply prints the same output as
#'\code{print.DiscreteFDR}, but also prints the p-value table.
#'
#'@return
#'\code{summary.DiscreteFDR} computes and returns a list that includes all the
#'data of an input \code{DiscreteFDR}, plus
#'\item{Table}{a \code{data.frame}, sorted by the raw p-values, that contains
#'             the indices, that raw p-values themselves, their respective
#'             critical values (if present), their adjusted p-values (if
#'             present) and a logical column to indicate rejection.}
#'
#'@template example
#'@examples
#'
#'DBH.sd.crit <- DBH(raw.pvalues, pCDFlist, direction = "sd", ret.crit.consts = TRUE)
#'summary(DBH.sd.crit)
#'
#'@rdname summary.DiscreteFDR
#'@export
## S3 method for class 'DiscreteFDR'
summary.DiscreteFDR <- function(object, ...){
  # include all data of x (DiscreteFDR)
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
  out <- c(out, list(Table = data.frame('Index' = i, 'P.value' = y)))
  if(exists('Critical.values', where = object)){
    out$Table <- data.frame(out$Table, 'Critical.value' = object$Critical.values)
  }
  if(exists('Adjusted', where = object)){
    out$Table <- data.frame(out$Table, 'Adjusted' = object$Adjusted[o])
  }
  out$Table <- data.frame(out$Table, 'Rejected' = r)
  
  # return output object
  class(out) <- "summary.DiscreteFDR" # basically a 'DiscreteFDR' object, but with a summary table (just like 'lm' and 'summary.lm' classes)
  return(out)
}

#'@rdname summary.DiscreteFDR
#'@export
## S3 method for class 'summary.DiscreteFDR'
print.summary.DiscreteFDR <- function(x, max = NULL, ...){
  # determine number of tests
  m <- length(x$Data$raw.pvalues)
  
  # print 'DiscreteFDR' part of the object
  print.DiscreteFDR(x)
  
  # rows to print: number of rejections + 5 (if not requested otherwise)
  max <- if(!is.null(max)) ncol(x$Table) * max else getOption("max.print")
  
  # print additional summary table
  print(x$Table, max = max, ...)
  
  cat("\n")
  invisible(x)
}