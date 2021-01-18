#' print.summaryAlphaPart.R
#'
#' Print method for objects of the class summaryAlphaPart.
#'
#' Print method for objects of the class \code{summaryAlphaPart} (result of
#' \code{summary(AlphaPart(...))}).
#'
#' @seealso
#' \code{\link[AlphaPart]{summary.AlphaPart}}
#'
#' @param x summaryAlphaPart, output object from \code{\link[AlphaPart]{summary.AlphaPart}} function.
#' @param ...  Arguments passed to other functions (not used at the moment).
#'
#' @example inst/examples/examples_summary.AlphaPart.R
#'
#' @useDynLib AlphaPart, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @export


print.summaryAlphaPart <- function(x, ...) {
  

  nI <- nrow(x[[1]])
  cat("\n\n Summary of partitions of breeding values \n")
  cat("   - paths: ",  x$info$nP, " (", paste(x$info$lP, collapse=", "), ")\n", sep="")
  cat("   - traits: ", x$info$nT, " (", paste(x$info$lT, collapse=", "), ")\n", sep="")
  if (length(x$info$warn) > 0) cat("   - warning: ", paste(x$info$warn, collapse="\n"), "\n", sep="")

  for (trt in x$info$lT) {
    cat("\n Trait:", trt, "\n\n")
    print(x[[trt]])
  }
  cat("\n")

}
  

