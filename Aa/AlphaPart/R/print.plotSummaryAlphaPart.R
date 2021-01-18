#' print.plotSummaryAlphaPart.R
#'
#' #' Print method for objects of the class \code{plotSummaryAlphaPart} (result of
#' \code{plot(summary(AlphaPart(...)))}).
#'
#' TODO: DETAILS
#'
#' @seealso
#' \code{\link[AlphaPart]{plot.summaryAlphaPart}}
#'
#' @param x plotSummaryAlphaPart, output object from
#' \code{\link[AlphaPart]{plot.summaryAlphaPart}} function
#' @param ask Logical, ask before printing another plot?
#' @param ...  Arguments passed to other functions (not used at the moment).
#'
#' @example inst/examples/examples_print.plotSummaryAlphaPart.R
#'
#' @useDynLib AlphaPart, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @export


print.plotSummaryAlphaPart <- function (x, ask=interactive(), ...) {

  k <- 1
  nT <- length(x)
  for (i in 1:nT) {
    print(x[[i]])
    if (ask) {
      if (k < nT) {
        msg <- paste("Press any key to print out the next plot (", k, "/", nT, ") ...\n", sep="")
        tmp <- readline(prompt=msg)
      }
    }
    k <- k + 1
  }
}




