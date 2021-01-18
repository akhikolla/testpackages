#' pedSetBase.R
#'
#' A function to set the base population in the pedigree.
#'
#' @details
#' Base population in the pedigree is set by removing rows for some individuals, while their
#' presence as parents is also removed.
#'
#' Arguments \code{down} and \code{na.rm} allow for repeated use of this function, i.e., with
#' \code{down=FALSE} and with \code{down=TRUE} (both in combination with \code{na.rm=TRUE}) in order to
#' propagate information over the pedigree until "convergence".
#'
#' This function can be very slow on large pedigrees with extensive missingness of birth years.
#'
#' @seealso
#' \code{\link[pedigree]{orderPed}} in \pkg{pedigree} package
#'
#' @param x data.frame , with (at least) the following columns: individual, father, and mother identification,
#' and year of birth; see arguments \code{colId},
#' \code{colFid}, \code{colMid}, and \code{colBY}
#' @param keep Logical, indicator that defines which individuals should stay in the
#' the pedigree; see details.
#' @param unknown Value used to represent unknown/missing identification
#' @param report Logical, report success.
#' @param colId Numeric or character, position or name of a column holding individual identification.
#' @param colFid Numeric or character, position or name of a column holding father identification.
#' @param colMid Numeric or character, position or name of a column holding mother identification.
#'
#' @example inst/examples/examples_pedSetBase.R
#'
#' @return Object \code{x} with removed rows for some individuals and their presence as parents.
#' If \code{report=TRUE} progress is printed on the screen.
#'
#' @useDynLib AlphaPart, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @export
#' @importFrom gdata unknownToNA
#' @importFrom gdata isUnknown


pedSetBase <- function (x, keep=NULL, unknown=NA, report=TRUE, colId=1, colFid=2,  colMid=3) {

  ## --- Setup ---

  if (is.null(keep)) stop

  test <- (length(colId) > 1 | length(colFid) > 1 | length(colMid) > 1)
  if (test) stop("arguments 'colId', 'colFid', and 'colMid' must be of length 1")

  n <- nrow(x)
  test <- (length(keep) != n)
  if (test) stop("the length of argument 'keep' must be the same as the number of entries (rows) in x")

  ## --- Action ---

  if (report) {
    cat("All individuals:", n, "\n")
    cat("Removing: ", sum(!keep, na.rm=TRUE), ", ",
                      round(sum(!keep, na.rm=TRUE) / length(keep) * 100), " %\n", sep="")
  }
  rem <- x[!keep, colId]
  ret <- x[keep, ]
  ret[ret[, colFid] %in% rem & !isUnknown(x=ret[, colFid], unknown=unknown), colFid] <- unknown
  ret[ret[, colMid] %in% rem & !isUnknown(x=ret[, colMid], unknown=unknown), colMid] <- unknown
  if (report) cat("Kept:", nrow(ret), "\n")

  ## --- Return ---

  ret


}
