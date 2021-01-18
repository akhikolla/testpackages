#' summary.AlphaPart.R
#'
#' A function to summarize AlphaPart object.
#'
#' Breedng values of individuals are often summarized, either by year of
#' birth or some other classification. Function \code{summary.AlphaPart} provides
#' a way to ease the computation of such summaries on partitions of breeding values.
#'
#' @seealso
#' \code{\link[AlphaPart]{AlphaPart}} for partitioning breeding values,
#' \code{\link[AlphaPart]{plot.summaryAlphaPart}} for plotting output of summary method
#'
#' @param object AlphaPart, output object from \code{\link[AlphaPart]{AlphaPart}} function.
#' @param by Character, the name of a column by which summary function FUN should
#' be applied; if \code{NULL} (default) summary is given for the whole table.
#' @param FUN Function, which function should be used in summary; function should
#' return single value per each level of by.
#' @param labelSum Character, label used for the overall breeding value.
#' @param subset Logical, perform summary only on a subset of \code{object} subsetted by
#' this argument.
#' @param sums Logical, link between \code{\link[AlphaPart]{AlphaPart}} and
#' \code{summary.AlphaPart()} (only for internal use!).
#' @param ...  Arguments passed to other functions (not used at the moment).
#'
#' @example inst/examples//examples_summary.AlphaPart.R
#'
#' @return An object of class \code{summaryAlphaPart}, which is a list of data frames
#' with summary statistics on breeding value partitions. For each trait there
#' a dataframe holds summary for the "whole/original" breeding value and its partitions.
#' In addition another list is added (named \code{info}) with the following
#' components holdinfg meta info:
#'   \item{path}{column name holding path information}
#'   \item{nP}{number of paths}
#'   \item{lP}{path labels}
#'   \item{nT}{number of traits}
#'   \item{lT}{trait labels}
#'   \item{by}{column name of variable by which summary was performed}
#'   \item{warn}{potential warning messages associated with this object}
#'   \item{labelSum}{column name of summary for "whole/original" breeding values}
#'
#' There is a handy plot method (\code{\link[AlphaPart]{plot.summaryAlphaPart}}) for output.
#'
#' @useDynLib AlphaPart, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @export


summary.AlphaPart <- function(object, by=NULL, FUN=mean, labelSum="Sum", subset=NULL, sums=FALSE, ...) {


  ## --- Setup ---

  groupSummary <- sums

  test <- !("AlphaPart" %in% class(object))
  if (test) {
    stop("'object' must be of a AlphaPart class")
  }

  if (groupSummary) by <- object$by

  if (!groupSummary) {
    test <- !is.null(by) && !(by %in% colnames(object[[1]]))
    if (test) {
      stop("argument 'by' must be NULL or one of column names in object")
    }
    test <- do.call(what=FUN, args=list(x=1:3, ...))
    if (length(test) > 1) {
      stop("function FUN must return a single value (scalar)")
    }
  }

  nC <- ncol(object[[1]]) ## number of columns
  nP <- object$info$nP    ## number of paths
  lP <- object$info$lP    ## names  of paths
  nT <- object$info$nT    ## number of traits
  lT <- object$info$lT    ## names  of traits
  ret <- vector(mode="list", length=nT+1)
  names(ret) <- c(lT, "info")

  ## Subset
  if (!is.null(subset)) {
    object[1:nT] <- lapply(object[1:nT], FUN=function(z) z[subset, ])
  }

  ret$info <- list(path=object$info$path, nP=nP, lP=lP, nT=nT, lT=lT, by=by, warn=object$info$warn, labelSum=labelSum)

  ## --- Compute ---

  z <- ifelse (groupSummary, 1, 2)
  for(i in 1:nT) { ## for each trait
     ## Setup
     cols <- c(lT[i], paste(lT[i], lP, sep="_"))
     paths <- cols
     paths[2:length(paths)] <- ret$info$lP
     paths[1] <- labelSum

     ## Summarize
     if (!groupSummary) {
       if (is.null(by)) {
           #pri length ne sme biti na.rm = TRUE
         tmp <- rep(1, times=nrow(object[[i]]))
         tmpM <- aggregate(x=object[[i]][, cols],    by=list(by=tmp),               FUN=FUN,  na.rm=TRUE)
         tmpN <- aggregate(x=object[[i]][, cols[1]], by=list(by=tmp),               FUN=length)
       } else {
         tmpM <- aggregate(x=object[[i]][, cols],    by=list(by=object[[i]][, by]), FUN=FUN, na.rm=TRUE)
         tmpN <- aggregate(x=object[[i]][, cols[1]], by=list(by=object[[i]][, by]), FUN=length)
       }
     } else {
       tmpN <- object$N[, c(1, i+1)]
       tmpM <- object[[i]][, -(1:2)]
       tmpM <- cbind(rowSums(tmpM), tmpM)
       tmpM <- tmpM / tmpN[, 2]
     }
     ## Add nice column names
     colnames(tmpN) <- c(by, "N")
     colnames(tmpM)[z:ncol(tmpM)] <- paths
     ## Combine FUN and number of records
     tmp <- cbind(tmpN, tmpM[, z:ncol(tmpM)])
     ## Store
     ret[[i]] <- tmp
  }

  ## --- Return ---

  class(ret) <- c("summaryAlphaPart", class(ret))
  ret


}