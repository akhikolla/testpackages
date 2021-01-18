#' AlphaPartSubset.R
#'
#' A function to choose the partition paths to keep.
#'
#' Displaying results of partitions for many paths is often confusing.
#' This function helps in selecting only paths of interest.
#' Unspecified paths are removed from the input object \code{x}.
#' Meta information is modified accordingly. Default setting does nothing.
#' @seealso
#' \code{\link[AlphaPart]{AlphaPart}} for the main method,
#' \code{\link[AlphaPart]{summary.AlphaPart}} for summary method that works on output of \code{AlphaPart},
#' \code{\link[AlphaPart]{AlphaPartSum}} for sum method.
#'
#' @param x AlphaPart or summaryAlphaPart, object from the \code{AlphaPart(...)} or \code{summary(AlphaPart(...), ...)} call.
#' @param paths Character, names of paths to be kept.
#'
#' @return An object of class \code{AlphaPart} or \code{summaryAlphaPart} with only some paths.
#' Meta information in slot "info" is modified as well.
#'
#' @example inst/examples/examples_AlphaPartSubset.R
#'


#' @useDynLib AlphaPart, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @export

AlphaPartSubset <- function(x, paths=NULL) {

  ## --- Setup ---

  test1 <- "AlphaPart"        %in% class(x)
  test2 <- "summaryAlphaPart" %in% class(x)
  if (!any(c(test1, test2))) stop("object 'x' must be of a 'AlphaPart' or 'summaryAlphaPart' class")

  ## Do nothing
  if (is.null(paths)) return(x)

  ## Keep only uniquely defined paths
  paths <- unique(paths)

  ## --- Action ---

  ## Create "identity" map for specified paths and call AlphaPartSum to ease the
  ##   work with a rather complex object structure ;)
  nP <- length(paths)
  map <- vector(mode="list", length=nP)
  for (i in 1:nP) {
    map[[i]] <- c(paths[i], paths[i])
  }
  ## Now add non-specified paths in the last mapping so AlphaPartSum will remove them
  map[[i]] <- c(map[[i]], x$info$lP[!(x$info$lP %in% paths)])
  ## Call AlphaPartSum
  AlphaPartSum(x=x, map=map, remove=TRUE, call="AlphaPartPathSubset")


}






