#' AlphaPartSum.R
#'
#' A function to sum partitions of several paths.
#'
#' @details
#' Sometimes partitions of particular paths are very small or we want to sum
#' paths that have some similarity. These actions are easy to achive manually
#' but this functions provides a way to do this consistently with the given
#' object \code{x}.
#'
#' Arguments \code{map} must be a list of vectors of length at least
#' two. Vectors of length one are skipped. The idea is that the first element
#' is the new or existing path into which we add up all the remaining specified
#' paths, say \code{list(c("A", "B"), c("X", "X", "Y"), c("Z", "X"))} would
#' imply A = B, X = X + Y, and Z = X = X + Y. Note that once X is changed its
#' changed value is used in further calculations. Specify different (new) names
#' for new targets if you want to avoid this.
#'
#' Be carefull with \code{remove=TRUE}, which is the default setting, as all
#' partitions defined after the first (target/new) partition in vector in list
#' will be removed, for example with \code{list(c("A", "B"), c("X", "X", "Y"),
#' c("Z", "X"))} partitions B and Y will be removed, while X will not be removed
#' as it is defined as a target/new partition.
#'
#' @seealso
#' \code{\link[AlphaPart]{AlphaPart}} for the main method,
#' \code{\link[AlphaPart]{summary.AlphaPart}} for summary method that works on output of \code{AlphaPart},
#' \code{\link[AlphaPart]{AlphaPartSubset}} for subset/keep method
#'
#' @param x summaryAlphaPart, object from the \code{AlphaPart(...)} or \code{summary(AlphaPart(...), ...)} call.
#' @param map List, a map of summing paths; see details and examples.
#' @param remove Logical, remove original paths or not.
#' @param zeroPath Logical, set called path to zero if it does not exist.
#' @param call character, for internal use with \code{AlphaPartSubset}).
#'
#' @example inst/examples/examples_AlphaPartSum.R
#'
#' @return An object of class \code{AlphaPart} or \code{summaryAlphaPart} with modified partitions.
#' Meta information in slot "info" is modified as well.
#'
#' @useDynLib AlphaPart, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @export

AlphaPartSum <- function (x, map=NULL, remove=TRUE, zeroPath=TRUE, call="AlphaPartSum") {
  ## --- Setup ---
  
  test1 <- "AlphaPart"        %in% class(x)
  test2 <- "summaryAlphaPart" %in% class(x)
  if (!any(c(test1, test2))) stop("object 'x' must be of a 'AlphaPart' or 'summaryAlphaPart' class")
  if (!is.list(map)) stop("object 'map' must be of a 'list' class")

  ## Initial number of columns
  nCOrig <- ifelse(test1, ncol(x[[1]]), ncol(x[[1]]))
  nPOrig <- x$info$nP 
  nCRem <- 0

  ## Any unknown path?
  mapP <- map
  mapT <- c()
  for (i in 1:length(mapP)) {
    ## Targets
    mapT <- c(mapT, mapP[[i]][1])
    ## Components
    mapP[[i]] <- mapP[[i]][2:length(mapP[[i]])]
  }
  mapM <- c()
  mapT <- c("", mapT) ## trick so that code bellow works with i=1
  for (i in 1:length(mapP)) {
    testE <- mapP[[i]] %in% x$info$lP ## path exists in the data?
    if (any(!testE)) {
      testT <- mapP[[i]] %in% mapT[1:i] ## path exists as a target defined up to now (at i)?
      if (any(!testT)) {
        mapM <- c(mapM, mapP[[i]][(!testE & !testT)]) 
      }
    }
  }
  if (length(mapM) > 0) {
    if (!zeroPath) {
      stop(paste("Unexisting path(s): ", paste(mapM, collapse=", ", sep=""), sep=""))
    }
  }

  ## --- Action ---

  ## Loop over traits
  for (t in 1:(length(x) - 1)) { ## t <- 1

    ## Sum up partitions given the map
    if (call == "AlphaPartSum") {
      if (zeroPath) {
        for (i in mapM) {
          x[[t]][, paste(x$info$lT[t], i, sep="_")] <- 0
        }
      }
      for (i in 1:length(map)) { ## i <- 1
        if (length(map[[i]]) > 1) {
          if (test1) { ## x comes from AlphaPart(...)
            if (length(map[[i]][2:length(map[[i]])]) > 1) { ## need this as rowSums() need an array of at least two dimensions
              x[[t]][, paste(x$info$lT[t], map[[i]][1], sep="_")] <- rowSums(x[[t]][, paste(x$info$lT[t], map[[i]][2:length(map[[i]])], sep="_")])
            } else {
              x[[t]][, paste(x$info$lT[t], map[[i]][1], sep="_")] <-         x[[t]][, paste(x$info$lT[t], map[[i]][2],                  sep="_")]
            }
          } else {    ## x comes from summary(AlphaPart(...), ...)
              if (length(map[[i]][2:length(map[[i]])]) > 1) {
                x[[t]][, map[[i]][1]] <- rowSums(x[[t]][, map[[i]][2:length(map[[i]])]])
              } else {
                x[[t]][, map[[i]][1]] <-         x[[t]][, map[[i]][2]]
              }            
            ## for (j in ...)
          } ## if (test1)
        } ## if (length...)
      } ## for (i in 1...)
    } ## if (call != "AlphaPartSum")

    ## Remove original partitions (we do this after we go through the whole map!)
    if (remove) {
      remY <- remN <- c()
      for (i in 1:length(map)) { ## i <- 1
        if (length(map[[i]]) > 1) {
          remN <- c(remN, map[[i]][1])
          remY <- c(remY, map[[i]][2:length(map[[i]])])
        } ## if (length...)
      } ## for (i in 1...)
      remY <- unique(remY)
      remN <- unique(remN)
      remY <- remY[!(remY %in% remN)]
      if (length(remY) > 0) {
        if (test1) { ## x comes from AlphaPart(...)
          for (i in remY) x[[t]][, paste(x$info$lT[t], i, sep="_")]     <- NULL          
        } else {    ## x comes from summary(AlphaPart(...), ...)
          for (i in remY) x[[t]][, i] <- NULL
        } ## if (test1)
      } ## if (length...)
    } ## if (remove)

  } ## for (t in ...)

  ## --- Fix meta info ---

  if (remove) nCRem <- length(remY)
  nC <- ifelse(test1, ncol(x[[1]]), ncol(x[[1]]))
  x$info$nP <- nC - (nCOrig - nPOrig)
  if (test1) {
    x$info$lP <- gsub(pattern=paste(x$info$lT[t], "_", sep=""), replacement="", x=colnames(x[[t]])[(nCOrig - nPOrig + 1):nC], fixed=TRUE)
  } else {
    x$info$lP <- colnames(x[[t]])[(nCOrig - nPOrig + 1):nC]
  }  
  x$info$warn <- c(x$info$warn, paste("Consistency of the overall sum of partitions might not be correct due to the previous use of '", call, "'", sep=""))

  ## --- Return ---

  x

}
