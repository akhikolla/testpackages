#' Combine multiple lattice plots layerwise
#' 
#' @description
#' This function combines multiple \strong{lattice} plot objects drawing each as 
#' a layer on top of the previous plots. Note that the global plot settings 
#' (e.g. xlim, ylim, ...) are taken from the first object. This is particularly 
#' useful when looping over large amounts of data using \code{\link{lapply}} 
#' (see examples).
#' 
#' @param trellis.list A \code{list} containing \strong{lattice} plot objects.
#' @param ... Additional arguments passed to \code{\link{as.layer}}.
#' 
#' @return
#' A single \strong{lattice} plot object.
#' 
#' @author
#' Tim Appelhans
#' 
#' @seealso
#' \code{\link{as.layer}}.
#' 
#' @examples
#' library(latticeExtra)
#' dat <- list(1:10,
#'             10:1,
#'             3:7,
#'             7:3)
#' 
#' plist <- lapply(seq(dat), function(i) {
#'  tmp <- xyplot(dat[[i]] ~ seq(dat[[i]]),
#'                type = "l", col = i)
#' })
#' 
#' p <- latticeCombineLayer(plist)
#' 
#' print(p)
#' 
#' @export latticeCombineLayer
#' @aliases latticeCombineLayer

latticeCombineLayer <- function(trellis.list, ...) {
  
  outLayer <- function(x, y, ...) {
    x + latticeExtra::as.layer(y, ...)
  }
  
  out <- Reduce(outLayer, trellis.list, ...)
  return(out)
}