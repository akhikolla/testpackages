#' Combine multiple lattice plots in a facetted grid (panels)
#' 
#' @description
#' This function combines multiple \strong{lattice} plot objects in a facetted
#' grid. Note that the global plot settings (e.g. xlim, ylim, ...) are taken 
#' from the first object though the user can specify whether \code{scales} 
#' should be identical or not. This is particularly useful when looping over 
#' large amounts of data using \code{\link{lapply}} (see examples).
#' 
#' @param trellis.list A \code{list} containing \strong{lattice} plot objects.
#' @param between Space between panels.
#' @param as.table If TRUE (default) drawing is top left to bottom right
#' @param ... Additional arguments passed to \code{\link{c.trellis}}.
#' 
#' @return
#' A single \strong{lattice} plot object.
#' 
#' @author
#' Tim Appelhans
#' 
#' @seealso
#' \code{\link{c.trellis}}.
#' 
#' @examples
#' #load data
#' #Use a probability map assuming high potential for city expansion is just 
#' #resulting from proximity to current urban area:
#' pred <- raster(system.file("extdata/probability.rst", package = "Orcs"))
#' 
#' #observed city growth between 1990 and 2006
#' obs <- raster(system.file("extdata/citygrowth.tif", package = "Orcs"))
#' 
#' #masking current urban area since these pixels have no potential for change
#' mask <- raster(system.file("extdata/citymask.tif", package = "Orcs"))
#' 
#' #create data list
#' dat <- list(pred, obs, mask)
#' 
#' #create list of lattice plots
#' plist <- lapply(seq(dat), function(i) {
#'   spplot(dat[[i]], scales = list(draw = TRUE))
#' })
#' 
#' \donttest{
#' #draw individually
#' plist[[1]]
#' plist[[2]]
#' plist[[3]]
#' }
#' 
#' #combine to grid, using c(1, 3) layout
#' p <- latticeCombineGrid(plist, layout = c(1, 3))
#' print(p)
#' 
#' @export latticeCombineGrid
#' @aliases latticeCombineGrid

latticeCombineGrid <- function(trellis.list,
                               between = list(y = 0.3, x = 0.3),
                               as.table = TRUE,
                               ...) {
  stopifnot(
    requireNamespace("lattice"),
    requireNamespace("latticeExtra")
  )
  
  ## combine plot objects
  outLayout <- function(x, y, ...) {
    update(c(x, y, ...), between = between, as.table = as.table)
  }
  
  out <- suppressWarnings(Reduce(outLayout, trellis.list))
  
  ## apply additional customizations
  dots = list(...)
  if (length(dots) > 0) {
    out <- do.call(update, args = append(list(object = out), dots))
  }
  
  return(out)
}