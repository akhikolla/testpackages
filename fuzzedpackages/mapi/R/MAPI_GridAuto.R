#' @title Function MAPI_GridAuto
#' @export
#' @description Wrapper that computes cell halfwidth for a given beta value, and then 
#'      builds a grid of hexagonal cells (call to \code{\link{MAPI_GridHexagonal}}). 
#'
#' @param samples a data.frame with names and geographical coordinates of samples. Column names must be: 'ind', 'x', 'y'.  
#'   Optional column 'errRad' with an error radius for sample locations (eg. GPS uncertainty). 
#'   Coordinates must be projected (not latitude/longitude).
#' @param crs coordinate reference system: integer with the EPSG code, or character with proj4string. 
#'   When using dummy coordinates (eg. simulation output) you may use EPSG:3857 (pseudo-Mercator) for example. 
#'   This allows computation but, of course, has no geographical meaning.
#' @param beta A value depending on sampling regularity: 0.5 for regular sampling, 0.25 for random sampling (Hengl, 2006).
#' @param buf optional. This parameter allows to expand or shrink the grid by a number of units in the 
#'   same reference system as the sample geographical coordinates (0 by default).
#'   
#' @details 
#' The halfwidth cell value used to build the grid is computed as 
#'   \eqn{h_w = \frac{\beta \sqrt{A/N}}{\sqrt{2.5980}}}{hw = p / sqrt(2.5980) with p = beta * sqrt(A/N)}, 
#'   where A is the study area (convex hull of sampling points) and N the number of samples. 
#' Parameter beta allows to respect the Nyquist-Shannon sampling theorem depending on sampling regularity 
#'   (call to \code{\link{MAPI_EstimateHalfwidth}}).
#' 
#' @return a spatial object of class 'sf' including the x and y coordinates of cell centers, cell geometry (polygons) and cell id (gid).
#' 
#' @examples
#' data("samples")
#' grid <- MAPI_GridAuto(samples, crs=3857, beta=0.5)
#'

MAPI_GridAuto <- function(samples, crs, beta=0.25, buf=0) {
	hw <- MAPI_EstimateHalfwidth(samples, crs=crs, beta=beta)
	grid <- MAPI_GridHexagonal(samples, crs=crs, hw=hw, buf=buf)
	return(grid)
}
