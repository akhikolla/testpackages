#' @title Function MAPI_GridHexagonal
#' @export
#' 
#' @description Build a grid of hexagonal cells according to samples coordinates and a given halfwidth 
#' cell value provided by users (can be computed using \code{\link{MAPI_EstimateHalfwidth}}).
#' 
#' @param samples a data.frame with names and geographical coordinates of samples. Column names must be: 'ind', 'x', 'y'.  
#'   Optional column 'errRad' with an error radius for sample locations (eg. GPS uncertainty). 
#'   Coordinates must be projected (not latitude/longitude).
#' @param hw Halfwidth : side length of hexagonal cells.
#' @param crs coordinate reference system: integer with the EPSG code, or character with proj4string. 
#'   When using dummy coordinates (eg. simulation output) you may use EPSG:3857 for example. 
#'   This allows computation but, of course, has no geographical meaning.
#' @param buf optional. This parameter allows to expand or shrink the grid by a number of units in 
#'   the same reference system as the sample geographical coordinates (0 by default).
#' @param shift boolean, use default FALSE value. DEPRECATED. 
#'   This parameter has been added for the sake of compatibility with previous PostgreSQL MAPI extension.
#' 
#' @return a spatial object of class 'sf' including the x and y coordinates of cell centers, cell geometry (polygons) and cell id (gid).
#' 
#' @examples
#' data("samples")
#' # Builds a grid of hexagonal cells according to samples coordinates (columns x and y) 
#' # using the EPSG:3857 projection and an halfwidth cell value of hw=250m.
#' grid <- MAPI_GridHexagonal(samples, crs=3857, hw=250)
#' 

MAPI_GridHexagonal <- function(samples, crs, hw, buf=0, shift=FALSE) {
  message("Building grid...")
  
  # Convert samples into spatial sf object
  samples2 <- sf::st_as_sf(samples, coords=c("x", "y"), crs=crs)
  
  # Computes study area (possibly modified by optionnal buffer)
  study.area <- sf::st_buffer(sf::st_convex_hull(sf::st_union(samples2$geometry)), buf)
  
  # Get the rectangular bounding box around study area
  dim <- sf::st_bbox(study.area)
  
  # Initialize variables for y
  if (shift) {
	y <- dim["ymin"] - hw
  } else {
	y <- dim["ymin"]
  }
  y_count <- 1
  
  # Initialize other variables
  rectGrid <- data.frame(x=numeric(),y=numeric())
  i <- 1
  vecteur <- list()
  
  # Start iteration on y axis
  while(y <= dim["ymax"]) {
    
    # Initialize x
    if (shift) {
		x <- dim["xmin"] - hw
	} else {
		x <- dim["xmin"]
	}
    x_count <- 1 
    
    while(x <= dim["xmax"]){
      
      # Shift y dependig wether row is odd or even
      y_offset <- (hw * sqrt(3)/2.0) * (0.5 - (x_count %% 2)) 
      y_adj <- y + y_offset
      
      centre <- c(x, y_adj)
      rectGrid[i,] <- centre
      
      # Computes hexagon summits
      x1 <- x - (0.5*hw)
      y1 <- y_adj + (sqrt(3)*0.5*hw)
      x2 <- x - hw
      y2 <- y_adj
      x3 <- x - (0.5*hw)
      y3 <- y_adj - (sqrt(3)*0.5*hw)
      x4 <- x + (0.5*hw)
      y4 <- y_adj - (sqrt(3)*0.5*hw)
      x5 <- x + hw 
      y5 <- y_adj
      x6 <- x + (0.5*hw)
      y6 <- y_adj + (sqrt(3)*0.5*hw)
      x7 <- x1
      y7 <- y1
      x_hexa <- c(x1,x2,x3,x4,x5,x6,x7)
      y_hexa <- c(y1,y2,y3,y4,y5,y6,y7)
      
      # Convert to sf polygon
      polygon <- sf::st_polygon(list(cbind(x_hexa,y_hexa)))
      vecteur[[i]] <- polygon
      
      # Increments x
      x <- x + hw * 3 / 2
      x_count <- x_count + 1
      i <- i+1
    }
    
    # Increments y
    y <- y + hw*sqrt(3)
    y_count <- y_count + 1
  }
  
  # Builds output table
  sf::st_geometry(rectGrid) <- sf::st_sfc( vecteur, crs=crs)
  rectGrid$gid <- as.integer(rownames(rectGrid))
  
  # Keep only cells intersecting with study area
  grid <- rectGrid[ c(sf::st_intersects(x=study.area, y=rectGrid)[[1]]) , ]
  
  # Work done.
  message(paste("...", nrow(grid), "cells"))
  return(grid)
}
