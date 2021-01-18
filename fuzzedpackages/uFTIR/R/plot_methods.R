#' @include classes.R
NULL

setClassUnion("plot.uFTIR", c("Tile", "SpectralPack", "SAM", "Smooth", "clipper"))

#' Plotting Objects
#' 
#' @description Plotting method for objects of class \code{\link[=Tile-class]{Tile}}, \code{\link[=SpectralPack-class]{SpectralPack}}, \code{\link[=SAM-class]{SAM}}, \code{\link[=Smooth-class]{Smooth}}, and \code{\link{clipper}}.
#' 
#' In objects of class SAM, always the clusters slot gets plotted.
#' 
#' @param x Tile, SpectralPack, SAM, Smooth, or clipper to plot. 
#' @param y Missing.
#' @param slice For objects of class SAM or Smooth, Which slice of the cube should be plotted?
#' @param FUN For objects of class SpectralPack and Tile, Which function should be used to collapse the cube to a matrix?
#' @param match_uFTIR The Agilent Microscope transposes and inverts the image (it inverts rows only -cols after transposition). Do you want the plot function to match this behaviour? Default FALSE.
#' @param ... Further arguments to \code{\link[raster]{plot}}
#'
#' @details 
#' The function requires the raster package to plot. It coerses the matrix to raster to do so.
#' 
#' @return
#' NULL
#' @export
#' @rdname plot_tile
#' @examples
#' # Tile objects:
#' x <- tile_read(base::system.file("extdata/tile.bsp", package = "uFTIR"))
#' plot(x)
#' # with arguments for raster::plot
#' plot(x, axes = FALSE, box = FALSE, legend = FALSE)
#'
#' # SpectralPack objects:
#' x <- tile_base_corr(x)
#' x <- wavealign(x, primpke)
#' plot(x) 
#' 
#' # SAM objects:
#' x <- tile_sam(x)
#' plot(x)
#' 
#' # Smooth objects
#' x <- smooth_sam(x, as.integer(length(primpke@clusternames)), window = 3, 1)
#' plot(x)
#' 
#' # clipper objects:
#' clip <- toClip(8,20,c(10,10))
#' polygon(clip@xycoords)
#' x <- clipper(x, clip@centre, clip@rad, 1)
#' plot(x)
setMethod("plot", c(x = "plot.uFTIR", y = "missing"),
          function(x, y, ...) plot_tile(x, ...))


#' @rdname plot_tile
#' @export
plot_tile <- function(x, slice = 1, FUN = sum, match_uFTIR = FALSE, ...){
  if("Tile" %in% class(x)){
    m <- x@Spectra
    m <- apply(m, c(1,2), FUN)
    #m <- apply(m, 1, rev)
  } else if("SpectralPack" %in% class(x)){
    m <- x@Readings@Spectra
    m <- apply(m, c(1,2), FUN)
    #m <- apply(m, 1, rev)
  } else if("SAM" %in% class(x)){
    m <- x@clusters[,,slice]
  } else if("Smooth" %in% class(x)){
    m <- x@smooth[,,slice]
  } else if("clipper" %in% class(x)){
    m <- x
    class(m) <- c("matrix")
  } else {
    warning(
      "Clase no reconocida\nMetodo implementado para Tile, SpectralPack, SAM y Smoothy"
    )
    return(NULL)
  }
  if(match_uFTIR){
    m <- apply(m, 1, rev)
  }
  m <- raster(m, xmn = 0, xmx = ncol(m), ymx = nrow(m), ymn = 0)
  raster::plot(m, ...)
  NULL
}

#' Highlight a selected substance
#' 
#' The function highlights a selected substance(s) -or cluster(s)- in a plot. It can add to an existing plot or create a new one. Currently implemented only for S3 objects returned by \code{\link{clipper}}.
#'
#' @param x "clipper" object to plot
#' @param dst_cluster cluster (or substance) to highlight. Numeric and character inputs are supported. If is char, then you should provide clusternames to coerse it to numeric.
#' @param new Do you want a new plot or the function should add the polygons to the current plot?
#' @param clusternames char vector given the list of analyzed substance clusters. Only relevant when dst_cluster is a character.
#' @param polygon_color color of the overploted polygons. Red by default.
#' @param match_uFTIR The Agilent Microscope transposes and inverts the image (it inverts rows only -cols after transposition). Do you want the plot function to match this behaviour? Default FALSE.
#' @param ... other arguments for \code{\link{plot_tile}}
#'
#' @return
#' NULL
#' @export
#' @examples
#' x <- tile_read(base::system.file("extdata/tile.bsp", package = "uFTIR"))
#' x <- tile_base_corr(x)
#' x <- wavealign(x, primpke)
#' x <- tile_sam(x)
#' x <- smooth_sam(x, as.integer(length(primpke@clusternames)), window = 3, 1)
#' clip <- toClip(8,20,c(10,10))
#' plot(x)
#' polygon(clip@xycoords)
#' x <- clipper(x, clip@centre, clip@rad, 1)
#' 
#' highlight_substance(x, 5)
highlight_substance <- function(x, dst_cluster, 
                                new = TRUE, clusternames = NULL, 
                                polygon_color = "red", 
                                match_uFTIR = FALSE, ...){
  
  if(!"clipper" %in% class(x)){
    warning("x should be a clipper object [returned by clipper function]")
    return(NULL)
  }
  
  if(is.character(dst_cluster)){
    if(!is.null(clusternames)){
      dst_cluster <- grep(dst_cluster, clusternames)
    } else {
      warning("To conver to dst_cluster from char to numeric I need the clusternames")
      return(NULL)
    }
  }
  
  hsub <- x
  hsub <- matrix(x %in% dst_cluster, nrow = nrow(hsub), ncol = ncol(hsub))
  hsub[is.na(x)] <- NA
  class(hsub) <- c("matrix")
  
  if(match_uFTIR){
    hsub <- apply(hsub, 1, rev)
    match_uFTIR <- FALSE
  }
  
  if(new){
    plot_tile(x, match_uFTIR = FALSE, ...)
  }
  raster::contour(raster::raster(hsub, 
                                 xmn = 0, xmx = ncol(hsub), ymn=0, ymx=nrow(hsub)),
                  add = TRUE,
                  col = polygon_color)
  NULL
}
