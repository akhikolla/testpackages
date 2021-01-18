#' @include classes.R
NULL

setClassUnion("asclipper", members = c("SAM", "Smooth"))

#' Visual clipper
#' 
#' @description Visual aid to find the right center and radius for the \code{\link{clipper}} function. It does not clip, but return an usable (plotable) object. 
#' 
#' Since the samples are placed under the uFTIR Microscope by hand, the cropping area to pass to \code{\link{clipper}} is not always the same and (usually) it has to be adjusted. To have a visual aid, you can use this function. It returns an S4 object of class \code{\link[=clipmask-class]{clipmask}} which can be over-plotted on top of a \code{\link[=SAM-class]{SAM}} or \code{\link[=Smooth-class]{Smooth}} object by calling \code{\link[graphics]{polygon}} and using the xycoords slot of the returned object. Using that process you can test manually different center points and radius for the clipping circle.
#'
#' @param rad The circle radius.
#' @param segments How many segments should the resulting polygon have? (because it is not a circle, 20 is ok).
#' @param centre The coordinates of the polygon centre (x,y).
#'
#' @return
#' An object of class \code{\link[=clipmask-class]{clipmask}}.
#' @export
#' 
#' @seealso 
#' \code{\link{clipper}}
#' @examples
#' toClip(1, 5, c(0,0))
#' 
#' x <- tile_read(base::system.file("extdata/tile.bsp", package = "uFTIR"))
#' x <- tile_base_corr(x)
#' x <- wavealign(x, primpke)
#' x <- tile_sam(x)
#' x <- smooth_sam(x, as.integer(length(primpke@clusternames)), window = 3, 1)
#' 
#' clip <- toClip(8,20,c(10,10))
#' polygon(clip@xycoords)
#' 
#' x <- clipper(x, clip@centre, clip@rad, 1)
#' 
toClip <- function(rad = 1, segments = 5, centre = c(0,0)){
  
  angle <- seq(0,2*pi, by = 2*pi/segments)
  x <- sin(angle) * rad + centre[1]
  y <- cos(angle) * rad + centre[2]
  
  new("clipmask",
      xycoords = cbind(x,y),
      rad = as.integer(trunc(rad)),
      centre = centre)
}

#' Clipper
#' 
#' A function to clip a \code{\link[=SAM-class]{SAM}} or \code{\link[=Smooth-class]{Smooth}} object object (also a matrix, but I recommed this only for debugging). ALthough you can use this function in a 'step by step' process it is otherwise directly called by \code{\link{summary_sam}} which uses a \code{\link[=clipmask-class]{clipmask}} object (returned by \code{\link{toClip}} as a clipping mask.
#'
#' @param tarjet object of class \code{\link[=SAM-class]{SAM}} or \code{\link[=Smooth-class]{Smooth}} (or matrix) to clip. 
#' @param centre The coordinates of the polygon centre (x,y).
#' @param rad The circle radius.
#' @param slice which slice of the \code{\link[=SAM-class]{SAM}} or \code{\link[=Smooth-class]{Smooth}} object should be clipped? (only one!)
#' 
#' @return
#' The tarjet object clipped as matrix (and S3 clipper -not for human consumption).
#' @export
#' @seealso 
#' \code{\link{toClip}} code{\link[=clipmask-class]{clipmask}} \code{\link{summary_sam}}.
#' @examples
#' x <- tile_read(base::system.file("extdata/tile.bsp", package = "uFTIR"))
#' x <- tile_base_corr(x)
#' x <- wavealign(x, primpke)
#' x <- tile_sam(x)
#' x <- smooth_sam(x, as.integer(length(primpke@clusternames)), window = 3, 1)
#' 
#' clip <- toClip(8,20,c(10,10))
#' plot(x)
#' polygon(clip@xycoords)
#' 
#' x <- clipper(x, clip@centre, clip@rad, 1)
#' 
clipper <- function(tarjet, centre = c(128,128), rad = 120, slice = 1){
  
  if(!class(tarjet) %in% c("SAM", "Smooth", "matrix")){
    stop("Invalid Class. Tarjet should be of class SAM, Smooth or matrix")
  }
  
  if(class(tarjet) == 'Smooth'){
    tarjet <- as.matrix(tarjet@smooth[,, slice])
    
  } else if(class(tarjet) == "SAM"){
    tarjet <- as.matrix(tarjet@clusters[,, slice])
    
  } else {
    warning("Tarjet's class is not Smooth\n Tarjet will be treated as matrix")
  }
  
  # as ploting the rows are inverted, then:
  centre[2] <- nrow(tarjet)-centre[2]
  
  g = expand.grid(1:ncol(tarjet), 1:nrow(tarjet)) #all the coordinates
  g$d2 = sqrt ((g$Var1-centre[1])^2 + (g$Var2-centre[2])^2) #distance to centre
  g$inside = g$d2<=rad #is the distance smaller than the radius
  
  tarjet[as.matrix(g[!g$inside,c("Var2","Var1")])] <- NA #NA to all values outside
  class(tarjet) <- c("clipper", "matrix")
  tarjet
}

#' As Clipper
#' 
#' Coerce a \code{\link[=SAM-class]{SAM}} or \code{\link[=Smooth-class]{Smooth}} object to clipper. A function useful to call \code{\link{highlight_substance}} when clipping is not necessary.
#'
#' @param object The \code{\link[=SAM-class]{SAM}} or \code{\link[=Smooth-class]{Smooth}} object to be coerced.
#' @param slice The slice to keep (a clipper is essentially a matrix).
#'
#' @return
#' and S3 clipper
#' 
#' @export
#' @examples
#' x <- tile_read(base::system.file("extdata/tile.bsp", package = "uFTIR"))
#' x <- tile_base_corr(x)
#' x <- wavealign(x, primpke)
#' x <- tile_sam(x)
#' x <- as.clipper(x)
as.clipper <- function(object, slice = 1){
  if(class(object) == "SAM"){
    x <- object@clusters[,,slice]
    x <- as.matrix(x)
  } else {
    x <- object@smooth[,,slice]
    x <- as.matrix(x)  
  }
  class(x) <- c("clipper", "matrix")  
  x
}
