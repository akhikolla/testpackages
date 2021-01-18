#' @include classes.R
NULL

setClassUnion("summary_sam", c("SAM", "Smooth", "clipper"))

#' Method to summarize \code{\link[=SAM-class]{SAM}} or \code{\link[=Smooth-class]{Smooth}} objects.
#' 
#' @description Summary method for SAM and Smooth objects.
#' 
#' @param object An object of class \code{\link[=SAM-class]{SAM}} or \code{\link[=Smooth-class]{Smooth}}. It accepts also object of S3 class "clipper" returned by the \code{\link{clipper}} function.
#' @param ... other parameters for summary_sam
#' @param mask a clipper mask returned by \code{\link{toClip}} function. It should be set to NULL To avoid clipping, the default.
#' @param clusternames The name of the clusters, it should match the clusternames list comprised in the \code{\link[=SpectralReference-class]{SpectralReference}} object passed to \code{\link{tile_sam}} or \code{\link{mosaic_sam}}. Mandatory if object is of class \code{\link[=SAM-class]{SAM}} when smooth = TRUE.
#' @param slice Which slice of object should be summarized?
#' @param window When smooth = TRUE the size of the window for \code{\link{smooth_sam}}.
#' @param smooth When object if of class \code{\link[=SAM-class]{SAM}}, should the function silently apply \code{\link{smooth_sam}} before summarizing?.
#' @param temporal if TRUE the raster and shape files from which the summary is extracted are written in the working directory.
#' 
#' @return 
#' A dataframe with three columns showing the cluster number, the area in pixels^2, and the cluster name as character (only if the parameeter clusternames is within the function call). Each row stand for a single particle. The function writes on disk two files "raster_out.tif" which is a raster of the object passed through the object paramenter, and a shape_out "ESRI Shapefile" folder with a shapefile holding a vectorized version of "raster_out.tif".
#' 
#' @export
#' @rdname summary_sam
#' @examples
#' x <- mosaic_info(base::system.file("extdata/mosaic.dmt", package = "uFTIR"))
#' mosaic_sam(x, primpke, n_cores = 1)
#' y <- mosaic_compose(x@path, clusterlist = primpke@clusterlist)
#' summary_sam(y, clusternames = primpke@clusternames, smooth = FALSE, temporal = TRUE)
setMethod("summary", c(object = "summary_sam"),
          function(object, ...) summary_sam(object, ...))

#' @export
#' @rdname summary_sam
summary_sam <- function(object, mask = NULL, clusternames = NULL, 
                        slice = 1, window = NULL, smooth = TRUE, temporal = FALSE){
  
  # Transform a SAM into a Smooth if the user did not.
  if("SAM" %in% class(object) & smooth == TRUE){
    if(is.null(window) | is.null(clusternames)){
      stop("When smooth = TRUE,\nyou need to provide clusternames & window size!", call. = F)
    }
    object <- smooth_sam(object, nclusters = as.integer(length(clusternames)),
                    window, nslices = slice)
  }
  
  if(!is.null(mask)){
    object <- clipper(object, c(mask@centre[1], mask@centre[2]), mask@rad, slice)
    class(object) <- c("matrix")
  }
  
  if(!is.matrix(object)){
    if("SAM" %in% class(object)){
      object <- as.matrix(object@clusters[,,slice])
    } else if("Smooth" %in% class(object)){
      object <- as.matrix(object@smooth[,,slice])
    }
  } else if("clipper" %in% class(object)){
    class(object) <- c("matrix")
  }
  
  src_raster <- raster::raster(object)
  src_raster@extent@xmax <- src_raster@ncols
  src_raster@extent@ymax <- src_raster@nrows
  # Is the proposed CRS the best choice?
  src_raster@crs <- CRS("+proj=utm +zone=11 +datum=WGS84 +units=m +no_defs +ellps=WGS84")
  
  if(temporal){
  raster::writeRaster(src_raster,
                      paste(tempdir(), "raster_out.tif", sep = .Platform$file.sep),
                      format = "GTiff", overwrite = TRUE)
  
  #require(GDALPolygonize)
  rgdal_polygonize(raster = paste(tempdir(), "raster_out.tif", sep = .Platform$file.sep), 
                                   folder = paste(tempdir(), "shape_out", sep = .Platform$file.sep),
                                   layer = "clusters", field = "cluster", overwrite = TRUE)
  
  tmp <- rgdal::readOGR(paste(tempdir(), "shape_out", sep = .Platform$file.sep), "clusters", verbose = F)
  
  } else {
    raster::writeRaster(src_raster,
                        "raster_out.tif",
                        format = "GTiff", overwrite = TRUE)
    
    #require(GDALPolygonize)
    rgdal_polygonize(raster = "raster_out.tif", 
                     folder = "shape_out",
                     layer = "clusters", field = "cluster", overwrite = TRUE)
    
    tmp <- rgdal::readOGR("shape_out", "clusters", verbose = F)
    
  }
  
  tmp_area <- raster::area(tmp)
  out <- cbind(tmp@data, tmp_area)
  if(!is.null(clusternames)){
    tmp_polymers <- clusternames[tmp@data[,1]]
    out <- cbind(out, tmp_polymers)
    colnames(out) <- c("cluster", "area", "clname") 
  } else {
    colnames(out) <- c("cluster", "area")  
  }
  out
}
