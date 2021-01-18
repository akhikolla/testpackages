#' @include classes.R
NULL

setClassUnion("togetprofile", members = c("SAM", "Smooth", "clipper"))

#' Get profile
#' 
#' @description To check whether the program is doing a good job or not (an the library is appropriate), it is usually a good idea to recover the expectra of all the pixels that match a given cluster or substance to plot the spectra and compare it to a standard (the library entrance for that cluster/substance for example).
#' 
#' This function allows you to retrieve the spectra of all pixels that match a given cluster/substance and returns a matrix that can then be plotted to do the checking. It also allows the user to plot the object and overlay the points that match a feature. The latter can be achieved also by using the \code{\link{highlight_substance}} function. However, that function is only implemented for S3 "clipper" objects, so this might be a way to overcome that limitation. 
#' 
#' get_profile_all will replace get_profile_sinfo some day. For now it does the same thing BUT it works only for clipper objects and you need to pass the fpa size and the wavenumber (waves) vector. It was designed to run in the cluster. It goes well with \code{\link{sam_write}}, \code{\link{sam_load}}.
#'
#' @param x Object of class \code{\link[=SAM-class]{sam}}, \code{\link[=Smooth-class]{Smooth}} or \code{\link{clipper}} that hold the processed tile/mosaic indicated in where.
#' @param where Original object processed. It might by of class \code{\link[=Tile-class]{Tile}} or \code{\link[=SpectralInfo-class]{SpectralInfo}}. However, in the second case, where should be a labeled list of two elements one the \code{\link[=SpectralInfo-class]{SpectralInfo}} object (labeled "info"), and the other the tarjet dmd file (as a character labeled "dmdfile"). For get_profile_all it is a character vector with the *.dmd files that you want to retrieve.
#' @param dst_cluster The expectra of the pixels that match which cluster/substance should be retrieved?
#' @param plotpol If true, the function will plot x and overlay countour lines highlighting dst_cluster. TRUE/FALSE.
#' @param plotpt If true, the function will plot x and draw points over the pixels that match dst_cluster. For mosaics, the points will be only in the chunk selected in where. TRUE/FALSE.
#' @param cluster If class(x) == "SAM", Should the substances or the clusters be used for analysis? TRUE means that the clusters are used, FALSE that the substances. If you set the argument to FALSE, dst_cluster and clusternames should be adjusted. Use the name of the substance and the substances list from the \code{\link[=SpectralReference-class]{SpectralReference}}.
#' @param slice If the object is of class \code{\link[=SAM-class]{sam}} or \code{\link[=Smooth-class]{Smooth}}, Which slice of them should be considered when retrieving the pixel locations?
#' @param clusternames if is.character(dst_cluster) you should provide a vector with the clusternames to coerse dst_cluster to numeric.
#' @param ... other parameters for plot.raster.
#' @param fpa The fpa size, used in get_profile_all
#' @param waves The wavenumbers, used in get_profile_all
#'
#' @return
#' matrix with the columns matching the wavenumbers of 'where'.
#' @export
#' @rdname get_profile
#' @examples
#' x <- mosaic_info(base::system.file("extdata/mosaic.dmt", package = "uFTIR"))
#' mosaic_sam(x, primpke, n_cores = 1)
#' y <- mosaic_compose(x@path, primpke@clusterlist)
#' y <- get_profile_sinfo(y,
#'                        where = list("info" = x, 
#'                                     "dmdfile" = "mosaic_0000_0000.dmd"),
#'                        5, FALSE, FALSE)
setGeneric("get_profile", function(x, where, ...){NULL})


#' @export
#' @rdname get_profile
setMethod("get_profile", c(x = "togetprofile", where = "Tile"),
          function(x, where, ...) get_profile_tile(x, where, ...))

#' @export
#' @rdname get_profile
setMethod("get_profile", c(x = "togetprofile", where = "list"),
          function(x, where, ...) get_profile_sinfo(x, where, ...))

#' @export
#' @rdname get_profile
setMethod("get_profile", c(x = "clipper", where = "character"),
          function(x, where, dst_cluster, fpa, waves){
            get_profile_all(x, where, dst_cluster, fpa, waves)
          })

#' @export
#' @rdname get_profile
get_profile_tile <- function(x, where, dst_cluster, plotpol = TRUE, 
                             plotpt = FALSE, cluster = TRUE, slice = 1, 
                             clusternames = NULL, ...){
  
  if("SAM" %in% class(x)){
    if(cluster) tarjet_sub <- x@clusters[,,slice]
    if(!cluster) tarjet_sub <- x@substances[,,slice]
    
  } else if("Smooth" %in% class(x)){
    tarjet_sub <- as.matrix(x@smooth[,,slice])
  } else {
    tarjet_sub <- x
  }
  
  if(is.character(dst_cluster)){
    if(!is.null(clusternames)){
      dst_cluster <- grep(dst_cluster, clusternames)
    } else {
      warning("To conver to dst_cluster from char to numeric I need the clusternames")
      return(NULL)
    }
  }
  
  if(plotpt | plotpol){
    buff_tarjet_sub <- tarjet_sub
    class(buff_tarjet_sub) <- c("clipper", "matrix")
  }
  
  if(length(class(tarjet_sub)) > 1 & "clipper" %in% class(tarjet_sub) ){
    class(tarjet_sub) <- c("matrix")
  }
  
  true_vector <- tarjet_sub == dst_cluster
  profile <- matrix(rep(NA), nrow = sum(true_vector), ncol = length(where@wavenumbers))
  a <- 1
  if(plotpt){
    xycords <- matrix(rep(NA), nrow = sum(true_vector), ncol = 2)
    b <- 1
  }
  for(i in 1:nrow(tarjet_sub)){
    for(j in 1:ncol(tarjet_sub)){
      if(true_vector[i, j]){
        profile[a,] <- where@Spectra[i,j,]
        if(plotpt){
          xycords[b,1] <- j
          xycords[b,2] <- i
          b <- b+1
        }
        a <- a+1
      }
    }
  }
    
  ## Plot if requested
  if(plotpol){
    highlight_substance(buff_tarjet_sub, dst_cluster, ...)
  }
  
  if(plotpt){
    plot(buff_tarjet_sub, ...)
    points(xycords[,2] ~ xycords[,1], cex = 0.5)
  }
  
  return(profile)
}

#' @export
#' @rdname get_profile
get_profile_sinfo <- function(x, where, dst_cluster, plotpol = TRUE, 
                              plotpt = FALSE, cluster = TRUE, slice = 1, 
                              clusternames = NULL, ...){
  
  if("SAM" %in% class(x)){
    if(cluster) tarjet_sub <- x@clusters[,,slice]
    if(!cluster) tarjet_sub <- x@substances[,,slice]
    
  } else if("Smooth" %in% class(x)){
    tarjet_sub <- as.matrix(x@smooth[,,slice])
  } else {
    tarjet_sub <- x
  }
  
  if(is.character(dst_cluster)){
    if(!is.null(clusternames)){
      dst_cluster <- grep(dst_cluster, clusternames)
    } else {
      warning("To conver to dst_cluster from char to numeric I need the clusternames")
      return(NULL)
    }
  }
  
  if(plotpt | plotpol){
    buff_tarjet_sub <- tarjet_sub
    class(buff_tarjet_sub) <- c("clipper", "matrix")
  }
  
  if( length(class(tarjet_sub)) > 1 & "clipper" %in% class(tarjet_sub) ){
    class(tarjet_sub) <- c("matrix")
  }
  
  true_vector <- tarjet_sub == dst_cluster 
  if(sum(true_vector) == 0){
    warning("No dst_cluster in x")
    return(0)
  }
  xycords <- matrix(rep(NA), nrow = sum(true_vector), ncol = 2)
  b <- 1
  for(i in 1:nrow(tarjet_sub)){
    for(j in 1:ncol(tarjet_sub)){
      if(true_vector[i, j]){
          xycords[b,1] <- j
          xycords[b,2] <- i
          b <- b+1
      }
    }
  }
  
  if(length(where) != 2 | sum(labels(where) != c("info", "dmdfile")) != 0 | class(where$info) != "SpectralInfo"){
    stop("where should be a labeled list c(\"info\", \"dmdfile\")")
  } 
  
  ## We need to fix the xycoords 
  ## Remember that Agilent file labeling is wrong
  ## I know it is super confusing
  file_size <- gsub("\\.dmd", "", 
                    gsub(gsub("\\.dmt", "", where$info@file), "", 
                         list.files(where$info@path, pattern = "\\.dmd", full.names = TRUE)))
  
  file_size <- unlist(strsplit(file_size, "_"))
  # get the original file numbers (Agilent)
  x_org <- as.numeric(file_size[seq(2, length(file_size), by = 3)])
  y_org <- as.numeric(file_size[seq(3, length(file_size), by = 3)])
  # get the correct version of the file numbers
  x_alt <- rep(seq(range(y_org)[1], range(y_org)[2]), each = max(x_org) +1)
  y_alt <- rep(seq(range(x_org)[1], range(x_org)[2]), times = max(y_org) +1)
  
  xy_alt <- cbind(x_alt, y_alt)
  # at which row is the file we want?
  loc <- grep(where$dmdfile, list.files(where$info@path, pattern = "\\.dmd"))
  
  # cropping the coordinates to the requested extent
  filter_row <- (where$info@fpasize * xy_alt[loc, 2]) < xycords[, 1] & 
    xycords[, 1] < ((where$info@fpasize * (xy_alt[loc, 2] + 1))+1)
  
  filter_col <- (where$info@fpasize * xy_alt[loc, 1]) < xycords[, 2] & 
    xycords[, 2] < ((where$info@fpasize * (xy_alt[loc, 1] + 1))+1)
  
  filter <- filter_row & filter_col
  
  xycords <- xycords[filter, ]
  
  if(nrow(xycords) == 0){
    #the cluster/sub is not in the chunk
    warning("No dst_cluster in x")
    return(0)
  }
 
  if(plotpt){
    buff_xycord <- xycords
  }
   
  #reshaping xycord to match the mosaic_chunk
  #the 0.00001 is to avoid 1/1 = 1 problem (the border condition)
  xycords[,1] <- xycords[,1] - trunc((max(xycords[,1])-0.00001) / where$info@fpasize) * where$info@fpasize
  xycords[,2] <- xycords[,2] - trunc((max(xycords[,2])-0.00001) / where$info@fpasize) * where$info@fpasize
  
  where <- mosaic_chunk(info = where$info, dmdfile = where$dmdfile)
  
  profile <- matrix(rep(NA), nrow = nrow(xycords), ncol = length(where@wavenumbers))
  for(i in 1:nrow(xycords)){
    profile[i, ]<- where@Spectra[xycords[i,2], xycords[i,1], ]
  }
  
  ## Plot if requested
  if(plotpol){
    highlight_substance(buff_tarjet_sub, dst_cluster, ...)
  }
  
  if(plotpt){
    plot(buff_tarjet_sub, ...)
    # The cols are rev...
    rev_cols <- (nrow(tarjet_sub)+1) - buff_xycord[,2]
    points(rev_cols ~ buff_xycord[,1], cex = 0.5)
  }
  
  return(profile)
}

###
# get_profile_all ----
###

#' @export
#' @rdname get_profile
get_profile_all <- function(x, where, dst_cluster, fpa, waves){
  
  true_vector <- x == dst_cluster 
  
  xycords <- matrix(rep(NA), nrow = sum(true_vector, na.rm=TRUE), ncol = 2)
  b <- 1
  for(i in 1:nrow(x)){
    for(j in 1:ncol(x)){
      if(!is.na(true_vector[i,j]) && true_vector[i, j]){
        xycords[b,1] <- j
        xycords[b,2] <- i
        b <- b+1
      }
    }
  }
  
  ## Get the Agilent possition of all single tiles
  fig_size <- list.files(pattern = "\\.dmd")
  fig_size <- strsplit(gsub("[^0-9_]+", "",fig_size), "_")
  fig_size <- t(sapply(fig_size, function(x){c(x[length(x)-1], x[length(x)])}))
  fig_size <- matrix(as.numeric(fig_size), nrow=nrow(fig_size),ncol=ncol(fig_size))
  
  # Get out version of the file numbers (which is correct)
  x_alt <- rep(seq(range(fig_size)[1], range(fig_size)[2]), each = max(fig_size) +1)
  y_alt <- rep(seq(range(fig_size)[1], range(fig_size)[2]), times = max(fig_size) +1)
  
  xy_alt <- cbind(x_alt, y_alt)
  
  b <- 1
  profile <- matrix(rep(NA), ncol =length(waves), nrow =nrow(xycords))
  for(i in where){
    loc <- grep(i, list.files(pattern = "\\.dmd"))
    
    # cropping the coordinates to the requested extent
    filter_row <- (fpa * xy_alt[loc, 2]) < xycords[, 1] & 
      xycords[, 1] < ((fpa * (xy_alt[loc, 2] + 1))+1)
    
    filter_col <- (fpa * xy_alt[loc, 1]) < xycords[, 2] & 
      xycords[, 2] < ((fpa * (xy_alt[loc, 1] + 1))+1)
    
    filter <- filter_row & filter_col
    
    chunk_xycords <- matrix(xycords[filter, ], ncol = 2)
    
    if(nrow(chunk_xycords) == 0){
      next
    }
    
    #reshaping xycord to match the mosaic_chunk
    #the 0.00001 is to avoid 1/1 = 1 problem (the border condition)
    chunk_xycords[,1] <- chunk_xycords[,1] - trunc((max(chunk_xycords[,1])-0.00001) / fpa) * fpa
    chunk_xycords[,2] <- chunk_xycords[,2] - trunc((max(chunk_xycords[,2])-0.00001) / fpa) * fpa
    
    target <- mosaic_chunk(dmdfile = i, fpa = fpa, wl = waves)
    
    for(j in 1:nrow(chunk_xycords)){
      profile[b, ]<- target@Spectra[chunk_xycords[j,2], chunk_xycords[j,1], ]
      b <- b+1
    }
    
  }
  return(profile)
}

