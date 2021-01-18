#' Mosaic compose
#' 
#' @description There are to ways to process a mosaic. One way -which should be the standard- is to do it internally in R by calling first \code{\link{mosaic_info}} and later \code{\link{mosaic_sam}}.
#' 
#' The function will write in the path specified in the \code{\link[=SpectralInfo-class]{SpectralInfo}} object returned by \code{\link{mosaic_info}} a series of binary files holding the SAM results. In other words, the function does not load the results to R automaticaly. You can call this function to load them back to R in a single object, which will be of class \code{\link[=SAM-class]{SAM}}.
#' 
#' @param path Where are the binary files? you can use the 'path' slot of the \code{\link[=SpectralInfo-class]{SpectralInfo}} instead of enter it manually.
#' @param clusterlist The clusterlist vector that you passed along with the \code{\link[=SpectralReference-class]{SpectralReference}} in the call to \code{\link{mosaic_sam}}.
#' @param nslices If you deal with large mosaics, you might want to load only a few of the sam matches. This argument allows you to define up to which match you want to load to R.
#' @param drop_raw If you are not interested in the angles you can set this argument to TRUE and avoid load them to R.
#'
#' @return
#' An object of class \code{\link[=SAM-class]{SAM}}.
#' 
#' @export
#' @seealso 
#' \code{\link{mosaic_sam}}
#' @examples
#' x <- mosaic_info(base::system.file("extdata/mosaic.dmt", package = "uFTIR"))
#' mosaic_sam(x, primpke, n_cores = 1)
#' y <- mosaic_compose(x@path, primpke@clusterlist)
mosaic_compose <- function(path = ".", clusterlist = NULL, nslices = NULL, drop_raw = FALSE){
  # get the list of images to process...
  # oldpath <- getwd()
  # setwd(path)

  imagefiles <- list.files(path = path, pattern = "*[0-9]{4}_[0-9]{4}\\.bin$",
                           full.names = TRUE)
  # for (i in 1:length(imagefiles)) {
  #   imagefiles[i] <- paste(getwd(), imagefiles[i], sep = .Platform$file.sep)
  # }
  # setwd(oldpath)
  
  if(is.null(nslices)){
    #the user wants all slices
    nslices <- -1
  } 
  
  # get the order of the images...
  # Do we need to change x,y to y,x? Check please!
  xy <- c()
  for(i in imagefiles){
    # Get the tile number for each file
    x <- unlist(strsplit(i, "_"))
    x <- x[grep('[0-9]{4}', x)]
    x <- suppressWarnings(
      # When the folders have numbered names an NA is returned... see next if clause
      as.numeric(gsub('\\.[[:alnum:]]*$', '', x))
    )
    
    # it was a bug: files with numbered names make problems
    # example_07082019/sample_1_0000_0000.bin
    if(length(x) > 2){
      x <- x[c(length(x)-1, length(x))] #keep only the last two numbers
    }
    
    if(length(xy) == 0){
      xy <- c(i, x)
    }else{
      xy <- rbind(xy, c(i, x))
    }
  }
  sam_files <- xy[1:(nrow(xy)/2), ]
  sub_files <- xy[(nrow(xy)/2+1):nrow(xy), ]
  
  # sam_files <- as.vector(sam_files[order(sam_files[,2]), 1])
  # sub_files <- as.vector(sub_files[order(sub_files[,2]), 1])
  sam_files <- as.vector(sam_files[order(sam_files[,3]), 1])
  sub_files <- as.vector(sub_files[order(sub_files[,3]), 1])
  
  xy_pos <- matrix(as.numeric(xy[ 1:(nrow(xy)/2), 2:3]), ncol = 2)
  
  ### BUG SOLVED:
  # Agilent Resolutions Pro software has serious issues when labeling the files.
  # They swap the dimmensions of the file. Then, when composing, the mosaic image is nonesensical. 
  # Ergo, we need to rename the files (virtually).
  
  # The files are named x_y so col_rows (by Agilent).
  # we need to check which vector repeat each number n times
  # and which repeat all numbers n times.
  if(isTRUE(sum(diff(grep(0, xy_pos[,1]), differences = 2)) == 0)){
    # Cond tests if x has each number repeated n times.
    alt.x <- rep(seq(min(xy_pos[,2]), max(xy_pos[,2])), each = max(xy_pos[,1])+1)
    alt.y <- rep(seq(min(xy_pos[,1]), max(xy_pos[,1])), times = max(xy_pos[,2])+1)
    xy_pos <- cbind(alt.x, alt.y)
  } else {
    # The first should be the only case.
    warning("TRUE location of tiles was not achieved!")
  }
  
  ### End of fix.
  
  out <- list()
  if(drop_raw){
    out$raw_sam <- array(NA, c(0,0,0))
  } else {
    out$raw_sam <- cmosaic_compose(sam_files, xy_pos, max(xy_pos[,1]), max(xy_pos[,2]), -1)  
  }
  out$match_sam <- cmosaic_compose(sub_files, xy_pos, max(xy_pos[,1]), max(xy_pos[,2]), nslices)
 
  if(!is.null(clusterlist)){
    out$clusters <- cmosaic_clusterfind(out$match_sam, as.vector(clusterlist))  
  }else{
    out$clusters <- array(NA, c(0,0,0))
  }
  
  # Reduce memory in the output
  mode(out$match_sam) <- 'integer'
  mode(out$clusters) <- 'integer'
  
  new("SAM",
      raw_sam = out$raw_sam,
      substances = out$match_sam,
      clusters = out$clusters
  )
}




