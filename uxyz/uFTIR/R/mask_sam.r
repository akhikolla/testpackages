#' Mask SAM
#' 
#' The spectral angle mapper algorithm always find a match, as it was designed to compare between to spectra. However, somethimes the match is far from good and the user could want to drop the poorly matched pixels. This function will create a mask for such a purpose.
#'
#' @param x an array. Tipically, the array at the raw_sam slot of a SAM object.
#' @param threshold cut off value above wich, the SAM relation found is deemed unsignificant
#'
#' @return
#' S3 matrix
#' @export
#' @examples
#' set.seed(4356)
#' x <- array(abs(rnorm(1000)), c(10, 10, 10))
#' x <- mask_sam(x, 0.1)
#' # trick to plot as clipper
#' class(x) <- c("clipper", "matrix")
#' plot(x, legend = FALSE)
mask_sam <- function(x, threshold = 1.1){
  out <- array(NA, dim = dim(x)[1:2])
  for(i in 1:nrow(x)){
    for(j in 1:ncol(x)){
      out[i,j] <- min(x[i,j,])
    }
  }
  # FALSE for NA values, that are over the threshold
  out[out > threshold] <- NA
  out <- !is.na(out)
  out
}

