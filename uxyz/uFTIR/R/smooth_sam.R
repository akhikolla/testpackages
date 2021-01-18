#' Smooth a Spectral Angle Mapper output
#'
#' @description The function makes a smooth version of the object returned by \code{\link{tile_sam}} or \code{\link{mosaic_compose}}. The function smooths out only the cluster slot, slicing it to a user-given range of slices.
#'
#' @param x Object of class \code{\link[=SAM-class]{SAM}}.
#' @param nclusters How many clusters did the spectral library had? This parameter refer back to the length of the clusterlist slot of the \code{\link[=SpectralReference-class]{SpectralReference}} object used when calling \code{\link{tile_sam}} or \code{\link{mosaic_sam}}. It should be an integer.
#' @param window The function smooths out the \code{\link[=SAM-class]{SAM}} objectc "x" using a moving window. You should decide which size (in pixels) this window should be.
#' @param nslices Integer. Starting from 1, up to which slice of the cube holded in the clusters slot of the \code{\link[=SAM-class]{SAM}} objectc "x" do you want the smooth_sam function to work? It will influence the output. See return.
#'
#' @return
#' An S4 object of class \code{\link[=Smooth-class]{Smooth}}. It has only one slot (smooth) that holds an array in which each slice is the smooth version of the corresponding slice of a \code{\link[=SAM-class]{SAM}} object.
#' 
#' The returned object has methods to plot, summarize, and extract profiles by calling the generics \code{\link[graphics]{plot}}, \code{\link[base]{summary}}, and \code{\link{get_profile}}. You might have to provide extra arguments. You can refer to the see also section to look for the object specific methods.
#'
#' @export
#' @seealso 
#' \code{\link{plot_tile}} \code{\link{summary_sam}} \code{\link{get_profile_tile}} \code{\link{get_profile_sinfo}} 
#'
#' @examples
#' a <- matrix(c(1,1,1,1,3,3,
#'               2,2,2,1,2,3,
#'               3,1,1,2,3,1,
#'               3,3,2,2,1,1,
#'               1,1,1,3,3,2,
#'               3,3,3,2,2,3), nrow = 6, byrow = TRUE)
#' a <- array(c(a,a), dim = c(6,6,2))
#' test_sam <- new(Class = "SAM",
#'                raw_sam = a,
#'                substances = a,
#'                clusters = a)
#' 
#' smooth_sam(test_sam, nclusters = 3, window = 3)
smooth_sam <- function(x, nclusters, window = 5, nslices = 1){
  nslices <- nslices -1
  out <- csmooth_sam(x@clusters, window, nclusters, nslices)
  #out
  new("Smooth",
      smooth = out)
}
