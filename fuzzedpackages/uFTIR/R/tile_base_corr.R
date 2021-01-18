#' Baseline correction
#'
#' A function to substract the minimun of a pixel (vector) to have a minimal signal = 0, avoiding negative readings. Currently the only preprocessing method implemented.
#'
#' @param data Object of class \code{\link[=Tile-class]{Tile}}.
#'
#' @return 
#' An object of class \code{\link[=Tile-class]{Tile}} congruent with data, with the spectra modified.
#' 
#' @export
#' @examples
#' x <- tile_read(base::system.file("extdata/tile.bsp", package = "uFTIR"))
#' x <- tile_base_corr(x)
tile_base_corr <- function(data){
  stopifnot(length(dim(data@Spectra)) == 3)

  # minimal signal = 0
  out <- apply(data@Spectra, c(1,2), function(x){x-min(x)})
  out <- aperm(out, c(2,3,1))
  data@Spectra <- out
  data
}
