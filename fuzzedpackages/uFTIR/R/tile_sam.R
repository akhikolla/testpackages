#' Spectral Angle Mapper
#' 
#' Performs the Spectral Angle Mapper to match the \code{\link[=SpectralReference-class]{SpectralReference}} with yout readings (as \code{\link[=Tile-class]{Tile}}).
#'
#' @param SpectralPack an object of class \code{\link[=SpectralPack-class]{SpectralPack}}.
#' @param derivative whether to apply the first (1) or second (2) derivative to preprocess the data before the algorithm is applied. Default NULL (It tells the program to not derivate).
#'
#' @return
#' An object of class \code{\link[=SAM-class]{SAM}}.
#' @export
#' @seealso 
#' For its application to mosaic images see \code{\link{mosaic_sam}}.
#' @examples
#' x <- tile_read(base::system.file("extdata/tile.bsp", package = "uFTIR"))
#' x <- tile_base_corr(x)
#' x <- wavealign(x, primpke)
#' x <- tile_sam(x)
tile_sam <- function(SpectralPack, derivative = NULL){
  stopifnot(class(SpectralPack) == "SpectralPack")
  
  if(is.null(derivative)){
    img <- SpectralPack@Readings@Spectra
    em <- SpectralPack@Reference@Spectra
    
  } else if(derivative == 1){
    img <- cderivate_cube(SpectralPack@Readings@Spectra, SpectralPack@Readings@wavenumbers)
    em <- cderivate_mat(SpectralPack@Reference@Spectra, SpectralPack@Reference@wavenumbers)
  
  } else if(derivative == 2){
    img <- cderivate_cube(SpectralPack@Readings@Spectra, SpectralPack@Readings@wavenumbers)
    em <- cderivate_mat(SpectralPack@Reference@Spectra, SpectralPack@Reference@wavenumbers)
    
    wave <- SpectralPack@Reference@wavenumbers
    wave_diff <- wave[-1] - diff(wave) / 2
    
    img <- cderivate_cube(img, wave_diff)
    em <- cderivate_mat(em, wave_diff)
  
  } else {
    stop("Invalid derivative!")
  }
  
  out <- ctile_sam(img, em)

  clusterlist <- SpectralPack@Reference@clusterlist
  clusters.out <- cmosaic_clusterfind(out$match_sam, as.vector(clusterlist))
  
  new("SAM",
      raw_sam = out$raw_sam,
      substances = out$match_sam,
      clusters = clusters.out
  )
}
