#' @include classes.R
NULL

#' Tile Preprocess
#' 
#' @description The function is a wrapper that allows the user to perform a user-define function to each pixel spectrum. It can be called internally in \code{\link{mosaic_sam}} using the FUN argument when processing mosaics.
#' 
#' As exemplified below, the function allows you to interact with other R packages that provide common features to analyze spectral data. Beware that you will need to adjust the wavenumbers manually if you are resampling... unless that 'data' is an object of class \code{\link[=SpectralPack-class]{SpectralPack}}.
#'
#' @param data An object of class \code{\link[=Tile-class]{Tile}} or \code{\link[=SpectralPack-class]{SpectralPack}}.
#' @param FUN A function to apply to the Spectra slot.
#'
#' @return
#' The same object but with the Spectra slot updated by FUN. \code{\link[=SpectralPack-class]{SpectralPack}} objects return the wavelengths adjusted to the new extent, numbered from 1 to n (original values are lost).
#' @export
#' @rdname preprocess
#' @examples
#' x <- tile_read(base::system.file("extdata/tile.bsp", package = "uFTIR"))
#' x <- preprocess(x, function(x){x+1})
#' 
#' # The function allows interacting with other R packages that provide common features
#' # for spectra analysis. For example, you can use the package "signal" to run
#' # a Savitzky-Golay filter.
#'
#' if(requireNamespace("signal", quietly = TRUE)){
#' 
#' # for Tile objects. NOTE that after the preprocess x@wavenumbers does not match dim(x@Spectra)[3]
#' x <- tile_read(base::system.file("extdata/tile.bsp", package = "uFTIR"))
#' x <- preprocess(x, function(x){signal::sgolayfilt(x)})
#' dim(x@Spectra)[3] == length(x@wavenumbers) # BEWARE!
#' 
#' # for SpectralPack objects
#' x <- tile_read(base::system.file("extdata/tile.bsp", package = "uFTIR"))
#' x <- tile_base_corr(x)
#' x <- wavealign(x, primpke)
#' preprocess(x, function(x){signal::sgolayfilt(x)})
#' 
#' # Here the problem with the wavenumbers is gone
#' dim(x@Readings@Spectra)[3] == length(x@Readings@wavenumbers)
#' dim(x@Reference@Spectra)[3] == length(x@Reference@wavenumbers)
#' length(x@Readings@wavenumbers) == length(x@Reference@wavenumbers)
#' 
#' }
#' 
setGeneric("preprocess", function(data, FUN){warning("Method defined only for Tiles...")})

#' @rdname preprocess
#' @export
setMethod("preprocess", c(data = "Tile", FUN = "function"),
          function(data, FUN) {
            out <- data@Spectra
            out <- apply(out, c(1,2), FUN)
            out <- aperm(out, c(2,3,1))
            data@Spectra <- out
            data
          })

#' @rdname preprocess
#' @export
setMethod("preprocess", c(data = "SpectralPack", FUN = "function"),
          function(data, FUN){
            out_0 <- data@Readings@Spectra
            out_1 <- data@Reference@Spectra
            
            out_0 <- apply(out_0, c(1,2), FUN)
            out_0 <- aperm(out_0, c(2,3,1))
            
            out_1 <- apply(out_1, 1, FUN)
            out_1 <- t(out_1)
            
            data@Readings@Spectra <- out_0
            data@Reference@Spectra <- out_1
            
            data@Readings@wavenumbers <- seq(1, dim(out_0)[3])
            data@Reference@wavenumbers <- seq(1, ncol(out_1))
            
            data
          })

