#' @include classes.R
NULL

#' @title Methods to resample the wavenumbers
#'
#' @description It clips and resamples the wavenumbers of data.x and data.y to a common extent. Typically, it uses the wavenumbers of data.x to resample data.y and then the function clips the largest base on the narrowest.
#'
#' @param data.x An object of class \code{\link[=SpectralReference-class]{SpectralReference}} or \code{\link[=Tile-class]{Tile}}.
#' @param data.y An object of class \code{\link[=SpectralReference-class]{SpectralReference}} or \code{\link[=Tile-class]{Tile}}. It should be the other one.
#'
#' @details
#' There are two methods defined:
#' \itemize{
#'   \item \bold{TileRead.wavealign} - This method should be the only one used by the user.
#'   \item \bold{MosaicGetRef.wavealign} - This method is called internally by \code{\link{mosaic_sam}}. The method is different since it expect a numeric vector in data.x, which is taken from the wavenumbers slot of a \code{\link[=SpectralInfo-class]{SpectralInfo}} object. The implication is that, for mosaics, the \code{\link[=SpectralReference-class]{SpectralReference}} is always resampled according to the wavenumbers of the mosaics. It was programed this way, since mosaics are not loaded into the memory (R) until they are processed by mosaic_sam.
#' }
#'
#' @return
#' \itemize{
#'   \item \bold{TileRead.wavealign} - An object of class \code{\link[=SpectralPack-class]{SpectralPack}}
#'   \item \bold{MosaicGetRef.wavealign} - The \code{\link[=SpectralReference-class]{SpectralReference}} object passed as data.y, resampled to data.x.
#' }
#' @export
#' @examples
#' x <- tile_read(base::system.file("extdata/tile.bsp", package = "uFTIR"))
#' x <- tile_base_corr(x)
#' x <- wavealign(x, primpke)
setGeneric("wavealign", function(data.x, data.y){NULL})

#' @rdname wavealign
#' @export
setMethod("wavealign", c(data.x = "numeric", data.y = "SpectralReference"),
          function(data.x, data.y) MosaicGetRef.wavealign(data.x, data.y))

#' @rdname wavealign
#' @export
setMethod("wavealign", c(data.x = "Tile", data.y = "SpectralReference"),
          function(data.x, data.y) TileRead.wavealign(data.x, data.y))

#' @rdname wavealign
#' @export
setMethod("wavealign", c(data.x = "SpectralReference", data.y = "Tile"),
          function(data.x, data.y) TileRead.wavealign(data.x, data.y))

#' @rdname wavealign
#' @export
MosaicGetRef.wavealign <- function(data.x, data.y){
  # This method for wavelengthAling
  # expect a numeric vector in data.x
  # that has the wavelength of a mosaic to be processed

  # Therefore, the method uses -invariably- the wavelength of the mosaic
  # to resample the SpectralReference object stored in the data.y param.
  min <- max(min(data.x), min(data.y@wavenumbers))
  max <- min(max(data.x), max(data.y@wavenumbers))

  # The reference spectra will be reshaped
  waveout <- data.x
  waveout <- waveout[waveout > min & waveout < max]

  data.y@Spectra <- interpol(data = data.y, wave = waveout)
  data.y@wavenumbers <- waveout

  data.y
}

#' @rdname wavealign
#' @export
TileRead.wavealign <- function(data.x, data.y){
  # This method should be the 'usual' call of the user.

  # Define the range
  min <- max(min(data.x@wavenumbers), min(data.y@wavenumbers))
  max <- min(max(data.x@wavenumbers), max(data.y@wavenumbers))

  # The spectra in data.y will be reshaped
  waveout <- data.x@wavenumbers
  waveout <- waveout[waveout > min & waveout < max]

  data.y@Spectra <- interpol(data = data.y, wave = waveout)
  data.y@wavenumbers <- waveout

  # Signal correction & resizing of data.x according to waveout
  if(length(dim(data.x@Spectra)) == 3){
    data.x@Spectra <- data.x@Spectra[ , ,data.x@wavenumbers %in% waveout]
  } else if(length(dim(data.x@Spectra)) == 2){
    data.x@Spectra <- data.x@Spectra[ ,data.x@wavenumbers %in% waveout]
  } else{
    print("No correction made to data.x")
  }

  data.x@wavenumbers <- waveout

  # Return objects as list
  # list(data.x = data.x, data.y = data.y)
  if(class(data.x) == "Tile"){
    out <-  new("SpectralPack",
                Readings = data.x,
                Reference = data.y)
  } else {
    out <- new("SpectralPack",
               Readings = data.y,
               Reference = data.x)
  }
  out
}

# Different methods to interpolate the wavelengths ----

setGeneric("interpol", function(data, wave){NULL})

setMethod("interpol", c(data = "SpectralReference", wave = "numeric"),
          function(data, wave) SpectralReference.interpol(data, wave)
          )

setMethod("interpol", c(data = "Tile", wave = "numeric"),
          function(data, wave) Tile.interpol(data, wave)
          )

SpectralReference.interpol <- function(data, wave){
  # Keep track of original row names
  outnames <- rownames(data@Spectra)
  # Interpolate [Linear]
  out <- apply(data@Spectra, 1, function(x){
    approx(x = data@wavenumbers, y = x, xout = wave)$y
  })
  # Tidy up
  out <- t(out)
  colnames(out) <- wave
  rownames(out) <- outnames
  out
}

Tile.interpol <- function(data, wave){
  stopifnot(length(dim(data@Spectra)) == 3)

  out <- apply(data@Spectra, c(1,2), function(x){
    approx(x = data@wavenumbers, y = x, xout = wave)$y
  })
  out <- aperm(out, c(2,3,1))
  out
}
