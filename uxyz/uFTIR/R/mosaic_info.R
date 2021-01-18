#' Mosaic Info
#' 
#' The function is equivalent to \code{\link{tile_read}}. It makes an \code{\link[=SpectralInfo-class]{SpectralInfo}} object using the .dmt file. It does not load the measured spectra, and it was written for mosaic images only. See the page of \code{\link[=SpectralInfo-class]{SpectralInfo}} class for a full description of the returned object.
#'
#' @param dmtfile Path to the .dmt file.
#'
#' @return
#' An object of class \code{\link[=SpectralInfo-class]{SpectralInfo}}.
#' @export
#' @seealso 
#' \code{\link[=SpectralInfo-class]{SpectralInfo}}
#' @examples
#' x <- mosaic_info(base::system.file("extdata/mosaic.dmt", package = "uFTIR"))
mosaic_info <- function(dmtfile){
  # Get the general information of the mosaic
  # size, date, wavelength, etc.
  
  # path.expand
  if(length(grep("^\\.", dmtfile)) == 1){
    dmtfile <- gsub("^\\.", getwd(), dmtfile)
    
  } else if(length(grep("^~", dmtfile)) == 1){
    dmtfile <- path.expand(dmtfile)
    
  } else if(.Platform$OS.type == "windows" & 
            length(grep("^[[:upper:]]:.*", dmtfile)) == 1){
    dmtfile <- dmtfile
    
  } else if(length(
    grep(paste("^[^", .Platform$file.sep, "]", sep = ""), dmtfile)) == 1){
    dmtfile <- paste(getwd(), dmtfile, sep = .Platform$file.sep)
    
  } 
  path <- gsub("/[^/]*.bsp$", "", dmtfile)
  path <- gsub("/[^/]*.dmt", "", dmtfile)
  
  # Get the wavenumbers and the number of measured points
  filesize <- file.info(dmtfile)$size
  fi <- file(dmtfile, open = "rb")

  seek(con = fi, where = 2228, origin = "start")
  startwavenumber <- readBin(con = fi, what = integer(), n = 1, endian = "little")

  seek(fi, 2236, "start")
  numberofpoints <- readBin(fi, integer(), 1, endian = "little")

  seek(fi, 2216, "start")
  wavenumberstep <- readBin(fi, "double", 1, endian = "little")
  
  wavenumbers <- 1:(numberofpoints + startwavenumber - 1)
  wavenumbers <- wavenumbers * wavenumberstep
  wavenumbers <- wavenumbers[startwavenumber:length(wavenumbers)]
  
  ##
  # Get the date form the dmt file
  oldlocale <- Sys.getlocale("LC_TIME")
  Sys.setlocale("LC_TIME", "C")
  
  seek(fi, 0, rw = "r")
  rdate <- readBin(fi, "character", filesize, endian = "little")
  
  suppressWarnings(
    dposition <- unlist(
      sapply(month.name, function(x) {grep(x, rdate)})
    )
  )
  
  rdate <- rdate[dposition]
  rdate <- as.POSIXct(rdate, format = c("%A, %B %d, %Y %H:%M:%S"))
  
  Sys.setlocale("LC_TIME", oldlocale)
  
  close(fi)

  # Get the FPA size from the tile 0000_0000.dmd file
  filesize <- file.info(gsub('.dmt$', '_0000_0000.dmd', dmtfile))$size
  fpasize <- sqrt(((filesize/4) - 255)/length(wavenumbers))

  # Return mosaic object
  new("SpectralInfo",
      #file = gsub(".dtm$", "", dmtfile),
      file = dmtfile,
      date = rdate,
      fpasize = fpasize,
      wavenumbers = wavenumbers,
      path = path
      )
}
