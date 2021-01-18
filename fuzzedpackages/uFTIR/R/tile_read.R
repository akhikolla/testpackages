#' Single tile read
#'
#' @description The function loads a single tile to R as an object of class \code{\link[=Tile-class]{Tile}}. 
#' 
#' The function is also capable to load the toy mosaic produced by the Agilent sofware. You can do so by passing as bspfile parameter the .bsp file comprised in a mosaic folder. However, beware! the toy mosaic IS NOT the complete mosaic, it is only a reduced version of it (i.e. you will lose resolution, and you will not longer know the pixel size). Use it with care.
#'
#' @param bspfile path to the .bsp file.
#' 
#' @export
#' @return
#' An S4 Object of class \code{\link[=Tile-class]{Tile}}.
#' 
#' @examples
#' x <- tile_read(base::system.file("extdata/tile.bsp", package = "uFTIR"))
tile_read <- function(bspfile){
  # path.expand
  if(length(grep("^\\.", bspfile)) == 1){
    bspfile <- gsub("^\\.", getwd(), bspfile)
    
  } else if(length(grep("^~", bspfile)) == 1){
    bspfile <- path.expand(bspfile)
    
  } else if(.Platform$OS.type == "windows" & 
            length(grep("^[[:upper:]]:.*", bspfile)) == 1){
    bspfile <- bspfile
    
  } else if(length(
    grep(paste("^[^", .Platform$file.sep, "]", sep = ""), bspfile)) == 1){
    bspfile <- paste(getwd(), bspfile, sep = .Platform$file.sep)
    
  } 
  path <- gsub("/[^/]*.bsp$", "", bspfile)
  name <- bspfile
  datfile <- gsub("\\.bsp$", "\\.dat", bspfile)
  
  # Retrieve data from de info file (.bsp)
  # Open connection
  filesize <- file.info(bspfile)$size
  bspfile <- file(bspfile, open = 'rb')
  
  # Wavenumbers
  seek(con = bspfile, where = 2228, origin = 'start')
  startwavenumber <- readBin(con = bspfile, what = integer(), n = 1, endian = 'little')
  
  seek(bspfile, 2236, 'start')
  numberofpoints <- readBin(bspfile, integer(), 1, endian = 'little')
  
  seek(bspfile, 2216, 'start')
  wavenumberstep <- readBin(bspfile, "double", 1, endian = 'little')
  
  wavenumbers <- 1:(numberofpoints+startwavenumber-1)
  wavenumbers <- wavenumbers * wavenumberstep
  wavenumbers <- wavenumbers[startwavenumber:length(wavenumbers)]
  
  # Date
  oldlocale <- Sys.getlocale("LC_TIME")
  Sys.setlocale("LC_TIME", "C") # Otherwise %A and %B don't work
  
  seek(bspfile, 0, rw = 'r')
  rdate <- readBin(bspfile, 'character', filesize, endian = 'little')
  suppressWarnings(
    dposition <- unlist(sapply(month.name, function(x){grep(x, rdate)}))
  )
  rdate <- rdate[dposition]
  rdate <- as.POSIXct(rdate, format = c('%A, %B %d, %Y %H:%M:%S'))
  
  Sys.setlocale("LC_TIME", oldlocale)
  
  close(bspfile)
  
  # Retrieve readings from data file (.dat)
  # datfile <- imagefiles[grep('.dat$', imagefiles)]
  filesize <- file.info(datfile)$size
  datfile <- file(datfile, open = 'rb')
  
  # We cannot read float32 in R.
  # We can handle the float value by setting the size argument to 4,
  # since a float has a size of 4 bytes (instead of the 8 of a double).
  sdata <- readBin(datfile, 'double', filesize, size=4, endian = 'little')
  close(datfile)
  
  # Get de FPA size [It should be 128x128 pixels in our case]
  fpasize <- sqrt(((filesize / 4) - 255) / length(wavenumbers))
  
  # Tidy up the data
  sdata <- sdata[256:length(sdata)]
  sdata <- array(sdata, dim = c(fpasize, fpasize, numberofpoints))
  
  # Rotate the image to match the spectrometer's output
  #sdata <- aperm(sdata, c(2, 1, 3))
  #sdata <- sdata[fpasize:1, ,]
  #sdata <- raster::brick(sdata)
  
  ##
  # Create the return object
  ##
  new("Tile",
      #file = imagefiles[grep('.bsp$', imagefiles)],
      file = name,
      date = rdate,
      fpasize = fpasize,
      Spectra = sdata,
      wavenumbers = wavenumbers,
      path = path
      )
}
