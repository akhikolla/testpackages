#' Mosaic chunk
#' 
#' The function read a single mosaic tile (.dmd file extension) using a \code{\link[=SpectralInfo-class]{SpectralInfo}} object as a guide to find the file. It returns an object of class \code{\link[=Tile-class]{Tile}}, which can be (pre/post)processed as if it were a single tile.
#'
#' @param info Optional. \code{\link[=SpectralInfo-class]{SpectralInfo}} object.
#' @param path Optional. Path to dmdfile. Overwrites the path on the info object. 
#' @param dmdfile Target *.dmd file to read.
#' @param fpa If info is missing, you should indicate the fpa size of the chunk.
#' @param wl  If info is missing, you should indicate the wavenumbers at which you read your data
#'
#' @return
#' A \code{\link[=Tile-class]{Tile}} object.
#' 
#' @export
#' @examples
#' x <- mosaic_info(base::system.file("extdata/mosaic.dmt", package = "uFTIR"))
#' y <- mosaic_chunk(info = x, dmdfile = "mosaic_0000_0000.dmd")
#' class(y)
mosaic_chunk <- function(info, path, dmdfile, fpa, wl){
  
  # re-writes dmdfile adding the path if info or path where provided.
  if(!missing(info) | !missing(path)){
    if(missing(info)){
      dmdfile <- paste(path, dmdfile, sep = .Platform$file.sep)
      
    } else if(missing(path)){
      dmdfile <- paste(info@path, dmdfile, sep = .Platform$file.sep)
      
    } else{
      # prioritize path over info.
      dmdfile <- paste(path, dmdfile, sep = .Platform$file.sep)
    }
  }
  
  # Get FPA size
  if(missing(info) & missing(fpa)){
    stop("You should provide the fpa size if you don't want to pass a SpectralInfo object through the info argument\n")
  }
  if(missing(fpa)){
    fpa <- info@fpasize
  }
  
  # Get wavenumbers
  if(missing(info) & missing(wl)){
    stop("You should provide the wavenumbers (as a vector) if you don't want to pass a SpectralInfo objecto through the info argument\n")
  }
  if(missing(wl)){
    wl <- info@wavenumbers
  }

  # Read the file  
  chunk <- mosaic_read_chunk(dmdfile, fpa, length(wl))
  
  file <- gsub("[[:alpha:]|[:punct:]|[:digit:]]*/","",dmdfile)
  path <- gsub(file, "", dmdfile)
  
  if(missing(info)){
    dates <- as.POSIXct("Monday, March 1, 1998 01:24:15", 
                       format = c("%A, %B %d, %Y %H:%M:%S"))
  }else{
    dates <- info@date
  }
                
  new("Tile",
      file = file,
      date = dates,
      fpasize = fpa,
      Spectra = chunk,
      wavenumbers = wl,
      path = path
  )
}