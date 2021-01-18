#' Perform SAM between to matrices
#'
#' @description Perform the spectral angle mapper algorithm between to matrices.
#' 
#' @param x matrix
#' @param y matrix
#' @param x.wl wavenumbers for x
#' @param y.wl wavenumbers for y
#'
#' @return
#' vector
#' @export
#'
matrix_sam <- function(x, y, x.wl, y.wl){
  
  # An internal function to align x and y
  .wavealign <- function(x, y, x.wl, y.wl){
    # two matrices
    min <- max(min(x.wl), min(y.wl))
    max <- min(max(x.wl), max(y.wl))
    
    # y will be reshaped
    waveout <- x.wl
    waveout <- waveout[waveout > min & waveout < max]
    
    #y <- interpol(y, waveout)
    
    y <- apply(y, 1, function(y){
      approx(y.wl, y, xout = waveout)$y
    })
    
    y <- t(y)
    
    # x will be clipped
    x <- x[x.wl %in% waveout]
    
    return(list(x=x,y=y))
  }
  
  xy_aligned <- .wavealign(x, y, x.wl, y.wl)
  x <- xy_aligned$x
  y <- xy_aligned$y
 
  # Fix is x as a single row matrix
  if(!"matrix" %in% class(x) & class(x) == "numeric"){
    x <- matrix(x, nrow=1)
  }
  
  # Calculate SAM
  out <- sam_internal(x, y)
  
  out
}


