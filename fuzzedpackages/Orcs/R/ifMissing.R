#' Take measures in case of nonexisting target files
#' 
#' @description 
#' If a target file already exists, it is simply being imported into R. However, 
#' if the specified target file does not exist, it is first created by a 
#' user-defined function and subsequently returned, thus rendering explicit 
#' calls to \code{\link{file.exists}} unnecessary.
#' 
#' @param ofl Target file name as \code{character}.
#' @param fun0 If 'ofl' exists, \code{function} to be applied to it (defaults to 
#' \code{\link[raster]{brick}}). 
#' @param fun1 If 'ofl' does not exist, \code{function} used to create it 
#' (defaults to \code{\link{writeRaster}}).
#' @param arg1 Argument in 'fun1' (as \code{character}) that corresponds to 
#' 'ofl', e.g. \code{filename} in \code{\link{writeRaster}} or \code{file} in 
#' \code{\link{write.table}}. If missing (default), the target file name passed 
#' to 'fun1' needs to be explicitly included via '...'.
#' @param ... Additional arguments passed to 'fun0,fun1'.
#' 
#' @return 
#' If 'ofl' has already existed, the contents of 'ofl' derived from 'fun0'; else 
#' the output resultant from 'fun1'.
#' 
#' @seealso \code{\link{file.exists}}, \code{\link{do.call}}.
#' 
#' @author Florian Detsch
#' 
#' @examples
#' # simply import existing file
#' logo <- system.file("external/rlogo.grd", package = "raster")
#' s <- ifMissing(logo) 
#' 
#' # create nonexisting file and import it afterwards
#' logo2 <- file.path(tempdir(), "rlogo.tif")
#' s2 <- ifMissing(logo2, arg1 = "filename", x = s, datatype = "INT1U")
#' 
#' # this also works with text files and more sophisticated custom functions
#' fun = function(x, file = "", ...) {
#'   write.csv(x, file, ...)
#'   read.csv(file)
#' }
#' 
#' data(iris)
#' ofl <- file.path(tempdir(), "iris.csv")
#' iris2 <- ifMissing(ofl, fun1 = fun, x = iris, file = ofl, quote = FALSE, row.names = FALSE)
#' 
#' @export ifMissing
#' @name ifMissing
ifMissing = function(ofl, fun0 = raster::brick, fun1 = raster::writeRaster, 
                     arg1, ...) {
  
  if (file.exists(ofl)) {
    do.call(fun0, args = list(ofl))
    
  } else {
    dots = list(...)
    
    if (!missing(arg1)) {
      ofl1 = list(ofl) 
      names(ofl1) = arg1
      dots = append(dots, ofl1)
    }
    
    do.call(fun1, args = dots)
  }
}