#' Load multiple packages 
#' 
#' @description
#' Load and attach multiple packages at once.
#' 
#' @param pkgs Packages to load as \code{character}. 
#' @param ... Additional arguments passed to \code{\link{library}}, except for 
#' 'character.only' which is set to \code{TRUE}.
#' 
#' @author 
#' Florian Detsch
#' 
#' @note 
#' Package startup messages are automatically disabled.
#' 
#' @seealso 
#' \code{\link{library}}.
#' 
#' @examples
#' loadPkgs(c("raster", "rgdal"))
#' 
#' @export loadPkgs
#' @name loadPkgs
loadPkgs <- function(pkgs, ...) {
  
  ## if 'character.only' has been specified, remove it from '...'
  dots = list(...)
  if ("character.only" %in% names(dots)) {
    dots = dots[-grep("character.only", names(dots))]
  }
  
  jnk = sapply(pkgs, function(x) {
    dots_sub = append(list(package = x, character.only = TRUE), dots)
    suppressPackageStartupMessages(
      do.call(library, args = dots_sub)
    )
  })

  return(invisible())  
}