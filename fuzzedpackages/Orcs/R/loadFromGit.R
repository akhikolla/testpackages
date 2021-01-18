#' Install and load a package from GitHub
#' 
#' @description
#' This function comprises multiple steps required to install and load a 
#' package directly from GitHub.
#' 
#' @param repo Repository address as \code{character}, defaults to 
#' "fdetsch/Orcs".
#' @param ... Additional arguments passed to 
#' \code{\link[remotes]{install_github}}.
#' 
#' @author 
#' Florian Detsch
#' 
#' @seealso
#' \code{\link[remotes]{install_github}}
#' 
#' @examples
#' \dontrun{
#' ## install 'Orcs' development version from GitHub
#' loadFromGit("fdetsch/Orcs", ref = "develop")
#' }
#' 
#' @export loadFromGit
#' @name loadFromGit
loadFromGit <- function(repo = "fdetsch/Orcs", ...) {
  ## install desired package
  remotes::install_github(repo, ...)
  
  ## load package
  ls_pkg <- strsplit(repo, "/")
  ch_pkg <- sapply(ls_pkg, "[[", 2)
  library(ch_pkg, character.only = TRUE)
  
  return(invisible())
}