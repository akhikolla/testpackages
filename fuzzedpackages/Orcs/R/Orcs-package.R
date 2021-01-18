#' Omnidirectional R Code Snippets
#'
#' The package provides a variety of functions which I regularly use during my 
#' everyday work.
#'
#' @name Orcs-package
#' @aliases orcspackage
#' @docType package
#' @title Omnidirectional R Code Snippets.
#' @author Florian Detsch \cr
#' \cr
#' \emph{Maintainer:} Florian Detsch \email{florian.detsch@@staff.uni-marburg.de}
#'
#' @import methods grid raster rgdal knitr lattice latticeExtra
#' @importFrom bookdown render_book
#' @importFrom grDevices rgb
#' @importFrom plotrix thigmophobe
#' @importFrom Rcpp sourceCpp
#' @importFrom remotes install_github
#' @importFrom sf st_as_sf
#' @importFrom sp proj4string
#' @importFrom stats coef complete.cases  median na.exclude sd
#' @rawNamespace useDynLib(Orcs, .registration = TRUE)
#' 
#' @keywords package
#'
NULL
#'
#' @docType data
#' @name KiLi
#' @title Bing Aerial Image of Kilimanjaro
#' @description Bing aerial image of Kilimanjaro downloaded from 
#' \href{https://www.openstreetmap.org/}{OpenStreetMap}.
#' @format A \code{"RasterStack-class"} with 3 bands (red, green, blue).
#' @details Copyright: OpenStreetMap contributors, see 
#' \url{https://www.openstreetmap.org/copyright}.
#' 
NULL