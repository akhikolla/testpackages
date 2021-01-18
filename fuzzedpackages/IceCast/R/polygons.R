#' Polygon of all regions
#'
#' All regions except the non-regional ocean converted into a single
#' \code{SpatialPolygons} object. Regions are modified from the region masks
#'provided by the National Snow and Ice Data
#' Center (NSIDC)
#' @docType data
#' @format \code{SpatialPolygons} object
#' @keywords datasets
#' @references  Region Mask: National Snow and Ice Data Center, 2017: Region
#'              mask for the northern hemisphere.
#'              \url{http://nsidc.org/data/polar-stereo/tools_masks.html}.
#' @examples
#' data(all_regions)
#' plot(all_regions)
"all_regions"


#' Polygon of the non-regional ocean
#'
#' The non-regional ocean converted into a single \code{SpatialPolygons}
#' object. The boundaries of the non-regional ocean were defined by modifying
#' the region masks provided by the National Snow and Ice Data  Center (NSIDC).
#' @docType data
#' @format \code{SpatialPolygons} object
#' @keywords datasets
#' @references  Region Mask:
#' National Snow and Ice Data Center, 2017: Region mask for the northern hemisphere.
#' \url{http://nsidc.org/data/polar-stereo/tools_masks.html}.
#' @examples
#' data(bg_water)
#' plot(bg_water)
"bg_water"


#' Polygon of land
#'
#' Land mask as a single \code{SpatialPolygons} object. The land mask was
#' obtained from the CM2.5 Forecast-oriented Low-Ocean Resolution (FLOR) model
#' produced by the National Oceanic and Atmospheric Administration’s Geophysical
#' Fluid Dynamics Laboratory converted to a Polar Stereographic grid.  (Vecchi
#' et al. 2014; Msadek et al. 2014). Weights for converting to a polar
#' stereograhic grid were obtained from the spherical
#' coordinate remapping and interpolation package (SCRIP) (Jones 1997).
#' @docType data
#' @format \code{SpatialPolygons} object
#' @keywords datasets
#' @references Vecchi, Gabriel A., et al.
#'             \href{http://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-14-00158.1}{"On the seasonal forecasting of regional tropical cyclone activity."} Journal of Climate
#'             27.21 (2014): 7994-8016.
#'
#'             Msadek, R., et al.
#'             \href{http://onlinelibrary.wiley.com/doi/10.1002/2014GL060799/full}{"Importance of initial conditions in seasonal predictions
#'             of Arctic sea ice extent."} Geophysical Research Letters
#'             41.14 (2014): 5208-5215.
#'
#'             Jones, P.W. "A user’s guide for SCRIP: A spherical coordinate
#'             remapping and interpolation package." Los Alamos National
#'             Laboratory, Los Alamos, NM (1997).
#'
#'             Region Mask: National Snow and Ice Data Center, 2017: Region mask
#'             for the northern hemisphere.
#'             \url{http://nsidc.org/data/polar-stereo/tools_masks.html}.
#'
#' @examples
#' data(land)
#' plot(land)
"land"

#' Spatial collection example
#'
#' Example of a \code{SpatialCollections} object that contains a
#' \code{SpatialPolygons} object and a \code{SpatialLines} object
#' @docType data
#' @format \code{SpatialCollections} object
#' @keywords datasets
#' @examples
#' data(SpatialCollEx)
#' plot(spatialCollEx)
#' plot(spatialCollEx@lineobj, col = "red", add = TRUE)
#' plot(spatialCollEx@polyobj, col = "blue", add = TRUE)
"spatialCollEx"


