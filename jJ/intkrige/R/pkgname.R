# See the DESCRIPTION FILE for the text that populates the description section
# in the RD file.
#' INTERVAL-VALUED KRIGING MODELS IN R
#' @aliases intkrige-package
#'
#' @section Functions:
#' The package contains several generic functions such as
#'  \code{\link[=plot,intsp,missing-method]{plot}},
#'  \code{\link[=print.intsp]{print}},
#'  and \code{\link[=summary.intsp]{summary}} to facilitate
#'  interval-valued analysis. In addition to these functions,
#'  the package also contains the following functions:
#'
#' \itemize{
#' \item \code{\link{intkrige}} Make predictions on
#'  interval-valued spatial data using interval-valued kriging.
#' \item \code{\link{interval<-}} Create an interval-valued spatial object.
#' \item \code{\link{intvariogram}} Simultaneously calculate the empirical variograms
#'  for an interval-valued spatial object.
#' \item \code{\link{fit.intvariogram}} Automatically fit theoretical variograms to
#'  empirical variograms obtained from interval-valued spatial objects.
#' \item \code{\link{intvCheck}} Visualize the variogram fits for interval-valued
#'  spatial objects.
#' \item \code{\link{dist_cpp}} A c++ function to calculate great circle
#'  or Euclidean distances.
#' }
#'
#' @section Data:
#' \itemize{
#' \item \link{utsnow} An interval-valued design ground snow load dataset for Utah.
#' \item \link{ohtemp} An interval-valued 30 year mean temperature dataset for the
#'  Ohio River Basin.
#' \item \link{ohMap} A \code{\link[sp:SpatialPolygons]{SpatialPolygons}}
#'   shapefile for the Ohio River Basin.
#' }
#'
#'
#' @author Brennan Bean \email{brennan.bean.20@@gmail.com}
#' @importFrom Rdpack reprompt
"_PACKAGE"


