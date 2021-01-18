# These classes are extensions of classes in the sp package
# Source code for the sp package can be obtained from
#
# https://github.com/edzer/sp
#
NULL

# The following link saved my bacon when trying to document s4 methods:
# - https://github.com/variani/pckdev/wiki/Documenting-with-roxygen2#s4-methods
#' @import utils
#' @import methods
#' @importFrom  raster raster
NULL

#' An interval extension of a SpatialPixelsDataFrame
#'
#' @slot interval A matrix of two columns representing the lower and upper
#'   endpoints of an interval.
#'
#' @export
intgrd <- setClass("intgrd",
                   contains = c("SpatialPixelsDataFrame"),
                   slots = c(interval = "matrix"))

#' An interval extension of a SpatialPointsDataFrame
#'
#' @slot interval A matrix of two columns representing the lower and upper
#'   endpoints of an interval.
#'
#' @export
intsp <- setClass("intsp",
                  contains = c("SpatialPointsDataFrame"),
                  slots = c(interval = "matrix"))



#=============================================================================
#' Function to extract the interval of an intsp or intgrd object
#' @param x An object of class intsp or intgrd.
#' @return A matrix containing the interval data.
#' @name interval
#' @rdname interval-methods
#' @exportMethod interval
setGeneric("interval",
           function(x)
             standardGeneric("interval"))

#' Function to reassign the contents of the interval slot
#' @param x An object of class intsp or intgrd.
#' @param value Either a character vector of length two specifying
#'   the column names which will occupy the interval slot. Or, a matrix
#'   of two columns to fill the slot.
#' @name interval<-
#' @rdname interval-methods-assign
#' @exportMethod interval<-
setGeneric("interval<-",
           function(x, value)
             standardGeneric("interval<-"))

# Generic to inititate the contents of the "intsp" interval slot.

#' Function to fit empirical variograms for an interval-valued spatial object
#'
#' @param x An object of class intsp or intgrd.
#' @param formulas A list of length two specifying the formulas related to the
#'   centers and radii respectively.
#' @param ... Additional arguments for sp::variogram().
#' @return An object of class 'intvariogram' containing empirical variograms
#'   for the center, radius, and center/radius interaction.
#' @name intvariogram
#' @rdname intvariogram-methods
#' @exportMethod intvariogram
setGeneric("intvariogram",
           function(x, formulas = list(center ~ 1, radius ~ 1), ...)
             standardGeneric("intvariogram"))
#=============================================================================

#============================================================================
#' Convert intgrd or intsp object back to a data frame
#'
#' @param x An object of class \code{intsp} or class \code{intgrd}.
#' @return An object of class \code{data.frame}.
#' @name as.data.frame
#' @rdname interval.as.data.frame-methods
NULL

#' Extract subset of an \code{intsp} or \code{intgrd} object
#'
#' @param x An object of class \code{intsp} or class \code{intgrd}
#'   from which to replace elements.
#' @param name Character vector corresponding to the name of the column that
#'   will be extracted or replaced.
#' @param i,j,... indices specifying elements to extract or replace. See
#'   generic function documentation for details.
#' @param drop The requested column that may be reassigned.
#' @param value The new data used to replaced drop data in the desired slot.
#' @return An object of class \code{intsp} or \code{intgrd}.
#' @name extract
#' @rdname extract-methods
NULL





