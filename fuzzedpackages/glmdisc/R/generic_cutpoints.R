#' @include allClasses.R
NULL

#' Obtaining the cutpoints and / or regroupments of a discretization.
#' @name cutpoints
#' @rdname cutpoints-method
#' @docType methods
#' @exportMethod cutpoints
#' @param object generic glmdisc object
#' @author Adrien Ehrhardt.
#' @description This defines the generic method "cutpoints" which will provide the cutpoints of a discretization scheme of S4 class \code{\link{glmdisc}}.
methods::setGeneric("cutpoints", function(object) attributes(object))
