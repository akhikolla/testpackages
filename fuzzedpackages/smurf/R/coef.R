###############################################
#
# Model coefficients
#
###############################################

#' @export
#' @title Coefficients of Estimated Model
#' 
#' @description Function to extract the coefficients of the estimated model. 
#'              \code{coefficients} is an \emph{alias} for it.
#' 
#' @param object An object of class '\code{\link[=glmsmurf-class]{glmsmurf}}', typically the result of a call to \code{\link{glmsmurf}} or \code{\link{glmsmurf.fit}}.
#' @param ... Additional arguments which are currently ignored.
#' 
#' @return A vector containing the coefficients of the estimated model in \code{object}.
#' 
#' @seealso \code{\link{coef_reest}}, \code{\link[stats]{coef}}, \code{\link{summary.glmsmurf}}, \code{\link{glmsmurf}}, \code{\link{glmsmurf-class}}
#' 
#' @examples ## See example(glmsmurf) for examples
#' 
coef.glmsmurf <- function(object, ...) {
  
  return(object$coefficients)
}

#' @export
#' @rdname coef.glmsmurf
coefficients.glmsmurf <- coef.glmsmurf


# Generic function since coef_reest does not exist yet for S3 objects
#' @export
#' @title Coefficients of Re-estimated Model
#' 
#' @description Function to extract the coefficients of the re-estimated model. 
#'              \code{coefficients_reest} is an \emph{alias} for it.
#' 
#' @param object An object for which the extraction of model coefficients is meaningful. 
#'               E.g. an object of class '\code{\link[=glmsmurf-class]{glmsmurf}}', typically the result of a call to \code{\link{glmsmurf}} or \code{\link{glmsmurf.fit}}.
#' @inheritParams coef.glmsmurf
#' 
#' @return A vector containing the coefficients of the re-estimated model in \code{object}, 
#'         when they are available, or, otherwise, the coefficients of the estimated model in \code{object} with a warning.
#' 
#' @seealso \code{\link{coef.glmsmurf}}, \code{\link[stats]{coef}}, \code{\link{summary.glmsmurf}}, \code{\link{glmsmurf}}, \code{\link{glmsmurf-class}}
#' 
#' 
#' @examples ## See example(glmsmurf) for examples
#' 
coef_reest <- function(object, ...) UseMethod("coef_reest", object)

#' @export
#' @rdname coef_reest
coef_reest.glmsmurf <- function(object, ...) {
  
  if (!exists("coefficients.reest", object)) {
    warning("Coefficients of the re-estimated model are not present in 'object', coefficients of the estimated model are returned.")
    return(object$coefficients)
    
  } else {
    
    return(object$coefficients.reest)
  }
  
}

# Alias for coef_reest
#' @export
#' @rdname coef_reest
coefficients_reest <- function(object, ...) UseMethod("coefficients_reest", object) 

#' @export
#' @rdname coef_reest
coefficients_reest.glmsmurf <- coef_reest.glmsmurf

