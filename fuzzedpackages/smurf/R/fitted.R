###############################################
#
# Fitted values
#
###############################################


#' @export
#' @title Fitted Values of Estimated Model
#' 
#' @description Function to extract the fitted values of the estimated model.
#' 
#' @param object An object of class '\code{\link[=glmsmurf-class]{glmsmurf}}', typically the result of a call to \code{\link{glmsmurf}} or \code{\link{glmsmurf.fit}}.
#' @param ... Additional arguments which are currently ignored.
#' 
#' @return A vector containing the fitted values of the estimated model in \code{object}.
#' 
#' @seealso \code{\link{fitted_reest}}, \code{\link[stats:fitted.values]{fitted}}, \code{\link{glmsmurf}}, \code{\link{glmsmurf-class}}
#' 
#' @examples ## See example(glmsmurf) for examples
#'     
fitted.glmsmurf <- function(object, ...) {
  
  return(object$fitted.values)
}


#' @export
#' @title Fitted Values of Re-estimated Model
#' 
#' @description Function to extract the fitted values of the re-estimated model.
#' 
#' @param object An object for which the extraction of fitted values is meaningful. 
#'               E.g. an object of class '\code{\link[=glmsmurf-class]{glmsmurf}}', typically the result of a call to \code{\link{glmsmurf}} or \code{\link{glmsmurf.fit}}.
#' @inheritParams fitted.glmsmurf
#' 
#' @return A vector containing the fitted values of the re-estimated model in \code{object}, 
#'         when they are available or, otherwise, the fitted values of the estimated model in \code{object} with a warning.
#'
#' @seealso \code{\link{fitted.glmsmurf}}, \code{\link[stats:fitted.values]{fitted}}, \code{\link{glmsmurf}}, \code{\link{glmsmurf-class}}
#
#' @examples ## See example(glmsmurf) for examples
#'  
fitted_reest <- function(object, ...) UseMethod("fitted_reest", object)
#fitted.values_reest <- function(object, ...) UseMethod("fitted_reest", object) # Alias for fitted_reest

#' @export
#' @rdname fitted_reest
fitted_reest.glmsmurf <- function(object, ...) {
  
  if (!exists("fitted.values.reest", object)) {
    warning("Fitted values of the re-estimated model are not present in 'object', fitted values of the estimated model in 'object' are used.")
    
    # Return fitted values of estimated model
    return(object$fitted.values)
    
  } else {
    # Return fitted values of re-estimated model
    return(object$fitted.values.reest)
  }
}
