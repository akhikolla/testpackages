###############################################
#
# Model deviance
#
###############################################


#' @export
#' @title Deviance of Estimated Model
#' 
#' @description Function to extract the deviance of the estimated model.
#' 
#' @param object An object of class '\code{\link[=glmsmurf-class]{glmsmurf}}', typically the result of a call to \code{\link{glmsmurf}} or \code{\link{glmsmurf.fit}}.
#' @param ... Additional arguments which are currently ignored.
#' 
#' @return The deviance of the estimated model in \code{object}.
#' 
#' @seealso \code{\link{deviance_reest}}, \code{\link[stats]{deviance}}, \code{\link{summary.glmsmurf}}, 
#'          \code{\link{glmsmurf}}, \code{\link{glmsmurf-class}}
#'          
#'
#' @examples ## See example(glmsmurf) for examples
#'          
deviance.glmsmurf <- function(object, ...) {
  
  return(object$deviance)
}


#' @export
#' @title Deviance of Re-estimated Model
#' 
#' @description Function to extract the deviance of the re-estimated model.
#' 
#' @param object An object for which the extraction of the deviance is meaningful. 
#'               E.g. an object of class '\code{\link[=glmsmurf-class]{glmsmurf}}', typically the result of a call to \code{\link{glmsmurf}} or \code{\link{glmsmurf.fit}}.
#' @inheritParams deviance.glmsmurf
#' 
#' @return The deviance of the re-estimated model in \code{object}, 
#'         when it is available or, otherwise, the deviance of the estimated model in \code{object} with a warning.
#' 
#' @seealso \code{\link{deviance.glmsmurf}}, \code{\link[stats]{deviance}}, \code{\link{summary.glmsmurf}}, 
#'          \code{\link{glmsmurf}}, \code{\link{glmsmurf-class}}
#' 
#' @examples ## See example(glmsmurf) for examples
#'    
deviance_reest <- function(object, ...) UseMethod("deviance_reest", object)


#' @export
#' @rdname deviance_reest
deviance_reest.glmsmurf <- function(object, ...) {
  
  if (!exists("deviance.reest", object)) {
    warning("Deviance of the re-estimated model is not present in 'object', deviance of the estimated model is used.")
    
    # Return deviance of estimated model
    return(object$deviance)
    
  } else {
    # Return deviance of re-estimated model
    return(object$deviance.reest)
  }
}
