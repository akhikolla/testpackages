###############################################
#
# Model residuals
#
###############################################


#' @export
#' @title Residuals of Estimated Model
#' 
#' @description Function to extract the residuals of the estimated model. 
#'              \code{resid} is an \emph{alias} for it.
#' 
#' @param object An object of class '\code{\link[=glmsmurf-class]{glmsmurf}}', typically the result of a call to \code{\link{glmsmurf}} or \code{\link{glmsmurf.fit}}.
#' @param type Type of residuals that should be returned. One of \code{"deviance"} (default), 
#'             \code{"pearson"}, \code{"working"}, \code{"response"} or \code{"partial"}.
#' @param ... Additional arguments which are currently ignored.
#' 
#' @return A vector containing the residuals of the estimated model in \code{object}.
#' 
#' @details See \code{\link[stats]{glm.summaries}} for an overview of the different types of residuals.
#' 
#' @seealso \code{\link{residuals_reest}}, \code{\link{residuals}}, \code{\link[stats]{glm.summaries}}, \code{\link{glmsmurf-class}}
#' 
#' @examples ## See example(glmsmurf) for examples
#'   
residuals.glmsmurf <- function(object, type = c("deviance", "pearson", "working", "response", "partial"), ...) {
  
  return(.residuals.glmsmurf.aux(object = object, type = type, reest = FALSE, ...))
}

# Alias for residuals.glmsmurf
#' @export
#' @rdname residuals.glmsmurf
resid.glmsmurf <- residuals.glmsmurf


#' @export
#' @title Residuals of Re-estimated Model
#' 
#' @description Function to extract the residuals of the re-estimated model. 
#'              \code{resid_reest} is an \emph{alias} for it.
#' 
#' @param object An object for which the extraction of model residuals is meaningful. 
#'               E.g. an object of class '\code{\link[=glmsmurf-class]{glmsmurf}}', typically the result of a call to \code{\link{glmsmurf}} or \code{\link{glmsmurf.fit}}.
#' @inheritParams residuals.glmsmurf
#' 
#' @return A vector containing the residuals of the re-estimated model in \code{object}
#'         when they are available, or, otherwise, the residuals of the estimated model in \code{object} with a warning.
#' 
#' @details See \code{\link[stats]{glm.summaries}} for an overview of the different types of residuals.
#' 
#' @seealso \code{\link{residuals.glmsmurf}}, \code{\link{residuals}}, \code{\link[stats]{glm.summaries}}, \code{\link{glmsmurf-class}}
#' 
#' @examples ## See example(glmsmurf) for examples
#'  
residuals_reest <- function(object, ...) UseMethod("residuals_reest", object)

#' @export
#' @param type Type of residuals that should be returned. One of \code{"deviance"} (default), 
#'             \code{"pearson"}, \code{"working"}, \code{"response"} or \code{"partial"}.
#' @rdname residuals_reest
residuals_reest.glmsmurf <- function(object, type = c("deviance", "pearson", "working", "response", "partial"), ...) {
  
  return(.residuals.glmsmurf.aux(object = object, type = type, reest = TRUE, ...))
}

# Alias for residuals_reest
#' @export
#' @rdname residuals_reest
resid_reest <- function(object, ...) UseMethod("residuals_reest", object)

#' @export
#' @rdname residuals_reest
resid_reest.glmsmurf <- residuals_reest.glmsmurf

# Auxiliary function for residuals of estimated model (reest=FALSE) or re-estimated model (reest=TRUE)
.residuals.glmsmurf.aux <- function(object, type = c("deviance", "pearson", "working", "response", "partial"), 
                                    reest = FALSE, ...) {
  
  # Get type of residuals
  type <- match.arg(type)
  
  if (!exists("residuals.reest", object) & reest) {
    warning("Residuals of the re-estimated model are not present in 'object', residuals of the estimated model in 'object' are used.")
    reest <- FALSE
  }
  
  # Response vector
  y <- object$y
  # Linear predictors
  if (reest) {
    eta <- object$linear.predictors.reest
  } else {
    eta <- object$linear.predictors
  }
  # Fitted values
  if (reest) {
    mu <- object$fitted.values.reest
  } else {
    mu <- object$fitted.values
  }
  # Weights
  weights <- object$weights
  
  # y - mu, see glm function
  # These are the residuals for type == "reponse"
  if (reest) {
    resid <- object$residuals.reest * object$family$mu.eta(eta)
  } else {
    resid <- object$residuals * object$family$mu.eta(eta)
  }
  
  
  if (type == "deviance") {
    
    if (is.null(y)) {
      # Compute y based on residuals
      y <- mu + resid
      # Round y if very close to integer
      ind.round <- which(abs(y - round(y)) < eps_num)
      y[ind.round] <- round(y)[ind.round]
    }
    
    if (object$df.residual > 0) {
      # sign(y - mu) * sqrt(dev_i)
      # Take pmax to avoid numeric issues
      resid <- sign(resid) * sqrt(pmax(object$family$dev.resids(y = y, mu = mu, wt = weights), 0))
      
    } else {
      resid <- rep(0, length(mu))
    }
    
  } else if (type == "pearson") {
    # (y - mu) * sqrt(weights) / sqrt(Var(mu))
    resid <- resid * sqrt(weights) / sqrt(object$family$variance(mu))
    
  } else if (type == "working") {
    # (y - mu) / (d mu / d eta(eta))
    if (reest) {
      resid <- object$residuals.reest
    } else {
      resid <- object$residuals
    }
    
  } else if (type == "partial") {
    # (y - mu) + terms
    if (reest) {
      resid <- object$residuals.reest + predict_reest.glmsmurf(object = object, type = "terms")
    } else {
      resid <- object$residuals + predict.glmsmurf(object = object, type = "terms")
    }
  }
  
  return(resid)
}
