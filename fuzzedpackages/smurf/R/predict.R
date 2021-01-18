###############################################
#
# Model predictions
#
###############################################


#' @export
#' @title Predictions Using Estimated Model
#' 
#' @description Function to obtain predictions using the estimated model.
#' 
#' @param object An object of class '\code{\link[=glmsmurf-class]{glmsmurf}}', typically the result of a call to \code{\link{glmsmurf}} or \code{\link{glmsmurf.fit}}.
#' @param newdata Optionally, a data frame containing the predictors used in the prediction. 
#'        This can only be used when \code{object} contains a formula. 
#'        When \code{newdata} is omitted, the predictions are based on the data used to fit the model in \code{object}.
#' @param newoffset Optionally, a vector containing a new offset to be used in the prediction.
#'                  When \code{newoffset} is omitted, the predictions use the offset which was used to fit the model in \code{object}.
#' @param type Type of prediction. The default is on the scale of the linear predictors (\code{"link"}).
#'             Another option is on the scale of the response variable (\code{"response"}). 
#'             For type \code{"terms"} a matrix containing the fitted values of each term in the model, on the linear predictor scale, is returned.
#' @param ... Additional arguments which are currently ignored.
#' 
#' @return A vector containing the predicted values using the estimated model in \code{object}.
#' 
#' @seealso \code{\link{predict_reest}}, \code{\link[stats]{predict.glm}}, \code{\link[stats]{predict}}, 
#'          \code{\link{glmsmurf}}, \code{\link{glmsmurf-class}}
#'          
#' @examples ## See example(glmsmurf) for examples
#'          
predict.glmsmurf <- function(object, newdata = NULL, newoffset = NULL, type = c("link", "response", "terms"), ...) {
  
  return(.predict.glmsmurf.aux(object = object, newdata = newdata, newoffset = newoffset, type = type, reest = FALSE, ...))
}

#' @export
#' @title Predictions Using Re-estimated Model
#' 
#' @description Function to obtain predictions using the re-estimated model.
#' 
#' @param object An object for which predictions are meaningful. 
#'               E.g. an object of class '\code{\link[=glmsmurf-class]{glmsmurf}}', typically the result of a call to \code{\link{glmsmurf}} or \code{\link{glmsmurf.fit}}.
#' @inheritParams predict.glmsmurf
#' 
#' @return A vector containing the predicted values using the re-estimated model in \code{object}, 
#'         when this is available, or, otherwise, the predicted values using the estimated model in \code{object} with a warning.
#' 
#' @seealso \code{\link{predict.glmsmurf}}, \code{\link[stats]{predict.glm}}, \code{\link[stats]{predict}}, 
#'          \code{\link{glmsmurf}}, \code{\link{glmsmurf-class}}
#'          
#' @examples ## See example(glmsmurf) for examples
#'           
predict_reest <- function(object, ...) UseMethod("predict_reest", object)

#' @export
#' @rdname predict_reest
predict_reest.glmsmurf <- function(object, newdata = NULL, newoffset = NULL, type = c("link", "response", "terms"), ...) {
  
  return(.predict.glmsmurf.aux(object = object, newdata = newdata, newoffset = newoffset, type = type, reest = TRUE, ...))
}


# Auxiliary function for predictions using estimated model (reest=FALSE) or re-estimated model (reest=TRUE)
.predict.glmsmurf.aux <- function(object, newdata = NULL, newoffset = NULL, 
                                type = c("link", "response", "terms"), reest = FALSE, ...) {
  
  # Get prediction type
  type <- match.arg(type)
  
  # Handle reest
  if (!exists("coefficients.reest", object) & reest) {
    warning("Coefficients of the re-estimated model are not present in 'object', coefficients of the estimated model are used.")
    reest <- FALSE
  }
  
  # Select coefficients to use (when type == "terms" or when newdata is provided)
  if (reest) {
    # Use re-estimated coefficients
    coefs <- as.vector(coef_reest(object))
    
  } else {
    # Use estimated coefficients
    coefs <- as.vector(coef(object))
  }
  
  # Set newoffset to NULL when missing
  if (missing(newoffset)) {
    newoffset <- NULL
  }
  
  # Input checking for new offset
  if (!is.null(newoffset)) {
    
    # Check if newoffset is numeric
    if (!is.numeric(newoffset)) {
      stop("'newoffset' must be a numeric vector or NULL.")
    }
    
    # Convert to vector if given as matrix
    if (!is.vector(newoffset)) {
      newoffset <- as.vector(newoffset)
    }
  }
  
  # Set newdata to NULL when missing
  if (missing(newdata)) {
    newdata <- NULL
  }
  
  # newdata is provided
  if (!is.null(newdata)) {
    
    if (is.null(object$formula)) {
      stop("'newdata' can only be used if a formula is present in 'object'.")
    }

    # Get terms of object and delete response
    Terms <- delete.response(object$terms)
    
    # Get original x-levels
    xlevels.orig <- object$xlevels
    xlevels <- xlevels.orig
    
    # Indices of predictors with a 2D Fused Lasso penalty
    ind.2dflasso <- which(object$pen.cov == "2dflasso")
    
    # Ignore x-levels of predictors with a 2D Fused Lasso penalty
    if (length(ind.2dflasso) > 0) {
      
      for (j in 1:length(ind.2dflasso)) {
        
        xlevels[[ind.2dflasso[j] - 1]] <- NULL
      }
    }
    
    # Make model frame for newdata based on terms and x-levels 
    # (ignore x-levels of predictors with a 2D Fused Lasso penalty)
    mf <- model.frame(formula = Terms, data = newdata, xlev = xlevels)

    # Treat predictors with a 2D Fused Lasso penalty separately.
    # Observations with at least one of the 1D levels equal to the 1D reference level
    # are changed to the 2D reference level
    if (length(ind.2dflasso) > 0) {
      
      for (j in 1:length(ind.2dflasso)) {
        
        # All levels that are present in original data after setting all levels, 
        # with at least one of the 1D levels equal to the 1D reference level, to the 2D reference level
        l <- xlevels.orig[[ind.2dflasso[j] - 1]]
        
        # Indices of all levels in interaction which are not present in l
        ind.ref <- which(!(levels(mf[, ind.2dflasso[j] - 1]) %in% l))
        
        
        for (i in 1:length(ind.ref)) {
          
          # Necessary but not sufficient check to see if at least one of the 1D levels 
          # is equal to the 1D reference level. 
          # Does not necessarily work when a dot is present in at least one of the 1D levels 
          # since we split on a dot for the interaction
          if (!any(strsplit(levels(mf[, ind.2dflasso[j] - 1])[ind.ref[i]], '[.]')[[1]] %in% 
                   strsplit(l[1], '[.]')[[1]])) {
            stop(paste0("Invalid level for the predictor '", attr(mf[, ind.2dflasso[j] - 1], "cov.name"), "."))
          }
        }

        # Change levels not in l to reference level to make sure that all levels with at least one 
        # of the 1D levels equal to the 1D reference level are now equal to the 2D reference level.
        # Note that the check above is not sufficient!
        levels(mf[, ind.2dflasso[j] - 1])[ind.ref] <- l[1]
      }
    }
    
    
    # Create sparse model matrix taking contrasts into account
    X.new <- sparse.model.matrix(object = Terms, data = mf, contrasts.arg = object$contrasts)

    if (ncol(X.new) != length(coefs)) {
      stop(paste0("Invalid 'newdata' provided."))
    }
    
    # Rename columns
    colnames(X.new) <- .rename_mm.cols(colnames(X.new))
    
    
    # Remove assign and contrasts attributes (if present)
    if ("assign" %in% names(attributes(X.new))) {
      attr(X.new, "assign") <- NULL
    }
    
    if ("contrasts" %in% names(attributes(X.new))) {
      attr(X.new, "contrasts") <- NULL
    }

    
    #################
    
    # Keep using old offset
    if (is.null(newoffset)) {
      newoffset <- object$offset
      
      # Make sure dimension of constant offset is correct
      if (length(unique(newoffset)) == 1L) {
        newoffset <- rep(unique(newoffset), nrow(X.new))
      }
    } 
    
    # Always check length of newoffset
    if (length(newoffset) != nrow(X.new)) {
      stop(paste0("'newoffset' must be a numeric vector of length ", nrow(X.new), "."))
    }
    
    
    #################
    
    if (type == "link") {
      # Linear predictors
      lp <- as.numeric(X.new %*% coefs + newoffset)
      names(lp) <- rownames(X.new)
      return(lp)
      
    } else if (type == "response") {
      # Fitted values
      fv <- object$family$linkinv(as.numeric(X.new %*% coefs + newoffset))
      names(fv) <- rownames(X.new)
      return(fv)
      
    } else if (type == "terms") {
      # Fitted values of each term
      return(sweep(X.new, 2L, coefs, "*"))
    }
    
  } else {
    # Use returned objects
    
    # Check length of newoffset
    if (!is.null(newoffset)) {
      if (length(newoffset) != length(object$offset)) {
        stop(paste0("'newoffset' must be a numeric vector of length ", length(object$offset), " or NULL."))
      }
    }

    
    if (type == "link") {
      # Linear predictors
      
      if (reest) {
        # Linear predictors after re-estimation
        lp <- object$linear.predictors.reest
        
      } else {
        # Linear predictors without re-estimation
        lp <- object$linear.predictors
      }
      
      if (is.null(newoffset)) {
        # No new offset given
        return(lp)
        
      } else {
        # Correction for new offset
        lp2 <- as.numeric(lp + (newoffset - object$offset))
        names(lp2) <- names(lp)
        return(lp2)
      }
      
    } else if (type == "response") {
      # Fitted values
      
      if (is.null(newoffset)) {
        # No new offset given
        
        if (reest) {
          # Fitted values after re-estimation
          return(object$fitted.values.reest)
          
        } else {
          # Fitted values without re-estimation
          return(object$fitted.values)
        }
        
      } else {
        
        if (reest) {
          # Linear predictors after re-estimation
          lp <- object$linear.predictors.reest
          
        } else {
          # Linear predictors without re-estimation
          lp <- object$linear.predictors
        }
        
        # Correction for new offset
        fv <- object$family$linkinv(as.numeric(lp + (newoffset - object$offset)))
        names(fv) <- names(lp)
        return(fv)
      }
      
    } else if (type == "terms") {
      
      if (is.null(object$X)) {
        stop("Terms cannot be predicted when 'object' does not contain 'X'. 
             Please provide the data in 'newdata' or use the glmsmurf function with option 'x = TRUE'.")
      }
      
      # Fitted values of each term
      return(sweep(object$X, 2L, coefs, "*"))
    }
  }
}
