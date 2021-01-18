###############################################
#
# Penalty function for formula
#
###############################################


#' @export
#' 
#' @title Define Individual Subpenalties for a Multi-Type Regularized GLM
#' 
#' @description Function used to define regularization terms in a \code{\link{glmsmurf}} model formula.
#' 
#' @param pred1 Name of the predictor used in the regularization term.
#' @param pred2 Either \code{NULL} (default) meaning that only one predictor is used in the regularization term, or the name of the second predictor that is used in a 2D Fused Lasso regularization term.
#' @param pen Type of penalty for this predictor, one of 
#' \itemize{
#'                  \item \code{"none"} (no penalty),
#'                  \item \code{"lasso"} (Lasso),
#'                  \item \code{"grouplasso"} (Group Lasso),
#'                  \item \code{"flasso"} (Fused Lasso), 
#'                  \item \code{"gflasso"} (Generalized Fused Lasso),
#'                  \item \code{"2dflasso"} (2D Fused Lasso),
#'                  \item \code{"ggflasso"} (Graph-Guided Fused Lasso).
#'               }  
#' Default is \code{"lasso"}.
#' @param refcat Reference level when \code{pred1} is a factor and \code{pen} is \code{"none"}, \code{"flasso"}, \code{"gflasso"}, or \code{"ggflasso"};
#'               otherwise \code{refcat} is ignored. Default is \code{NULL} which means that the first level of \code{pred1} is used as the reference level (if \code{refcat} is not ignored).
#' @param group Group to which the predictor belongs, only used for a Group Lasso penalty. 
#'              Default is \code{NULL} which means that predictor does not belong to a group.
#' @details Predictors with no penalty, a Lasso penalty or a Group Lasso penalty should be numeric or a factor which can be non-numeric.  
#' Predictors with a Fused Lasso, Generalized Fused Lasso, Graph-Guided Fused Lasso or 2D Fused Lasso penalty should be given as a factor 
#' which can also be non-numeric. When a predictor is given as a factor, there cannot be any unused levels.
#' 
#' For a predictor with a Fused Lasso penalty, the levels should be ordered from smallest to largest. 
#' The first level will be the reference level, but this can be changed using the \code{refcat} argument.
#' 
#' When \code{lambda * lambda1 > 0} or \code{lambda * lambda2 > 0} in \code{\link{glmsmurf}}, no reference level is used
#' for the Fused Lasso, Generalized Fused Lasso and Graph-Guided Fused Lasso penalties, and \code{refcat} will hence be ignored.
#' 
#' If \code{pred2} is different from \code{NULL}, \code{pen} should be set to \code{"2dflasso"}, and vice versa.
#' Note that there cannot be any unused levels in the interaction between \code{pred1} and \code{pred2}.
#' 
#' When adding an interaction between \code{pred1} and \code{pred2} with a 2D Fused Lasso penalty, the 1D effects
#' should also be present in the model and the reference categories for the 1D predictors need to be the respective first levels. 
#' The reference level for the 2D predictor will then be the 2D level where it least one of the 1D components is equal to the 1D reference levels. 
#' It is also allowed to add binned factors, of predictors
#' that are included in the model, in the interaction. They should have the original predictor name + '.binned' as predictor names.
#' For example: the original predictors 'age' and 'power' are included in the model and
#' the interaction of 'age.binned' and 'power.binned' can also be present in the model formula.
#' 
#' An overview of the different penalty types and their usage can be found in the package vignette.
#' 
#' @seealso \code{\link{glmsmurf}}
#' 
#' @example /inst/examples/Rent_example3.R
p <- function(pred1, pred2 = NULL, pen = "lasso", refcat = NULL, group = NULL) {
  
  
  #######################
  # Penalty type
  
  pen.names <- c("none", "lasso", "grouplasso", "flasso", "gflasso", "2dflasso", "ggflasso")
  
  # Check if valid penalty name
  if (!(pen %in% pen.names)) {
    stop("Invalid penalty name.")
  }
  
  
  #######################
  # Check input
 
  # Check if factor for certain penalty types
  if (pen %in% c("flasso", "gflasso", "ggflasso") & (!is.factor(pred1)) ) {
    
    stop(paste0("Predictor '", deparse(substitute(pred1)),"' needs to be given as a factor."))
  }
  
  # Check if numeric or factor for no penalty, Lasso and Group Lasso
  if (pen %in% c("none", "lasso", "grouplasso") & (!is.factor(pred1) & !is.numeric(pred1)) ) {
    
    stop(paste0("Predictor '", deparse(substitute(pred1)),"' needs to be given as a factor or numeric."))
  }
  
  
  #######################
  # 2dflasso 
  
  # Check if type is "2dflasso" when second predictor is given
  if (!is.null(pred2) & pen != "2dflasso") {
    stop("Two predictors are given, penalty type should be '2dflasso'.")
  }
  
  # Check if second predictor is given when type is "2dflasso"
  if (is.null(pred2) & pen == "2dflasso") {
    stop("No second predictor is given for penalty type '2dflasso'.")
  }
  
  # Check for factors if 2dflasso
  if (pen == "2dflasso" & (!is.factor(pred1) | !is.factor(pred2)) ) {
    
    stop(paste0("Predictors '", deparse(substitute(pred1)),"' and '", 
                deparse(substitute(pred2)),"' in the interaction need to be given as factors."))
  }
  
  if (pen == "2dflasso") {
    
    # Change object 
    x <- interaction(pred1, pred2, lex.order = TRUE)
    # Combine names of two covariates together with "interaction"
    cov.name <- paste0("inter.", deparse(substitute(pred1)), ".", deparse(substitute(pred2)))
    cov.names <- list(deparse(substitute(pred1)), deparse(substitute(pred2)))
    
  } else {
    
    # Name of covariate
    cov.name <- deparse(substitute(pred1))
    cov.names <- cov.name
  }
  
  #######################
  # group 
  
  # Check if group is properly provided
  if (!is.null(group)) {
    
    if (!is.numeric(group) | length(group) != 1) {
      stop("'group' needs to be a strictly positive integer of length 1 or NULL.")
    }
    
    if (is.nan(group) | is.infinite(group)) {
      stop("'group' needs to be a strictly positive integer or NULL.")
    }
    
    if (group <= 0 | group %% 1 !=0) {
      stop("'group' needs to be a strictly positive integer or NULL.")
    }
    # Convert to integer
    group <- as.integer(group)
    
  } else {
    group <- 0L
  }
  
  
  #######################
  # refcat 
  
  if (pen != "2dflasso") {
    
    # Keep object
    x <- pred1

    # Transform to ordinary factor when of class "ordered"
    if (is.ordered(x)) {
      
      class(x) <- class(x)[!class(x) == "ordered"]
      warning(paste0("Predictor '", deparse(substitute(pred1)),"' is transformed from an ordered factor to a factor."))
    }
    
    if (!is.factor(pred1) | pen %in% c("lasso", "grouplasso")) {
      # Set reference category to zero when the predictor is not a factor, or the penalty type is 
      # Lasso or Group Lasso
      refcat <- 0L
      
    } else {
      
      if (is.null(refcat)) {
        # First level will be used as reference level
        refcat <- 1L
        
      } else {
      
        if (length(refcat) != 1) {
          stop(paste0("The 'refcat' argument for the predictor '", substitute(pred1), "' should be NULL or have length 1."))
        }
        
        # Index of provided reference level
        ind.refcat <- which(levels(pred1) == as.character(refcat))
        
        # Give error if provided level is not one of the levels of pred1
        if (length(ind.refcat) == 0) {
          stop(paste0("'", as.character(refcat), "' is not a level of the predictor '", substitute(pred1), "'."))
        } 
        
        # Return index of provided reference level  
        refcat <- ind.refcat
      }
    }
    
  } else {
    
    if (!is.null(refcat)) {
      warning("The 'refcat' argument is ignored for 'pen = \"2dflasso\"'.")
    }
    refcat <- 1L
  } 
  
  #######################
  
  # Add penalty and covariate name as attribute
  attributes(x) <- c(attributes(x), list(penalty = pen, cov.name = cov.name, cov.names = cov.names, 
                                         refcat = refcat, group = group))
  
  return(x)
}
