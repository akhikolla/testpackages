###############################################
#
# Main fitting function
#
###############################################


#' @export 
#' @title Fit a Multi-Type Regularized GLM Using the SMuRF Algorithm
#' 
#' @description SMuRF algorithm to fit a generalized linear model (GLM) with multiple types of predictors via regularized maximum likelihood.
#' \code{glmsmurf.fit} contains the fitting function for a given design matrix.
#' 
#' @param formula A \code{\link[stats]{formula}} object describing the model to be fitted. 
#' Penalties are specified using the \code{\link{p}} function. For \code{glmsmurf.fit} this is an optional argument which is only used when penalty weights are computed using a generalized additive model (GAM).
#' @param data A data frame containing the model response and predictors for \code{n} observations.
#' @param family A \code{\link[stats]{family}} object specifying the error distribution and link function for the model. 
#' @param weights An optional vector of prior weights to use in the likelihood. It should be a numeric vector of length \code{n} (the number of observations),
#'  or \code{NULL}. When \code{NULL} or nothing is given, equal prior weights (all ones) will be used.
#' @param start A vector containing the starting values for the coefficients. It should either be a numeric vector 
#' of length \code{p+1} (with \code{p} the number of parameters excluding the intercept) or \code{NULL}. In the latter case, the link function applied to the weighted average of the response vector is used 
#' as starting value for the intercept and zero for the other coefficients.
#' @param offset A vector containing the offset for the model. It should be a vector of size \code{n} or NULL (no offset). 
#' Offset(s) specified using the \code{formula} object will be ignored!
#' @param lambda Either the penalty parameter, a positive number; or a string describing the method and measure used to select the penalty parameter:
#'               \itemize{
#'                  \item \code{"is.aic"} (in-sample; Akaike Information Criterion (AIC)),
#'                  \item \code{"is.bic"} (in-sample; Bayesian Information Criterion (BIC)),
#'                  \item \code{"is.gcv"} (in-sample; Generalized Cross-Validation (GCV) score),
#'                  \item \code{"oos.dev"} (out-of-sample; deviance), 
#'                  \item \code{"oos.mse"} (out-of-sample; Mean Squared Error (MSE)),
#'                  \item \code{"oos.dss"} (out-of-sample; Dawid-Sebastiani Score (DSS)),
#'                  \item \code{"cv.dev"} (cross-validation (CV); deviance),
#'                  \item \code{"cv.mse"} (CV; MSE),
#'                  \item \code{"cv.dss"} (CV; DSS),
#'                  \item \code{"cv1se.dev"} (CV with one standard error (SE) rule; deviance),
#'                  \item \code{"cv1se.mse"} (CV with one SE rule; MSE),
#'                  \item \code{"cv1se.dss"} (CV with one SE rule; DSS).
#'               }                 
#'               E.g. \code{"is.aic"} indicates in-sample selection of lambda with the AIC as measure.
#'               When \code{lambda} is missing or \code{NULL}, it will be selected using cross-validation with the one standard error rule and the deviance as measure (\code{"cv1se.dev"}).
#' @param lambda1 The penalty parameter for the \eqn{L_1}-penalty in Sparse (Generalized) Fused Lasso or Sparse Graph-Guided Fused Lasso is \eqn{\lambda \times \lambda_1}. 
#'                 A positive numeric with default 0 (no extra  \eqn{L_1}-penalty).
#' @param lambda2 The penalty parameter for the \eqn{L_2}-penalty in Group (Generalized) Fused Lasso or Group Graph-Guided Fused Lasso is \eqn{\lambda \times \lambda_2}. A positive numeric with default 0 (no extra \eqn{L_2}-penalty).
#' @param pen.weights Either a string describing the method to compute the penalty weights: 
#'  \itemize{
#'                  \item \code{"eq"} (default; equal penalty weights),
#'                  \item \code{"stand"} (standardization penalty weights),
#'                  \item \code{"glm"} (adaptive GLM penalty weights),
#'                  \item \code{"glm.stand"} (stand. ad. GLM penalty weights), 
#'                  \item \code{"gam"} (ad. GAM penalty weights),
#'                  \item \code{"gam.stand"} (stand. ad. GAM penalty weights);
#'               }                 
#'                or a list with the penalty weight vector per predictor. This list should have length equal to the number of predictors and predictor names as element names.                
#' @param adj.matrix A named list containing the adjacency matrices (a.k.a. neighbor matrices) for each of the predictors with a Graph-Guided Fused Lasso penalty. 
#'                The list elements should have the names of the corresponding predictors. If only one predictor has a Graph-Guided Fused Lasso penalty, 
#'                it is also possible to only give the adjacency matrix itself (not in a list).
#' @param standardize Logical indicating if predictors with a Lasso or Group Lasso penalty are standardized, default is \code{TRUE}.
#'                   The returned coefficients are always on the original (i.e. non-standardized) scale.
#' @param x.return Logical indicating if the used model matrix should be returned in the output object, default is \code{FALSE}.
#' @param y.return Logical indicating if the used response vector should be returned in the output object, default is \code{TRUE}.
#' @param pen.weights.return Logical indicating if the list of the used penalty weight vector per predictor should be returned in the output object, default is \code{FALSE}.
#' @param control A list of parameters used in the fitting process. This is passed to \code{\link{glmsmurf.control}}.
#' @return An object of class '\code{glmsmurf}' is returned. See \code{\link{glmsmurf-class}} for more details about this class and its generic functions.
#'         
#' @details See the package vignette for more details and a complete description of a use case.
#' 
#' As a user, it is important to take the following into acocunt:
#' \itemize{
#'  \item The estimated coefficients are rounded to 7 digits.
#'  \item The cross-validation folds are not deterministic. The validation sample for selecting lambda out-of-sample is determined at random when no indices are provided 
#'  in 'validation.index' in the control object argument. In these cases, the selected value of lambda is hence not deterministic. 
#'  When selecting lambda in-sample, or out-of-sample when indices are provided in 'validation.index' in the control object argument, the selected value of lambda is deterministic.
#'  \item The \code{glmsmurf} function can handle many use cases and is preferred for general use.
#'  The \code{glmsmurf.fit} function requires a more thorough understanding of the package internals and should hence be used with care!
#' }
#' 
#'          
#' @seealso \code{\link{glmsmurf-class}}, \code{\link{glmsmurf.control}}, \code{\link{p}}, \code{\link[stats]{glm}}
#' 
#' @references Devriendt, S., Antonio, K., Reynkens, T. and Verbelen, R. (2018). "Sparse Regression with Multi-type Regularized Feature Modeling." \emph{arXiv:1810.03136}.
#' 
#' Hastie, T., Tibshirani, R., and Wainwright, M. (2015). \emph{Statistical Learning with Sparsity: The Lasso and Generalizations}. CRC Press.
#' @example /inst/examples/Rent_example1.R
glmsmurf <- function(formula, family, data, weights, start, offset, lambda, lambda1 = 0, lambda2 = 0, 
                     pen.weights, adj.matrix, standardize = TRUE, control = list(), 
                     x.return = FALSE, y.return = TRUE, pen.weights.return = FALSE) {
  
  # Match call
  call <- match.call()

  ###########
  # Handle family object
  
  # If family given by name, convert to function with that name
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  
  # If family given as function, convert to object
  if (is.function(family)) {
    family <- family()
  }
  
  # Throw error if invalid family
  if (is.null(family$family)) {
    stop("Invalid family.")
  }

  
  #########
  
  # Terms
  terms <- terms(formula, specials = "p", data = data)
  # Model frame, add "p" as special function
  mf <- model.frame(formula = terms, data = data, drop.unused.levels = TRUE)
  
  # Phoney model frame without dropped levels to get penalty attributes (which can be dropped 
  # with drop.unused.levels = TRUE)
  mf_nodrop <- model.frame(formula = terms, data = data, drop.unused.levels = FALSE)
  
  if (ncol(mf) != ncol(mf_nodrop)) {
    stop("Something went wrong with the model frame.")
  }
  
  # Check if offset in formula
  if (!is.null(model.offset(mf))) {
    stop("No offset(s) can be given in the formula. Please use the 'offset' argument to specify the offset.")
  }
  
  
  #########
  # Prepare response vector (y)
  
  # Extract response
  y <- model.response(mf, type = "numeric")
  
  # Convert factor to numeric
  if (is.factor(y)) {
    y <- as.numeric(levels(y))[y]
  }
  
  
  #########
  # Prepare pen.cov, n.par.cov, group.cov and refcat.cov
  
  # Logical indicating if intercept is present
  inter <- attr(attr(mf, "terms"), "intercept")
  
  # -1 since y is not counted
  # Penalty type
  pen.cov <- vector("list", ncol(mf) + (inter - 1))
  # Number of parameters to estimate per covariate
  n.par.cov <- vector("list", ncol(mf) + (inter - 1))
  # Group of the covariate, where 0 means no group
  group.cov <- vector("list", ncol(mf) + (inter - 1))
  # Index of reference level (in original order of levels), where 0 means no reference category
  refcat.cov <- vector("list", ncol(mf) + (inter - 1))
  # Covariate names
  cov.names <- character(ncol(mf) + (inter - 1))
  # List of contrasts
  contrasts <- list()
  # List of levels
  xlevels <- list()
  
  # Intercept is present
  if (inter) {
    pen.cov[[1]] <- "none"
    n.par.cov[[1]] <- 1
    cov.names[1] <- "intercept"
    group.cov[1] <- 0
    refcat.cov[1] <- 0
    
  } else {
    stop("An intercept needs to present in the formula, please do not use '-1' in the formula.")
  }
  
  # lambda used to test for reference category in flasso, gflasso and ggflasso
  lambda.test <- .check_input_lambda(lambda)
  lambda.test <- ifelse(is.character(lambda.test), 1, lambda.test)

  # Check input for lambda1 and lambda2
  .check_input_lambda12(lambda1 = lambda1, lambda2 = lambda2)

  for (j in 2:ncol(mf_nodrop)) {
    
    ####
    # Covariate name
    if (is.null(attr(mf_nodrop[, j], "cov.name"))) {
      # No covariate name given, use column name of mf_nodrop
      cov.names[[j + (inter - 1)]] <- colnames(mf_nodrop)[j]
      
    } else {
      # Get covariate name
      cov.names[[j + (inter - 1)]] <- attr(mf_nodrop[, j], "cov.name")
    }
    
    
    ####
    # Penalty type
    if (is.null(attr(mf_nodrop[, j], "penalty"))) {
      # No penalty type given, set type to "none"
      
      pen.cov[[j + (inter - 1)]] <- "none"
      group.cov[[j + (inter - 1)]] <- 0
      refcat.cov[[j + (inter - 1)]] <- ifelse(is.factor(mf[, j]), 1, 0)
      
    } else {
      
      # Get penalty type
      pen.cov[[j + (inter - 1)]] <- attr(mf_nodrop[, j], "penalty")
      
      # Get group
      group.cov[[j + (inter - 1)]] <- attr(mf_nodrop[, j], "group")
      
      # Get index of reference category
      refcat.cov[[j + (inter - 1)]] <- attr(mf_nodrop[, j], "refcat")
    }
    
    
    # Check for unused factor levels
    if (is.factor(mf_nodrop[, j])) {
      if (!isTRUE(all.equal(levels(mf_nodrop[, j]), levels(mf[, j])))) {
        stop(paste0("No unused levels can be present in the predictor '", cov.names[j + (inter - 1)],"'.",
                    " Please remove the unused factor levels."))
      }
    }
    
    ####
    # Handle different penalty types
    if (pen.cov[[j + (inter - 1)]] == "none") {
      
      # Number of parameters to estimate
      n.par.cov[[j + (inter - 1)]] <- ifelse(is.factor(mf[, j]), length(levels(mf[, j])) - 1, 1)
      
      # Relevel to make chosen category the reference category
      if (refcat.cov[[j + (inter - 1)]] > 1) {
        mf[, j] <- relevel(mf[, j], ref = levels(mf[, j])[refcat.cov[[j + (inter - 1)]]])
      }
      
    } else if (pen.cov[[j + (inter - 1)]] %in% c("lasso", "grouplasso")) {
      # No reference category for Lasso and Group Lasso penalized predictors, 
      # this is imposed through contrast.args
      
      # Only use contrasts if a factor
      if (is.factor(mf[, j])) {
        name <- colnames(mf)[j]
        # Set contrasts to FALSE for predictor => no reference category
        # Note that the list element has the same name as the predictor
        contrasts[[name]] <- contrasts(mf[, j], contrasts = FALSE)
      }
      
      # Number of parameters to estimate, no reference category
      n.par.cov[[j + (inter - 1)]] <- ifelse(is.factor(mf[, j]), length(levels(mf[, j])), 1)
      
    } else if (pen.cov[[j + (inter - 1)]] %in% c("flasso", "gflasso", "ggflasso")) {
      
      if (lambda.test * lambda1 > 0 | lambda.test * lambda2 > 0) {
        # No reference category, this is imposed through contrast.args
        
        name <- colnames(mf)[j]
        # Set contrasts to FALSE for predictor => no reference category
        # Note that the list element has the same name as the predictor
        contrasts[[name]] <- contrasts(mf[, j], contrasts = FALSE)
        
        # Number of parameters to estimate
        n.par.cov[[j + (inter - 1)]] <- length(levels(mf[, j]))
        
      } else {
        
        # Relevel to make chosen category the reference category
        if (refcat.cov[[j + (inter - 1)]] > 1) {
          mf[, j] <- relevel(mf[, j], ref = levels(mf[, j])[refcat.cov[[j + (inter - 1)]]])
        }
        
        # Number of parameters to estimate
        n.par.cov[[j + (inter - 1)]] <- length(levels(mf[, j])) - 1
      }
      
    } else if (pen.cov[[j + (inter - 1)]] == "2dflasso") {
      
      # Check if 1D effects (or non-binned 1D effects) are also included
      error2d <- paste0(" should also be included in the model ",
                        "since this predictor (or a binned version) is present in an interaction. ",
                        "Please mind that 1D effects need to be included in the formula before 2D effects!")
      
      
      if (!(gsub(".binned*", "", attr(mf_nodrop[, j], "cov.names")[[1]]) %in% cov.names[(inter - 1) + (1:j)])) {
        stop(paste0("Variabele ", attr(mf_nodrop[, j], "cov.names")[[1]], error2d))
      }
      
      if (!(gsub(".binned*", "", attr(mf_nodrop[, j], "cov.names")[[2]]) %in% cov.names[(inter - 1) + (1:j)])) {
        stop(paste0("Variabele ", attr(mf_nodrop[, j], "cov.names")[[2]], error2d))
      }
      
      # Levels of pred1
      l1 <- levels(data[, attr(mf_nodrop[, j], "cov.names")[[1]]])
      # Levels of pred2
      l2 <- levels(data[, attr(mf_nodrop[, j], "cov.names")[[2]]])

      # Check that reference category is not changed for both 1D predictors (or their non-binned versions)
      n1 <- gsub(".binned*", "", attr(mf_nodrop[, j], "cov.names")[[1]])
      if (refcat.cov[[which(cov.names == n1)]] > 1) {
        stop(paste0("The reference category for the predictor '", n1,
                    " 'cannot be changed as it (or its binned version) is included in a 2D effect."))
      }

      n2 <- gsub(".binned*", "", attr(mf_nodrop[, j], "cov.names")[[2]])
      if (refcat.cov[[which(cov.names == n2)]] > 1) {
        stop(paste0("The reference category for the predictor '", n2,
                    " 'cannot be changed as it (or its binned version) is included in a 2D effect."))
      }

      
      # Change all observations with at least one of the 1D levels equal to the 1D reference level (l1[1] and l2[1])
      # to the 2D reference level (l1[1].l2[1]). This is done by renaming the corresponding levels
      l <- c(sapply(l2, function(x) paste0(l1[1], ".", x)), sapply(l1, function(x) paste0(x, ".", l2[1])))
      levels(mf[, j])[levels(mf[, j]) %in% l] <- paste0(l1[1], ".", l2[1])
      
      # Number of parameters that needs to be estimated for this covariate
      # and the number of levels for the 1D effects excluding the reference category
      n.par.cov[[j + (inter - 1)]] <- c(length(levels(mf[, j])) - 1, length(l1) - 1, length(l2) - 1)
      
      # Check for missing levels
      if (n.par.cov[[j + (inter - 1)]][1] != n.par.cov[[j + (inter - 1)]][2] * n.par.cov[[j + (inter - 1)]][3]) {
        stop(paste0("There are some missing levels in the interaction of '", attr(mf_nodrop[, j], "cov.names")[[1]],
                    "' and '", attr(mf_nodrop[, j], "cov.names")[[2]], "'.",
                    " Please make sure that there is at least one observation in each interaction category."))
      }
      
    } else {
      stop("Invalid penalty type.")
    }
    
    if (is.factor(mf[, j])) {
      # Save levels (needs to be done after possible releveling!)
      xlevels[[colnames(mf)[j]]] <- levels(mf[, j])
    }
    
  }
  
  # contrasts can be an empty list and should then be set to zero
  if (length(contrasts) == 0) {
    contrasts <- NULL
  }
  
  # Rename list elements
  names(pen.cov) <- names(n.par.cov) <- names(group.cov) <- names(refcat.cov) <- cov.names
  
  if (anyDuplicated(cov.names) > 0) {
    stop("A predictor can only be included once in the formula (except in a 2D Fused Lasso).")
  }
  
  #########
  # Prepare model matrix (X)
  
  # Make sparse model matrix
  X <- sparse.model.matrix(object = attr(mf, "terms"), data = mf, contrasts.arg = contrasts)
  # Rename columns
  colnames(X) <- .rename_mm.cols(colnames(X))
  
  
  #########
  # Call fitting function
  fit <- glmsmurf.fit(y = y, X = X, family = family, weights = weights, start = start, offset = offset,
                      n.par.cov = n.par.cov, pen.cov = pen.cov, group.cov = group.cov, refcat.cov = refcat.cov,
                      lambda = lambda, lambda1 = lambda1, lambda2 = lambda2,
                      pen.weights = pen.weights, adj.matrix = adj.matrix, standardize = standardize, formula = formula, data = data,
                      x.return = x.return, y.return = y.return, pen.weights.return = pen.weights.return, 
                      control = control)
  
  
  # Add call
  fit$call <- call
  
  # Add formula
  fit$formula <- formula

  # Add terms
  fit$terms <- terms
  
  # Add contrasts
  fit$contrasts <- contrasts
  
  # Add x-levels
  fit$xlevels <- xlevels
  
  return(fit)
}  




# Rename column names of model matrix
#
# col.names: Column names of matrix
.rename_mm.cols <- function(col.names) {
  
  for (j in 1:length(col.names)) {
    
    # Boolean for interaction
    inter_bool <- (length(grep("\"2dflasso\"", col.names[j])) > 0)
    
    # Change "(Intercept)" to Intercept
    col.names[j] <- gsub("\\(Intercept\\)", "Intercept", col.names[j])

    # String without everything within p(), i.e. factor label
    factor_label <- gsub("^p\\s*\\([^\\)]+\\)", "", col.names[j])
    
    # Everything within p() (including "p(" and ")")
    tmp <- regexpr("^p\\s*\\([^\\)]+\\)", col.names[j])
    
    # Only if match, otherwise keep column name
    # No match happens when a variable is added without a penalty (e.g. intercept)
    if (tmp > 0) {
      
      # Everything within p() (including "p(" and ")")
      col.names[j] <- regmatches(col.names[j], tmp)
      
      # Remove "p(" in beginning of string
      col.names[j] <- gsub("^p\\s*\\(+", "", col.names[j])
      
      # Remove 'pred1 = " "' with optional spaces
      col.names[j] <- gsub("\\s*pred1\\s*=\\s*", "", col.names[j])
      
      # Remove 'pred2 = " "' with optional spaces
      col.names[j] <- gsub("\\s*pred2\\s*=\\s*", "", col.names[j])
      
      # Remove ', pen = " "' with optional spaces and any character string or numeric (or space) for pen
      col.names[j] <- gsub(",*\\s*pen\\s*=\\s*\"[[:digit:][:alpha:]]*\"", "", col.names[j])
      
      # Remove ', group = ' with optional spaces and any or numeric or NULL (or space) for group
      col.names[j] <- gsub(",*\\s*group\\s*=\\s*[[:digit:][:alpha:]]*", "", col.names[j])
      
      # Remove ', refcat = ' with optional spaces and optional quotes and any or numeric or NULL (or space) for refcat
      col.names[j] <- gsub(",*\\s*refcat\\s*=\\s*\"*[[:digit:][:alpha:]]*\"*", "", col.names[j])
      
      # Remove ")" at the end
      col.names[j] <- gsub("\\)$", "", col.names[j])
      
      if (inter_bool) {
        
        # Remove "\"2dflasso\"" if still present (which happens when no "pen=" is used)
        col.names[j] <- gsub(",\\s*\"2dflasso\"", "", col.names[j])
        
        # Combine two predictor names
        col.names[j] <- gsub(",\\s*", ".", col.names[j])
        
        # Add "inter."
        col.names[j] <- paste0("inter.", col.names[j])
        
      } else {
        
        # Remove optional commas
        col.names[j] <- gsub(",\\s*", "", col.names[j])
      }
      
      # Combine with factor label since this was removed
      col.names[j] <- paste0(col.names[j], factor_label) 
    }
     
  }

  return(col.names)
}
