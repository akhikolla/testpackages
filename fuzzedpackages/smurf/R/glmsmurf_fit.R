###############################################
#
# Fitting function with given model matrix
#
###############################################


#' @export 
#' 
#' @rdname glmsmurf
#' 
#' @param X Only for \code{glmsmurf.fit}: the design matrix including ones for the intercept. A \code{n} by \code{(p+1)} matrix which can
#'           be of numeric matrix class (\code{\link[methods:StructureClasses]{matrix-class}}) or of class Matrix (\code{\link[Matrix]{Matrix-class}}) including sparse matrix class (\code{\link[Matrix]{dgCMatrix-class}}). 
#' @param y Only for \code{glmsmurf.fit}: the response vector, a numeric vector of size \code{n}. 
#' @param pen.cov Only for \code{glmsmurf.fit}: a list with the penalty type per predictor (covariate). A named list of strings with predictor names as element names. 
#'        Possible types: \code{"none"} (no penalty, e.g. for intercept), \code{"lasso"} (Lasso), \code{"grouplasso"} (Group Lasso), 
#'        \code{"flasso"} (Fused Lasso), \code{"gflasso"} (Generalized Fused Lasso), 
#'        \code{"2dflasso"} (2D Fused Lasso) or \code{"ggflasso"} (Graph-Guided Fused Lasso).  
#' @param n.par.cov Only for \code{glmsmurf.fit}: a list with the number of parameters to estimate per predictor (covariate). A named list of strictly positive integers with predictor names as element names.
#' @param group.cov Only for \code{glmsmurf.fit}: a list with the group of each predictor (covariate) which is only used for the Group Lasso penalty. 
#' A named list of positive integers with predictor names as element names where 0 means no group. 
#' @param refcat.cov Only for \code{glmsmurf.fit}: a list with the number of the reference category in the original order of the levels of each predictor (covariate).
#' When the predictor is not a factor or no reference category is present, it is equal to 0. This number will only be taken into account for a Fused Lasso, Generalized Fused Lasso or Graph-Guided Fused Lasso penalty 
#' when a reference category is present. 
glmsmurf.fit <- function(X, y, weights, start, offset, family, pen.cov, n.par.cov, group.cov, refcat.cov,
                         lambda, lambda1 = 0, lambda2 = 0, pen.weights, adj.matrix, standardize = TRUE, control = list(), 
                         formula = NULL, data = NULL, x.return = FALSE, y.return = FALSE, pen.weights.return = FALSE) {
  
  
  ######################################
  # Input checking
  
  
  ###########
  # y and X
  
  if (!is.numeric(y)) {
    stop("'y' should be a numeric vector.")
  }
  
  # Convert to vector if given as matrix
  if (!is.vector(y)) {
    y <- as.vector(y)
  }
  
  # Sample size
  n <- length(y)
  
  # Check if X is of class "Matrix"
  if (!(class(X)[1] %in% c("Matrix", "dgeMatrix", "dgCMatrix"))) {
    
    # Check if X is numeric
    if (!is.numeric(X)) {
      stop("'X' should be a numeric matrix or of class 'Matrix'.")
      
    } 
  } 
  
  if (nrow(X) != n) {
    stop(paste0("The number of rows of 'X' should be equal to the length of 'y': ", n, "."))
  }
  
  # Check if intercept is present. This is done by checking if at least one column contains only ones.
  # We only need to look at columns with sum nrow(X), which speeds up calculations and saves memory
  inter <- any(apply(X[, which(colSums(X) == nrow(X)), drop = FALSE], 
                     2L, function(x) all(x == 1)))
  if (!inter) {
    stop("An intercept needs to be included in the model.")
  }
  # Number of coefficients (excluding intercept)
  p <- ncol(X) - inter
  
  
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
  
  
  ###########
  # GLM weights
  
  weights <- .check_input_weights(weights, n = n) 
  
  
  ###########
  # Control parameters
  
  control <- do.call("glmsmurf.control", control)

  
  ###########
  # Starting value
  
  start <- .check_input_start(start, y = y, family = family, weights = weights, n = n, p = p, inter = inter)

  
  ###########
  # Offset
  
  offset <- .check_input_offset(offset, n = n)
 
  
  ###########
  # n.par.cov and pen.cov
  
  # Number of covariates
  n.cov <- length(pen.cov)
  
  # Check input for pen.cov
  .check_input_pen.cov(pen.cov)
  
  # Check input for n.par.cov
  n.par.cov.list <- .check_input_n.par.cov(n.par.cov, pen.cov = pen.cov)
  n.par.cov <- n.par.cov.list$n.par.cov
  n.par.cov.2dflasso <- n.par.cov.list$n.par.cov.2dflasso

  
  ###########
  # group.cov
  
  group.cov <- .check_input_group.cov(group.cov, pen.cov = pen.cov)
  
  
  ###########
  # refcat.cov
  
  .check_input_refcat.cov(refcat.cov, n.cov = n.cov)

  
  ###########
  # lambda
  
  lambda <- .check_input_lambda(lambda)
  
  # Indicator for out-of-sample selection of lambda
  oos <- lambda %in% c("oos.dev", "oos.mse", "oos.dss")
  
  # Validation data for out-of-sample selection
  valdata <- NULL
  
  if (oos) {
    # Out-of-sample selection of lambda
    
    # Indices of rows corresponding to validation sample
    validation.index <- control$validation.index
    
    if (is.null(validation.index)) {
      
      # No indices are provided, select them randomly based on oos.prop
      validation.index <- sample(1:n, size = n * control$oos.prop, replace = FALSE)
      
    } else {
      
      if (length(validation.index) >= n) {
        stop(paste0("'validation.index' should have length at most ", n, "."))
      }
      
    }
    
    # Indices of rows corresponding to training sample
    training.index <- (1:n)[which(!((1:n) %in% validation.index))]

    if (length(training.index) < 1) {
      stop("At least one observation should be present in the training data.")
    }
    
    if (length(validation.index) < 1) {
      stop("At least one observation should be present in the validation data.")
    }
    
    # Validation sample, save in valdata
    valdata$X.oos <- X[validation.index, , drop = FALSE]
    valdata$y.oos <- y[validation.index]
    valdata$weights.oos <- weights[validation.index]
    valdata$offset.oos <- offset[validation.index]
    
    # Check if all levels are present in training sample
    if (any(colSums(abs(X[training.index, , drop = FALSE])) < eps_num & colSums(abs(X)) > eps_num)) {
      stop(paste0("Some levels are missing in the training sample. Please provide different indices in 'validation.index' in the control object",
                  " or use a different (smaller) value for 'oos.prop' in the control object."))
    }
    
    # Training sample, reselect X, etc.
    X <- X[training.index, , drop = FALSE]
    y <- y[training.index]
    weights <- weights[training.index]
    offset <- offset[training.index]
    n <- length(training.index)
  }
  
 
  ###########
  # lambda1 and lambda2
  
  # Check input for lambda1 and lambda2
  .check_input_lambda12(lambda1 = lambda1, lambda2 = lambda2)
  
  
  # Change to list with lambda1 per covariate
  lambda1.orig <- lambda1
  lambda1 <- n.par.cov
  
  # Change to list with lambda2 per covariate
  lambda2.orig <- lambda2
  lambda2 <- n.par.cov
  
  for (j in 1:length(n.par.cov)) {
    
    # lambda1 and lambda2 are only needed for these penalty types
    if (pen.cov[[j]] %in% c("flasso", "gflasso", "ggflasso")) {
      lambda1[[j]] <- lambda1.orig
      lambda2[[j]] <- lambda2.orig
      
    } else {
      # Set lambda1 and lambda2 to 0 as they are not needed for this penalty type
      lambda1[[j]] <- lambda2[[j]] <- 0
    }
  }
  
  # lambda used to test for reference category
  lambda.test <- ifelse(is.character(lambda), 1, lambda)
  
  # Indicator for reference category for flasso, gflasso and ggflasso
  refcat <- (lambda.test * lambda1.orig == 0 & lambda.test * lambda2.orig == 0)

  
  ###########
  # formula
  
  # Check if valid formula when not NULL
  if (!is.null(formula)) {
    
    if (!inherits(formula, "formula")) {
      stop("'formula' needs to be a formula object.")
    }
  }
  
  
  ###########
  # Step size
  
  # Set step size to 0.1 times sample size if NULL
  if (is.null(control$step)) {
    control$step <- 0.1 * n
  }
  
  
  ###########
  # adj.matrix
  
  # Split column names of X based on n.par.cov, only used for ggflasso
  col.names.split <- split(colnames(X), rep(1:n.cov, n.par.cov))
  
  # Check input for adjacency matrix
  adj.matrix <- .check_input_adj(adj.matrix, pen.cov = pen.cov, n.par.cov = n.par.cov, refcat.cov = refcat.cov, 
                                 refcat = refcat, col.names.split = col.names.split)
  
  
  ######################################
  # Standardize
  
  # Indices of predictors with Lasso or Group Lasso penalty
  ind.stand <- which(pen.cov %in% c("lasso", "grouplasso"))
  
  # Standardize X matrix for predictors with Lasso or Group Lasso penalty (if standardize = TRUE)
  list.stand <- .X.stand(X = X, standardize = standardize, ind.stand = ind.stand, n.par.cov = n.par.cov, weights = weights)

  # Possibly standardized X matrix
  X <- list.stand$X
  # Set attribute indicating if the predictors for Lasso and Group Lasso are standardized
  attr(X, "standardize") <- standardize
  # Add X.means and X.sds as attributes of X
  attr(X, "X.means") <- list.stand$X.means
  attr(X, "X.sds") <- list.stand$X.sds

  
  ######################################
  # Penalty matrices
  
  pen.mat <- pen.cov
  pen.mat.transform <- pen.cov
  
  # Compute penalty matrices based on penalty type
  for (j in 1:n.cov) {
    pen.mat[[j]] <- switch(pen.cov[[j]],
                           none = as.matrix(0),
                           lasso = .pen.mat.lasso(n.par.cov[[j]]),
                           grouplasso = .pen.mat.grouplasso(n.par.cov[[j]]),
                           # Provide number of level that is reference category if present and 0 otherwise
                           flasso = .pen.mat.flasso(n.par.cov[[j]], refcat = ifelse(refcat, refcat.cov[[j]], refcat)),
                           # No change needed if reference category is not the first level since order does not matter!
                           gflasso = .pen.mat.gflasso(n.par.cov[[j]], refcat = refcat),
                           "2dflasso" = if(length(n.par.cov.2dflasso[[j]]) > 1) {
                             # Second and third element are number of levels from 1D predictors minus 1
                             .pen.mat.2dflasso(n.par.cov.2dflasso[[j]][2], n.par.cov.2dflasso[[j]][3])
                           } else {
                             .pen.mat.2dflasso(sqrt(n.par.cov.2dflasso[[j]]), sqrt(n.par.cov.2dflasso[[j]]))
                           },
                           ggflasso = .pen.mat.ggflasso(adj.matrix = get(names(pen.cov)[j], adj.matrix), 
                                                        # Level names contain covariate name, remove this first
                                                        lev.names = gsub(names(pen.cov)[j], "", unlist(col.names.split[j])), 
                                                        refcat = ifelse(refcat, refcat.cov[[j]], refcat))
    )
  }
  
  
  ######################################
  # Penalty weights
  
  
  # Check input
  pen.weights <- .check_input_pen_weights(pen.weights = pen.weights, pen.cov = pen.cov, pen.mat = pen.mat, n.cov = n.cov)
  
  
  # Compute penalty weights if needed (when character describing the method to compute them)
  if (is.character(pen.weights)) {
    
    L <- .compute_pen_weights(pen.weights = pen.weights, pen.cov = pen.cov, n.par.cov = n.par.cov, group.cov = group.cov,
                              pen.mat = pen.mat, data = data, X = X, y = y, weights = weights, offset = offset, family = family, formula = formula,
                              standardize = standardize, ind.stand = ind.stand, refcat = refcat, n = n,
                              lambda1 = lambda1, lambda1.orig = lambda1.orig, lambda2 = lambda2, lambda2.orig = lambda2.orig,
                              control = control)
    
    # Get computed penalty weights
    pen.weights <- L$pen.weights
    # Get list with lambda1 multiplied with penalty weights (for Lasso penalty) per predictor 
    lambda1 <- L$lambda1
    # Get list with lambda2 multiplied with penalty weights (for Group Lasso penalty) per predictor 
    lambda2 <- L$lambda2
    
  } else {
    
    # Show warning when lambda1.orig > 0 or lambda2.orig > 0 and penalty weights are inputted manually by the user
    if (lambda1.orig > 0 | lambda2.orig > 0) {
      warning("When lambda1 > 0 or lambda2 > 0, it is more reliable to let 'glmsmurf' or 'glmsmurf.fit' compute the penalty weights instead of inputting them manually via 'pen.weights'. ")
    }
  }
  
  
  
  # Add original lambda1 as final element (not done earlier to avoid problems with weights)
  lambda1$lambda1.orig <- lambda1.orig
  # Add original lambda2 as final element
  lambda2$lambda2.orig <- lambda2.orig
  
  
  ###########  
  # Multiply penalty matrices with penalty weights
  
  for (j in 1:n.cov) {
    pen.mat[[j]] <- pen.mat[[j]] * pen.weights[[j]]
    
    # Compute inverse of pen.mat for ordinal predictors
    if (pen.cov[[j]] == "flasso") {
      pen.mat.transform[[j]] <- t(upper.tri(pen.mat[[j]], diag=TRUE) * (0 * pen.mat[[j]] + 1) / pen.weights[[j]])
      
    } else {
      pen.mat.transform[[j]] <- 0
    }
    
  }
  
  # Compute eigenvalue decomposition of t(pen.mat[[j]]) %*% pen.mat[[j]] for all penalty types
  # except "none", "lasso" and "grouplasso"
  pen.mat.aux <- .pen.mat.eig(pen.mat = pen.mat, pen.cov = pen.cov)
  
  
  ######################################
  # Selection of lambda using cross-validation or in-sample
  
  lambda.orig <- lambda
  
  if (is.character(lambda.orig)) {
    
    # Print info
    if (control$print) {
      print("Selecting lambda")
    }
    
    lambda.select.list <- .select.lambda(X = X, y = y, weights = weights, start = start, offset = offset, family = family, 
                                         pen.cov = pen.cov, n.par.cov = n.par.cov, group.cov = group.cov, refcat.cov = refcat.cov, 
                                         lambda = lambda, lambda1 = lambda1, lambda2 = lambda2, 
                                         pen.mat = pen.mat, pen.mat.aux = pen.mat.aux, pen.mat.transform = pen.mat.transform, 
                                         valdata = valdata, control = control)
    
    # Optimal value of lambda
    lambda <- lambda.select.list$lambda
    
    # Print info
    if (control$print) {
      print("Lambda selected")
    }
  }
  
  ######################################
  # Run algorithm with given value of lambda or selected value of lambda
  
  # Print info
  if (control$print) {
    print(paste0("Running algorithm with inputted or selected value of lambda: ", round(lambda, 6)))
  }
  
  # Call auxiliary function to fit model where lambda is either the input value or 
  # the selected value based on cross-validation or in-sample selection.
  fit <- .glmsmurf.fit.internal(X = X, y = y, weights = weights, start = start, offset = offset, family = family,
                                pen.cov = pen.cov, n.par.cov = n.par.cov, group.cov = group.cov, refcat.cov = refcat.cov,
                                lambda = lambda, lambda1 = lambda1, lambda2 = lambda2, pen.mat = pen.mat, pen.mat.aux = pen.mat.aux, 
                                x.return = x.return, y.return = y.return, control = control, full.output = TRUE)
  
  # Print info
  if (control$print) {
    print("Algorithm finished")
  }
  
  if (pen.weights.return) {
    # Add pen.weights
    fit$pen.weights <- pen.weights
  }
  
  # Add n.par.cov
  fit$n.par.cov <- n.par.cov
  
  # Add pen.cov
  fit$pen.cov <- pen.cov
  
  # Add group.cov
  group.cov <- as.list(group.cov)
  names(group.cov) <- names(pen.cov)
  fit$group.cov <- group.cov
  
  # Add refcat.cov
  fit$refcat.cov <- refcat.cov

  # Add control
  fit$control <- control
  
  # Add results from cross-validation if performed
  if (is.character(lambda.orig)) {
    fit$lambda.method <- lambda.orig
    fit$lambda.vector <- lambda.select.list$lambda.vector
    fit$lambda.measures <- lambda.select.list$lambda.measures
    fit$lambda.coefficients <- lambda.select.list$lambda.coef
  }
  
  # Make fit object of class glmsmurf which (partially) extends the classes glm and lm
  class(fit) <- c("glmsmurf", "glm", "lm", class(fit))
  return(fit)
}
