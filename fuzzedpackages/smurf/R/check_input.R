###############################################
#
# Functions for input checking
#
###############################################

# Note that input checking for pen.weights is in Penalty_weights.R


# Check input for weights
#
# weights: Vector of prior weights
# n: Sample size
.check_input_weights <- function(weights, n) {
  
  if (missing(weights)) {
    weights <- NULL
  }
  
  # If weights are not given or NULL, set to default (all ones)
  if (is.null(weights)) {
    weights <- rep(1, n)
    
  } else {
    
    # Check if weights are numeric
    if (!is.numeric(weights)) {
      stop("'weights' must be a numeric vector or NULL.")
    }
    
    # Check if any is NaN or Inf
    if (any(is.nan(weights)) | any(is.infinite(weights))) {
      stop("'weights' must be a vector of finite numbers.")
    }
    
    # Check if all weights are non-negative
    if (any(weights < 0)) {
      stop("All weights need to be non-negative.")
    }
    
    # Check if weights have same length as y
    if (length(weights) != n) {
      stop(paste0("'weights' must be a numeric vector of length ", n, " or NULL."))
    }
    
    # Convert to vector if given as matrix
    if (!is.vector(weights)) {
      weights <- as.vector(weights)
    }
  }
  
  return(weights)
}


# Check input for start
#
# start: vector of starting values for the coefficients
# y: Response vector
# weights: Vector of prior weights
# family: Family object
# n: Sample size
# p: Number of coefficients excluding intercept
# inter: Logical indicating if an intercept is present
.check_input_start <- function(start, y, weights, family, n, p, inter) {
  
  if (missing(start)) {
    start <- NULL
  }
  
  # If start is not given or NULL, use (linkfun(mean(y * weights) / mean(weights)),0,...,0)
  if (is.null(start)) {
    
    start <- c(family$linkfun(weighted.mean(y, w = weights)), rep(0, p))
    
  } else {
    
    # Check if start is numeric
    if (!is.numeric(start)) {
      stop("'start' must be a numeric vector or NULL.")
    }
    
    # Check if any is NaN or Inf
    if (any(is.nan(start)) | any(is.infinite(start))) {
      stop("'start' must be a vector of finite numbers.")
    }
    
    # Check if start has length p+inter
    if (length(start) != p + inter) {
      stop(paste0("'start' must be a numeric vector of length ", p + inter, " or NULL."))
    }
    
    # Convert to vector if given as matrix
    if (!is.vector(start)) {
      start <- as.vector(start)
    }
  }
  
  return(start)
}


# Check input for offset
#
# offset: Vector containing the offset for the model
# n: Sample size
.check_input_offset <- function(offset, n) {
  
  if (missing(offset)) {
    offset <- NULL
  }
  
  # If offset is not given or NULL, set to default (all zeros)
  if (is.null(offset)) {
    offset <- rep(0, n)
    
  } else {
    # Check if offset is numeric
    if (!is.numeric(offset)) {
      stop("'offset' must be a numeric vector or NULL.")
    }
    
    # Check if any is NaN or Inf
    if (any(is.nan(offset)) | any(is.infinite(offset))) {
      stop("'offset' must be a vector of finite numbers.")
    }
    
    # Check if offset has same length as y
    if (length(offset) != n) {
      stop(paste0("'offset' must be a numeric vector of length ", n, " or NULL."))
    }
    
    # Convert to vector if given as matrix
    if (!is.vector(offset)) {
      offset <- as.vector(offset)
    }
  
  }
  
  return(offset)
}


# Check input for pen.cov
#
# pen.cov: List with penalty type per predictor
.check_input_pen.cov <- function(pen.cov) {
  
  # Check if pen.cov is a character list
  if (!is.list(pen.cov)) {
    stop("'pen.cov' should be a list.")
  }
  
  if (any(!sapply(pen.cov, is.character))) {
    stop("'pen.cov' should be a list of characters.")
  }
  
  # Check if correct types of penalties
  pen.cov.vec <- c("none", "lasso", "grouplasso", "flasso", "gflasso", "2dflasso", "ggflasso")
  if (any(!(pen.cov %in% pen.cov.vec))) {
    stop(paste0("'pen.cov' should be a list of characters with possible elements: ", 
                paste(sprintf("'%s'", pen.cov.vec[-length(pen.cov.vec)]), collapse=", "), 
                " or ", sprintf("'%s'", pen.cov.vec[length(pen.cov.vec)]), "."))
  }
  
}


# Check input for n.par.cov
#
# n.par.cov: List with number of parameters to estimate per predictor
# pen.cov: List with penalty type per predictor
.check_input_n.par.cov <- function(n.par.cov, pen.cov) {
  
  # Check if n.par.cov is a numeric list of integers
  if (!is.list(n.par.cov)) {
    stop("'n.par.cov' should be a list.")
  }
  
  if (any(!sapply(n.par.cov, is.numeric))) {
    stop("'n.par.cov' should be a list of positive integers.")
  }
  
  if (any(sapply(n.par.cov, function(x) any(x <= 0 |  x %% 1 != 0) ))) {
    stop("'n.par.cov' should be a list of positive integers.")
  }
  
  # Check if n.par.cov and pen.cov have the same length
  if (length(n.par.cov) != length(pen.cov) ) {
    stop("'n.par.cov' and 'pen.cov' should be lists of equal length.")
  }
  
  # Index of predictors with 2D fused Lasso penalty
  ind2d <- which(pen.cov == "2dflasso") 
  
  if (length(ind2d) > 0) {
    # Remove extra parts of n.par.cov for 2dflasso
    
    # Save n.par.cov with extra parts as n.par.cov.2dflasso
    n.par.cov.2dflasso <- n.par.cov
    
    for (j in ind2d) {
      
      # More than 1 element
      if(length(n.par.cov[[j]]) != 1) {
        
        if(length(n.par.cov[[j]]) != 3) {
          stop(paste0("Element ", j ," of n.par.cov needs to have length 1 or 3."))
        }
        
        # Keep only first element (and hence remove extra parts)
        n.par.cov[[j]] <- n.par.cov[[j]][1]
        
        # Check if first element is product of second and third element
        if (n.par.cov[[j]][1] != n.par.cov.2dflasso[[j]][2] * n.par.cov.2dflasso[[j]][3]) {
          stop(paste0("The first element of n.par.cov[[", j, "]] should be the product of the second ",
                      "and the third element of n.par.cov[[", j, "]]."))
        }
      }
    }
    
  } else {
    n.par.cov.2dflasso <- NULL
  }
  
  # Index of predictors without 2D fused Lasso penalty
  ind.non.2d <- which(pen.cov != "2dflasso") 
  
  if (length(ind.non.2d) > 0) {
    
    for (j in ind.non.2d) {
      
      if(length(n.par.cov[[j]]) != 1) {
        stop(paste0("Element ", j ," of n.par.cov needs to have length 1."))
      }
    }
  }
  
  return(list(n.par.cov = n.par.cov, n.par.cov.2dflasso = n.par.cov.2dflasso))
}


# Check input for group.cov
#
# group.cov: List with group of each predictor which is only used for the Group Lasso penalty, 0 means no group
# pen.cov: List with penalty type per predictor
.check_input_group.cov <- function(group.cov, pen.cov) {
  
  # Check if group.cov is a numeric list of integers
  if (!is.list(group.cov)) {
    stop("'group.cov' should be a list.")
  }
  
  if (any(!sapply(group.cov, is.numeric))) {
    stop("'group.cov' should be a list of positive integers.")
  }
  
  if (any(sapply(group.cov, function(x) x < 0 |  x %% 1 != 0 ))) {
    stop("'group.cov' should be a list of positive integers.")
  }
  
  # Check if group.cov and pen.cov have the same length
  if (length(group.cov) != length(pen.cov)) {
    stop("'group.cov' and 'pen.cov' should be lists of equal length.")
  }
  
  # Set group.cov to 0 when no Group Lasso penalty
  group.cov[pen.cov != "grouplasso"] <- 0L
  
  # Unlist group.cov
  group.cov <- unlist(group.cov)
  
  return(group.cov)
}


# Check input for refcat.cov
#
# refcat.cov: List with number of the reference category in the original order of the levels of each predictor.
# When the predictor is not a factor or no reference category is present, it is equal to 0. This number will only be taken into account for a Fused Lasso, Generalized Fused Lasso or Graph-Guided Fused Lasso penalty 
# when a reference category is present.
# n.cov: Number of covariates (including intercept), i.e. p+1
.check_input_refcat.cov <- function(refcat.cov, n.cov) {
  
  # Check if refcat.cov is a numeric list of integers
  if (!is.list(refcat.cov)) {
    stop("'refcat.cov' should be a list.")
  }
  
  if (any(!sapply(refcat.cov, is.numeric))) {
    stop("'refcat.cov' should be a list of positive integers.")
  }
  
  if (any(sapply(refcat.cov, function(x) x < 0 |  x %% 1 != 0 ))) {
    stop("'refcat.cov' should be a list of positive integers.")
  }
  
  # Check if refcat.cov and pen.cov have the same length
  if (length(refcat.cov) != n.cov) {
    stop("'refcat.cov' and 'pen.cov' should be lists of equal length.")
  }
}


# Check input for lambda
#
# lambda: penalty parameter
.check_input_lambda <- function(lambda) {
  
  # If lambda is not given, set to NULL
  if (missing(lambda)) {
    lambda <- NULL
  }
  
  # If lambda is NULL (or not given), set to default
  if (is.null(lambda)) {
    lambda <- "cv1se.dev"
    
  } else if (is.character(lambda)) {
    
    # Check if lambda is selected using one of the supported methods
    lambda_char_vec <- c("cv.dev", "cv.mse", "cv.dss",
                         "cv1se.dev", "cv1se.mse", "cv1se.dss",
                         "is.aic", "is.bic", "is.gcv",
                         "oos.dev", "oos.mse", "oos.dss")
    
    if (!(lambda %in% lambda_char_vec)) {
      stop(paste0("'lambda' must be a numeric, NULL, or one of ",  
                  paste(sprintf("'%s'", lambda_char_vec[-length(lambda_char_vec)]), collapse=", "), 
                  " or ", sprintf("'%s'", lambda_char_vec[length(lambda_char_vec)]), "."))
    }
    
  } else {
    
    # Check if lambda is numeric
    if (!is.numeric(lambda)) {
      stop("'lambda' must be numeric.")
    }
    
    # Check if lambda has length 1
    if (length(lambda) != 1) {
      stop("'lambda' must be a numeric of length 1.")
    }
    
    # Check if NA, NaN or Inf
    if (is.na(lambda) | is.nan(lambda) | is.infinite(lambda)) {
      stop("'lambda' must be a finite number.")
    }
    
    # Check if lambda is non-negative
    if (lambda < 0) {
      stop("'lambda' must be non-negative.")
    }
  }
  
  return(lambda)
}


# Check input for lambda1 and lambda2
#
# lambda1: The penalty parameter for the L_1-penalty in Sparse (Generalized) Fused Lasso or Sparse Graph-Guided Fused Lasso is lambda*lambda1
# lambda2: The penalty parameter for the L_2-penalty in Group (Generalized) Fused Lasso or Group Graph-Guided Fused Lasso is lambda*lambda2
.check_input_lambda12 <- function(lambda1, lambda2) {
  
  ###########
  # lambda1 and lambda2
  
  # Check if lambda1 is numeric
  if (!is.numeric(lambda1)) {
    stop("'lambda1' must be numeric.")
  }
  
  # Check if lambda1 has length 1
  if (length(lambda1) != 1) {
    stop("'lambda1' must be a numeric of length 1.")
  }
  
  # Check if NA, NaN or Inf
  if (is.na(lambda1) | is.nan(lambda1) | is.infinite(lambda1)) {
    stop("'lambda1' must be a finite number.")
  }
  
  # Check if lambda1 is non-negative
  if (lambda1 < 0) {
    stop("'lambda1' must be non-negative.")
  }
  
  
  # Check if lambda2 is numeric
  if (!is.numeric(lambda2)) {
    stop("'lambda2' must be numeric.")
  }
  
  # Check if lambda2 has length 1
  if (length(lambda2) != 1) {
    stop("'lambda2' must be a numeric of length 1.")
  }
  
  # Check if NA, NaN or Inf
  if (is.na(lambda2) | is.nan(lambda2) | is.infinite(lambda2)) {
    stop("'lambda2' must be a finite number.")
  }
  
  # Check if lambda2 is non-negative
  if (lambda2 < 0) {
    stop("'lambda2' must be non-negative.")
  }
}


# Check input for adjacency matrix
#
# adj.matrix: Named list containing the adjacency matrices (a.k.a. neighbor matrices) for each of the predictors with a Graph-Guided Fused Lasso penalty 
# pen.cov: List with penalty type per predictor
# n.par.cov: List with number of parameters to estimate per predictor
# refcat.cov: List with number of the reference category in the original order of the levels of each predictor.
# When the predictor is not a factor or no reference category is present, it is equal to 0. This number will only be taken into account for a Fused Lasso, Generalized Fused Lasso or Graph-Guided Fused Lasso penalty 
# when a reference category is present. 
# refcat: Indicator for reference category for Fused Lasso, Generalized Fused Lasso and Graph-Guided Fused Lasso
# col.names.split: Vector with column names of model matrix X split using n.par.cov 
.check_input_adj <- function(adj.matrix, pen.cov, n.par.cov, refcat.cov, refcat, col.names.split) {
  
  # Default is empty list
  if (missing(adj.matrix)) {
    adj.matrix <- list()
  }
  
  # Indices of predictors with ggflasso penalty
  ind.adj <- which(pen.cov == "ggflasso")
  
  # Transform to list if adjacency matrix is not given as a list and only one predictor with ggflasso penalty
  if (length(ind.adj) == 1 & !is.list(adj.matrix)) {
    adj.matrix <- list(adj.matrix)
    names(adj.matrix)[1] <- names(pen.cov)[ind.adj]
  }
  
  # Check if adj.matrix is a list
  if (!is.list(adj.matrix)) {
    stop("'adj.matrix' needs to be a named list or NULL.")
  }
  
  # Number of list elements
  adj.length <- length(adj.matrix)
  

  
  if (length(ind.adj) > 0) {
    
    if (length(ind.adj) != adj.length) {
      stop(paste0("The number of elements of 'adj.matrix' needs to be the same as the number of ",
                  "predictors that is penalized using the Graph-Guided Fused Lasso, i.e. ", length(ind.adj), "."))
    }
    
    for (i in 1:adj.length) {
      
      # Check if list element of adj.matrix have correct names
      if (is.null(names(adj.matrix)[i])) {
        stop("Element ", i, " of 'adj.matrix' needs to have name '", names(pen.cov)[ind.adj[i]], "'.")
      }
      if (names(adj.matrix)[i] != names(pen.cov)[ind.adj[i]]) {
        stop("Element ", i, " of 'adj.matrix' needs to have name '", names(pen.cov)[ind.adj[i]], "'.")
      }
      
      # Check if list element of adj.matrix has correct class
      if (!(class(adj.matrix[[i]])[1] %in% c("matrix", "Matrix", "dgeMatrix", "dgCMatrix"))) {
        stop("A numeric matrix or element of class Matrix was expected in element ", i, " of 'adj.matrix'.")
      }
      
      # Check if list element of adj.matrix has correct number of rows (number of levels of the corresponding predictor)
      if ((n.par.cov[[ind.adj[i]]] + refcat) != nrow(adj.matrix[[i]])) {
        stop("Element ", i, " of 'adj.matrix' needs to be a square matrix with ", (n.par.cov[[ind.adj[i]]] + refcat), 
             " rows (the number of levels of the corresponding predictor).")
      }
      
      # Level names (without reference level if refcat=TRUE) of predictor
      lev.names <- gsub(names(pen.cov)[ind.adj[i]], "", unlist(col.names.split[ind.adj[i]]))
      # Row names of adjacency matrix
      spat.names <- rownames(adj.matrix[[i]])
      
      if (refcat) {
        # Remove reference category
        spat.names <- spat.names[-refcat.cov[[ind.adj[i]]]]
      }
      
      if (!isTRUE(all.equal(as.character(lev.names), spat.names))) {
        stop("The rownames of element ", i, " of 'adj.matrix' are not the same as the level names of the corresponding predictor.",
             " Note that the order of the names is also important.")
      }
    }
  }
  
  return(adj.matrix)
}

