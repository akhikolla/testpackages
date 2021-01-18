######################################################################
#
# Functions for input checking and computation of penalty weights
#
######################################################################


# Check input for pen.weights
#
# pen.weights: List with penalty weight vector per predictor
# pen.cov: List with penalty type per predictor
# pen.mat: List with penalty matrix per predictor
# n.cov: Number of covariates (including intercept), i.e. p+1
.check_input_pen_weights <- function(pen.weights, pen.cov, pen.mat, n.cov) {
  
  # If pen.weights are not given or NULL, determine them using a method
  if (missing(pen.weights)) {
    pen.weights <- NULL
  }
  
  # If pen.weights are not given or NULL, set to default
  if (is.null(pen.weights)) {
    pen.weights <- "eq"
    
  } else if (is.character(pen.weights)) {
    
    # Check if pen.weights is one of the supported methods
    pen.weights_char_vec <- c("eq", "stand", "glm", "glm.stand", "gam", "gam.stand")
    if (!(pen.weights %in% pen.weights_char_vec)) {
      stop(paste0("'pen.weights' must be a list of numerics, NULL, or one of ",  
                  paste(sprintf("'%s'", pen.weights_char_vec[-length(pen.weights_char_vec)]), collapse=", "), 
                  " or ", sprintf("'%s'", pen.weights_char_vec[length(pen.weights_char_vec)]), "."))
    }
    
  } else {
    
    # Check if pen.weights is a list of numeric vectors
    if (!is.list(pen.weights) | !all(sapply(pen.weights, is.numeric))) {
      stop("'pen.weights' must be a list of numeric vectors.")
    }
    
    if (any(sapply(pen.weights, function(x) any(is.na(x) | is.nan(x) | is.infinite(x))))) {
      stop("'pen.weights' must be a list of finite numeric vectors.")
    }
    
    # Check if number of elements of pen.weights is correct
    if (length(pen.weights) != n.cov) {
      stop(paste0("'pen.weights' must be a list of length ", n.cov, "."))
    }
    
    # Check if elements of pen.weights have correct length
    ind <- which(sapply(pen.weights, length) != sapply(pen.mat, nrow))
    
    # One element of pen.weights has wrong length
    if (length(ind) == 1) {
      stop(paste0("Element ", sprintf("'%s'", ifelse(is.null(names(ind)), ind, names(ind))), 
                  " of 'pen.weights' has the wrong length."))
      
      # Multiple elements of pen.weights have wrong length  
    } else if (length(ind) > 1) {
      
      # Give predictor names if not empty and number otherwise
      if (is.null(names(ind))) {
        tmp <- paste(as.character(ind), collapse=", ")
      } else {
        tmp <- paste(sprintf("'%s'", names(ind)), collapse=", ")
      }
      stop(paste0("Elements ", tmp, " of 'pen.weights' have the wrong length."))
      
    }
    
    
    
    # Check if elements of pen.weights have correct names
    ind <- which(names(pen.weights) != names(pen.cov))
    
    if (length(ind) > 0) {
      names(ind) <- names(pen.weights)[ind]
    }
    
    # One element of pen.weights has wrong name
    if (length(ind) == 1) {
      stop(paste0("Element ", sprintf("'%s'", ifelse(is.null(names(ind)), ind, names(ind))), 
                  " of 'pen.weights' should have name '", names(pen.cov)[ind], "'."))
      
      # Multiple elements of pen.weights have wrong names  
    } else if (length(ind) > 1) {
      
      # Give predictor names if not empty and number otherwise
      if (is.null(names(ind))) {
        tmp <- paste(as.character(ind), collapse=", ")
      } else {
        tmp <- paste(sprintf("'%s'", names(ind)), collapse=", ")
      }
      tmp2 <- paste(sprintf("'%s'", names(pen.cov)[ind]), collapse=", ")
      stop(paste0("Elements ", tmp, " of 'pen.weights' should have names ", tmp2, ", respectively."))
    }
  }
  
  return(pen.weights)
}


# Compute penalty weights if not given as input (check happens in fit function)
#
# pen.weights: List with penalty weight vector per predictor
# pen.cov: List with penalty type per predictor
# n.par.cov: List with number of parameters to estimate per predictor
# group.cov: List with group of each predictor which is only used for the Group Lasso penalty, 0 means no group
# pen.mat: List with penalty matrix per predictor
# data: Data frame containing the model response and predictors
# X: The design matrix including ones for the intercept
# y: Response vector
# weights: Vector of prior weights
# offset: Vector containing the offset for the model
# family: Family object
# standardize: Logical indicating if the columns of X corresponding to Lasso or Group Lasso penalties 
# ind.stand: Indices of predictors with a Lasso or a Group Lasso penalty  
# refcat: Indicator for reference category for Fused Lasso, Generalized Fused Lasso and Graph-Guided Fused Lasso
# n: Sample size
# lambda1: List with lambda1 per predictor. The penalty parameter for the L_1-penalty in Sparse (Generalized) Fused Lasso or Sparse Graph-Guided Fused Lasso is lambda*lambda1
# lambda2: List with lambda2 per predictor. The penalty parameter for the L_2-penalty in Group (Generalized) Fused Lasso or Group Graph-Guided Fused Lasso is lambda*lambda2
# lambda1.orig: Input value for lambda1
# lambda2.orig: Input value for lambda2
.compute_pen_weights <- function(pen.weights, pen.cov, n.par.cov, group.cov, pen.mat, 
                                 data, X, y, weights, offset, family, formula, standardize, ind.stand, refcat, n,
                                 lambda1, lambda2, lambda1.orig, lambda2.orig, control) {
  
  # Print info
  if (control$print) {
    print("Computing penalty weights")
  }
  
  # Save type of penalty weights
  pen.weights.type <- pen.weights
  
  # Logical indicating if standardization penalty weights are used
  pen.weights.stand <- (pen.weights.type %in% c("stand", "glm.stand", "gam.stand")) 
  
  pen.weights <- n.par.cov
  
  # Zero weight for non-penalized predictors, hence no need to compute weights
  # if only penalties of these types
  if (!all(pen.cov == "none")) {
    
    # Equal weights or only standardization weights
    if (pen.weights.type %in% c("eq", "stand")) {
      
      for (j in 1:length(pen.weights)) {
        pen.weights[[j]] <- rep(1, nrow(pen.mat[[j]]))
      }
      
      
    } else {
      
      if (pen.weights.type %in% c("glm", "glm.stand")) {
        
        # Check if full rank: this is not the case when Lasso or Group Lasso penalties are present or 
        # when Fused Lasso, Generalized Fused Lasso or Graph-Guided Fused Lasso penalties are present 
        # without reference category (i.e. lambda1>0 or lambda2>0)
        fullRank <- !(any(pen.cov %in% c("lasso", "grouplasso")) | 
                        (!refcat & any(pen.cov %in% c("flasso", "gflasso", "ggflasso"))))
        
        if (fullRank) {
          
          # Compute GLM coefficients used for penalty weights
          glm.fit <- speedglm.wfit(y = y, X = X, intercept = TRUE, family = family, weights = weights, offset = offset,
                                   start = NULL, eigendec = FALSE, 
                                   sparse = ifelse(class(X)[1] %in% c("dgCMatrix"), TRUE, NULL), 
                                   trace = FALSE)
          
          beta.weights <- glm.fit$coefficients
          
        } else {
          
          
          family_glmnet <- family
          
          # Use named families for glmnet whenever possible as this is faster
          if ((family$family == "binomial" & family$link == "logit") | 
              (family$family == "gaussian" & family$link == "identity") |
              (family$family == "poisson" & family$link == "log")) {
            family_glmnet <- family$family
          }
          
          
          # Compute GLM coefficients using GLM with small ridge penalty (alpha=0) to avoid problems with multicollinearity
          # Do not include column for intercept
          glm.fit <- glmnet(y = y, x = X[, -1L], family = family_glmnet, weights = weights, offset = offset,
                            nlambda = 2, alpha = 0, intercept = TRUE, standardize = FALSE)
          # Use coefficients obtained with smallest lambda and add fitted intercept
          beta.weights <- c(glm.fit$a0[2L], glm.fit$beta[, 2L])
          names(beta.weights)[1] <- colnames(X)[1]
        }
        
      } else {
        
        if (is.null(formula)) {
          stop("GAM penalty weights are only computed if a formula is given in 'glmsmurf'.")
        } 
        
        if (standardize & length(ind.stand) > 0) {
          # Note that the GAM model does not take standardization for Lasso and Group Lasso into account!
          warning("GAM penalty weights can be unreliable when 'standardize = TRUE', please use GLM penalty weights or equal penalty weights if possible.")
        }
        
        
        # Compute GAM coefficients used for penalty weights
        beta.weights <- .gam.coef.pen.weights(formula = formula, family = family, data = data, weights = weights, 
                                              offset = offset, n.par.cov = n.par.cov, refcat = refcat)
        
      }
      
      # Split coefficient vector
      beta.weights.split <- split(beta.weights, rep(1:length(n.par.cov), n.par.cov))
      
      # Convert to correct penalty weights
      for (j in 1:length(pen.weights)) {
        
        if (pen.cov[[j]] == "grouplasso") {
          
          if (group.cov[[j]] == 0) {
            # Only for predictors in group 0 (i.e. no group)
            pen.weights[[j]] <- as.numeric(1 / sqrt(sum((pen.mat[[j]] %*% beta.weights.split[[j]])^2)))
            
          } else {
            # Predictors not in group 0 are treated later
            pen.weights[[j]] <- 0
          }
          
        } else {
          
          pen.weights[[j]] <- as.numeric(1 / abs(pen.mat[[j]] %*% beta.weights.split[[j]]))
          
          # lambda1 and lambda2 are only needed for these penalty types
          if (pen.cov[[j]] %in% c("flasso", "gflasso", "ggflasso")) {
            # Multiply lambda1 with weights: inverse of absolute value of coefficients
            lambda1[[j]] <- lambda1[[j]] * as.numeric(1 / abs(beta.weights.split[[j]]))
            # Multiply lambda2 with weights: inverse of 2-norm of coefficients
            lambda2[[j]] <- lambda2[[j]] * as.numeric(1 / sqrt(sum(beta.weights.split[[j]]^2)))
          }
        }
      }
      
      #########
      # Correction for Group Lasso (with group != 0)
      
      # Unique group numbers (excluding zero)
      groups.unique.nz <- unique(group.cov[group.cov != 0])
      
      if (length(groups.unique.nz) > 0) {
        
        # Loop over groups (except group 0)
        for (j in 1:length(groups.unique.nz)) {
          # Indices of predictors in the group
          ind.group <- which(group.cov == groups.unique.nz[j])
          
          swn <- 0
          for (l in ind.group) {
            # Compute contribution of predictor to squared weighted norm (swn) and add to swn
            swn <- swn + sum((pen.mat[[l]] %*% beta.weights.split[[l]])^2)
          }
          
          # Convert to correct penalty weights
          pen.weights[ind.group] <- as.numeric(1 / sqrt(swn))
        }
      }
      
      #########
      
    } 
    
  } else {
    
    # No need for standardization since zero weight for non-penalized predictors
    pen.weights.stand <- FALSE
  }
  
  # Zero weight for non-penalized covariates
  pen.weights[pen.cov == "none"] <- 0
  
  if (!all(pen.cov == "none")) {
    
    tmp <- unlist(pen.weights)
    # Add penalty weights in lambda1 to tmp, divide by lambda1.orig to only consider the penalty weights (and not the weighted lambda1)
    if (lambda1.orig > 0) {
      tmp <- c(tmp, unlist(lambda1[pen.cov %in% c("flasso", "gflasso", "ggflasso")]) / lambda1.orig) 
    }
    
    # Add penalty weights in lambda2 to tmp, divide by lambda2.orig to only consider the penalty weights (and not the weighted lambda2)
    if (lambda2.orig > 0) {
      tmp <- c(tmp, unlist(lambda2[pen.cov %in% c("flasso", "gflasso", "ggflasso")]) / lambda2.orig)
    }
    
    # Transform penalty weights
    pen.weights <- Map("*", pen.weights, (length(tmp) - 1) / sum(tmp))
    # Transform lambda1 and lambda2
    lambda1 <- Map("*", lambda1, (length(tmp) - 1) / sum(tmp))
    lambda2 <- Map("*", lambda2, (length(tmp) - 1) / sum(tmp))
  }
  
  # Standardization penalty weights (Bondell and Reich, 2009)
  if (pen.weights.stand) {
    
    STAND.pen.weights <- pen.weights
    
    # Cumulative sum of number of levels per covariate (excluding reference category if this is present)
    ind.cov <- cumsum(n.par.cov)
    
    for (j in 1:length(n.par.cov)) {
      
      # Compute standardization factors, not for "none", "lasso" or "grouplasso"
      if (!(pen.cov[[j]] %in% c("none", "lasso", "grouplasso"))) {
        
        # Adjustment factor relative to Fused Lasso
        pen.count.adj <- ifelse(pen.cov[[j]] == "flasso", 1, n.par.cov[[j]] / nrow(pen.mat[[j]]))
        
        # Auxiliary matrix to compute sums of counts for all subsequent categories ("flasso"),
        # all possible combinations of categories ("gflasso") or other graph structure ("2dflasso" and "ggflasso")
        
        if (!refcat & pen.cov[[j]] != "2dflasso") {
          
          # Compute number of observations per category for given covariate and convert to matrix
          n.level.var <- as.matrix(colSums(X[, ind.cov[j-1L] + (1:n.par.cov[[j]]), drop = FALSE]))
          
          level.pen.mat <- abs(pen.mat[[j]])
          
        } else {
          
          # Compute number of observations per category for given covariate
          # All except reference category
          cs1 <- colSums(X[, ind.cov[j-1L] + (1:n.par.cov[[j]]), drop = FALSE])
          # Also include count for reference category
          n.level.var <- matrix(c(n - sum(cs1), cs1), ncol = 1) 
          
          # Also add zero column for reference category
          level.pen.mat <- cbind(rep(0, nrow(pen.mat[[j]])), abs(pen.mat[[j]]))
          
          # If only 1 category in a row (and hence rowsum 1), compare with reference category
          level.pen.mat[which(rowSums(level.pen.mat) == 1), 1L] <- 1
        }
        
        
        # Compute standardization factor
        STAND.pen.weights[[j]] <- as.numeric(sqrt((level.pen.mat %*% n.level.var) / n)) * pen.count.adj
        
      } else {
        
        # No standardization factor needed since no standardization for this type of penalty
        STAND.pen.weights[[j]] <- rep(1, length(pen.weights[[j]]))
      }
      
    }
    
    pen.weights <- Map("*", pen.weights, STAND.pen.weights)
    
    tmp <- unlist(pen.weights)
    # Add penalty weights in lambda1 to tmp, divide by lambda1.orig to only consider the penalty weights (and not the weighted lambda1)
    if (lambda1.orig > 0) {
      tmp <- c(tmp, unlist(lambda1[pen.cov %in% c("flasso", "gflasso", "ggflasso")]) / lambda1.orig) 
    }
    
    # Add penalty weights in lambda2 to tmp, divide by lambda2.orig to only consider the penalty weights (and not the weighted lambda2)
    if (lambda2.orig > 0) {
      tmp <- c(tmp, unlist(lambda2[pen.cov %in% c("flasso", "gflasso", "ggflasso")]) / lambda2.orig)
    }
    
    # Transform penalty weights
    pen.weights <- Map("*", pen.weights, (length(tmp) - 1) / sum(tmp))
    # Transform lambda1 and lambda2
    lambda1 <- Map("*", lambda1, (length(tmp) - 1) / sum(tmp))
    lambda2 <- Map("*", lambda2, (length(tmp) - 1) / sum(tmp))
  }
  
  # Print info
  if (control$print) {
    print("Penalty weights computed")
  }
  
  # Return penalty weights (list), lambda1 (list) and lambda2 (list)
  return(list(pen.weights = pen.weights, lambda1 = lambda1, lambda2 = lambda2))
}