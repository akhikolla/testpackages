###############################################
#
# Internal fitting function
#
###############################################


# Actual fitting function with given penalty matrices (with penalty weights) and value for lambda.
# The SMuRF algorithm itself is contained in the .glmsmurf.algorithm. 
# .glmsmurf.fit.internal is a wrapper around this function and handles re-estimation and output.
#
# X: The design matrix including ones for the intercept
# y: Response vector
# weights: Vector of prior weights
# start: vector of starting values for the coefficients
# offset: Vector containing the offset for the model
# family: Family object
# pen.cov: List with penalty type per predictor
# n.par.cov: List with number of parameters to estimate per predictor
# group.cov: List with group of each predictor which is only used for the Group Lasso penalty, 0 means no group
# refcat.cov: List with number of the reference category in the original order of the levels of each predictor.
# When the predictor is not a factor or no reference category is present, it is equal to 0. This number will only be taken into account for a Fused Lasso, Generalized Fused Lasso or Graph-Guided Fused Lasso penalty 
# when a reference category is present. 
# lambda: Penalty parameter
# lambda1: List with lambda1 multiplied with penalty weights (for Lasso penalty) per predictor. The penalty parameter for the L_1-penalty in Sparse (Generalized) Fused Lasso or Sparse Graph-Guided Fused Lasso is lambda*lambda1
# lambda2: List with lambda2 multiplied with penalty weights (for Group Lasso penalty) per predictor. The penalty parameter for the L_2-penalty in Group (Generalized) Fused Lasso or Group Graph-Guided Fused Lasso is lambda*lambda2
# pen.mat: List with (weighted) penalty matrix per predictor
# pen.mat.aux: List with eigenvector matrix and eigenvalue vector of (weighted) penalty matrix per predictor
# x.return: Logical indicating if the used model matrix should be returned in the output object
# y.return: Logical indicating if the used response vector should be returned in the output object
# control: List of parameters used in the fitting process
# full.output: A logical indicating if all output is returned. Set to TRUE unless when using this function in the selection of lambda, see Lambda_select.R
.glmsmurf.fit.internal <- function(X, y, weights, start, offset, family, pen.cov, n.par.cov, group.cov, refcat.cov,
                                   lambda, lambda1, lambda2, pen.mat, pen.mat.aux, x.return, y.return, control, full.output = TRUE) {
  
  ###########
  # Objects from control

  # Numeric tolerance value for stopping criterion
  epsilon <- control$epsilon
  # Maximum number of iterations of the SMuRF algorithm
  maxiter <- control$maxiter
  # Initial step size
  step <- control$step
  # Parameter for backtracking the step size
  tau <- control$tau
  # Logical indicating if the obtained (reduced) model is re-estimated using a standard GLM
  reest <- control$reest
  # Number of cores used when computing the proximal operators
  po.ncores <- control$po.ncores
  # Logical indicating if intermediate results need to be printed
  print <- control$print
  
  
  ###########
  # Objects from family
  
  # Inverse of the link function
  linkinv <- family$linkinv
  
  # Deviance residuals as a function of (y, mu, wt)
  dev.resids <- family$dev.resids
  
  # Scaled log-likelihood function
  #
  # y: Response vector
  # n: Sample size
  # mu: Fitted mean values
  # wt: Vector of prior weights
  .scaled.ll <- function(y, n, mu, wt) {
    return(-0.5 * sum(dev.resids(y = y, mu = mu, wt = wt)) / sum(wt != 0))
  }
  
  
  ###########
  # Miscellaneous
  
  # Sample size
  n <- length(y)
  
  # Number of covariates
  n.cov <- length(pen.cov)
  
  # Check if lambda1 is a list
  if (!is.list(lambda1)) {
    
    # Change to list with lambda1 per covariate
    lambda1.orig <- lambda1
    lambda1 <- n.par.cov
    
    for (j in 1:n.cov) {
      lambda1[[j]] <- lambda1.orig
    }
    
    # Add original lambda1 as final element
    lambda1$lambda1.orig <- lambda1.orig
    
  } else {
    # Get original lambda1
    lambda1.orig <- lambda1$lambda1.orig
  }
  
  # Check if lambda2 is a list
  if (!is.list(lambda2)) {
    
    # Change to list with lambda2 per covariate
    lambda2.orig <- lambda2
    lambda2 <- n.par.cov
    
    for (j in 1:n.cov) {
      lambda2[[j]] <- lambda2.orig
    }
    
    # Add original lambda2 as final element
    lambda2$lambda2.orig <- lambda2.orig
    
  } else {
    # Get original lambda2
    lambda2.orig <- lambda2$lambda2.orig
  }
  
  ###########
  # Run SMuRF algorithm
  
  res <- .glmsmurf.algorithm(X = X, y = y, weights = weights, start = start, offset = offset, family = family,
                             pen.cov = pen.cov, n.par.cov = n.par.cov, group.cov = group.cov,
                             lambda = lambda, lambda1 = lambda1, lambda2 = lambda2,
                             .scaled.ll = .scaled.ll, pen.mat = pen.mat, pen.mat.aux = pen.mat.aux, 
                             epsilon = epsilon, maxiter = maxiter, step = step, tau = tau, po.ncores = po.ncores, print = print) 
  
  # Coefficient estimates
  beta.new <- res$beta.new
  # Convergence indicator
  conv <- res$conv
  # Number of iterations performed
  iter <- res$iter
  # Final step size
  step <- res$step
  
  
  # Round estimates for beta to 7 digits
  beta.new <- round(beta.new, 7)
  
  # Degrees of freedom: number of non-zero unique beta's
  edf <- length(unique(beta.new)[unique(beta.new) != 0])
  
  ###########
  # Model re-estimation
  
  if (reest) {
    
    # Print info
    if (print) {
      print("Re-estimation of fitted model")
    }
    
    ######
    # Indices of predictors with penalty different from 2dflasso
    ind <- which(pen.cov %in% c("lasso", "grouplasso", "flasso", "gflasso", "ggflasso"))
    
    if (length(ind) > 0) {
      
      # Split beta.new vector based on predictors
      beta.new.split <- split(beta.new, rep(1:n.cov, n.par.cov))
      nms <- names(beta.new)
      
      # Warning indicator
      warn <- 0
      
      for (j in ind) {
        
        if (sum(beta.new.split[[j]] != 0) == n.par.cov[[j]] & n.par.cov[[j]] > 1) {
   
          # Correction is not necessary for these penalty types and reference class present
          if (!(pen.cov[[j]] %in% c("flasso", "gflasso", "ggflasso") & 
              (lambda * lambda1.orig == 0 & lambda * lambda2.orig == 0))) {
            
            # If all coefficients for the predictor are non-zero, and there is more than one coefficient,
            # problems arise due to multicollinearity. This is solved by setting the first coefficient
            # (and all other coefficients which are equal to this one) to 0 to make the first (or the chosen) category the reference category
            ind.ref <- ifelse(refcat.cov[[j]] == 0, 1L, refcat.cov[[j]])
            beta.new.split[[j]][beta.new.split[[j]] == beta.new.split[[j]][ind.ref]] <- 0
            
            # Change warning indicator 
            warn <- 1
          }
        }
      }
      
      if (warn > 0) {
        warning("Reference classes are introduced for some predictor(s) in the re-estimated model.")
      }
      # Transform back to beta vector and use original names
      # Use new vector such that output for beta.new is not altered!
      beta.new2 <- unlist(beta.new.split)
      names(beta.new2) <- nms
      
    } else {
      beta.new2 <- beta.new
    }
    
    ######
    
    # Get unique beta's that are non-zero
    beta.reduced <- unique(beta.new2)[unique(beta.new2) != 0]
    # Transformation matrix which contains ones when a covariate has a non-zero estimate for beta
    X.transform <- matrix(0, nrow = ncol(X), ncol = length(beta.reduced))
    for (i in 1:length(beta.reduced)) {
      par <- beta.reduced[i]
      X.transform[which(beta.new2 == par), i] <- 1
    }
    
    if (length(beta.reduced) > 1) {
      # Combine covariates into the groups and convert to sparse matrix
      model.matrix.reest <- Matrix(X %*% X.transform, sparse = TRUE)
      # Re-estimate model using groups, use speedglm.wfit 
      glm.reest <- speedglm.wfit(y = y, X = model.matrix.reest, intercept = TRUE, family = family, 
                                 start = NULL, offset = offset, weights = weights,
                                 sparse = TRUE, trace = FALSE)
    } else {
      # Only intercept is present, use glm.fit function from "stats" to avoid problems
      # Combine covariates into the groups
      model.matrix.reest <- as.matrix(X %*% X.transform)
      # Re-estimate model using groups, use glm.fit 
      glm.reest <- glm.fit(y = y, x = model.matrix.reest, intercept = TRUE, family = family, 
                           start = NULL, offset = offset, weights = weights)
    }
    
    # Re-estimated coefficients per group
    beta.reest <- glm.reest$coefficients
    beta.reest.long <- rep(0, length(beta.new2))
    # Make new coefficient vector with re-estimated coefficients per predictor
    for (i in 1:length(beta.reduced)) {
      par <- beta.reduced[i]
      beta.reest.long[which(beta.new2 == par)] <- beta.reest[i]
    }
    
    # Round beta.reest.long to 7 digits
    beta.reest.long <- round(beta.reest.long, 7)
    
    # Degrees of freedom of re-estimated model
    edf.reest <- length(unique(beta.reest.long)[unique(beta.reest.long) != 0])
  }

  ###########
  # Transform X and beta and beta.reest.long to original scale
  
  # Logical indicating if the predictors for Lasso and Group Lasso are standardized
  if (!is.null(attr(X, "standardize"))) {
    standardize <- attr(X, "standardize")
    
  } else {
    warning("'standardize is not an attribute of X.")
    standardize <- FALSE
  }
  
  # Indices of predictors with Lasso or Group Lasso penalty
  ind.stand <- which(pen.cov %in% c("lasso", "grouplasso"))
  
  if (standardize & length(ind.stand) > 0) {
    
    # Column means of X which are set to 0
    # when the penalty type is not Lasso or Group Lasso
    X.means <- attr(X, "X.means")
    # Standard deviations of columns of X which are set to 1
    # when the penalty type is not Lasso or Group Lasso
    X.sds <- attr(X, "X.sds")
    
    for (j in ind.stand) {
      # Loop over predictors with Lasso or Group Lasso penalty
      
      # Index of first coefficient corresponding to this predictor
      ind.start <- ifelse(j == 1, 1L, sum(unlist(n.par.cov[1:(j-1L)])) + 1L)
      # Index of last coefficient corresponding to this predictor
      ind.end <- sum(unlist(n.par.cov[1:j]))
      # Indices
      ind.s.e <- ind.start:ind.end
      
      # Transform columns of X back to original scale
      if (is.matrix(X)) {
        # Standard matrix type
        
        X[, ind.s.e] <- sweep(sweep(X[, ind.s.e, drop = FALSE], 2L, X.sds[ind.s.e], "*"), 
                              2L, X.means[ind.s.e], "+")
        
      } else {
        # Sparse matrix type
        # Handle differently to avoid serious performance drop
        
        # Column names of X
        X.cnames <- colnames(X)
        
        if (ind.start == 1) {
          # ind.start is first column
          X <- cbind(sweep(sweep(X[, ind.s.e, drop = FALSE], 2L, X.sds[ind.s.e], "*"), 
                           2L, X.means[ind.s.e], "+"), 
                     X[, -ind.s.e])
          
        } else if (ind.end == ncol(X)) {
          # ind.end is last column
          X <- cbind(X[, -ind.s.e],
                     sweep(sweep(X[, ind.s.e, drop = FALSE], 2L, X.sds[ind.s.e], "*"), 
                           2L, X.means[ind.s.e], "+"))
          
        } else {
          
          X <- cbind(X[, 1:(ind.start-1L)], 
                     sweep(sweep(X[, ind.s.e, drop = FALSE], 2L, X.sds[ind.s.e], "*"), 
                           2L, X.means[ind.s.e], "+"), 
                     X[, (ind.end+1L):ncol(X)])
        }
        # Correct column names
        colnames(X) <- X.cnames 
      }
      
        
      # Correct beta.new[ind.s.e]
      beta.new[ind.s.e] <- beta.new[ind.s.e] / X.sds[ind.s.e]
      
      # Correct intercept
      beta.new[1] <- beta.new[1] - sum(X.means[ind.s.e] * beta.new[ind.s.e])

      
      if (reest) {
        # Correct beta.reest.long[ind.s.e]
        beta.reest.long[ind.s.e] <- beta.reest.long[ind.s.e] / X.sds[ind.s.e]
        
        # Correct re-estimated intercept
        beta.reest.long[1] <- beta.reest.long[1] - sum(X.means[ind.s.e] * beta.reest.long[ind.s.e])
      }
    }  
    
    if (reest) {
      # Recompute model.matrix.reest
      model.matrix.reest <- Matrix(X %*% X.transform, sparse = TRUE)
    }
  }
  
  ##############
  # Output for estimated model
  
  # Set names of beta.new
  names(beta.new) <- colnames(X)
  
  if (!full.output) {
    
    # Output with only coefficients and rank (used for selection of lambda)
    output.list <- list(coefficients = beta.new, rank = edf)
    
  } else {

  
    # Linear predictors
    eta <- as.numeric(X %*% beta.new + offset)
    # Fitted values
    mu <- linkinv(eta)
    
    # Minus scaled log-likelihood of estimated model
    fbeta <- -.scaled.ll(y = y, n = n, mu = mu, wt = weights)
    # Total penalty of estimated model
    gbeta <- .pen.tot(beta = beta.new, pen.cov = pen.cov, n.par.cov = n.par.cov, group.cov = group.cov, pen.mat.cov = pen.mat, 
                      lambda = lambda, lambda1 = lambda1, lambda2 = lambda2)
    
    # Deviance
    deviance <- sum(dev.resids(y = y, mu = mu, wt = weights))
    # Generalised cross-validation score
    GCVscore <- .GCV.score(n = n, weights = weights, deviance = deviance, edf = edf)
    # Akaike Information Criterion
    AICscore <- family$aic(y = y, n = n, mu = mu, wt = weights, dev = deviance) + 2 * edf
    # Bayesian Information Criterion
    BICscore <- AICscore + (log(sum(weights != 0)) - 2) * edf
    
    # Set names of eta, mu, weights and offset
    names(eta) <- names(y)
    names(mu) <- names(y)
    names(weights) <- names(y)
    names(offset) <- names(y)
    
    # Residuals (as in GLM)
    residuals <- .resid(y = y, eta = eta, mu = mu, family = family)
    
    # Residual degrees of freedom
    df.residual <- .df.residual(n = n, weights = weights, edf = edf)
    
    # Compute null deviance
    
    # Check if intercept is present. This is done by checking if at least one column contains only ones.
    # We only need to look at columns with sum nrow(X), which speeds up calculations and saves memory
    inter <- any(apply(X[, which(colSums(X) == nrow(X)), drop = FALSE], 
                       2, function(x) all(x == 1)))
    if (inter) {
      # Intercept is present

      # Intercept-only model, use fitted intercept as starting value
      fit.null <- tryCatch(glm(y ~ 1, family = family, weights = weights, offset = offset, start = beta.new[1]), 
                           error = function(e) list(deviance = NaN))
      # Null deviance
      null.deviance <- fit.null$deviance
      
      # Residual degrees of freedom for null model: number of observations (excluding those with zero weight) - intercept
      df.null <- (n - sum(weights == 0)) - 1
      
      
    } else {
      # No intercept is present, compute null deviance using offset
      null.deviance <- sum(dev.resids(y = y, mu = linkinv(offset), wt = weights))
      # Residual degrees of freedom for null model: number of observations (excluding those with zero weight)
      df.null <- n - sum(weights == 0)
    } 
    
    
    # Create full output list
    output.list <- list(coefficients = beta.new, residuals = residuals, fitted.values = mu, rank = edf,
                        family = family, linear.predictors = eta, deviance = deviance, 
                        aic = AICscore, bic = BICscore, gcv = GCVscore, null.deviance = null.deviance, 
                        df.residual = df.residual, df.null = df.null,
                        obj.fun = fbeta + gbeta, weights = weights, offset = offset,
                        lambda = lambda, lambda1 = lambda1.orig, lambda2 = lambda2.orig,
                        iter = iter, converged = conv, final.stepsize = step)
  }

  # Return model matrix if required
  if (x.return) {
    output.list$X <- X
  }
  
  # Return response vector if required
  if (y.return) {
    output.list$y <- y
  }
  
  ##############
  # Output for re-estimated model
  
  
  if (reest) {
    
    # Set names of beta.reest.long
    names(beta.reest.long) <- colnames(X)
    
    if (!full.output) {
      
      # Add only re-estimated coefficients and rank to output
      output.list <- c(output.list, list(coefficients.reest = beta.reest.long, rank.reest = edf.reest)) 
                       
    } else {
    
      # Re-estimated mu-values
      eta.reest <- as.numeric(X %*% beta.reest.long + offset)
      mu.reest <- linkinv(eta.reest)
      # Minus scaled log-likelihood of re-estimated model
      fbeta.reest <- -.scaled.ll(y = y, n = n, mu = mu.reest, wt = weights)
      # Total penalty of re-estimated model
      gbeta.reest <- .pen.tot(beta = beta.reest.long, pen.cov = pen.cov, n.par.cov = n.par.cov, group.cov = group.cov, pen.mat.cov = pen.mat, 
                              lambda = lambda, lambda1 = lambda1, lambda2 = lambda2)
      
      # Deviance of re-estimated model
      deviance.reest <- sum(dev.resids(y = y, mu = mu.reest, wt = weights))
      
      # Generalised cross-validation score of re-estimated model
      GCVscore.reest <- .GCV.score(n = n, weights = weights, deviance = deviance.reest, edf = edf.reest)
      # Akaike Information Criterion of re-estimated model
      AICscore.reest <- family$aic(y = y, n = n, mu = mu.reest, wt = weights, dev = deviance.reest) + 2 * edf.reest
      # Bayesian Information Criterion of re-estimated model
      BICscore.reest <- AICscore.reest + (log(sum(weights != 0)) - 2) * edf.reest

      
      # Set names of eta.reest and mu.reest
      names(eta.reest) <- names(y)
      names(mu.reest) <- names(y)
      
      # Return residuals.reest in output
      
      # Residuals (as in GLM)
      residuals.reest <- .resid(y = y, eta = eta.reest, mu = mu.reest, family = family)
      
      # Residual degrees of freedom (for re-estimated model) 
      df.residual.reest <- .df.residual(n = n, weights = weights, edf = edf.reest)
      
      # Add to output list
      output.list <- c(output.list, list(glm.reest = glm.reest, 
                                         coefficients.reest = beta.reest.long, residuals.reest = residuals.reest,
                                         fitted.values.reest = mu.reest, rank.reest = edf.reest, 
                                         linear.predictors.reest = eta.reest, deviance.reest = deviance.reest, 
                                         aic.reest = AICscore.reest, bic.reest = BICscore.reest, gcv.reest = GCVscore.reest, 
                                         df.residual.reest = df.residual.reest, obj.fun.reest = fbeta.reest + gbeta.reest))
    }
    # Return model matrix used for re-estimation if required
    # Always on original scale!
    if (x.return) {
      output.list$X.reest <- model.matrix.reest
    }
  }
  
  
  # Return output list
  return(output.list)
  
}

# Compute generalized cross-validation score
# 
# n: Sample size
# weights: Vector of prior weights
# deviance: Deviance of the model
# edf: Estimated degrees-of-freedom
.GCV.score <- function(n, weights, deviance, edf) {
  
  n2 <- n - sum(weights == 0)
  return(deviance / n2 / (1 - edf / n2) ^ 2)
}


# Compute residuals (as in GLM)
#
# y: Response vector
# eta: Fitted linear predictors (link scale)
# mu: Fitted mean values
# family: Family object
.resid <- function(y, eta, mu, family) {
  
  residuals <- (y - mu) / family$mu.eta(eta)
  names(residuals) <- names(y)
  
  return(residuals)
}


# Compute residual degrees of freedom
#
# n: Sample size
# weights: Vector of prior weights
# edf: Estimated degrees-of-freedom
.df.residual <- function(n, weights, edf) {
  
  # Number of observations (excluding those with zero weight) - model df
  return((n - sum(weights == 0)) - edf)
}