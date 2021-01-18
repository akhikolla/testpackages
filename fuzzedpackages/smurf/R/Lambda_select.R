###############################################
#
# Selection of lambda
#
###############################################


# Select lambda via cross-validation, in-sample or out-of-sample
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
# pen.mat.transform: List with inverse of (weighted) penalty matrix for ordinal predictors (Fused Lasso penalty)
# valdata: Validation data when selecting lambda out-of-sample
# control: List of parameters used in the fitting process
.select.lambda <- function(X, y, weights, start, offset, family, pen.cov, n.par.cov, group.cov, refcat.cov,
                           lambda, lambda1, lambda2, pen.mat, pen.mat.aux, pen.mat.transform, valdata, control) {
  
  
  # Number of folds
  # In-sample selection of lambda using measures is the same as cross-validation with 1 fold
  # For out-of-sample selection we only have one fold, validation data is in 'valdata'
  n.folds <- ifelse(!(lambda %in% c("cv.dev", "cv.mse", "cv.dss", "cv1se.dev", "cv1se.mse", "cv1se.dss")),
                    1, control$k)
  
  n <- length(y)
  
  if (n.folds > n) {
    stop("The number of cross-validation folds cannot be larger than the number of observations.")
  }
  
  # Indicator for out-of-sample selection
  oos <- lambda %in% c("oos.dev", "oos.mse", "oos.dss")
  
  if (oos & is.null(valdata)) {
    stop("'valdata' cannot be equal to NULL when lambda is selected out-of-sample.")
  }
  
  
  ##################################
  # Allocation of stratified cross-validation folds
  
  # Vector of results
  lambda.results <- vector("list", n.folds)
  
  if (n.folds == 1) {
    # Only one fold
    cv.fold.allocation <- rep(1, n)
    
  } else {
    # Multiple folds
    
    # Matrix with indices of observations sorted based on y (from largest to smallest), and fold number
    tmp <- cbind(sort(y, decreasing = TRUE, index.return = TRUE)$ix, rep(0, n))
    # Index
    ind.strat <- 0
    
    # Loop over all urns of size n.folds. Note that the last urn can have less than n.folds observations.
    for (i in 1:ceiling(n / n.folds)) {
      
      if (i != ceiling(n / n.folds)) {
        # Sample n.folds times from 1 to n.folds without replacement
        tmp[ind.strat + (1:n.folds), 2L] <- sample(1:n.folds, size = n.folds, replace = FALSE)
        
      } else {
        # Treat last urn separately as possibly less than n.folds observations.
        tmp[(ind.strat + 1):n, 2L] <- sample(1:(n - ind.strat), size = n - ind.strat, replace = FALSE)
      }
      
      # Increase index with number of observations in urn
      ind.strat <- ind.strat + n.folds
    }
    
    cv.fold.allocation <- numeric(n)
    
    # Actual folds per observation (this is the original sample, not sorted based on y!
    cv.fold.allocation[tmp[, 1L]] <- tmp[, 2L]
  }
  
  # Indices of predictors with Lasso or Group Lasso penalty
  ind.stand <- which(pen.cov %in% c("lasso", "grouplasso"))
  
  ##################################
  
  if (is.null(control$lambda.vector)) {
    # No lambda vector is provided as input: get values for lambda.min and lambda.max from input
    # or compute them
    
    ##################################
    # lambda.max
    
    if (is.null(control$lambda.max)) {
      # Compute lambda.max
      
      if (control$print) {
        print("Computing maximum value of lambda")
      }
      
      # Return NULL if error
      tmp <- tryCatch(.max.lambda(X = X, y = y, weights = weights, start = start, offset = offset, family = family,
                                  pen.cov = pen.cov, n.par.cov = n.par.cov, pen.mat = pen.mat, 
                                  pen.mat.transform = pen.mat.transform), 
                      error = function(e) NULL)
      
      
      if (is.null(tmp)) {
        stop("An error occured when determining the maximum value of lambda. Please set it manually using 'lambda.max' in the control object.")
        
      } else {
        lambda.max <- max(tmp)
      }
      
      
      if (lambda.max <= 0) {
        stop("The computed value of 'lambda.max' is not strictly positive,",
             " please provide a suitable value for 'lambda.max' using the 'control' argument.")
      }
      
      if (control$print) {
        print(paste0("Maximum value of lambda computed: ", round(lambda.max, 4)))
      }
      
    } else {
      # Use lambda.max from control object
      lambda.max <- control$lambda.max
    }
    
    
    ##################################
    # lambda.min
    
    if (is.null(control$lambda.min)) {
      # Compute lambda.min based on lambda.max
      lambda.min <- lambda.max * 1e-4
      
    } else {
      # Use lambda.min from control object
      lambda.min <- control$lambda.min
    }
    
    if (lambda.max <= lambda.min) {
      stop("The computed value of 'lambda.max' is not strictly larger than 'lambda.min',",
           " please provide suitable values for 'lambda.max' and 'lambda.min' using the 'control' argument.")
    }
    
    
    ##################################
    # lambda.vector
    
    # Values of lambda to consider
    lambda.vector <- exp(seq(log(lambda.max), log(lambda.min), length.out = control$lambda.length))
    
  } else {
    # lambda vector is provided as input
    
    lambda.vector <- control$lambda.vector
    
    # Check if decreasing vector
    if (is.unsorted(rev(lambda.vector))) {
      warning("'lambda.vector' is not a decreasing sequence which is less efficient as we make use of warm starts.")
    }
  }
  
  
  ##################################
  # In-sample, out-of-sample or cross-validation samples
  
  if (control$print) {
    # Print how lambda is selected
    
    if (n.folds == 1 & !oos) {
      print("In-sample selection of lambda")
      
    } else if (n.folds == 1 & oos) {
      print("Out-of-sample selection of lambda")
      
    } else {
      print("Selection of lambda using cross-validation")
    }
  }
  
  # Auxiliary function for fold k in cross-validation,
  # or training data (k=1) for in-sample or out-of-sample selection
  .select.lambda.aux <- function(k) {
    
    # Make sure that reest is always equal to lambda.reest here
    control$reest <- control$lambda.reest
    
    if (n.folds > 1) {
      
      if (control$print) {
        
        print(paste0("Cross-validation fold: ", k))
      }
      
      # Cross-validation training sample
      y.train <- y[cv.fold.allocation != k]
      X.train <- X[cv.fold.allocation != k, ]
      weights.train <- weights[cv.fold.allocation != k]
      offset.train <- offset[cv.fold.allocation != k]

      
      # Standardize X.train matrix for predictors with Lasso or Group Lasso penalty (if standardize = TRUE)
      list.stand <- .X.stand(X = X.train, standardize = attr(X, "standardize"), ind.stand = ind.stand, 
                             n.par.cov = n.par.cov, weights = weights.train)
      # Possibly standardized X.train matrix
      X.train <- list.stand$X

            
      # Check if all levels are present in cross-validation training sample
      if (any(colSums(abs(X.train)) < eps_num & colSums(abs(X)) > eps_num)) {
        stop(paste0("Some levels are missing in cross-validation sample ", k, ".", 
                    " Please use a different (larger) value for 'k' in the control object."))
      }
      
      # Set attributes (about standardization)
      attr(X.train, "standardize") <- attr(X, "standardize")
      attr(X.train, "X.means") <- attr(X, "X.means") + list.stand$X.means * attr(X, "X.sds")
      attr(X.train, "X.sds") <- attr(X, "X.sds") * list.stand$X.sds
      
    } else {
      
      # Training sample is original sample if n.folds=1
      y.train <- y
      X.train <- X
      weights.train <- weights
      offset.train <- offset
    }
    
    coef_est <- matrix(0, nrow = length(lambda.vector), ncol = ncol(X))
    colnames(coef_est) <- colnames(X)
    edf <- rep(0, length(lambda.vector))
    
    # Starting value for algorithm
    start.train <- start
    
    for (i in 1:length(lambda.vector)) {
      
      if (control$print) {
        print(paste0(" Index: ", i))
        print(paste0(" lambda = ", round(lambda.vector[i], 6)))
      }
      
      # Run fitting algorithm
      fit.train <- .glmsmurf.fit.internal(X = X.train, y = y.train, weights = weights.train, start = start.train, 
                                          offset = offset.train, family = family, pen.cov = pen.cov, 
                                          n.par.cov = n.par.cov, group.cov = group.cov, refcat.cov = refcat.cov, 
                                          lambda = lambda.vector[i], lambda1 = lambda1, lambda2 = lambda2, 
                                          pen.mat = pen.mat, pen.mat.aux = pen.mat.aux, 
                                          x.return = FALSE, y.return = FALSE, control = control, full.output = FALSE)
      
      # Save (re-)estimated coefficients and degrees of freedom
      if (control$lambda.reest) {
        coef_est[i, ] <- fit.train$coefficients.reest
        edf[i] <- fit.train$rank.reest
        
      } else {
        coef_est[i, ] <- fit.train$coefficients
        edf[i] <- fit.train$rank
      }
      
      # Use warm starts, i.e. use estimates for beta of this value of lambda as starting value
      # for next value of lambda
      start.train <- fit.train$coefficients
    }
    
    # Return list with estimated coefficients and degrees of freedom
    return(list(coef_est = coef_est, edf = edf))
  }
  
  if (control$ncores == 1 | n.folds == 1) {
    # Run code on one thread
    
    # Apply function over all values for k (i.e. folds)
    lambda.results <- lapply(1:n.folds, .select.lambda.aux)
    
  } else {
    # Run code in parallel
    
    # Create clusters, use forking on Unix platforms
    # Use minimum of number of cores and number of folds
    cl.lambda <- makeCluster(min(control$ncores, n.folds), 
                         type = ifelse(.Platform$OS.type == "unix", "FORK", "PSOCK"), outfile = "") 
    
    # Export predictors to use in cluster
    clusterExport(cl.lambda, varlist = c("n.folds", "y", "X", "family", "weights", "start", 
                                         "offset", "n.par.cov", "pen.cov", "ind.stand", "group.cov", "refcat.cov",
                                         "lambda.vector", "lambda1", "lambda2",
                                         "pen.mat", "pen.mat.aux", "control", "cv.fold.allocation"), 
                  envir = environment())
    
    # Apply function in parallel over all values for k (i.e. folds)
    lambda.results <- parLapply(cl.lambda, 1:n.folds, .select.lambda.aux)
    
    # Stop cluster on exit of function
    on.exit(stopCluster(cl.lambda))
  }
  
  
  
  ##################################
  # Compute validation scores
  
  # Number of lambdas
  nlambda <- length(lambda.vector)
  
  if (attr(X, "standardize") & length(ind.stand) > 0) {
    
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
    }
  }
  
  
  # Create matrices for measures
  
  dev.val <- matrix(0, nrow = nlambda, ncol = n.folds)
  # Change row names to values of lambda
  rownames(dev.val) <- round(lambda.vector, 4)
  # Change column names
  if (n.folds == 1 & !oos) {
    colnames(dev.val) <- "In-sample"
    
  } else if (n.folds == 1 & oos) {
    colnames(dev.val) <- "Out-of-sample"
    
  } else {
    colnames(dev.val) <- paste("Fold", 1:n.folds)
  }
  
  mse.val <- dss.val <- dev.val
  aic.val <- bic.val <- gcv.val <- dev.val
  
  for (k in 1:n.folds) {
    
    if (oos) {
      # Out-of-sample calculations
      y.val <- valdata$y.oos
      X.val <- valdata$X.oos
      weights.val <- valdata$weights.oos
      offset.val <- valdata$offset.oos
      
    } else {
      
      # Out-of-sample calculations for cross-validation if n.folds > 1, otherwise in-sample calculations
      y.val <- y[cv.fold.allocation == k]
      X.val <- X[cv.fold.allocation == k, ]
      weights.val <- weights[cv.fold.allocation == k]
      offset.val <- offset[cv.fold.allocation == k]
    }
    n.val <- length(y.val)
    
    
    for (i in 1:nlambda) {
      
      edf <- lambda.results[[k]]$edf[i]
      
      # mu-values
      mu.val <- family$linkinv(as.numeric(X.val %*% lambda.results[[k]]$coef_est[i, ] + offset.val))
      
      # Deviance
      dev.val[i, k] <- sum(family$dev.resids(y = y.val, mu = mu.val, wt = weights.val))
      
      # Mean Squared Error
      mse.val[i, k] <- weighted.mean((y.val - mu.val)^2, w = weights.val)
      # Dawid-Sebastiani Score (Gneiting and Raftery, 2007)
      sigma <- sqrt(family$variance(mu.val))
      dss.val[i, k] <- weighted.mean(((y.val - mu.val) / sigma)^2 + 2 * log(sigma), w = weights.val)
      
      # AIC
      aic.val[i, k] <- family$aic(y = y.val, n = n.val, mu = mu.val, wt = weights.val, dev = dev.val[i, k]) + 2 * edf
      # BIC
      bic.val[i, k] <- aic.val[i, k] + (log(sum(weights.val != 0)) - 2) * edf
      
      # GCV score
      gcv.val[i, k] <- .GCV.score(n = length(weights.val), weights = weights.val, deviance = dev.val[i, k], edf = edf)
    }
    
  }
  
  ##################################
  # Select optimal value for lambda 
  
  if (lambda %in% c("cv.dev", "oos.dev")) {
    
    lambda.opt <- lambda.vector[.which.min.last(rowMeans(dev.val))]
    
  } else if (lambda %in% c("cv.mse", "oos.mse")) {
    
    lambda.opt <- lambda.vector[.which.min.last(rowMeans(mse.val))]
    
  } else if (lambda %in% c("cv.dss", "oos.dss")) {
    
    lambda.opt <- lambda.vector[.which.min.last(rowMeans(dss.val))]
    
  } else if (lambda == "cv1se.dev") {
    
    lambda.opt <- .cv1se(x = dev.val, lambda.vector = lambda.vector, n.folds = n.folds)
    
  } else if (lambda == "cv1se.mse") {
    
    lambda.opt <- .cv1se(x = mse.val, lambda.vector = lambda.vector, n.folds = n.folds)
    
  } else if (lambda == "cv1se.dss") {
    
    lambda.opt <- .cv1se(x = dss.val, lambda.vector = lambda.vector, n.folds = n.folds)
    
  } else if (lambda == "is.aic") {
    
    lambda.opt <- lambda.vector[.which.min.last(rowMeans(aic.val))]
    
  } else if (lambda == "is.bic") {
    
    lambda.opt <- lambda.vector[.which.min.last(rowMeans(bic.val))]
    
  } else if (lambda == "is.gcv") {
    
    lambda.opt <- lambda.vector[.which.min.last(rowMeans(gcv.val))]
    
  } else {
    stop("Invalid method to select lambda.")
  }
  
  
  
  # Measures
  
  
  if (lambda %in% c("is.aic", "is.bic", "is.gcv")) {
    # In-sample
    lambda.measures <- list(aic = aic.val, 
                            bic = bic.val, 
                            gcv = gcv.val)
    
  } else {
    # Cross-validation or out-of-sample
    lambda.measures <- list(dev = dev.val, 
                            mse = mse.val,
                            dss = dss.val)
  }
  
  
  
  
  # Coefficients for considered values of lambda
  if (n.folds > 1) {
    lambda.coef <- NULL
    
  } else {
    lambda.coef <- lambda.results[[1]]$coef_est
    rownames(lambda.coef) <- rownames(lambda.measures[[1]])
  }
  
  # Return optimal value of lambda, used values of lambda, validation scores and coefficients
  return(list(lambda = lambda.opt, 
              lambda.vector = lambda.vector, 
              lambda.measures = lambda.measures,
              lambda.coef = lambda.coef))
}

# Return last index of minimum of x
.which.min.last <- function(x) {
  return(max(which(x == min(x, na.rm = TRUE))))
}


# Perform selection of lambda using cross-validation with
# one standard error rule
#
# x: Matrix with error measures per value of lambda (rows) and fold (columns)
# lambda.vector: Vector of considered values of lambda
# n.folds: Number of folds
.cv1se <- function(x, lambda.vector, n.folds) {
  
  # Index of lambda with minimum average error measure
  ind0 <- .which.min.last(rowMeans(x))
  
  # Minimum average error measure
  cv.mean.opt0 <- mean(x[ind0, ])
  
  # Standard error of error measures for value of lambda that minimizes average error measure
  cv.se.opt0 <- sd(x[ind0, ]) / sqrt(n.folds)
  
  # Maximum value of lambda such that average error measure is still
  # smaller than minimum average error measure plus 1 standard error
  lambda.opt <- max(lambda.vector[rowMeans(x) <= cv.mean.opt0 + cv.se.opt0])
  
  return(lambda.opt)
}
