###############################################
#
# Standardization function
#
###############################################


# Standardize columns of X matrix corresponding to Lasso or Group Lasso penalties
#
# X: Model matrix
# standardize: Logical indicating if the columns of X corresponding to Lasso or Group Lasso penalties 
# need to be standardized
# ind.stand: Indices of predictors with a Lasso or a Group Lasso penalty  
# n.par.cov: List with number of parameters to estimate per predictor
# weights: Vector of prior weights
.X.stand <- function(X, standardize, ind.stand, n.par.cov, weights) {
  
  # Column means of X which are set to 0
  # when the penalty type is not Lasso or Group Lasso
  X.means <- rep(0, ncol(X))
  # Standard deviations of columns of X which are set to 1
  # when the penalty type is not Lasso or Group Lasso
  X.sds <- rep(1, ncol(X))
  
  if (standardize & length(ind.stand) > 0) {
    
    # Check if the weights sum to one
    if (abs(sum(weights) - 1) < eps_num) {
      # Nothing added in denominator which results in biased weighted standard deviation
      wsd <- 0
      warning("Weights sum to one, biased weighted standard deviation is used in standardization.")
      
    } else {
      # -1 added in denominator which results in unbiased weighted standard deviation
      wsd <- 1
    }
    
    for (j in ind.stand) {
      # Loop over predictors with Lasso or Group Lasso penalty
      
      # Index of first coefficient corresponding to this predictor
      ind.start <- ifelse(j == 1, 1L, sum(unlist(n.par.cov[1:(j-1L)])) + 1L)
      # Index of last coefficient corresponding to this predictor
      ind.end <- sum(unlist(n.par.cov[1:j]))
      # Indices
      ind.s.e <- ind.start:ind.end
      
      # Compute (weighted) mean and (weighted) sd of columns ind.s.e
      X.means[ind.s.e] <- colSums(sweep(X[, ind.s.e, drop = FALSE], 1L, weights, "*")) / sum(weights)
      X.sds[ind.s.e] <- sqrt(colSums(sweep(sweep(X[, ind.s.e, drop = FALSE], 2L, X.means[ind.s.e], "-") ^ 2,
                                           1L, weights, "*")) / (sum(weights) - wsd))
      
      #       # Slower for-loop
      #       for (l in ind.s.e) {
      #         # Compute (weighted) mean and (weighted) sd of column l
      #         X.means[l] <- weighted.mean(X[, l], w = weights)
      #         X.sds[l] <- sqrt(sum(weights * (X[, l] - X.means[l]) ^ 2) / (sum(weights) - wsd))
      #       }
      
      # Standardize columns in ind.s.e
      if (is.matrix(X)) {
        # Standard matrix type
        
        X[, ind.s.e] <- sweep(sweep(X[, ind.s.e, drop = FALSE], 2L, X.means[ind.s.e], "-"), 
                              2L, X.sds[ind.s.e], "/")
        
      } else {
        # Sparse matrix type
        # Handle differently to avoid serious performance drop
        
        # Column names of X
        X.cnames <- colnames(X)
        
        if (ind.start == 1) {
          # ind.start is first column
          X <- cbind(sweep(sweep(X[, ind.s.e, drop = FALSE], 2L, X.means[ind.s.e], "-"), 
                           2L, X.sds[ind.s.e], "/"), 
                     X[, -ind.s.e])
          
        } else if (ind.end == ncol(X)) {
          # ind.end is last column
          X <- cbind(X[, -ind.s.e],
                     sweep(sweep(X[, ind.s.e, drop = FALSE], 2L, X.means[ind.s.e], "-"), 
                           2L, X.sds[ind.s.e], "/"))
          
        } else {
          
          X <- cbind(X[, 1:(ind.start-1L)], 
                     sweep(sweep(X[, ind.s.e, drop = FALSE], 2L, X.means[ind.s.e], "-"), 
                           2L, X.sds[ind.s.e], "/"), 
                     X[, (ind.end+1L):ncol(X)])
        }
        # Correct column names
        colnames(X) <- X.cnames 
      }
    } 
  }
  
  return(list(X = X, X.means = X.means, X.sds = X.sds))
}
