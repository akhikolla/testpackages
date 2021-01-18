###############################################
#
# Maximum value of lambda
#
###############################################


# Compute maximum value of lambda per covariate and return maximum over all covariates
#
# X: The design matrix including ones for the intercept
# y: Response vector
# weights: Vector of prior weights
# start: vector of starting values for the coefficients
# offset: Vector containing the offset for the model
# family: Family object
# pen.cov: List with penalty type per predictor
# n.par.cov: List with number of parameters to estimate per predictor
# pen.mat: List with (weighted) penalty matrix per predictor
# pen.mat.transform: List with inverse of (weighted) penalty matrix for ordinal predictors (Fused Lasso penalty)
.max.lambda <- function(X, y, weights, start, offset, family, pen.cov, n.par.cov, pen.mat, pen.mat.transform) {
  
  # Number of covariates
  n.cov <- length(n.par.cov)
  
  ###########
  # Initialisation
  
  # Starting value for mu (using starting value for beta)
  eta <- as.numeric(X %*% start + offset)
  mu <- family$linkinv(eta)
  
  ###########
  # Gradient update
  grad <- as.numeric((weights * (mu - y) / family$variance(mu) * family$mu.eta(eta)) %*% X / sum(weights != 0))
  # Note that this is the same as
  #grad <- as.numeric(Matrix::t(X) %*% as.vector(weights * (mu - y) / family$variance(mu) * family$mu.eta(eta)) / sum(weights != 0))
  
  # Compute gradient update
  beta.tilde <- start - grad
  
  # Split beta.tilde into groups defined by covariates
  beta.tilde.split <- split(beta.tilde, rep(1:n.cov, n.par.cov))
  
  lambda.max <- numeric(length(n.par.cov))
  
  for (j in 1:length(beta.tilde.split)) {
    
    if (pen.cov[[j]] == "lasso") {
      
      # Based on Lasso problem with y = beta, X = diag and penalty = lambda * sum(pen.weights * |beta|)
      lambda.max[j] <- max(abs(beta.tilde.split[[j]]) / diag(pen.mat[[j]]))
      
    } else if (pen.cov[[j]] == "grouplasso") {
      
      # Based on Group Lasso problem with y = beta, X = diag and penalty = lambda * pen.weights * ||beta||_2
      lambda.max[j] <- sqrt(sum((beta.tilde.split[[j]] / diag(pen.mat[[j]]))^2))
      
    } else if (pen.cov[[j]] == "flasso") {
      
      # Based on dual problem, Lasso with y = beta, X = pen.mat.transform and penalty = lambda * sum(|beta|)
      lambda.max[j] <- max(abs(t(pen.mat.transform[[j]]) %*% beta.tilde.split[[j]]))
      
    } else if (pen.cov[[j]] %in% c("gflasso", "2dflasso", "ggflasso")) {
      # Minimal lambda such that all elements of u are zero is equal to
      # maximum norm of minimal x (in sense of maximum norm) such that pen.mat^t*x=beta
      if (nrow(pen.mat[[j]]) > ncol(pen.mat[[j]])) {
        lambda.max[j] <- .min.maxnorm(A = t(pen.mat[[j]]), b = beta.tilde.split[[j]])
        
      } else {
        warning(paste0("'lambda.max' cannot be determined for predictor '", names(pen.cov)[j], "'."))
      }
    }
  }
  
  return(lambda.max)
}


# Compute maximum norm of minimal x (in sense of maximum norm) such that A*x = b
# where A is a matrix of full rank, and length(b) < length(x)
# Earle et al. (2017, J Optim Theory Appl)
.min.maxnorm <- function(A, b) {
  
  # Solution to A*x = b that minimises maximum norm of x
  x0 <- .PPF(A, b)
  
  # Compute generalised inverse of A
  A.ginv <- ginv(A)
  A.kernel <- round(diag(nrow = nrow(A.ginv)) - A.ginv %*% A, 10)
  
  # Find minimum of ||x0 + A.kernel * x||_infty
  # Note that all solutions of A*x = b are of the form x0 + A.kernel * x
  opt <- optim(par = 0*x0, fn = .maxnorm, method = "BFGS",
               control = list(maxit = 1e6, reltol = 1e-10), x0 = x0, A.kernel = A.kernel)
  
  return(opt$value)
}

# Maximum norm of x0 + A.kernel * z
.maxnorm <- function(z, x0, A.kernel) {
  return(max(abs(x0 + A.kernel %*% z)))
}


# Primal Path Following (PPF) method (for underdetermined systems), phase 1 
# Earle et al. (2017, J Optim Theory Appl)
.PPF <- function(A, b) {
  
  # Initialisation, m<n for underdetermined system
  m <- nrow(A)
  n <- ncol(A)
  if (n <= m) {
    stop("Matrix 'A' is not underdetermined.")
  }
  
  # Initial starting point
  xk <- ginv(A) %*% b
  dk <- 0 * xk
  ak <- 0
  
  # Initial starting index
  Active.Set <- which(abs(xk) == max(abs(xk)))
  
  # Perform step 1 of algorithm to obtain maximal active set
  while (length(Active.Set) < n-m+1) {
    
    # Number of non-active components
    Non.Active.Set <- (1:length(xk))[-Active.Set]
    n.ASc <- n - length(Active.Set)
    
    # Split into active and non-active components
    A.AS <- matrix(A[, Active.Set], ncol = length(Active.Set))
    A.ASc <- matrix(A[, -Active.Set], ncol = n.ASc)
    xk.AS <- xk[Active.Set]
    #xk.aSc <- xk[-Active.Set]
    
    if (rankMatrix(A.ASc) < m) {
      break
    }
    
    # Calculate descent direction
    dk.AS <- -sign(xk.AS)
    dk.ASc <- ginv(A.ASc) %*% A.AS %*% sign(xk.AS)
    
    dk[Active.Set] <- dk.AS
    dk[-Active.Set] <- dk.ASc
    
    # Calculate step size
    a.possible.set <- rep(-1, n)
    xi <- xk[Active.Set[1]]
    di <- dk[Active.Set[1]]
    
    for (j in Non.Active.Set) {
      a.one <- (xi - xk[j]) / (dk[j] - di)
      a.two <- (-xi - xk[j]) / (dk[j] - (-di))
      a.cand <- c(a.one, a.two)
      
      a.possible.set[j] <- ifelse(a.one < 0 & a.two < 0, -1, min(a.cand[which(a.cand > 0)]))
    }
    ak <- min(a.possible.set[which(a.possible.set > 0)])
    
    xk <- xk + ak * dk
    Active.Set <- which(abs(xk) == max(abs(xk)))
    
  }
  
  return(xk)
  
}