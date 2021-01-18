###############################################
#
# SMuRF algorithm
#
###############################################


# SMuRF algorithm; re-estimation and output is handled in the .glmsmurf.fit.internal function
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
# lambda: Penalty parameter
# lambda1: List with lambda1 multiplied with penalty weights (for Lasso penalty) per predictor. The penalty parameter for the L_1-penalty in Sparse (Generalized) Fused Lasso or Sparse Graph-Guided Fused Lasso is lambda*lambda1
# lambda2: List with lambda2 multiplied with penalty weights (for Group Lasso penalty) per predictor. The penalty parameter for the L_2-penalty in Group (Generalized) Fused Lasso or Group Graph-Guided Fused Lasso is lambda*lambda2
# .scaled.ll: The scaled log-likelihood function
# pen.mat: List with (weighted) penalty matrix per predictor
# pen.mat.aux: List with eigenvector matrix and eigenvalue vector of (weighted) penalty matrix per predictor
# epsilon: Numeric tolerance value for stopping criterion
# maxiter: Maximum number of iterations of the SMuRF algorithm
# step: Initial step size
# tau: Parameter for backtracking the step size
# po.ncores: Number of cores used when computing the proximal operators
# print: A logical indicating if intermediate results need to be printed
.glmsmurf.algorithm <- function(X, y, weights, start, offset, family, pen.cov, n.par.cov, group.cov,
                                lambda, lambda1, lambda2, .scaled.ll,
                                pen.mat, pen.mat.aux, epsilon, maxiter, step, tau, po.ncores, print) {
  
  
  # Inverse of the link function
  linkinv <- family$linkinv
  
  # Sample size
  n <- length(y)
  
  ###########
  # Initialisation
  
  # Initialise values of beta for use later on
  beta.old <- start
  beta.new <- start

  
  # Initialise counter for while-loop
  iter <- 0L
  
  # Initialise alpha's for acceleration
  alpha.old <- 0 # This value is overwritten before it will be used
  alpha.new <- 1
  # Initialise theta for acceleration
  theta <- beta.new
  # Count the subsequent number of restarts (to avoid infinite loops)
  subs.restart <- 0L
  # Lower bound for step size
  step.low <- 1e-14
  
  # Initialise old objective function
  obj.beta.old <- 0
  # Starting value for mu (using starting value for beta)
  mu.cand <- linkinv(as.numeric(X %*% start + offset))
  
  # f function: minus scaled log-likelihood function
  f.beta.cand <- -.scaled.ll(y = y, n = n, mu = mu.cand, wt = weights)
  # Set f.beta.cand to infinity if numerical problems occur
  if (is.nan(f.beta.cand)) f.beta.cand <- Inf
  
  # g function: total penalty
  g.beta.cand <- .pen.tot(beta = start, pen.cov = pen.cov, n.par.cov = n.par.cov, group.cov = group.cov, 
                          pen.mat.cov = pen.mat, lambda = lambda, lambda1 = lambda1, lambda2 = lambda2)
  
  # Initialise new objective function
  obj.beta.new <- f.beta.cand + g.beta.cand
  
  ###########
  # While-loop
  
  # Update estimates for beta if stopping criterion not met
  while ((abs((obj.beta.old - obj.beta.new) / obj.beta.old) > epsilon | iter==0) & iter < maxiter) {
    
    # Update old objective function
    obj.beta.old <- obj.beta.new
    # Update old estimates for beta
    beta.old <- beta.new
    
    ###########
    # Updates for values depending on theta
    
    # Compute eta and mu based on estimates for theta
    eta.theta <- as.numeric(X %*% theta + offset)
    mu.theta <- linkinv(eta.theta)
    
    # Compute f: minus scaled log-likelihood using estimates for theta
    f.theta <- -.scaled.ll(y = y, n = n, mu = mu.theta, wt = weights)
    # Set f.theta to infinity if numerical problems occur
    if (is.nan(f.theta)) f.theta <- Inf
    
    
    ###########
    # Gradient update
    
    # Compute gradient
    grad <- as.numeric((weights * (mu.theta - y) / family$variance(mu.theta) * family$mu.eta(eta.theta)) %*% X / sum(weights != 0))
    # Note that this is the same as
    # grad <- as.numeric(Matrix::t(X) %*% as.vector(weights * (mu.theta - y) / family$variance(mu.theta) * family$mu.eta(eta.theta)) / sum(weights != 0))
    
    # Compute gradient update
    beta.tilde <- theta - step * grad
    
    
    ###########
    # Proximal operators
    
    # Compute proximal operators per covariate to update beta
    beta.cand <- .PO(beta.tilde = beta.tilde, beta.old = beta.old, pen.cov = pen.cov, n.par.cov = n.par.cov, group.cov = group.cov, pen.mat.cov = pen.mat, pen.mat.cov.aux = pen.mat.aux, 
                     lambda = lambda, lambda1 = lambda1, lambda2 = lambda2, step = step, po.ncores = po.ncores)
    
    # Compute mu based on estimates for beta
    mu.cand <- linkinv(as.numeric(X %*% beta.cand + offset))
    
    # Compute f: minus scaled log-likelihood using estimates for beta
    f.beta.cand <- -.scaled.ll(y = y, n = n, mu = mu.cand, wt = weights)
    # Set f.beta.cand to infinity if numerical problems occur
    if (is.nan(f.beta.cand)) f.beta.cand <- Inf
    
    # Compute g: total penalty
    g.beta.cand <- .pen.tot(beta = beta.cand, pen.cov = pen.cov, n.par.cov = n.par.cov, group.cov = group.cov, pen.mat.cov = pen.mat, 
                            lambda = lambda, lambda1 = lambda1, lambda2 = lambda2)
    
    # Objective function
    obj.beta.cand <- f.beta.cand + g.beta.cand
    
    ###########
    # Backtracking step size
    
    # Auxiliary function
    h.theta <- f.theta + grad %*% (beta.cand - theta) + 1/(2*step) * sum((beta.cand - theta)^2) + g.beta.cand
    
    # Reduce step size until objective function in beta.cand is below auxiliary function in beta.cand
    # or too small step size
    while (obj.beta.cand > h.theta & step >= step.low) {
      
      # Reduce step size
      step <- step * tau
      
      ###########
      # Gradient update
      
      beta.tilde <- theta - step * grad
      
      ###########
      # Proximal operators
      
      # Compute proximal operators per covariate to update beta
      beta.cand <- .PO(beta.tilde = beta.tilde, beta.old = beta.old, pen.cov = pen.cov, n.par.cov = n.par.cov, group.cov = group.cov, pen.mat.cov = pen.mat, pen.mat.cov.aux = pen.mat.aux, 
                       lambda = lambda, lambda1 = lambda1, lambda2 = lambda2, step = step, po.ncores = po.ncores)
      
      # Compute mu based on updated estimates for beta
      mu.cand <- linkinv(as.numeric(X %*% beta.cand + offset))
      
      # Compute f: minus scaled log-likelihood using updated estimates for beta
      f.beta.cand <- -.scaled.ll(y = y, n = n, mu = mu.cand, wt = weights)
      # Set f.beta.cand to infinity if numerical problems occur
      if (is.nan(f.beta.cand)) f.beta.cand <- Inf
      
      # Compute g: total penalty
      g.beta.cand <- .pen.tot(beta = beta.cand, pen.cov = pen.cov, n.par.cov = n.par.cov, group.cov = group.cov, pen.mat.cov = pen.mat, 
                              lambda = lambda, lambda1 = lambda1, lambda2 = lambda2)
      
      # Objective function
      obj.beta.cand <- f.beta.cand + g.beta.cand
      
      # Auxiliary function
      h.theta <- f.theta + grad %*% (beta.cand - theta) + 1/(2*step) * sum((beta.cand - theta)^2) + g.beta.cand
    }
    
    
    # Obtained new estimate for beta
    beta.new <- beta.cand
    # Obtained new objective function
    obj.beta.new <- obj.beta.cand
    
    ###########
    # Acceleration updates
    
    # Update values for alpha
    alpha.old <- alpha.new 
    alpha.new <- (1 + sqrt(1 + 4 * alpha.old^2)) / 2
    
    # If objective function increases in this iteration:
    # keep old estimate for beta and reset values for alpha 
    if (obj.beta.new > obj.beta.old * (1 + epsilon)) {
      
      beta.new <- beta.old
      obj.beta.new <- obj.beta.old
      obj.beta.old <- 0
      
      # This value is overwritten before it will be used since beta.new = beta.old (see theta below)
      alpha.old <- 0
      # This value is used in next iteration
      alpha.new <- 1
      
      # Add 1 to subsequent restart counter
      subs.restart <- subs.restart + 1L
      
    } else { 
      # Reset subsequent restart counter
      subs.restart <- 0L 
    }
    
    # Extra if-loop to evade possible subsequent restarting loop 
    # (e.g. can happen when accuracy of starting value is too high)
    if (subs.restart > 1) {
      warning("Two subsequent restarts were performed, while-loop is ended.")
      break
    }
    
    # Update value for theta
    # Note that theta = beta.new if objective function increased in this iteration since beta.new = beta.old then (see if-loop above)
    theta <- beta.new + (alpha.old - 1) / alpha.new  * (beta.new - beta.old)
    
    # Increase counter
    iter <- iter + 1L
    
    # Print info every 100 iterations
    if (print & (iter %% 100 == 0)) {
      print(paste("   Iteration:", iter))
      print(paste("   Objective function:", obj.beta.new)) 
      print(paste("   Step size:", step)) 
    }
  }
  
  # Convergence indicator
  conv <- numeric(0)
  # Maximum number of iterations reached
  if (iter == maxiter) conv <- 1L
  # Two subsequent restarts
  if (subs.restart > 1) conv <- c(conv, 2L) 
  # Low step size
  if (step < step.low) {
    conv <- c(conv, 3L) 
    # Issue warning
    warning(paste0("The step size is below ", step.low, " and is no longer reduced. ",
                   "It might be better to start the algorithm with a larger step size."))
  }
  # Succesful convergence
  if (length(conv) == 0) conv <- 0L
  
  # Return estimated coefficients, number of performed iterations, final step size and convergence indicator
  return(list(beta.new = beta.new, iter = iter, step = step, conv = conv))
}