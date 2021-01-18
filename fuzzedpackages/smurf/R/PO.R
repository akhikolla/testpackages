###############################################
#
# Proximal operators
#
###############################################


# Compute proximal operators based on penalty types
#
# beta.tilde: Vector with coefficients after gradient update
# beta.old: Vector with old coefficient estimates
# pen.cov: List with penalty type per predictor
# n.par.cov: List with number of parameters to estimate per predictor
# group.cov: List with group of each predictor which is only used for the Group Lasso penalty, 0 means no group
# pen.mat.cov: List with (weighted) penalty matrix per predictor
# pen.mat.cov.aux: List with eigenvector matrix and eigenvalue vector of (weighted) penalty matrix per predictor
# lambda: Penalty parameter
# lambda1: List with lambda1 multiplied with penalty weights (for Lasso penalty) per predictor. The penalty parameter for the L_1-penalty in Sparse (Generalized) Fused Lasso or Sparse Graph-Guided Fused Lasso is lambda*lambda1
# lambda2: List with lambda2 multiplied with penalty weights (for Group Lasso penalty) per predictor. The penalty parameter for the L_2-penalty in Group (Generalized) Fused Lasso or Group Graph-Guided Fused Lasso is lambda*lambda2
# step: Current step size
# po.ncores: Number of cores used when computing the proximal operators
.PO <- function(beta.tilde, beta.old, pen.cov, n.par.cov, group.cov, pen.mat.cov, pen.mat.cov.aux, 
                lambda, lambda1, lambda2, step, po.ncores) {
  
  n.cov <- length(pen.cov)
  # Split beta.tilde into groups defined by covariates
  beta.tilde.split <- split(beta.tilde, rep(1:n.cov, n.par.cov))
  # Split beta.old into groups defined by covariates to use as starting values
  beta.old.split <- split(beta.old, rep(1:n.cov, n.par.cov))
  
  ###
  # Vector of norms for Group Lasso
  norms.group <- group.cov
  # Compute norms per predictor (for Group Lasso penalty)
  for (j in which(pen.cov == "grouplasso")) {
    norms.group[j] <- sqrt(sum(beta.tilde.split[[j]]^2))
  }
  # Unique group numbers (excluding zero)
  groups.unique.nz <- unique(group.cov[group.cov != 0])
  
  # Compute norms per group (except group zero)
  if (length(groups.unique.nz) > 0) {
    
    for (j in 1:length(groups.unique.nz)) {
      # Indices of predictors in the group
      ind.group <- which(group.cov == groups.unique.nz[j])
      # Combine norms in a group by summing the squared norms per predictor
      # For every predictor in the group, the vector element is now equal to the norm of the group!
      norms.group[ind.group] <- sqrt(sum(norms.group[ind.group]^2))
    }
  }
  ###
  
  # Output list
  prox.op.list <- beta.tilde.split
  
  # Auxiliary function for lapply or parLapply
  .PO.aux <- function(j) {
    
    if (pen.cov[[j]] == "none") {
      return(beta.tilde.split[[j]])
      
    } else if (pen.cov[[j]] == "lasso") {
      return(.PO.Lasso(beta.tilde = beta.tilde.split[[j]],
                       slambda = diag(pen.mat.cov[[j]]) * lambda * step))
      
    } else if (pen.cov[[j]] == "grouplasso") {
      # Note that the norm of all coefficients in the group is used.
      # For group 0 (i.e. no group), only the norm of the coefficients for predictor j is used.
      return(.PO.GroupLasso(beta.tilde = beta.tilde.split[[j]],
                            slambda = diag(pen.mat.cov[[j]]) * lambda * step,
                            norm = norms.group[j])) 
      
    } else if (pen.cov[[j]] %in% c("flasso", "gflasso", "2dflasso", "ggflasso")) {
      
      # Avoid numerical problems with zero eigenvalues, use old ADMM then (fast = FALSE)
      return(admm_po_cpp(beta_tilde = as.numeric(beta.tilde.split[[j]]),
                         slambda = lambda * step, 
                         # lambda1[[j]] and lambda2[[j]] are 0 for 2dflasso
                         # lambda1[[j]] can be a vector, lambda2[[j]] is always a single number
                         lambda1 = lambda1[[j]], lambda2 = lambda2[[j]],
                         penmat = pen.mat.cov[[j]], Q = pen.mat.cov.aux[[j]]$Q, eigval = pen.mat.cov.aux[[j]]$eigval, 
                         fast = all(abs(pen.mat.cov.aux[[j]]$eigval) >= eps_num), maxiter = 1e4, rho = 1, 
                         beta_old = beta.old.split[[j]]))
      
    } else {
      stop("Invalid penalty type in '.PO'.")
    }
  }      
  
  
  if (po.ncores == 1) {
    # Run code on one thread
    
    # Apply function over all values for j
    prox.op.list <- lapply(1:n.cov, .PO.aux)
    
  } else {
    # Run code in parallel
    
    # Create clusters, use forking on Unix platforms
    cl <- makeCluster(po.ncores, type = ifelse(.Platform$OS.type == "unix", "FORK", "PSOCK")) 
    
    # Export predictors to use in cluster
    clusterExport(cl, varlist = c("beta.tilde.split", "lambda", "lambda1", "lambda2", "step", "pen.mat.cov", 
                                  "pen.mat.cov.aux", "beta.old.split", "norms.group"), 
                  envir = environment())
    
    # Apply function in parallel over all values for j
    prox.op.list <- parLapply(cl, 1:n.cov, .PO.aux)
    
    # Stop cluster on exit of function
    on.exit(stopCluster(cl))
  }
  
  # Return vector with proximal operators
  return(unlist(prox.op.list))
  
}



# Proximal operator for Lasso
#
# beta.tilde: Vector with coefficients after gradient update
# slambda: penalty weights * penalty parameter (lambda) * step size
.PO.Lasso <- function(beta.tilde, slambda) {
  
  return(sign(beta.tilde) * pmax(abs(beta.tilde) - slambda, 0))
}


# Proximal operator for Group Lasso
#
# beta.tilde: Vector with coefficients after gradient update
# slambda: penalty weights * penalty parameter (lambda) * step size
# norm: Norm of all coefficients (elements of beta.tilde) in the group
.PO.GroupLasso <- function(beta.tilde, slambda, norm) {
  
  if (norm <= 0) {
    # Avoid problems if negative norm
    po <- 0 * beta.tilde
    
  } else {
    po <- beta.tilde * pmax(1 - slambda / norm, 0)
  }
  
  return(po)
}
