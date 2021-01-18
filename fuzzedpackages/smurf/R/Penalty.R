###############################################
#
# Total penalty
#
###############################################


# Compute total penalty based on penalty types
#
# beta: Coefficient vector
# pen.cov: List with penalty type per predictor
# n.par.cov: List with number of parameters to estimate per predictor
# group.cov: List with group of each predictor which is only used for the Group Lasso penalty, 0 means no group
# pen.mat.cov: List with (weighted) penalty matrix per predictor
# lambda: Penalty parameter
# lambda1: The penalty parameter for the L_1-penalty in Sparse (Generalized) Fused Lasso or Sparse Graph-Guided Fused Lasso is lambda*lambda1
# lambda2: The penalty parameter for the L_2-penalty in Group (Generalized) Fused Lasso or Group Graph-Guided Fused Lasso is lambda*lambda2
.pen.tot <- function(beta, pen.cov, n.par.cov, group.cov, pen.mat.cov, lambda, lambda1, lambda2) {
  
  pen.tot <- 0
  n.cov <- length(pen.cov)
  # Split beta per covariate
  beta.split <- split(beta, rep(1:n.cov, n.par.cov))
  
  for (j in 1:n.cov) {
    beta <- beta.split[[j]]
    
    if (pen.cov[[j]] == "none") {
      pen <- 0
      
    } else if (pen.cov[[j]] == "lasso") {
      pen <- lambda * sum(abs(pen.mat.cov[[j]] %*% beta))

    } else if (pen.cov[[j]] == "grouplasso") {
      
      if (group.cov[[j]] == 0) {
        # Only for predictors in group 0 (i.e. no group)
        pen <- lambda * sqrt(sum((pen.mat.cov[[j]] %*% beta)^2))
        
      } else {
        # Predictors not in group 0 are treated later
        pen <- 0
      }

    } else if (pen.cov[[j]] %in% c("flasso", "gflasso", "2dflasso", "ggflasso")) {
      # lambda1[[j]] and lambda2[[j]] are 0 for 2dflasso
      # lambda1[[j]] can be a vector, lambda2[[j]] is always a single number
      pen <- lambda * (sum(abs(lambda1[[j]] * beta)) + lambda2[[j]] * sqrt(sum(beta^2)) + 
                       sum(abs(pen.mat.cov[[j]] %*% beta)))
      
    } else {
      stop("Invalid penalty type in '.pen.tot'.")
    }
    
    pen.tot <- pen.tot + pen
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
        swn <- swn + sum((pen.mat.cov[[l]] %*% beta.split[[l]])^2)
      }
      # Multiply weighted norm for group with lambda and add to total penalty
      pen.tot <- pen.tot + lambda * sqrt(swn)
    }
  }
  
  #########
  
  # Return total penalty
  return(pen.tot)
}
