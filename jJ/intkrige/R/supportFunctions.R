# ============================================================================
# Wrapper to the Newton - Rhapson algorithm that returns a two column matrix
# of interval values predictions at the locations of interest.
# Inputs:
# - lam - vector of lambdas from previous iteration
# - pCovC  - pairwise Covariance matrix for centers
# - pCovR  - pairwise Covariance matrix for Radii
# - pCovCR - pairwise Covariance matrix for Center Radius interaction
# - lCovC  - vector of covariances (centers) at location of interest
# - lCovR  - vector of covariances (radius) at location of interest
# - lCovCR - vector of covariances (center/radius) at location of interest
# - measurements - A two column matrix holding the centers and radii at the
#                  measurement locations.
# - eta    - growth/shink parmameter for penalty criteria.
# - ctrend - When defined, used to detrend the centers in the simple
#            kriging calculation. When trend = -9999, ordinary kriging
#            is used instead.
# - A      - Vector of three generalized L2 distance weights.
# - thresh - when |lambda_i| < (1/n*thresh), value is set to 0.
# - tolq   - convergence tolerance for NR step
# - maxq   - max number of iterations for NR step
# - tolp   - convergence tolerance for penalty criteria
# - maxp   - max number of iterations for penalty criteria
# - r      - initial (and somewhat large) penalty parmater that ensures
#            intitial solutions resides in feasible region. Used for
#            ordinary kriging only.
# Outputs:
# - predicts - A two column matrix representing the centers and radii at each
#            - resective location of interest.
# ============================================================================
nrShell_R = function(pCovC, pCovR, pCovCR, lCovC, lCovR, lCovCR,
                     measurements, eta, ctrend, A, thresh, tolq, maxq,
                     tolp, maxp, r, fast, weights, cores){

  # Make predictions at each prediction location.
  predictionC <- predictionR <- predVar <- vector("numeric", ncol(lCovC))

  if(weights){
    updt2 <- vector("list", ncol(lCovC))
  }

  # Calculate the effective magnitude of 0. Lambda values within this
  # magnitude are set equal to 0.
  len = nrow(measurements)
  threshold = abs(1/(len*thresh))

  # Pre-allocate a matrix to store the intkrige outputs depending on whether
  # or not weights are requested.
  if(weights){
    tcol <- 4 + nrow(pCovC)
  }else{
    tcol <- 4
  }

  # If cores > 1, run in parallel. If cores = 1, do not load parallel packages
  # to save time.
  if(cores > 1){
    # Create an register the clusters.
    `%dopar%` <- foreach::`%dopar%`
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)

    # If this function fails, we need to make sure to close the parallel connection.
    predictionCR <- try({foreach::foreach(i = 1:ncol(lCovC),
                                     .combine = rbind,
                                     #.packages = "intkrige",
                                     .export = c("optimShell",
                                                 "newRap_2_R",
                                                 "newRap_long_R",
                                                 "newRap_R",
                                                 "nrstep_2_R",
                                                 "nrStep_long_R",
                                                 "nrStep_R")) %dopar%
                                       optimShell(pCovC, pCovR,
                                                  pCovCR, lCovC[, i],
                                                  lCovR[, i], lCovCR[, i],
                                                  measurements, eta,
                                                  ctrend, A,
                                                  threshold, tolq,
                                                  maxq, tolp,
                                                  maxp, r,
                                                  fast, weights)})

    # Close the parallel connection even if the function fails.
    if(inherits(predictionCR, "try-error")){
      parallel::stopCluster(cl) # Close the parallel connection.
      stop("Kriging prediction failed to execute in parallel. Try
           sequential version instead.")
    }

    predictionCR <- unname(predictionCR)

    parallel::stopCluster(cl) # Close the parallel connection.

  }else{ # If not in parallel, compute the usual way.

    predictionCR <- matrix(0, nrow = ncol(lCovC), ncol = tcol)
    for(i in 1:ncol(lCovC)){
      predictionCR[i, ] <- optimShell(pCovC, pCovR, pCovCR, lCovC[, i], lCovR[, i],
                                      lCovCR[, i], measurements, eta, ctrend, A,
                                      threshold, tolq, maxq, tolp, maxp, r,
                                      fast, weights)
    }
  }

  # Return predictions for the center and radius, as well as the list of
  # lambdas used for each prediction.
  return(predictionCR)
}

# ============================================================================
# This function was added to accomodate the parallelization. It is essentially
# a switch function that determines which optimization technique is most
# appropriate, decides whether or not to complete the task in parallel,
# and completes the optimization.
# The arguments are identical to the arguments passed to the newRap
# functions. The weights simply determines whether or not to append
# the vector of prediction weights for each iteration.
# Returns:
# - predictionCR: a vector that includes predictions for center and radius,
#                 the prediction variance, and possibly a sub-vector
# ============================================================================
optimShell <- function(pCovC, pCovR, pCovCR, lCovC, lCovR, lCovCR,
                       measurements, eta, ctrend, A, threshold, tolq, maxq,
                       tolp, maxp, r, fast, weights){
  # Use naive simple kriging for intitial guess
  lambda.0 <- solve(pCovC, lCovC)

  # Don't allow initial guess to be too close to 0.
  # This should help to avoid singularities.
  lambda.0[lambda.0 < threshold & lambda.0 >= 0] <- threshold
  lambda.0[lambda.0 > -threshold & lambda.0 < 0] <- -threshold

  if(is.null(ctrend)){
    lambda.0 <- abs(lambda.0)

    updt <- newRap_2_R(lambda.0, pCovC, pCovR, pCovCR,
                       lCovC, lCovR, lCovCR, maxq, tolq,
                       threshold, eta, maxp, tolp, r, A)

    # Extract the warning flag and subset the lambdas accordingly.
    is.warn = updt[1]
    updt = updt[-1]

    predictionC <- t(measurements[, 1]) %*% updt
  }else{
    if(fast){
      updt <- newRap_R(lambda.0, pCovC, pCovR, pCovCR,
                       lCovC, lCovR, lCovCR, maxq, tolq,
                       threshold, eta, maxp, tolp, A)
    }else{
      updt <- newRap_long_R(lambda.0, pCovC, pCovR, pCovCR,
                            lCovC, lCovR, lCovCR, maxq, tolq,
                            threshold, eta, maxp, tolp, A, r)
    }


    # Extract the warning flag and subset the lambdas accordingly.
    is.warn = updt[1]
    updt = updt[-1]

    predictionC <- (t(measurements[, 1] - ctrend) %*% updt) + ctrend
  }

  # Use of abs() ensures that radius multiplier is always positive.
  predictionR <- t(measurements[, 2]) %*% abs(updt)

  # Compute the prediction variance. Remember that the diagonals of the
  # covariance matrix represent the variance at length 0. These are added so
  # that the prediction variance remains positive.
  lamlam <- updt %*% t(updt)
  lamlam2 <- updt %*% t(abs(updt))

  predVar <-
    A[1]*( (t(as.vector(lamlam)) %*% as.vector(pCovC)) -
             2*(t(updt) %*% lCovC) + pCovC[1, 1] ) +
    A[2]*( (t(as.vector(abs(lamlam))) %*% as.vector(pCovR)) -
             2*(t(abs(updt)) %*% lCovR) + pCovR[1, 1] ) +
    A[3]*( (t(as.vector(lamlam2)) %*% as.vector(pCovCR)) -
             (t(abs(updt) + updt) %*% lCovCR) + pCovCR[1, 1] )

  predictionCR <- c(predictionC, predictionR, predVar, is.warn)

  # Include the weights if requested.
  # If desired, store the lambda values.
  if(weights){
    predictionCR <- c(predictionCR, updt)
  }

  return(predictionCR)
}

# ============================================================================
### Function that carries out the newton raphson algorithm
### for simple kriging
# Inputs:
# lam - current values of the weights (lambdas)
# pCovC  - pairwise covariances for the interval centers (C^C(u_i - u_j))
# pCovR  - pairwise covariances for radii
# pCovCR - pairwise cross covariance between center and radius
# lCovC  - location of interest covariances for centers
# lCovR  - location of interest covariances for radii
# lCovCR - location of interest cross covariances for center/radius
# maxiter- for a set penalty term, maximum number of iterations allowed
#          for convergence.
# tolq    - for a set penalty, tolerance of max difference in lambda vector
#          required to suspend the algorithm.
# threshold - values with magnitude less than threshold are set to 0.
# eta    - growth parameter of penalty term (eta > 1).
# maxp   - maximum number of iterations allowed for penalty shrinkage.
# tolp   - convergence tolerance for sum(lam) = 1. Note that the algorithm
#          will not allow any lambda term to be negative, thus satisfying
#          these inequality constraints.
# A      - vector of length 3 with weights for generalized L2 distance.
# Outputs:
# lamFinal - A (hopefully) convergent vector of lambdas.
# ============================================================================
newRap_R = function(lam, pCovC, pCovR, pCovCR, lCovC, lCovR,
                    lCovCR, maxiter, tolq, threshold, eta, maxp,
                    tolp, A){ #, reset = TRUE
  # Determine the length of the original set of lambdas
  len <- length(lam)

  # Include an indices for the lambdas, to keep track of which ones are nonzero.
  tlam <- matrix(c(1:len, lam), ncol = 2, byrow = FALSE)

  # Make copies of original parameters, the dimensions of this will possibly
  # change throughout.
  tpCovC <- pCovC
  tpCovR <- pCovR
  tpCovCR <- pCovCR
  tlCovC <- lCovC
  tlCovR <- lCovR
  tlCovCR <- lCovCR

  # Parameters required for the penalty update.
  adj <- 0
  lastTry <- FALSE

  # Keep track of non-convergent solutions.
  is.warn = 0

  # Run the newton-rhapson algorithm.
  for(j in 1:(maxp+1)){ # +1 - the first iteration is run without penalty.
    for(i in 1:maxiter){
      lamup <- nrStep_R(tlam[,2], tpCovC, tpCovR, tpCovCR,
                        tlCovC, tlCovR, tlCovCR,
                        pen = adj, A = A)

      # Check the difference in the new predictions versus the old.
      tdiff <- abs(tlam[,2] - lamup)

      # Update the values of the lambdas
      tlam[, 2] <- lamup

      # Now, remove all lambdas below the threshold.
      # This allows for an expanding 0 threshold.
      include = which(abs(tlam[, 2]) >= threshold)
      tlam = tlam[include, ]

      # If only one variable is recommended for inclusion, end the algorithm by
      # assigning a value of 1 to the remaining lambda and zero to all others.
      if(length(include) < 2){
        lamFinal <- vector("numeric", len)

        if(length(include) < 1){
          warning("No viable lambas exist, returning prediction = 0...")
          is.warn = 1
          return(c(is.warn, lamFinal))
        }

        # If exactly one lambda remains, set this equal to 0, all others equal
        # to 0, and return.
        lamFinal[tlam[1]] <- 1
        return(c(is.warn, lamFinal))
      }

      # In addition, reduce the size of the hessian and gradient.
      tpCovC <- tpCovC[include, include]
      tpCovR <- tpCovR[include, include]
      tpCovCR <- tpCovCR[include, include]
      tlCovC <- tlCovC[include]
      tlCovR <- tlCovR[include]
      tlCovCR <- tlCovCR[include]

      maxD <- max(tdiff)

      if(maxD <= tolq){break}

    } # End the inner-loop

    # Break the outer loop when we have a feasible solution
    penalty <- sum(abs(tlam[, 2])) - 1

    # Look at magnitude of penalty, not raw values.
    # End the search if we have a convergent solution that meets the
    # criteria. If not, try one more iteration of the penalty to see
    # if things improve.
    if(abs(penalty) <= tolp){
      if(maxD <= tolq || lastTry){
        if(maxD > tolq){
          warning( paste("Feasible solution obtained from a non-convergent",
                         "optimization step with max(tdiff) =", maxD ) )
          is.warn = 1
        }
        break
      }else{
        lastTry <- TRUE
      }
    }

    if(j == maxp){
      warning( paste("Convergent, feasible solution not obtained with c =",
                     adj, "and sum(abs(lambda))=", penalty + 1,
                     "and max(tdiff) = ", maxD) )
      is.warn = 1
    }

    # Set penalty to 1 after first try, then grow thereafter.
    if(j == 1){
      # Adjust the penalty function to be roughly the same scale as the
      # minimization equation. Do this by ensuring that the first
      # non-zero penalty is on the order of 1.
      adj <- 1/(penalty^2)

    }else{
      adj = adj/eta
    }
  }

  # Return final results, with warning in mind.
  lamFinal <- vector("numeric", len)
  lamFinal[tlam[,1]] = tlam[,2]

  return(c(is.warn, lamFinal))
}

# ============================================================================
### Function that carries out the newton raphson algorithm
### for simple kriging. Includes penalty term that discourages 0-valued lambda
# Inputs:
# lam - current values of the weights (lambdas)
# pCovC  - pairwise covariances for the interval centers (C^C(u_i - u_j))
# pCovR  - pairwise covariances for radii
# pCovCR - pairwise cross covariance between center and radius
# lCovC  - location of interest covariances for centers
# lCovR  - location of interest covariances for radii
# lCovCR - location of interest cross covariances for center/radius
# maxiter- for a set penalty term, maximum number of iterations allowed
#          for convergence.
# tolq    - for a set penalty, tolerance of max difference in lambda vector
#          required to suspend the algorithm.
# threshold - values with magnitude less than threshold are set to 0.
# eta    - growth parameter of penalty term (eta > 1).
# maxp   - maximum number of iterations allowed for penalty shrinkage.
# tolp   - convergence tolerance for sum(lam) = 1. Note that the algorithm
#          will not allow any lambda term to be negative, thus satisfying
#          these inequality constraints.
# A      - vector of length 3 with weights for generalized L2 distance.
# r      - initial penalty parameter.
# Outputs:
# lamFinal - A (hopefully) convergent vector of lambdas.
# ============================================================================
newRap_long_R = function(lam, pCovC, pCovR, pCovCR, lCovC, lCovR,
                         lCovCR, maxiter, tolq, threshold, eta, maxp,
                         tolp, A, r){ #, reset = TRUE

  # Determine the length of the original set of lambdas
  len <- length(lam)

  # Include an indices for the lambdas, to keep track of which ones are nonzero
  tlam <- matrix(c(1:len, lam), ncol = 2, byrow = FALSE)

  # Make copies of original parameters, the dimensions of this will possibly
  # change throughout.
  tpCovC <- pCovC
  tpCovR <- pCovR
  tpCovCR <- pCovCR
  tlCovC <- lCovC
  tlCovR <- lCovR
  tlCovCR <- lCovCR

  lastTry <- FALSE
  is.warn = 0
  # Run the newton-rhapson algorithm.
  for(j in 1:maxp){

    for(i in 1:maxiter){
      lamup <- nrStep_long_R(tlam[,2], tpCovC, tpCovR, tpCovCR,
                             tlCovC, tlCovR, tlCovCR,
                             pen = r, A = A, len)

      # Check the difference in the new predictions versus the old.
      tdiff <- abs(tlam[,2] - lamup)

      # Update the values of the lambdas
      tlam[, 2] <- lamup

      # Now, remove all lambdas below the threshold.
      # This allows for an expanding 0 threshold.
      include = which(abs(tlam[, 2]) >= threshold)
      tlam = tlam[include, ]

      # If only one variable is recommended for inclusion, end the
      # algorithm by assigning a value of 1 to the remaining lambda and
      # zero to all others.
      if(length(include) < 2){
        lamFinal <- vector("numeric", len)

        if(length(include) < 1){
          warning("No viable lambas exist, returning prediction = 0...")
          is.warn = 1
          return(c(is.warn, lamFinal))
        }

        # If exactly one lambda remains, set this equal to 0, all others equal
        # to 0, and return.
        lamFinal[tlam[1]] <- 1
        return(c(is.warn, lamFinal))
      }

      # In addition, reduce the size of the hessian and gradient.
      tpCovC <- tpCovC[include, include]
      tpCovR <- tpCovR[include, include]
      tpCovCR <- tpCovCR[include, include]
      tlCovC <- tlCovC[include]
      tlCovR <- tlCovR[include]
      tlCovCR <- tlCovCR[include]

      maxD <- max(tdiff)
      if(maxD <= tolq){break}

    } # End the inner-loop.

    # Break the outer loop when we have a feasible solution
    penalty <- sum(abs(tlam[, 2])) - 1

    # Look at magnitude of penalty, not raw values.
    # End the search if we have a convergent solution that meets the
    # criteria. If not, try one more iteration of the penalty to see
    # if things improve.
    if(abs(penalty) <= tolp){
      if(maxD <= tolq || lastTry){
        if(maxD > tolq){
          warning( paste("Feasible solution obtained from a non-convergent",
                         "optimization step with max(tdiff) =", maxD ) )
          is.warn = 1
        }
        break
      }else{
        lastTry <- TRUE
      }
    }

    # Shrink the penalty term and try again.
    r = r*eta

  } # End the outer loop.

  # Warn the user if the iteration stopped before a feasible solution
  # was obtained.
  if(j == maxp){
    warning( paste("Convergent, feasible solution not obtained with c =", r,
                   "and sum(abs(lambda))=", penalty + 1, "and max(tdiff) = ",
                   maxD) )
    is.warn = 1
  }

  # Return final results, with warning in mind.
  lamFinal <- vector("numeric", len)
  lamFinal[tlam[,1]] = tlam[,2]
  return(c(is.warn, lamFinal))
}


# ============================================================================
### Function that carries out the newton raphson algorithm
### for ordinary kriging
# Inputs:
# lam - current values of the weights (lambdas)
# pCovC  - pairwise covariances for the interval centers (C^C(u_i - u_j))
# pCovR  - pairwise covariances for radii
# pCovCR - pairwise cross covariance between center and radius
# lCovC  - location of interest covariances for centers
# lCovR  - location of interest covariances for radii
# lCovCR - location of interest cross covariances for center/radius
# maxiter- for a set penalty term, maximum number of iterations allowed
#          for convergence.
# tolq    - for a set penalty, tolerance of max difference in lambda vector
#          required to suspend the algorithm.
# threshold - values with magnitude less than threshold are set to 0.
# eta    - shink parameter of penalty term (eta < 1).
# maxp   - maximum number of iterations allowed for penalty shrinkage.
# tolp   - convergence tolerance for sum(lam) = 1. Note that the algorithm
#          will not allow any lambda term to be negative, thus satisfying
#          these inequality constraints.
# r      - penalty parameter for constrained optimization
# A      - vector of length 3 with weights for generalized L2 distance.
# Outputs:
# lamFinal - A (hopefully) convergent vector of lambdas.
# ============================================================================
newRap_2_R = function(lam, pCovC, pCovR, pCovCR,
                      lCovC, lCovR, lCovCR, maxiter,
                      tolq, threshold, eta, maxp, tolp, r, A){

  # Throw error if given values not in feasible region.
  if(min(lam) < 0){stop("Initial lambda must be within feasible region
                        (i.e. all positive values of lambda).")}

  if(eta >= 1){stop("Shrink parameter eta must be less than one.")}

  # Determine the length of the original set of lambdas
  len <- length(lam)

  # Include an indices for the lambdas, to keep track of which ones are nonzero.
  tlam <- matrix(c(1:len, lam), ncol = 2, byrow = FALSE)

  # Make copies of original parameters, the dimensions of this will possibly
  # change throughout.
  tpCovC <- pCovC
  tpCovR <- pCovR
  tpCovCR <- pCovCR
  tlCovC <- lCovC
  tlCovR <- lCovR
  tlCovCR <- lCovCR

  # Parameters required for the penalty update.
  pen <- 0
  adj <- 0

  done <- FALSE
  lastTry <- FALSE
  is.warn = 0

  for(j in 1:maxp){
    for(i in 1:maxiter){

      lamup <- nrstep_2_R(tlam[,2], tpCovC, tpCovR, tpCovCR,
                          tlCovC, tlCovR, tlCovCR,
                          r = r, A = A, threshold = threshold,
                          len = len)

      # Break the inner loop. Resetting the done variable breaks the
      # outer loop also.
      if(min(lamup) < -threshold){
        done <- TRUE
        break
      }

      # Check the difference in the new predictions versus the old.
      tdiff <- abs(tlam[,2] - lamup)

      # Update the values of the lambdas
      tlam[, 2] <- lamup

      # Now, remove all lambdas at or below the threshold.
      include = which(abs(tlam[, 2]) >= threshold)
      tlam = tlam[include, ]

      # If only one variable is recommended for inclusion, end the algorithm by
      # assigning a value of 1 to the remaining lambda and zero to all others.
      if(length(include) < 2){
        lamFinal <- vector("numeric", len)

        if(length(include) < 1){
          warning("No viable lambas exist, returning prediction = 0...")
          is.warn = 1
          return(c(is.warn, lamFinal))
        }

        # If exactly one lambda remains, set this equal to 1, all others equal
        # to 0, and return.
        lamFinal[tlam[, 1]] <- 1
        return(c(is.warn, lamFinal))
      }

      # In addition, reduce the size of the hessian and gradient.
      tpCovC <- tpCovC[include, include]
      tpCovR <- tpCovR[include, include]
      tpCovCR <- tpCovCR[include, include]
      tlCovC <- tlCovC[include]
      tlCovR <- tlCovR[include]
      tlCovCR <- tlCovCR[include]

      # If we are within the tolerance, break the loop
      maxD <- max(tdiff)
      if(maxD <= tolq){break}

    } # End the inner for-loop

    if(done){
      warning(paste("Left feasible region with r = ", r,
                    ". Solution not feasible.", sep = ""))
      is.warn = 1
      break
    } # Break the loop if we have left the feasible region

    # End the search if the equality condition is satisfied.
    penalty <- sum(abs(tlam[, 2])) - 1

    # Look at magnitude of penalty, not raw values.
    # End the search if we have a convergent solution that meets the
    # criteria. If not, try one more iteration of the penalty to see
    # if things improve.
    if(abs(penalty) <= tolp){
      if(maxD <= tolq || lastTry){
        if(maxD > tolq){
          warning( paste("Feasible solution obtained from a non-convergent",
                         "optimization step with max(tdiff) =", maxD ) )
          is.warn = 1
        }
        break
      }else{
        lastTry <- TRUE
      }
    }

    if(j == maxp){
      warning( paste("Convergent, feasible solution not obtained with c =", r,
                     "and sum(abs(lambda))=", penalty + 1, "and max(tdiff) = ",
                     maxD) )
      is.warn = 1
    }

    r = r*eta # Reduce the size of r penalty parameter
  } # End the outer loop.


  if(j < 2){
    warning("Barrier method converged in a single iteration, confirm results
            are valid...")
    is.warn = 1
  }

  # Return final results, with warning in mind.
  lamFinal <- vector("numeric", len)
  lamFinal[tlam[, 1]] = tlam[, 2]

  return(c(is.warn, lamFinal))
} # End the function

# ============================================================================
### Function to update the value of lambda according to
# the Newton-Raphson Method (simple kriging)
# Inputs:
# lam - current values of the weights (lambdas)
# pCovC  - pairwise covariances for the interval centers (C^C(u_i - u_j))
# pCovR  - pairwise covariances for radii
# pCovCR - pairwise cross covariance between center and radius
# lCovC  - location of interest covariances for centers
# lCovR  - location of interest covariances for radii
# lCovCR - location of interest cross covariances for center/radius
# pen    - penalty parameter for constrained optimization
# A      - vector of length 3 with weights for generalized L2 distance.
# Outputs:
# lamUp  - and updated version of the lambda vector.
# Notes:
# Implementation of the penalty assumes the square difference of the
# contraint equation is being used.
# ============================================================================
nrStep_R = function(lam, pCovC, pCovR, pCovCR,
                    lCovC, lCovR, lCovCR, pen, A){
  l = length(lam)

  cGrad = (matrix(lam, nrow = 1) %*% pCovC) - lCovC # No need to store entries.

  rGrad = (matrix(abs(lam), nrow = 1) %*% pCovR) -
    lCovR

  calc.cr <- matrix(lam, nrow = 1) %*% pCovCR

  crGrad = (matrix(abs(lam), nrow = 1) %*% pCovCR) +
    sign(lam)*calc.cr -
    (sign(lam) + 1)*lCovCR

  penalty <- pen*(sum(abs(lam)) - 1)

  grad = 2*(A[1]*cGrad + sign(lam)*(A[2]*rGrad + penalty) + A[3]*crGrad)

  grad = matrix(grad, ncol = 1)

  ## Define the Hessian Matrix
  # First, define the matrix of sign functions to multiply with the pwCovR
  lamSign = sign(lam) %*% matrix(sign(lam), nrow = 1)

  lamSign.cr <- outer(sign(lam), sign(lam), "+")

  # Use Quadratic Approximation
  hessianR <- lamSign*(pCovR) + diag(((as.vector(rGrad)/abs(lam))))

  hessianCR <- lamSign.cr*pCovCR + diag((as.vector(calc.cr - lCovCR)/abs(lam)))

  penalty.hessian <- lamSign*pen + diag((penalty/abs(lam)))

  hessian = 2*(A[1]*pCovC + A[2]*hessianR + A[3]*hessianCR + penalty.hessian)
  b <- (hessian %*% lam) - grad

  lamUp = solve(hessian, b)

  return(lamUp)
}

# ============================================================================
### Function to update the value of lambda according to
# the Newton-Raphson Method (simple kriging). Includes a penalty parameter that
# discourages 0-valued lambdas.
# Inputs:
# lam - current values of the weights (lambdas)
# pCovC  - pairwise covariances for the interval centers (C^C(u_i - u_j))
# pCovR  - pairwise covariances for radii
# pCovCR - pairwise cross covariance between center and radius
# lCovC  - location of interest covariances for centers
# lCovR  - location of interest covariances for radii
# lCovCR - location of interest cross covariances for center/radius
# pen    - penalty parameter for constrained optimization
# A      - vector of length 3 with weights for generalized L2 distance.
# threshold - The magnitude representing "effective 0"
# Outputs:
# lamUp  - and updated version of the lambda vector.
# Notes:
# Implementation of the penalty assumes the square difference of the
# contraint equation is being used.
# ============================================================================
nrStep_long_R = function(lam, pCovC, pCovR, pCovCR,
                         lCovC, lCovR, lCovCR, pen, A, len){
  l = length(lam)

  cGrad = (matrix(lam, nrow = 1) %*% pCovC) - lCovC # No need to store entries.

  rGrad = (matrix(abs(lam), nrow = 1) %*% pCovR) -
    lCovR

  calc.cr <- matrix(lam, nrow = 1) %*% pCovCR

  crGrad = (matrix(abs(lam), nrow = 1) %*% pCovCR) +
    sign(lam)*calc.cr -
    (sign(lam) + 1)*lCovCR

  penalty <- (sum(abs(lam)) - 1)/pen

  grad = 2*( A[1]*cGrad + sign(lam)*(A[2]*rGrad + penalty) +
               A[3]*crGrad - (pen/(lam*len*len)) )

  grad = matrix(grad, ncol = 1)

  ## Define the Hessian Matrix
  # First, define the matrix of sign functions to multiply with the pwCovR
  lamSign = sign(lam) %*% matrix(sign(lam), nrow = 1)

  lamSign.cr <- outer(sign(lam), sign(lam), "+")

  # Use Quadratic Approximation
  hessian = 2*( A[1]*pCovC + A[2]*lamSign*pCovR +
                  A[3]*lamSign.cr*pCovCR + lamSign/pen ) +
    2*diag( pen/(lam*lam*len*len) )

  b <- (hessian %*% lam) - grad

  lamUp = base::solve(hessian, b)

  return(lamUp)
}

# ============================================================================
### Function to update the value of lambda according to
# the Newton-Raphson Method (ordinary kriging)
# Inputs:
# lam - current values of the weights (lambdas)
# pCovC  - pairwise covariances for the interval centers (C^C(u_i - u_j))
# pCovR  - pairwise covariances for radii
# pCovCR - pairwise cross covariance between center and radius
# lCovC  - location of interest covariances for centers
# lCovR  - location of interest covariances for radii
# lCovCR - location of interest cross covariances for center/radius
# r      - penalty parameter for constrained optimization
# A      - vector of length 3 with weights for generalized L2 distance.
# threshold - Threshold for which lambdas with magnitude less than are set to 0.
#             Calculated within newRap_2.
# len       - Original length of lambda vector, calculated in newRap_2.
# Outputs:
# lamUp  - and updated version of the lambda vector.
# Notes:
# Implementation of the penalty assumes the square difference of the
# contraint equation is being used.
# ============================================================================
nrstep_2_R = function(lam, pCovC, pCovR, pCovCR, lCovC, lCovR, lCovCR, r,
                      A, threshold, len){

  # First, calculate the gradient
  t.grad <- 2*((t(lam) %*% (A[1]*pCovC + A[2]*pCovR + 2*A[3]*pCovCR)) -
                 (A[1]*lCovC + A[2]*lCovR + 2*A[3]*lCovCR))
  pen.grad <-  -( r/( (lam + threshold)
  ) ) + (2*(sum(lam)-1))/(r*len)

  grad <- matrix(t.grad + pen.grad, ncol = 1)

  t.hessian <- 2*(A[1]*pCovC + A[2]*pCovR + 2*A[3]*pCovCR)
  pen.hessian <- diag(r/( (lam + threshold)^2 ) ) + (2/(r*len))

  hessian <- t.hessian + pen.hessian

  b <- hessian %*% lam - grad

  lamUp = solve(hessian, b)

  return(lamUp)
}
