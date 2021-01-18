# Functions needed for AMUSE and SOBI ladles

# function to measure the bootstrap variation of the eigenvectors
fi <- function(EVboot, EVdata, rank) {
  fni <- numeric(rank)
  for (ii in 1:rank) {
    fni[ii] <- det(crossprod(EVdata[, 1:ii], EVboot[, 1:ii]))
  }
  1 - abs(fni)
}

# Mfunction for AMUSE
MAmuse <- function(X, k) {
  n <- nrow(X)
  prep <- .Call("PREPBSS", X, n, PACKAGE = "tsBSS") #calling the function PREPBSS
  Y <- prep$Y 
  
  M <- crossprod(Y[1:(n - k), ], Y[(k + 1):n, ])/(n - k)
  M.sym <- (M + t(M))/2
  crossprod(M.sym)
}

# Bootstrapping function for AMUSE
AMUSEbootLADLE <- function(X, EVdata, tau, rank) {
  Mboot <- MAmuse(X, k = tau)
  EVboot <- .Call("EIGEN", Mboot, PACKAGE = "tsBSS")$vectors
  fi(EVboot, EVdata, rank)
}

# Mfunction for SOBI
MSobi <- function(X, k_set) {
  n <- nrow(X)
  p <- ncol(X)
  prep <- .Call("PREPBSS", X, n, PACKAGE = "tsBSS") #calling the function PREPBSS
  Y <- prep$Y 

  M_array <- array(0, dim = c(p, p, length(k_set)))
  for (t in 1:length(k_set)) {
    M_array[ , , t] <- crossprod(Y[1:(n - k_set[t]), ], Y[(k_set[t] + 1):n, ])/(n - k_set[t])
    M_array[ , , t] <- (M_array[ , , t] + t(M_array[, , t]))/2
  }
  M_array
}

# Bootstrapping function for SOBI
SOBIbootLADLE <- function(X, EVdata, tau, rank, maxiter, eps) {
  Mboot <- MSobi(X, k_set = tau)
  frjdboot <- JADE::frjd.int(Mboot, maxiter = maxiter, eps = eps)
  Dfrjd <- diag(apply(frjdboot$D^2, 1:2, sum))
  EVboot <- frjdboot$V[ , order(Dfrjd, decreasing = TRUE)]
  fi(EVboot, EVdata, rank)
}
