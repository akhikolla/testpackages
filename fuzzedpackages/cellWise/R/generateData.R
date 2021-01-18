

generateCorMat <- function(d, corrType = "ALYZ", CN = 100, seed = NULL) {
  
  if (!corrType %in% c("A09", "ALYZ")) {
    stop("corrType should be either \"A09\" or \"ALYZ\"")
  }
  
  if (is.numeric(seed)) {
    set.seed(seed)
  }
  
  
  if (corrType == "ALYZ") {
    #generate random correlation matrix with given condition number
    l <- sort(runif(d - 2, 1, CN))
    L <- diag(c(1, l, CN))
    Y <- matrix(rnorm(d * d), ncol = d)
    U <- eigen(t(Y) %*% Y)$vectors 
    Sig_0 <- U %*% L %*% t(U)
    
    Sig_temp <- Sig_0
    while (abs(kappa(cov2cor(Sig_temp), exact = TRUE) - CN) > 1e-4) { 
      R_0_temp <- cov2cor(Sig_temp)
      U_0_temp <- eigen(R_0_temp)$vectors[, d:1]
      L_0_temp <- diag(eigen(R_0_temp)$values[d:1])
      L_0_temp[d,d] <- L_0_temp[1,1] * CN
      Sig_temp <- U_0_temp %*% L_0_temp %*% t(U_0_temp)
    }
    R_0 <- cov2cor(Sig_temp) #the correlation matrix
  } else if (corrType == "A09") { 
    columns <- matrix(data = (1:d), nrow = d, ncol = d, byrow = TRUE)
    rows    <- matrix(data = (1:d), nrow = d, ncol = d, byrow = FALSE)
    R_0     <- matrix(-0.9, nrow = d, ncol = d)
    R_0     <- R_0^(abs(columns - rows)) # A09 correlation matrix
  } 
  return(R_0)
}



generateData = function(n, d,
                        mu, Sigma,
                        perout, gamma,
                        outlierType = "casewise", 
                        seed = NULL) {
  
  # outlierType:
  #   casewise: eplacement value: gamma*smallest eigenvector where smallest eigenv has md 1
  #   cellwisePlain: cellwise naive
  #   cellwiseStructured: cellwise smallest eigenvector
  #  both
  
  if (!outlierType %in% c("casewise", "cellwisePlain", "cellwiseStructured",
                         "both")) {
    stop("outlierType should be one of \"casewise\", \"cellwisePlain\",
    \"cellwiseStructured\" or \"both\"")
  }
  
  if (is.numeric(seed)) {
    set.seed(seed)
  }
  
  # generating random data
  X = MASS::mvrnorm(n, mu = mu, Sigma = Sigma)
  indcells <- c()
  indrows  <- c()
  #adding contamination
  if (perout > 0) {
    if (outlierType == "casewise") {
      replacement <- eigen(Sigma)$vectors[, d] /
        sqrt(mahalanobis(eigen(Sigma)$vectors[, d], mu, Sigma))
      replacement <- gamma * replacement * sqrt(d) * d
      ind         <- 1:(floor(perout * n)) 
      X[ind, ]    <- matrix(replacement, nrow = length(ind), ncol = d, byrow = TRUE)
      indrows <- ind 
    } else if (outlierType == "cellwisePlain") {
      ind <- replicate(d, sample(1:n, perout * n, replace = FALSE))
      ind <- as.vector(t(t(ind) + n * (0:(d - 1))))
      X[ind] <- gamma
      indcells <- ind
    } else if (outlierType == "cellwiseStructured") {
      ind <- replicate(d, sample(1:n, perout * n, replace = FALSE))
      ind <- as.vector(t(t(ind) + n * (0:(d - 1))))
      W   <- array(0, dim(X)); W[ind] <- 1
      for (i in 1:n) {
        continds <- which(W[i, ] == 1)
        if (length(continds) > 0) {
          eigen_out <- eigen(Sigma[continds, continds])$vectors
          replacement <- eigen_out[, length(continds)] /
            sqrt(mahalanobis(eigen_out[, length(continds)], mu[continds], Sigma[continds, continds]))
          X[i, continds] <- replacement * gamma * sqrt(length(continds))
          indcells <- c(indcells, i + (continds - 1) * n)
        }
      }
    } else  if (outlierType == "both") {
      replacement <- eigen(Sigma)$vectors[, d] /
        sqrt(mahalanobis(eigen(Sigma)$vectors[, d], mu, Sigma))
      replacement <- gamma * replacement * d * sqrt(d)
      ind         <- 1:(floor(perout / 2 * n))
      X[ind, ]    <- matrix(replacement, nrow = length(ind), ncol = d, byrow = TRUE)
      indrows <- ind
      
      # Now add cells
      startind <- (floor(perout / 2 * n) + 1) # first row which is not a casewise outliers
      ind <- replicate(d, sample(startind:n, ceiling(perout / 2 * n),
                                 replace = FALSE))
      ind <- as.vector(t(t(ind) + n * (0:(d - 1))))
      W   <- array(0, dim(X)); W[ind] <- 1
      for (i in startind:n) {
        continds <- which(W[i, ] == 1)
        if (length(continds) > 0) {
          
          eigen_out <- eigen(Sigma[continds, continds])$vectors
          replacement <- eigen_out[, length(continds)] /
            sqrt(mahalanobis(eigen_out[, length(continds)], mu[continds], Sigma[continds, continds]))
          X[i, continds] <- replacement * gamma * sqrt(length(continds))
          indcells <- c(indcells,  i + (continds - 1) * n)
        }
      }
    }  
  }
  return(list(X = X, indcells = indcells,
              indrows = indrows))
}
