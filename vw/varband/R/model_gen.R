#' Generate random samples.
#'
#' Generate \code{n} random samples from multivariate Gaussian distribution N(0, (L^TL)^{-1})
#'
#' @param L p-dimensional inverse Cholesky factor of true covariance matrix.
#' @param n number of samples to generate.
#' @return returns a n-by-p matrix with each row a random sample generated.
#' @examples
#' set.seed(123)
#' true <- varband_gen(p = 50, block = 5)
#' x <- sample_gen(L = true, n = 100)
#' @export
sample_gen<-function(L, n){
  #### This function generates normal random vectors
  # Y ~ N(0, (L^T L)^{-1}) = N(0, L^{-1}L^{-T})
  # Sigma = L^{-1}L^{-T}
  p <- nrow(L)
  # Z~N(0,1_p)
  # Z<-mvrnorm(n, rep(0,p), diag(rep(1,p)))
  Z <- matrix(rnorm(n * p), n, p)
  # X[i,]%*%L^{-T} ~ N (0, L^{-1}L^{-T}) = N(0, Sigma)
  X <- t(solve(L, t(Z)))
  return(X)
}

#' Generate an autoregressive model.
#'
#' Generate lower triangular matrix with strict bandwidth. See, e.g., Model 1 in the paper.
#'
#' @param p the dimension of L
#' @param phi_vec a K-dimensional vector for off-diagonal values
#' @return a p-by-p strictly banded lower triangular matrix
#' @examples
#' true_ar <- ar_gen(p = 50, phi = c(0.5, -0.4, 0.1))
#' @export
#' @import stats
ar_gen<-function(p, phi_vec){
  set.seed(123)
  ii <- toeplitz(1:p)
  K <- length(phi_vec) + 1
  L <- ii <= K

  for (k in seq(2, K)) L[ii == k] <- phi_vec[k - 1]
  diag(L) <- rep(0, p)
  L <- L - upper.tri(L) * L

  L <- diag(rep(1, p)) - L
  # Note that Omega = L^T L
  # NEWLY ADDED
  L <- diag(1/runif(p, 2, 5)) %*% L
  return(L)
}

#' Generate a model with variable bandwidth.
#'
#' Generate lower triangular matrix with variable bandwidth. See, e.g., Model 2 and 3 in the paper.
#' @param p the dimension of L
#' @param block the number of block diagonal structures in the resulting model, assumed to divide p
#' @return a p-by-p lower triangular matrix with variable bandwidth
#' @examples
#' set.seed(123)
#' # small block size (big number of blocks)
#' true_small <- varband_gen(p = 50, block = 10)
#' # large block size (small number of blocks)
#' true_large <- varband_gen(p = 50, block = 2)
#' @export
varband_gen <- function(p, block = 10){
  set.seed(123)
  L <- matrix(0, p, p)
  stopifnot(p%%block == 0)
  block_size <- p / block
  for(k in seq(block)){
    L[((k-1) * block_size + 1):(k * block_size),
      ((k-1)*block_size + 1):(k * block_size)] <- block_gen(block_size)
  }
  diag(L) <- 1
  L <- diag(1/runif(p, 2, 5)) %*% L
  return(L)
}

block_gen <- function(block_size){
  smallL <- matrix(0, block_size, block_size)
  for(i in seq(2, block_size)){
    flag <- rbinom(1, 1, 0.5)
    if (flag == 1){
      Ji <- sample(i - 1, 1)
      sign <- ((runif(i - 1, 0, 1) > 0.5) - 1/2 ) * 2
      smallL[i, 1 : (i - 1)] <- sign * runif(i - 1, 0.1, 0.4)
      smallL[i, 1 : Ji] <- 0
    }
    else
      smallL[i, 1:(i - 1)] <- 0
  }
  return(smallL)
}

#' Generate a model with block-diagonal structure
#'
#' @param p the dimension of L
#' @return a p-by-p lower triangular matrix with block-diagonal structure from p/4-th row to 3p/4-th row
#' @examples
#' set.seed(123)
#' true_L_block_diag <- block_diag_gen(p = 50)
#' @export
block_diag_gen <- function(p){
  L <- matrix(0, p, p)
  set.seed(123)
  for(i in seq(p / 4 + 2, p * 3 / 4)){
    sign <- ((runif(i - p / 4 - 1, 0, 1) > 0.5) - 1/2 ) * 2
    L[i, (p / 4 + 1) : (i - 1)] <- sign * runif(i - p / 4 - 1, 0.1, 0.2)
  }
  diag(L) <- 1
  L <- diag(1/runif(p, 2, 5)) %*% L
  return(L)
}
