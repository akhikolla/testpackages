#' @title Simulate data
#' @description the function that generates the simulation data set
#' @param n number of observations
#' @param p number of variables
#' @param a numveric vector of length k+1 that contains the correlations between the variables in each block (the last block contains uncorrelated variables); by default is (0.9, 0.5, 0)
#' @param bLength the number of correlated variables in the first k blocks
#' @param SD  numveric vector of length k+1 that contains the standard deviation of the variables in each block (the last block contains uncorrelated variables); by default is (10, 5, 2)
#' @param eps proportion of outliers, default is 0
#' @param seed random seed used to simulate the data
#' @param eta parameter that contols the outlyingness, default is 25
#' @param setting type of outliers: \code{setting}="1" generates the outliers which are outlying in the first two variables in the second block; \code{setting}="2" generates
#' score outliers; \code{setting}="3" generates the orthogonal outliers which are easy to detect (the setting used in Hubert, et al (2016)); default is "3"
#' @param vc controls the direction of the score outliers within the PC subspace, default is NULL
#' @return a list with components \cr
#' \item{data}{generated data matrix}
#' \item{ind}{row indices of outliers}
#' \item{R}{Correlation matrix of the data}
#' \item{Sigma}{Covariance matrix of the data}
dataSim <- function(n = 200, p = 20,  bLength = 4, a = c(0.9, 0.5, 0), SD = c(10, 5, 2),
                    eps = 0, eta = 25, setting = "3", seed = 123, vc = NULL){

  if (p < 8) {
    stop("Dimension needs to be at least 8")
  }
  R <- genMat(p = p, a = a, bLength = bLength)
  k <- length(a) - 1
  if (k < 2) {
    stop("a should have length at least 3.")
  }
  if (length(SD) != length(a)) {
    stop("SD have the same length as a")
  }
  sdx <- rep(0, p)
  for (i in 1:k) {
    sdx[(i * bLength) - (0:(bLength - 1))] <- rep(SD[i],
                                                  bLength)
  }
  sdx[(k * bLength + 1):p] <- rep(SD[k + 1], p - k * bLength)
  SDx <- diag(sdx[1:(k*bLength)])
  Sigma <- diag(p)
  Sigma[(k*bLength+1):p,(k*bLength+1):p] <- diag(SD[k+1]^2,(p-k*bLength),(p-k*bLength))
  Sigma[1:(k*bLength),1:(k*bLength)] <- SDx %*% R[1:(k*bLength),1:(k*bLength)] %*% SDx

  X <- matrix(NA,n,p)
  set.seed(seed)
  X[,1:(k*bLength)] <- mvtnorm::rmvnorm(n = n, mean = rep(0,(k*bLength)), sigma = Sigma[1:(k*bLength),1:(k*bLength)])
  X[,(k*bLength+1):p] <- matrix(rnorm((n*(p-(k*bLength))),mean=0,sd=SD[k + 1]),n,(p-(k*bLength)))
  X <- X + matrix(rnorm(n*p),n,p)


  if (eps >= 1) {
    eps <- eps/100
  }
  n_out <- floor(eps * n)

  if (eps > 0){
    index <- sample(1:n, n_out, replace = FALSE)
    if (setting == "1"){
      X[index,(bLength+1)] <- rnorm(n_out,mean = -eta,sd= 1)
      X[index,(bLength+2)] <- rnorm(n_out,mean = -eta,sd= 1)
    }
    if (setting == "2"){
      v <- matrix(0, nrow = p,ncol = k)
      v[1:(k*bLength),1:k] <- svd(R[1:(k*bLength),1:(k*bLength)])$v[,1:k]
      if (k == 2){
        vc <- c(-1,1)
      }else{
        if (is.null(vc)){
          print("This vector needs to be specified by the user to generate the score outliers when k is larger than 2!")
        }
      }

      X[index,1:(k*bLength)] <- mvtnorm::rmvnorm(n_out, mean = eta * v[1:(k*bLength),] %*% vc , sigma = Sigma[1:(k*bLength),1:(k*bLength)])
      X[index,(k*bLength+1):p] <- matrix(rnorm((n_out*(p-k*bLength)),mean = 0, sd = 2),nrow=n_out,ncol=(p-k*bLength))
      X[index,] <- X[index,] + matrix(rnorm(n_out*p),nrow = n_out,ncol = p)
    }
    if (setting == "3"){
      mu_out <- eta * c(c(0, -4, 4, 2, 0, 4, -4, 2), rep(c(3, -3), length.out = p - 8))
      sigma2 <- sqrt(20)
      for (j in 1:p){
        X[index,j] <- rnorm(eps*n,mean=mu_out[j],sd=sigma2)
      }
    }
  }else{
    index = 0
  }
  return(list(data = X, ind = index, R = R, Sigma = Sigma))
}

genMat <- function (a = c(0.9, 0.5, 0), p = 10, bLength = 4){
  k <- length(a) - 1
  if (p <= bLength * k) {
    stop("Dimension needs to be at least bLength*length(a)")
  }
  if (any(abs(a) > 1)) {
    stop("Off-diagonal elements need to be smaller than 1 in absolute value!")
  }
  X <- diag(p)
  for (t in 0:(k - 1)) {
    for (i in (t * bLength + 1):((t + 1) * bLength)) {
      for (j in (t * bLength + 1):((t + 1) * bLength)) {
        if (j != i) {
          X[i, j] = a[t + 1]
        }
      }
    }
  }
  if (a[k + 1] != 0) {
    for (i in (k * bLength + 1):p) {
      for (j in (k * bLength + 1):p) {
        if (j != i) {
          X[i, j] = a[k + 1]
        }
      }
    }
  }
  return(X)
}
