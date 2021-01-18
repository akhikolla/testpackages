#' Simulate a stationary Vector Autoregressive (VAR) time series.
#'
#' @param T An integer giving the number of timepoints.
#' @param A A d x d transition matrix.
#' @param Sigma A d x d innovation covariance matrix.
#' @keywords var simulate
#'
#' @examples
#'
#' theta    <- diag(c(.7,.8,.9,.6,.7,.9))
#' data     <- t(var_sim(100, theta, diag(.1,6)))
#' datalag  <- embed(data, 2)
#' b        <- datalag[,1:6]
#' A        <- datalag[,7:12]
#' A_est    <- fista_sparse(A, b, 1, theta, niter = 10, backtrack = TRUE)$out.x
#' var_forecast(t(b), 2, A_est)
#'
#' @export
var_sim = function(T, A, Sigma){
  k    <- dim(A)[1]
  p    <- dim(A)[2]/k
  burn <- 500
  inno <- MASS::mvrnorm(n=T+burn, rep(0, k), Sigma)
  init <- MASS::mvrnorm(n=p, rep(0, k), Sigma)
  init <- matrix(init, nrow=p)
	j    <- 1
	id   <- seq(from= j+p-1, to = j, by=-1)
  Y    <-  matrix(0, (T+burn), k)
  for(r in 1:(T+burn)){
  	Y[r,] <-  A %*% as.vector(t(init[id,])) + inno[r,]
     init <- rbind(init[-1,], Y[r,])
  }
  return(t(Y[-(1:burn),]))
}
