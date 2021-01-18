THRESH <- 1e-10
# #' Refit a row problem with given support
# #'
# #' @param r row index
# #' @param ind the index set of support
# #' @param S p-by-p sample covariance matrix
# #' @param delta a nonnegative tuning parameter for the ridge penalty for numerical stability, default is 0.01
refit_row <- function(r, ind, S, delta = 0.01){
  # Refitting by plain MLE with
  # bandwidth k = r-1-J given our estimator
  p <- ncol(S)
  res <- rep(0, r)
  if(is.null(ind)){
    res[r] <- 1/sqrt(S[r, r])
  }
  else{
    # if J < r-1
    # ind <- (J+1):(r-1)
    # If S[ind,ind] is not invertible
    # add a little ridge penalty to that
    tmpVec <- solve((S[ind, ind] +
                       delta * diag(rep(1, length(ind)))),
                    S[ind, r])
    res[r] <- 1/sqrt(S[r, r] + delta -
                       crossprod(S[r, ind], tmpVec))
    res[ind] <- -tmpVec * res[r]
  }
  return(res)
}

# #' Refit the estimate of lower triangular matrix L with given support
# #'
# #' @param S p-by-p sample covariance matrix
# #' @param mat p-by-p estimate of lower triangular matrix L
# #' @param delta a nonnegative tuning parameter for the ridge penalty for numerical stability, default is 0.01
refit_matrix <- function(S, mat, delta = 0.01){
  p <- ncol(S)
  refit <- matrix(0, p, p)
  refit[1, 1] <- 1/sqrt(S[1, 1])
  for(r in seq(2, p)){
    ind <- c()
    for(j in seq(1,r-1)){
      if(abs(mat[r, j]) >= THRESH){
        ind <- c(ind, j)
      }
    }
    refit[r, 1:r] <- refit_row(r = r, ind = ind,
                               S = S, delta = delta)
  }
  return(refit)
}

# #' Refit a path of estimates of lower triangular matrix L with given support
# #'
# #' @param  S: p-by-p sample covariance matrix
# #' @param  path: a list of p-by-p estimate of lower triangular matrix L along a path of tuning parameters
# #' @param delta: a nonnegative tuning parameter for the ridge penalty for numerical stability, default is 0.01
refit_path <- function(S, path, delta = 0.01){
  p <- dim(path)[1]
  nlam <- dim(path)[3]
  refit <- array(NA, c(p, p, nlam))
  for(i in seq(nlam)){
    refit[, , i] <- refit_matrix(S, path[, , i], delta)
  }
  return(refit)
}
