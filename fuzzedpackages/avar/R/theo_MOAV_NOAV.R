#' @title Non-stationary Maximal-overlapping Allan Variance
#' @description
#' Calculation of the theoretical Maximal-overlapping Allan variance for constant-mean non-stationary time series data.
#' @export
#' @param n        An \code{integer} indicating the length of each vector of consecutive observations considered for the average.
#' @param covmat   A \code{matrix} indicating the T-by-T covariance matrix of the time series with length T.
#' @return A \code{field <numeric>} that is the theoretical Maximal-overlapping Allan variance for constant-mean non-stationary time series data.
#' @details
#' This calculation of Maximal-overlapping Allan variance is based on the definition on "A Study of the Allan Variance for Constant-Mean Non-Stationary Processes" by Xu et al., 2017, IEEE Signal Processing Letters, 24(8): 1257–1260.
#' Here n is an integer larger than 1 and smaller than \eqn{floor\left(log_2 \left(dim\left(covmat\right)[1]\right)\right)-1}{floor(log2(dim(T)[1]))-1}.
#' @author Haotian Xu
#' @examples
#' \donttest{
#' set.seed(999)
#' Xt = arima.sim(n = 100, list(ar = 0.3))
#' avar(Xt, type = "to")
#'
#' a = matrix(rep(0, 1000^2), nrow = 1000)
#' for (i in 1:1000){
#'   a[,i] = seq(from = 1 - i, length.out = 1000)
#' }
#' a.diag = diag(a)
#' a[upper.tri(a,diag=TRUE)] = 0
#' a = a + t(a) + diag(a.diag)
#' covmat = 0.3^a
#' sapply(1:8, function(y){MOAV(2^y, covmat)})
#' }

MOAV = function(n, covmat){
  # k: starting position of filter
  # n: length of each block
  # covmat: T by T covariance matrix
  # T: number of obversations
  mat = matrix(rep(0, n^2), nrow = n)
  T = dim(covmat)[1]
  for (k in 1:(T-2*n+1)){
    v1 = seq.int(from = k, to = k+n-1)
    v2 = seq.int(from = k+n, to = k+2*n-1)

    sigma_mat1 = covmat[v1, v1]
    sigma_mat2 = covmat[v2, v2]
    gamma_mat = covmat[v1, v2]

    mat = mat + sigma_mat1 + sigma_mat2 - 2*gamma_mat
  }
  sum(mat)/(n^2*2*(T-2*n+1))
}


#' @title Non-stationary Non-overlapping Allan Variance
#' @description
#' Calculation of the theoretical Non-overlapping Allan variance for constant-mean non-stationary time series data.
#' @export
#' @param n        An \code{integer} indicating the length of each vector of consecutive observations considered for the average.
#' @param covmat   A \code{matrix} indicating the T-by-T covariance matrix of the time series with length T.
#' @return A \code{field <numeric>} that is the theoretical Non-overlapping Allan variance for constant-mean non-stationary time series data.
#' @details
#' This calculation of Non-overlapping Allan variance is based on the definition on "A Study of the Allan Variance for Constant-Mean Non-Stationary Processes" by Xu et al., 2017, IEEE Signal Processing Letters, 24(8): 1257–1260.
#' Here n is an integer larger than 1 and smaller than \eqn{floor\left(log_2 \left(dim\left(covmat\right)[1]\right)\right)-1}{floor(log2(dim(T)[1]))-1}.
#' @author Haotian Xu
#' @examples
#' \donttest{
#' set.seed(999)
#' Xt = arima.sim(n = 100, list(ar = 0.3))
#' avar(Xt, type = "to")
#'
#' a = matrix(rep(0, 1000^2), nrow = 1000)
#' for (i in 1:1000){
#'   a[,i] = seq(from = 1 - i, length.out = 1000)
#' }
#' a.diag = diag(a)
#' a[upper.tri(a,diag=TRUE)] = 0
#' a = a + t(a) + diag(a.diag)
#' covmat = 0.3^a
#' sapply(1:8, function(y){NOAV(2^y, covmat)})
#' }
#'

NOAV = function(n, covmat){
  # k: index of blocks
  # n: length of each block
  # covmat: T by T covariance matrix
  # T: number of obversations
  mat = matrix(rep(0, n^2), nrow = n)
  m = floor(dim(covmat)[1]/(2*n))
  for (k in 1:m){
    v1 = seq.int(from = (2*k-1)*n+1, to = 2*k*n)
    v2 = seq.int(from = (2*k-2)*n+1, to = (2*k-1)*n)

    sigma_mat1 = covmat[v1, v1]
    sigma_mat2 = covmat[v2, v2]
    gamma_mat = covmat[v1, v2]

    mat = mat + sigma_mat1 + sigma_mat2 - 2*gamma_mat
  }
  sum(mat)/(n^2*2*m)
}
