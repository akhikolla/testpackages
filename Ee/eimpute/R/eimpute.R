#' @title Efficiently impute missing values for a large scale matrix
#'
#' @description Fit a low-rank matrix approximation to a matrix with missing values.
#' The algorithm iterates like EM: filling the missing values with the current guess,
#' and then approximating the complete matrix via truncated SVD.
#'
#' @aliases eimpute
#'
#' @param x an \eqn{m} by \eqn{n} matrix with \code{NA}s.
#' @param r the rank of low-rank matrix for approximating \code{x}
#' @param svd.method a character string indicating the truncated SVD method.
#' If \code{svd.method = "rsvd"}, a randomized SVD is used,
#' else if \code{svd.method = "tsvd"}, standard truncated SVD is used.
#' Any unambiguous substring can be given. Default \code{svd.method = "tsvd"}.
#' @param thresh convergence threshold, measured as the relative change in the Frobenius norm between two successive estimates.
#' @param maxit	 maximal number of iterations.
#' @param override logical value indicating whether the observed elements in \code{x} should be overwritten by its low-rank approximation.
#' @param control a list of parameters that control details of standard procedure, See \link{biscale.control}.
#' @param ... arguments to be used to form the default control argument if it is not supplied directly.
#'
#'
#' @return A list containing the following components
#' \item{\code{x.imp}}{the matrix after completion.}
#' \item{\code{rmse}}{the relative mean square error of matrix completion, i.e., training error.}
#' \item{\code{iter.count}}{the number of iterations.}
#'
#' @references Rahul Mazumder, Trevor Hastie and Rob Tibshirani (2010) Spectral Regularization Algorithms for Learning Large Incomplete Matrices, Journal of Machine Learning Research 11, 2287-2322
#' @references Nathan Halko, Per-Gunnar Martinsson, Joel A. Tropp (2011) Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate Matrix Decompositions, Siam Review Vol. 53, num. 2, pp. 217-288
#'
#' @rdname eimpute
#' @useDynLib eimpute, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
#' @examples
#' ################# Quick Start #################
#' m <- 100
#' n <- 100
#' r <- 10
#' x_na <- incomplete.generator(m, n, r)
#' head(x_na[, 1:6])
#' x_impute <- eimpute(x_na, r)
#' head(x_impute[["x.imp"]][, 1:6])
#' x_impute[["rmse"]]
eimpute <- function(x, r, svd.method = c("tsvd", "rsvd"),
                    thresh = 1e-05, maxit = 100, override = FALSE, control = list(...), ...){
  m <- nrow(x)
  n <- ncol(x)

  x_sd <- biscale(x, control = control)

  x_train <- x_sd[[1]]
  # x_train <- x
  ind_ob <- which(!is.na(x_train), arr.ind = TRUE)
  x_ob <- x_train[ind_ob]
  ind <- ind_ob - 1

  svdm <- match.arg(svd.method)
  type <- ifelse(svdm == "tsvd", 1, 2)


  Z_temp <- kkt_fix(ind, x_ob, m, n, r, maxit, thresh, type)
  Z.fit <- Z_temp[[1]] * (x_sd[[4]] %*% t(x_sd[[5]])) + matrix(rep(x_sd[[2]], n), nrow = m) + t(matrix(rep(x_sd[[3]], m), nrow = n))


  if (!override) {
    Z.fit[ind_ob] <- x_ob
  }

  list(x.imp = Z.fit, rmse = Z_temp[[3]], iter.count = Z_temp[[2]])
}


