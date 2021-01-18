#' @title Hard Graphical Thresholding Algorithm
#' @description Estimates sparse inverse covariance matrix.
#'
#' @param x There are 2 options: (1) \code{x} is an \eqn{n} by \eqn{p} data matrix; (2) a \eqn{p} by \eqn{p} sample covariance matrix. The program automatically identifies the input matrix by checking the symmetry. (\eqn{n} is the sample size and \eqn{p} is the dimension.)
#' @param size A non-negative integer for determining the model size, i.e., the number of non-zero off-diagonal entries in the upper-triangular precision matrix,
#' which is also the number of edges in the graph. \code{size} must range from 0 to \eqn{(p^2 - p) / 2}.
#' @param active.entry Pre-determined non-zero off-diagonal entries positions of the precision matrix. Default: \code{active.entry = NULL}.
#' @param bcd.opt A list of options that control details of the block coordinate descent algorithm.
#'
#' @return A list with following components:
#' \item{\code{Omega}}{Estimated inverse covariance matrix.}
#' \item{\code{active.entry}}{The position of the non-zero off-diagonal entries of \code{Omega} in the upper-triangular part.}
#'
#' @details
#' Hard Graphical Thresholding (HGT) algorithm proceeds by thresholding the sample correlation matrix and
#' estimating the inverse covariance matrix with block coordinate descent algorithm.
#' HGT algorithm could recover the inverse covariance matrix given model size or given active entries.
#' When active entries are given directly, model fitting is the so-called covariance selection.
#'
#' @note
#' Either \code{size} or \code{active.entry} should be specified when function \code{hgt} is called.
#' If both arguments are given, \code{size} would be omitted and the inverse covariance matrix would be estimated based on the given \code{active.entry}.
#'
#' If arguments \code{active.entry} is specified, only one of the entries in symmetric positions should be given.
#'
#' @references Luo, Shikai, Rui Song, and Daniela Witten (2014). Sure Screening for Gaussian Graphical Models. arXiv preprint arXiv:1407.7819. URL https://arxiv.org/abs/1407.7819.
#' @references Dempster, A.P. (1972). Covariance Selection. Biometrics, 28(1), 157-175. doi:10.2307/2528966
#'
#' @export
#'
#' @examples
#' library(gif)
#'
#' data("ar1")
#' p <- 100
#' non_zero_num <- sum(ar1[["Omega"]] != 0) - p
#' res <- hgt(ar1[["x"]], size = non_zero_num / 2)
hgt <- function(x, size, active.entry = NULL, bcd.opt = list("max.iter" = 10, "eps" = 1e-3)) {
  stopifnot(!anyNA(x))
  stopifnot(!missing(size) | !is.null(active.entry))

  p <- ncol(x)

  if(!isSymmetric(x)) {
    s <- var(x)
  } else {
    s <- x
  }

  max.iter <- bcd.opt[["max.iter"]]
  eps <- bcd.opt[["eps"]]
  stopifnot(round(max.iter, 0) == max.iter)
  stopifnot(max.iter > 0)
  stopifnot(eps > 0)

  if(!is.null(active.entry)) {
    if(ncol(active.entry) != 2) {
      stop("Arguments active.entry is invalid!")
    }
    if(any(active.entry[, 1] == active.entry[, 2])) {
      stop("Digonal entries should not be included in arguments active.entry!")
    }
    active_index <- as.vector(active.entry)
    if(any(active_index > p) | any(active_index < 1)) {
      stop("Arguments active.entry is invalid since subscript out of bounds!")
    }
    double_active_entry <- rbind(active.entry, active.entry[, 2:1])
    if(any(duplicated(double_active_entry))) {
      stop("Symmetric positions of off-diagonal entries should be considered once in arguments active.entry!")
    }

    act_set <- active.entry - 1
    Omega <- bcd(s, act_set, iter_max = max.iter, eps = eps)
  } else {
    if (size > (p ^ 2 - p) / 2) {
      stop("Arguments size is too large!")
    }
    if (size < 0) {
      stop("Arguments size should be greater than zero!")
    }

    S_init <- cov2cor(s)
    diag(S_init) <- 0
    S_vec <- as.vector(abs(S_init))
    lambda <- sort(S_vec, decreasing = TRUE)[2 * size + 1]
    S_init[abs(S_init) <= lambda] <- 0
    non_zero_index <- which(as.matrix(S_init) != 0, arr.ind = TRUE)
    active.entry <- non_zero_index[which(non_zero_index[,1] < non_zero_index[,2]),]
    act_set <- active.entry - 1
    Omega <- bcd(s, act_set, iter_max = max.iter, eps = eps)
  }

  return(list(Omega = Omega, active.entry = active.entry))
}

#' @title Soft Graphical Thresholding Algorithm
#' @description Estimates a sparse inverse covariance matrix using the closed form solution of graphical lasso under acyclic graph structure.
#'
#' @inheritParams hgt
#' @param lambda The regularization parameter for graphical lasso.
#'
#' @return A list with following components:
#' \item{\code{Omega}}{Estimated inverse covariance matrix.}
#' \item{\code{active.entry}}{The position of the non-zero entries of \code{Omega}.}
#' \item{\code{is.acyclic}}{The boolean flag of whether the detected graph structure is acyclic or not.}
#'
#' @details
#' Soft Graphical Thresholding (SGT) algorithm proceeds by thresholding the sample covariance matrix and
#' estimating the inverse covariance matrix with a closed-form formula.
#' If the graph structure detected by the thresholding procedure is acyclic,
#' then the estimation is equivalent to the solution of graphical lasso.
#'
#' @note
#' Either \code{lambda} or \code{size} should specified when function \code{sgt} is called.
#' If both arguments are given, only \code{lambda} would be considered.
#'
#' @references Fattahi, Salar, and Somayeh Sojoudi. Graphical Lasso and Thresholding: Equivalence and Closed-form Solutions. Journal of Machine Learning Research 20.10 (2019): 1-44. doi: 10.5555/3322706.3322716
#'
#' @export
#'
#' @examples
#' library(gif)
#'
#' data("ar1")
#' res <- sgt(ar1[["x"]], lambda = 0.01)
sgt <- function(x, lambda, size = NULL) {
  stopifnot(!missing(lambda) | !is.null(size))

  stopifnot(!anyNA(x))
  if(!isSymmetric(x)) {
    s <- cov(x)
  } else {
    s <- x
  }
  p <- ncol(x)

  if(!missing(lambda)) {
    S_init <- s
    diag(S_init) <- 0
  } else {
    if (size > (p ^ 2 - p) / 2) {
      stop("Arguments size is too large!")
    }
    if (size < 0) {
      stop("Arguments size should be greater than zero!")
    }

    S_init <- s
    diag(S_init) <- 0
    S_vec <- as.vector(abs(S_init))
    lambda <- sort(S_vec, decreasing = TRUE)[2 * size + 1]
  }

  S_init[abs(S_init) <= lambda] <- 0
  non_zero_index <- which(as.matrix(S_init) != 0, arr.ind = TRUE)
  active.entry <- non_zero_index[which(non_zero_index[,1] < non_zero_index[,2]),]
  act_set <- active.entry - 1
  res <- soft_GT(s, lambda, act_set)
  Omega <- res[["Omega"]]
  is.acyclic <- res[["is.acyclic"]]

  return(list(Omega = Omega, active.entry = active.entry, is.acyclic = is.acyclic))
}
