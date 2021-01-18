#' Solve main optimization problem along a path of lambda
#'
#' Compute the varband estimates along a path of tuning parameter values.
#'
#' @param S The sample covariance matrix
#' @param w Logical. Should we use weighted version of the penalty or not? If \code{TRUE}, we use general weight. If \code{FALSE}, use unweighted penalty. Default is \code{FALSE}.
#' @param lasso Logical. Should we use l1 penalty instead of hierarchical group lasso penalty? Note that by using l1 penalty, we lose the banded structure in the resulting estimate. And when using l1 penalty, the becomes CSCS (Convex Sparse Cholesky Selection) introduced in Khare et al. (2016). Default value for \code{lasso} is \code{FALSE}.
#' @param lamlist A list of non-negative tuning parameters \code{lambda}.
#' @param nlam If lamlist is not provided, create a lamlist with length \code{node}. Default is 60.
#' @param flmin if lamlist is not provided, create a lamlist with ratio of the smallest and largest lambda in the list. Default is 0.01.
#' @return A list object containing \describe{
#' \item{path: }{A array of dim (\code{p}, \code{p}, \code{nlam}) of estimates of L}
#' \item{lamlist: }{a grid values of tuning parameters}
#' }
#' @examples
#' set.seed(123)
#' n <- 50
#' true <- varband_gen(p = 50, block = 5)
#' x <- sample_gen(L = true, n = n)
#' S <- crossprod(scale(x, center = TRUE, scale = FALSE))/n
#' path_res <- varband_path(S = S, w = FALSE, nlam = 40, flmin = 0.03)
#' @export
#'
#' @seealso \code{\link{varband}} \code{\link{varband_cv}}
varband_path <- function(S, w = FALSE, lasso = FALSE, lamlist = NULL, nlam = 60, flmin = 0.01){
  p <- ncol(S)
  stopifnot(p == nrow(S))

  if (is.null(lamlist)) {
    lam_max <- lammax(S = S)
    lamlist <- pathGen(nlam = nlam, lam_max = lam_max,
                       flmin = flmin, S = S)
  } else {
    nlam <- length(lamlist)
  }
  result<- array(NA, c(p, p, nlam))

  for (i in seq(nlam)) {
    if(i==1){
      #      cat(i)
      result[, , i] <- diag(1/sqrt(diag(S)))
    }
    else
    {
      #      cat(i)
      result[, , i] <- varband(S = S, lambda = lamlist[i],
                               init = result[, , i-1], w = w, lasso = lasso)
    }
  }
  return(list(path = result, lamlist = lamlist))
}

lammax <- function(S){
  #### This function calculates the max value in the tuning parameter list
  # such that the estimator L_{\lambda} is a diagonal matrix
  # NOTE: this is not necessarily true, but generally
  # a upper bound of the value we are looking for.

  # Args:
  #     S: the p-by-p sample covariance matrix

  p <- ncol(S)
  sighat <- rep(NA, p-1)
  for (r in seq(2, p)){
    sighat[r-1] <- max(abs(S[(1:(r-1)), r]))/sqrt(S[r, r])
  }
  2 * max(sighat)
}

pathGen <- function(nlam, lam_max, flmin, S){
  # Generate a path of lambda, with
  # nlam/2 decreasing exponentially
  # nlam/2 decreasing linearly
  # lam_max <- lammax(S)
  lamlist_lin <- lam_max * exp(seq(0, log(flmin), length = nlam/2))
  lamlist_exp <- seq(lam_max - 1e-8, lam_max*flmin - 1e-8, length.out = nlam/2)
  return(sort(unique(c(lamlist_lin, lamlist_exp)), decreasing = T))
}
