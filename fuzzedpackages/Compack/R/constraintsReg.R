#' @title
#' Fit a linearly constrained linear regression model with group lasso regularization.
#'
#' @description
#' Fit a linearly constrained regression model with group lasso regularization.
#'
#' @usage
#' cglasso(y, Z, Zc = NULL, k, W = rep(1, times = p), intercept = TRUE,
#'         A =  kronecker(matrix(1, ncol = p), diag(k)), b = rep(0, times = k),
#'         u = 1, mu_ratio = 1.01,
#'         lam = NULL, nlam = 100,lambda.factor = ifelse(n < p1, 0.05, 0.001),
#'         dfmax = p, pfmax = min(dfmax * 1.5, p), tol = 1e-8,
#'         outer_maxiter = 1e+6, outer_eps = 1e-8,
#'         inner_maxiter = 1e+4, inner_eps = 1e-8)
#'
#' @param y respones vector with length \eqn{n}.
#'
#' @param Z design matrix of dimension \eqn{n \times p1}{n*p1}.
#'
#' @param k the group size in \eqn{Z}. The number of groups is \eqn{p = p1 / k }.
#'
#' @param Zc design matrix for unpenalized variables. Default value is NULL.
#'
#' @param A,b linear equalities of the form \eqn{A\beta_{p1} = b}, where \eqn{b} is a vector with length \eqn{k}, and \eqn{A} is
#'            a \eqn{k \times p1}{k*p1} matrix.
#'            Default values: \eqn{b} is a vector of 0's and
#'            \cr
#'            \code{A = kronecker(}\code{matrix(1, ncol = p), diag(k))}.
#'
#' @param W a vector in length p (the total number of groups), or a matrix with dimension \code{p1*p1}.
#'          Default value is rep(1, times = p).
#'          \itemize{
#'          \item a vector of penalization weights for the groups of coefficients. A zero weight implies no shrinkage.
#'          \item a diagonal matrix with positive diagonal elements.
#'          }
#'
#' @param lambda.factor the factor for getting the minimal lambda in \code{lam} sequence,
#'                      where \code{min(lam)} = \code{lambda.factor} * \code{max(lam)}.
#'                      \code{max(lam)} is the smallest value of \code{lam} for which all penalized group are 0's.
#'                      If \eqn{n >= p1}, the default is \code{0.001}. If \eqn{n < p1}, the default is \code{0.05}.
#'
#' @param nlam the length of the \code{lam} sequence. Default is 100. No effect if \code{lam} is
#'             provided.
#'
#' @param lam a user supplied lambda sequence.
#'            If \code{lam} is provided as a scaler and \code{nlam}\eqn{>1}, \code{lam} sequence is created starting from
#'            \code{lam}. To run a single value of \code{lam}, set \code{nlam}\eqn{=1}.
#'            The program will sort user-defined \code{lambda} sequence in decreasing order.
#'
#' @param dfmax limit the maximum number of groups in the model. Useful for handling very large \eqn{p}, if a partial path is desired.
#'              Default is \eqn{p}.
#'
#' @param pfmax limit the maximum number of groups ever to be nonzero. For example once a group enters the model along the path,
#'              no matter how many times it re-enters the model through the path, it will be counted only once.
#'              Default is \code{min(dfmax*1.5, p)}.
#
#' @param u the inital value of the penalty parameter of the augmented Lagrange method adopted in the outer loop. Default value is 1.
#'
#' @param mu_ratio the increasing ratio of the penalty parameter \code{u}. Default value is 1.01.
#'                 Inital values for scaled Lagrange multipliers are set as 0's.
#'                 If \code{mu_ratio} < 1, the program automatically set
#'                 the initial penalty parameter \code{u} as 0
#'                 and \code{outer_maxiter} as 1, indicating
#'                 that there is no linear constraint.
#'
#' @param tol tolerance for coefficient to be considered as non-zero.
#'            Once the convergence criterion is satisfied, for each element \eqn{\beta_j} in coefficient vector
#'            \eqn{\beta}, \eqn{\beta_j = 0} if \eqn{\beta_j < tol}.
#'
#             tolerance for coefficient vector for each group to be considered as none zero's. For example, considering
#             coefficient vector \eqn{\beta_j} for group j, if \eqn{max(abs(\beta_j))} < \code{tol}, set \eqn{\beta_j} as 0's.
#             Default value is 1e-8.
#'
#' @param outer_maxiter,outer_eps \code{outer_maxiter} is the maximun number of loops allowed for the augmented Lagrange method;
#'                                and \code{outer_eps} is the corresponding convergence tolerance.
#'
#' @param inner_maxiter,inner_eps \code{inner_maxiter} is the maximum number of loops allowed for blockwise-GMD;
#'                                and \code{inner_eps} is the corresponding convergence tolerance.
#'
#' @param intercept Boolean, specifying whether to include an intercept.
#'                  Default is TRUE.
#'
#' @return A list of
#' \item{beta}{a matrix of coefficients.}
#' \item{lam}{the sequence of lambda values.}
#' \item{df}{a vector, the number of nonzero groups in estimated coefficients for \code{Z} at each value of lambda.}
#' \item{npass}{total number of iteration.}
#' \item{error}{a vector of error flag.}
#'
## usethis namespace: start
#' @useDynLib Compack, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
#'
#' @export


cglasso <- function(y, Z, Zc = NULL, k,
                    W = rep(1, times = p),
                    intercept = TRUE,
                    A =  kronecker(matrix(1, ncol = p), diag(k)), b = rep(0, times = k),
                    u = 1, mu_ratio = 1.01,
                    lam = NULL, nlam = 100,lambda.factor = ifelse(n < p1, 0.05, 0.001),
                    dfmax = p, pfmax = min(dfmax * 1.5, p),
                    tol = 1e-8,
                    outer_maxiter = 1e+6, outer_eps = 1e-8,
                    inner_maxiter = 1e+4, inner_eps = 1e-8) {

  y <- drop(y)
  Z <- as.matrix(Z)
  n <- length(y)
  p1 <- dim(Z)[2]
  p <- p1 / k
  group.index <- matrix(1:p1, nrow = k)



  inter <- as.integer(intercept)
  Zc <- as.matrix(cbind(Zc, rep(inter, n)))
  if(inter == 1) {
    Zc_proj <- solve(crossprod(Zc)) %*% t(Zc)
  } else {
    if(dim(Zc)[2] == 1) Zc_proj <- t(Zc) else Zc_proj <- rbind(solve(crossprod(Zc[, -ncol(Zc)])) %*% t(Zc[, -ncol(Zc)]), rep(inter, n))
  }
  beta_c <- as.vector(Zc_proj %*% y)
  beta.ini <- c(rep(0, times = p1), beta_c)

  if( is.vector(W) ){
    if(length(W) != p) stop("W should be a vector of length p")
    W_inver <- diag(p1)
    pf <- W
  } else if( is.matrix(W) ) {
    if(any(dim(W) - c(p1, p1) != 0)) stop('Wrong dimension of Weights matrix')
    pf <- rep(1, times = p)
    W_inver <- W
  }

  Z <- Z %*% W_inver
  A <- A %*% W_inver



  # if(is.null(lam)) {
  #   lam0 <- lam.ini(Z = Z, y = y - Zc %*% beta_c, ix = group.index[1, ], iy = group.index[k ,], pf = pf)
  #   lam <- exp(seq(from=log(lam0), to=log(lambda.factor * lam0),length=nlam))
  # } else if(length(lam) == 1) {
  #   lam <- exp(seq(from=log(lam), to=log(lambda.factor * lam),length=nlam))
  # }

  if(is.null(lam)) {
    lam0 <- lam.ini(Z = Z, y = y - Zc %*% beta_c, ix = group.index[1, ], iy = group.index[k ,], pf = pf)
    lam <- exp(seq(from=log(lam0), to=log(lambda.factor * lam0),length=nlam))
  } else {
    if(length(lam) == 1 && nlam > 1) {
      lam <- exp(seq(from=log(lam), to=log(lambda.factor * lam),length=nlam))
    } else {
      lam <- sort(lam, decreasing = TRUE)
    }
  }

  pfmax <- as.integer(pfmax)
  dfmax <- as.integer(dfmax)
  inner_maxiter <- as.integer(inner_maxiter)
  outer_maxiter <- as.integer(outer_maxiter)
  inner_eps <- as.double(inner_eps)
  outer_eps <- as.double(outer_eps)
  tol <- as.double(tol)
  u <- as.double(u)
  mu_ratio <- as.double(mu_ratio)
  #estar <- as.double(estar)

  output <- ALM_GMD(y = y, Z = Z, Zc = Zc, Zc_proj = Zc_proj, beta = beta.ini, lambda = lam,
                    pf = pf, dfmax = dfmax, pfmax = pfmax, A = A, b = b,
                    group_index = group.index, u_ini = u, mu_ratio = mu_ratio,
                    inner_eps = inner_eps, outer_eps = outer_eps,
                    inner_maxiter = inner_maxiter, outer_maxiter = outer_maxiter, tol = tol
                    #, estar = estar
  )

  output$beta[1:p1, ] <-  W_inver %*% output$beta[1:p1, ]
  return(output)

}



#'
#' @title
#' Fit a linearly constrained linear regression model with lasso regularization.
#'
#' @description
#' Fit a linearly constrained linear model with lasso regularization.
#'
#' @usage
#' classo(y, Z, Zc = NULL, intercept = TRUE, pf = rep(1, times = p),
#'        lam = NULL, nlam = 100,lambda.factor = ifelse(n < p, 0.05, 0.001),
#'        dfmax = p, pfmax = min(dfmax * 1.5, p),
#'        u = 1, mu_ratio = 1.01, tol = 1e-10,
#'        outer_maxiter = 3e+08, outer_eps = 1e-8,
#'        inner_maxiter = 1e+6, inner_eps = 1e-8,
#'        A = rep(1, times = p), b = 0, beta.ini)
#'
#' @param y a response vector with length n.
#'
#' @param Z a design matrix, with dimension \eqn{n \times p}{n*p}.
#'
#' @param Zc design matrix for unpenalized variables. Default value is NULL.
#'
#' @param pf penalty factor, a vector of length p. Zero implies no shrinkage. Default value for each entry is 1.
#'
#' @param A,b linear equalities of the form \eqn{A\beta_p = b}, where \eqn{b} is a scaler,
#'            and \eqn{A} is a row-vector of length \code{p}. Default values: \eqn{b} is 0 and \code{A = matrix(1, ncol = p)}.
#'
#' @param lambda.factor the factor for getting the minimal lambda in the \code{lam} sequence,
#'                      where \code{min(lam)} = \code{lambda.factor} * \code{max(lam)}.
#'                      \code{max(lam)} is the smallest value of \code{lam} for which all penalized coefficients become zero.
#'                      If \eqn{n >= p}, the default is \code{0.001}. If \eqn{n < p}, the default is \code{0.05}.
#'
#' @param mu_ratio the increasing ratio, with value at least 1, for \code{u}. Default value is 1.01.
#'                 Inital values for scaled Lagrange multipliers are set as 0.
#'                 If \code{mu_ratio} < 1, the program automatically set \code{u} as 0 and \code{outer_maxiter} as 1, indicating
#'                 that there is no linear constraint.
#'
#' @param tol tolerance for the estimated coefficients to be considered as non-zero, i.e., if \eqn{abs(\beta_j)} < \code{tol}, set \eqn{\beta_j} as 0.
#'            Default value is 1e-10.
#'
#' @param outer_maxiter,outer_eps \code{outer_maxiter} is the maximum number of loops allowed for the augmented Lanrange method;
#'                                and \code{outer_eps} is the corresponding convergence tolerance.
#'
#' @param inner_maxiter,inner_eps \code{inner_maxiter} is the maximum number of loops allowed for the coordinate descent;
#'                                and \code{inner_eps} is the corresponding convergence tolerance.
#'
#' @param beta.ini inital value of the coefficients. Can be unspecified.
#'
#'
#' @inheritParams cglasso
#'
#'
#' @return A list of
#' \item{beta}{a matrix of coefficients.}
#' \item{lam}{the sequence of lambda values.}
#' \item{df}{a vector, the number of nonzero coefficients for \code{Z} at each value of lambda.}
#' \item{npass}{total number of iteration.}
#' \item{error}{a vector of error flag.}
#'
#' @export


classo <- function(y, Z, Zc = NULL, intercept = TRUE,
                   pf = rep(1, times = p),
                   lam = NULL, nlam = 100,lambda.factor = ifelse(n < p, 0.05, 0.001),
                   dfmax = p, pfmax = min(dfmax * 1.5, p),
                   u = 1, mu_ratio = 1.01, tol = 1e-10,
                   outer_maxiter = 3e+08, outer_eps = 1e-8,
                   inner_maxiter = 1e+6, inner_eps = 1e-8,
                   A = rep(1, times = p), b = 0,
                   beta.ini) {

  this.call <- match.call()
  y <- drop(y)
  Z <- as.matrix(Z)
  n <- length(y)
  p <- dim(Z)[2]
  if(length(A) != p) stop("Length of vector A in Ax = b is wrong")
  inter <- as.integer(intercept)
  Znames <- colnames(Z)
  if (is.null(Znames)) Znames <- paste0("Z", seq(p))

  Zc <- as.matrix(cbind(Zc, rep(inter, n)))
  if(inter == 1) {
    Zc_proj <- solve(crossprod(Zc)) %*% t(Zc)
  } else {
    if(dim(Zc)[2] == 1) Zc_proj <- t(Zc) else Zc_proj <- rbind(solve(crossprod(Zc[, -ncol(Zc)])) %*% t(Zc[, -ncol(Zc)]), rep(inter, n))
  }

  if(missing(beta.ini)) {
    beta_c <- as.vector(Zc_proj %*% y)
    beta.ini <- c(rep(0, times = p), beta_c)
  }

  m <- dim(Zc)[2]
  if(m > 1) {
    Zcnames <- colnames(Z)
    if (is.null(Zcnames)) Zcnames <- paste0("Zc", seq(m-1))
  } else {
    Zcnames <- NULL
  }

  if(is.null(lam)) {
    lam0 <- t(Z) %*% (y - Zc %*% beta_c) / n
    lam0 <- max(abs(lam0) / pf)
    lam <- exp(seq(from=log(lam0), to=log(lambda.factor * lam0),length=nlam))
  }

  pfmax <- as.integer(pfmax)
  dfmax <- as.integer(dfmax)
  inner_maxiter <- as.integer(inner_maxiter)
  outer_maxiter <- as.integer(outer_maxiter)
  inner_eps <- as.double(inner_eps)
  outer_eps <- as.double(outer_eps)
  tol <- as.double(tol)
  u <- as.double(u)
  mu_ratio <- as.double(mu_ratio)
  b <- as.double(b)
  #estar <- as.double(estar)

  output <- ALM_CD(y = y, Z = Z, Zc = Zc, Zc_proj = Zc_proj, beta = beta.ini,
                   lambda = lam, pf = pf, dfmax = dfmax, pfmax = pfmax,
                   inner_eps = inner_eps, outer_eps = outer_eps,
                   inner_maxiter = inner_maxiter, outer_maxiter = outer_maxiter,
                   u_ini = u, mu_ratio = mu_ratio, tol = tol,
                   A = A, b = b

  )

  output$call <- this.call
  class(output) <- "compCL"
  output$dim <- dim(output$beta)
  output$lam <- drop(output$lam)
  dimnames(output$beta) <- list(c(Znames, Zcnames, "Intercept"), paste0("L", seq(along = output$lam)) )
  return(output)

}


