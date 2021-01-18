#### This file contains all utility functions
#### Not necessarily useful in the package

#' Plot the sparsity pattern of a square matrix
#'
#' Black, white and gray stand for positive, zero and negative respectively
#'
#' @param Mat A matrix to plot.
#' @param main A plot title.
#'
#' @examples
#' set.seed(123)
#' p <- 50
#' n <- 50
#' phi <- 0.4
#' true <- varband_gen(p = p, block = 5)
#' matimage(true)
#' @export
#' @import graphics
matimage <- function(Mat, main = NULL){
  tmppar <- par(pty = "s")
  image(sign(t(apply(Mat, 2, rev))), axes = FALSE, col = c("gray50","white","black"), main = main)
  par(tmppar)
}

KL <- function(Omega_est, Omega_true){
  # Kullback-Leibler divergence for gaussian likelihood
  # see (14) in Nested Lasso (Levina et,al 2008)
  res <- determinant(x = Omega_est, logarithm = TRUE)$modulus
  res <- -res + determinant(x = Omega_true, logarithm = TRUE)$modulus
  res <- sum(diag(solve(Omega_true, Omega_est))) + res - nrow(Omega_est)
  return(as.numeric(res))
}

Irrep <- function(L){
  # Calculate the irrepresentability condition of true L
  p <- ncol(L)
  Sigma <- solve(crossprod(L))
  res <- c()
  for (r in seq(2, p)){
    # For each row
    Jr <- 0
    for (j in seq(r-1)){
      if(L[r, j] != 0){
        Jr <- j-1
        break
      }
    }
    # Got Ir
    if(Jr > 0){
      Ir <- ((Jr+1):(r-1))
      tmp <- rep(NA, Jr)
      for (ll in seq(Jr)){
        tmp[ll] <- norm(as.matrix(Sigma[ll, Ir]%*%solve(Sigma[Ir, Ir])), "1")
      }
      cat("r = ", r, "  ", "Jr = ", length(tmp), "  ", "Kr = ", length(Ir), fill = TRUE)
      res <- c(res, max(tmp))
    }
  }
  return(max(res))
}

#################################################################
####  The following functions are the first R version of main
####  functions used in the package
####  They have been replaced by their c++ version
#################################################################
# #' Evaluate the proximal operator of the hierarchical group lasso with general weights
# #'
# #' This function solves (7) in the paper for general weight w
# #' by solving its dual by performing Newton's method on at most
# #' r-1 univariate functions.
# #' See Algorithm 1 and Theorem 1 in the online supplemenatry.
# #'
# #' @param y An r-dimensional vector.
# #' @param tau lambda/rho
elliproj_w_R <- function(y, tau){
  # This function performs the ellipsoid projection
  # See supplementary material
  r <- length(y)

  nu <- rep(NA, r-1)
  # pp is the z vector in the paper
  pp <- y
  for(l in seq(r-1)){
    # ww[m] = w_{lm}
    ww <- 1 / (seq(l,1)^2)
    # check if it lies in the ellipsoid
    if (sum((pp[1:l]/ww)^2) <= tau^2){
      nu[l] <- 0
      pp[1:l] <- 0
    }
    else{
      # project onto the elliposid
      f <- function(nu)
        1-tau/sqrt(sum((pp[1:l]/(ww+nu/ww))^2))
      # lower and upper bound of root
      nu.u <- sqrt(sum((ww*pp[1:l])^2))/tau
      nu.l <- max(nu.u-ww[l]^2, 0)
      # find root
      if(abs(f(nu.u)) < SMALL)
        nu[l] <- nu.u
      else if (abs(f(nu.l)) < SMALL)
        nu[l] <- nu.l
      else
        nu[l] <- uniroot(f = f, interval = c(nu.l, nu.u), tol = 1e-15)$root
      pp[1:l] <- pp[1:l] * nu[l] / ( ww^2 + nu[l] )
    }
  }
  return(pp)
}

# #' Evaluate the proximal operator of the hierarchical group lasso with simple weights
# #'
# #' This function solves (7) in the paper for unweighted version(w = 1)
# #' by solving its dual by performing Newton's method on at most
# #' r-1 univariate functions.
# #' See Algorithm 2 in the paper
# #'
# #' @param y An r-dimensional vector.
# #' @param tau lambda/rho
elliproj_u_R <- function(y, tau){
  # This function performs the ellipsoid projection
  # of the unweighted estimator, which is very easy
  # See algorithm 2 in the paper
  r <- length(y)
  # pp is the z vector in the paper
  pp <- y
  for(l in seq(r-1)){
    tmpnorm <- sqrt(sum(pp[1:l]^2))
    if(tmpnorm <= tau)
      pp[1:l] <- 0
    else
      pp[1:l] <- (1-tau/tmpnorm)*pp[1:l]
  }
  return(pp)
}

SMALL <- 1e-15
# #' Compute the varband estimate for a fixed tuning parameter value.
# #'
# #' Solves the main optimization problem in Yu & Bien (2016):
# #' \deqn{min_L -2 \sum_{r=1}^p L_{rr} + tr(SLL^T) + lam * \sum_{r=2}^p P_r(L_{r.})}{min_L -2 sum_{r=1}^p L_{rr} + tr(SLL^T) + lam * sum_{r=2}^p P_r(L_{r.})}
# #' where \deqn{P_r(L_{r.}) = \sum_{\ell = 2}^{r-1} \left(\sum_{m=1}^\ell w_{\ell m}^2 L_{rm}^2\right)^{1/2}}{P_r(L_r.) = sum_{l=2}^{r-1} (sum_m=1^l w^2_lm L^2_rm)^{1/2}}
# #'
# #' The function decomposes into p independent row problems,
# #' each of which is solved by an ADMM algorithm.
# #' see paper for more explanation.
# #' @param S The sample covariance matrix
# #' @param lambda Non-negative tuning parameter. Controls sparsity level.
# #' @param w Logical. Should we use weighted version of the penalty or not? If \code{TRUE}, we use general weight. If \code{FALSE}, use unweighted penalty. Default is \code{FALSE}.
# #' @param init Initial estimate of L. Default is a closed-form diagonal estimate of L.
# #' @return Returns the variable banding estimate of L, where L^TL = Omega.
# #'
# #' @seealso \code{\link{varband_path}} \code{\link{varband_cv}}
# #'
varband_R <- function(S, lambda, init = NULL, w = FALSE){
  p <- ncol(S)
  # check that S is square
  stopifnot(p == nrow(S))
  if (is.null(init))
    init <- diag(1/sqrt(diag(S)))

  L <- matrix(0,p,p)
  L[1,1] <- 1/(sqrt(S[1,1]))

  # for the second to the p-th row
  for(r in seq(2,p)){
    L[r,] <- c(rowadmm_R(S = S[1:r,1:r],
                         init_row = init[r,1:r],
                         lambda = lambda, w = w),
               rep(0,p-r))
  }
  return(L)
}

# #' Compute one row of varband estimate for a fixed tuning parameter
# #'
# #'This function solve the following r-th row estimation problem \deqn{min_{beta_r>0} -2 log beta_r + 1/n ||X beta||^2 + lambda P(beta)}
# #' using an ADMM algorithm with changing rho.
# #'
# #' See algorithm 1 in the paper.
# #'
# #' @param S An r-by-r submatrix of sample covariance matrix.
# #' @param init_row The initial estimate of the row.
# #' @param lambda Non-negative tuning parameter. Controls sparsity level.
# #' @param w Logical. Should we use weighted version of the penalty or not? If \code{TRUE}, we use general weight. If \code{FALSE}, use unweighted penalty. Default is \code{FALSE}.
# #' @param tol Tolerance for convergence.
# #' @param itermax Maximum number of iterations of ADMM to perform.
rowadmm_R <- function(S, init_row, lambda, w = FALSE, tol = 1e-4, itermax = 1e5){
  # This function solve the following row estimation problem
  # \min_{\beta_r>0} -2 log \beta_r + 1/n ||X\beta||^2
  # + \lambda P(\beta)
  # using an ADMM algorithm with changing rho

  r <- ncol(S)
  stopifnot(r == nrow(S))
  # could use a lower tolerance for unweighted version
  if (!w)
    tol <- 1e-8

  # Default parameter in ADMM
  tolabs <- tol
  tolrel <- tol
  # Changing rho
  rho <- 2
  mu <- 10
  inc <- 2
  dec <- 2

  # Initialize the result
  beta <- init_row
  gamma <- init_row

  # dual variable
  u <- rep(0,r)

  S_inv <- inverse_update_R(S = as.matrix(S[-r,-r]),
                            r = r, rho = rho)
  for (i in seq(itermax)) {
    # Primal&Dual Updates
    beta_new <- close_update_R(S = S, S_inv = S_inv,
                               r = r, rho = rho, u = u,
                               gamma = gamma)
    if (w)
      gamma_new <- elliproj_w_R(y = beta_new + u/rho,
                                tau = lambda/rho)
    else
      gamma_new <- elliproj_u_R(y = beta_new + u/rho,
                                tau = lambda/rho)

    u <- u + rho*(beta_new - gamma_new)

    # check convergence See pp 22 Boyd(2011)
    # primal residual
    pres <- sqrt(crossprod(beta_new - gamma_new))
    # dual residual
    dres <- rho*sqrt(crossprod(gamma_new - gamma))
    # primal tolerance
    peps <- tolabs*sqrt(r) +
      tolrel*max(sqrt(crossprod(beta_new)),
                 sqrt(crossprod(gamma_new)))
    # dual tolerance
    deps <- tolabs*sqrt(r) + tolrel*sqrt(crossprod(u))

    if(pres <= peps & dres <= deps)
    {
      #cat("ADMM converges after",i,"iterations",fill=TRUE)
      break
    }
    else{
      # if not, update estimates and rho
      beta <- beta_new
      gamma <- gamma_new

      # Update rho if needed and corresponding S_inv
      if(pres > mu*dres){
        rho <- rho*inc
        S_inv <- inverse_update_R(S = as.matrix(S[-r,-r]),
                                  r = r, rho = rho)
      }
      else if(dres > mu*pres){
        rho <- rho/dec
        S_inv <- inverse_update_R(S = as.matrix(S[-r,-r]),
                                  r = r, rho = rho)
      }
    }
  }
  if(i==itermax)
    cat("ADMM fails to converge",fill=TRUE)
  return (gamma_new)
}

inverse_update_R <- function(S, r, rho){
  if(r==1){
    S_inv <- as.matrix(1/(2*S[1,1] + rho))
    return(S_inv)
  }
  else{
    S_inv <- 2*S
    diag(S_inv) <- diag(S_inv) + rho
    S_inv <- solve(S_inv)
    return(S_inv)
  }
}

# #' Close-form update of beta in Algorithm 1
# #'
# #' This function solves (6) in the paper with a closed form solution
# #'
# #' @param S An r-by-r submatrix of sample covariance matrix.
# #' @param S_inv inverse of (2S_{-r,-r} + rho I)
# #' @param r row index
# #' @param rho parameter rho used in ADMM
# #' @param u dual variable in ADMM
# #' @param gamma priaml variable in ADMM
close_update_R <- function(S, S_inv, r, rho, u, gamma){
  # This performs the closed updated for beta
  # see Section 3 in the paper
  vec.tmp <- S_inv%*%S[-r,r]
  A <- 4*crossprod(vec.tmp,S[-r,r]) - 2*S[r,r] - rho
  B <- 2*crossprod(vec.tmp, (u[-r] - rho*gamma[-r])) - u[r] + rho*gamma[r]
  res <- rep(0,r)
  # beta_r
  res[r] <- (-sqrt(B^2 - 8 * A) - B)/(2 * A)
  # beta_[-r]
  res[-r] <- -2*res[r]*vec.tmp - S_inv%*%(u[-r] - rho*gamma[-r])
  return(res)
}

