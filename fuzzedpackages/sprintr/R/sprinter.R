#' Sparse Reluctant Interaction Modeling
#'
#' This is the main function that fits interaction models with a path of tuning parameters (for Step 3).
#'
#' @param x An \code{n} by \code{p} design matrix of main effects. Each row is an observation of \code{p} main effects.
#' @param y A response vector of size \code{n}.
#' @param num_keep Number of candidate interactions to keep in Step 2. If \code{num_keep} is not specified (as default), it will be set to \code{[n / log n]}.
#' @param square Indicator of whether squared effects should be fitted in Step 1. Default to be FALSE.
#' @param lambda A user specified list of tuning parameter. Default to be NULL, and the program will compute its own \code{lambda} path based on \code{nlam} and \code{lam_min_ratio}.
#' @param nlam The number of \code{lambda} values. Default value is \code{100}.
#' @param lam_min_ratio The ratio of the smallest and the largest values in \code{lambda}. The largest value in \code{lambda} is usually the smallest value for which all coefficients are set to zero. Default to be \code{1e-2} in the \code{n} < \code{p} setting.
#' @return An object of S3 class "\code{sprinter}".
#'  \describe{
#'   \item{\code{n}}{The sample size.}
#'   \item{\code{p}}{The number of main effects.}
#'   \item{\code{a0}}{Estimate of intercept.}
#'   \item{\code{coef}}{Estimate of regression coefficients.}
#'   \item{\code{idx}}{Indices of all main effects and interactions in Step 3.}
#'   \item{\code{fitted}}{Fitted response value. It is a \code{n}-by-\code{nlam} matrix, with each column representing a fitted response vector for a value of lambda.}
#'   \item{\code{lambda}}{The sequence of \code{lambda} values used.}
#'   \item{\code{call}}{Function call.}
#'  }
#' @seealso
#'   \code{\link{cv.sprinter}}
#' @examples
#' set.seed(123)
#' n <- 100
#' p <- 200
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x[, 1] - 2 * x[, 2] + 3 * x[, 1] * x[, 3] - 4 * x[, 4] * x[, 5] + rnorm(n)
#' mod <- sprinter(x = x, y = y)
#'
#' @import glmnet
#' @export
sprinter <- function(x, y, num_keep = NULL, square = FALSE, lambda = NULL, nlam = 100, lam_min_ratio = ifelse(nrow(x) < ncol(x), 0.01, 1e-04)){
  n <- nrow(x)
  p <- ncol(x)
  # by default, num_keep is set to be [n / log(n)]
  if(is.null(num_keep)){
    num_keep = ceiling(n / log(n))
  }
  else{
    num_keep = min(ceiling(num_keep), ifelse(square, p*(p - 1) / 2, p * (p - 1)/ 2 + p))
  }

  # we always standardize the design matrix to get main effects
  # squared effects and interactions are built upon standardized main effects
  x <- myscale(x)
  # xm is the (standardized) design matrix of main effects
  col_mean <- attr(x = x, which = "scaled:center")
  col_sd <- attr(x = x, which = "scaled:scale")
  xm <- x
  mean_y <- mean(y)

  # The First Step
  # run CV-lasso on
  # (1) main effects M (square == FALSE)
  # (2) main effects M + squared main effects M^2 (square == TRUE)
  # return the fitted value of response
  if(square){
    x_sq <- myscale(x^2)

    col_mean <- c(col_mean, attr(x_sq, which = "scaled:center"))
    col_sd <- c(col_sd, attr(x_sq, which = "scaled:scale"))
    x <- cbind(x, x_sq)
  }

  fit <- glmnet::cv.glmnet(x = x, y = y - mean_y,
                           lambda = get_lambda(x = x, y = y - mean_y),
                           intercept = FALSE,
                           standardize = FALSE)
  # residual
  fitted_first <- mean_y + x %*% as.numeric(fit$glmnet.fit$beta[, which.min(fit$cvm)])
  r <- y - fitted_first

  #cat("step 1 finished!", fill = TRUE)

  # The Second Step
  # find num_keep higher order terms from
  # (1) squared main effects M^2 + Interaction effects I
  #     (square == FALSE)
  # (2) Interaction effects I (square == TRUE)
  # with largest absolute correlation with the residuals r from first step
  # return the selected variables set B
  idx <- screen_cpp(x = xm, y = r, num_keep = num_keep, square = square)
  #cat("step 2 finished!", fill = TRUE)

  # construct design matrix of selected interactions
  design <- myscale(xm[, idx[, 1]] * xm[, idx[, 2]])

  col_mean <- c(col_mean,
                attr(design, which = "scaled:center"))
  col_sd <- c(col_sd,
              attr(design, which = "scaled:scale"))
  # the total design matrix
  design <- cbind(x, design)
  if(square){
    idx <- rbind(cbind(rep(0, p), seq(p)), cbind(seq(p), seq(p)), idx)
  }
  else{
    idx <- rbind(cbind(rep(0, p), seq(p)), idx)
  }
  idx <- as.matrix(idx)
  colnames(idx) <- c("index_1", "index_2")
  # idx has two columns
  # which are the j,k indices of nonzero elements
  # main effect index is of form (0, k)
  # squared effect index is of form (k, k)
  # interaction effect index is of form (j, k) for j < k

  if(is.null(lambda)){
    lambda <- get_lambda(x = design, y = y - mean_y)
  }
  nlam <- length(lambda)

  # The Third Step:
  #     run lasso of response y on A + B
  #     corresponding to the best lambda
  fit <- glmnet::glmnet(x = design, y = y - mean_y,
                        lambda = lambda,
                        intercept = FALSE,
                        standardize = FALSE)
  coef <- fit$beta
  # drop the names of the matrix object returned by glmnet
  fitted <- mean_y + design %*% coef
  # fitted_main <- mean_y + design[, head(seq(ncol(design)), -num_keep)] %*% coef[head(seq(ncol(design)), -num_keep), ]
  # scale estimates back to the original scale of x
  coef <- coef / col_sd
  a0 <- as.numeric(mean_y - crossprod(col_mean, as.matrix(coef)))
  #cat("step 3 finished!", fill = TRUE)

  out <- list(n = n,
              p = p,
              a0 = a0,
              coef = coef,
              idx = idx,
              #fitted_first = fitted_first,
              #fitted_main = fitted_main,
              fitted = fitted,
              lambda = lambda,
              call = match.call())
  class(out) <- "sprinter"
  return(out)
}
