
#' @title Split Regularized Regression algorithm with a sparsity and diversity penalty.
#' @param x Design matrix.
#' @param y Response vector.
#' @param num_lambdas_sparsity Length of the grid of sparsity penalties.
#' @param num_lambdas_diversity Length of the grid of diversity penalties.
#' @param alpha Elastic Net tuning constant: the value must be between 0 and 1. Default is 1 (Lasso).
#' @param num_models Number of models to build.
#' @param tolerance Tolerance parameter to stop the iterations while cycling over the models.
#' @param max_iter Maximum number of iterations before stopping the iterations while cycling over the models.
#' @param num_folds Number of folds for cross-validating.
#' @param num_threads Number of threads used for parallel computation over the folds.

#' @return An object of class cv.SplitReg, a list with entries
#' \item{betas}{Coefficients computed over the path of penalties for sparsity; the penalty for diversity is fixed at the optimal value.}
#' \item{intercepts}{Intercepts for each of the models along the path of penalties for sparsity.}
#' \item{index_opt}{Index of the optimal penalty parameter for sparsity.}
#' \item{lambda_sparsity_opt}{Optimal penalty parameter for sparsity.}
#' \item{lambda_diversity_opt}{Optimal penalty parameter for diversity.}
#' \item{lambdas_sparsity}{Grid of sparsity parameters.}
#' \item{lambdas_diversity}{Grid of diversity parameters.}
#' \item{cv_mse_opt}{Optimal CV MSE.}
#' \item{call}{The matched call.}

#' 
#' @description
#' Computes a split regularized regression estimator. The sparsity and diversity penalty
#' parameters are chosen automatically.
#' 
#' @details
#' Computes a split regularized regression estimator with \code{num_models} (\eqn{G}) models, defined as the linear models
#' \eqn{\boldsymbol{\beta}^{1},\dots, \boldsymbol{\beta}^{G}} that minimize
#' \deqn{\sum\limits_{g=1}^{G}\left( \frac{1}{2n}\Vert \mathbf{y}-\mathbf{X} \boldsymbol{\beta}^{g}\Vert^{2} 
#' +\lambda_{S}\left( \frac{(1-\alpha)}{2}\Vert \boldsymbol{\beta}^{g}\Vert_{2}^{2}+\alpha \Vert \boldsymbol{
#' \beta \Vert_1}\right)+\frac{\lambda_{D}}{2}\sum\limits_{h\neq g}\sum_{j=1}^{p}\vert \beta_{j}^{h}\beta_{j}^{g}\vert \right),}
#' over grids for the penalty parameters \eqn{\lambda_{S}} and \eqn{\lambda_{D}} that are built automatically.
#' Larger values of \eqn{\lambda_{S}} encourage more sparsity within the models and larger values of \eqn{\lambda_{D}} encourage more diversity
#' among them. 
#' If \eqn{\lambda_{D}=0}, then all of the models are equal to the Elastic Net regularized
#' least squares estimator with penalty parameter \eqn{\lambda_{S}}. Optimal penalty parameters are found by
#' \code{num_folds} cross-validation, where the prediction of the ensemble is formed by simple averaging.
#' The predictors and the response are standardized to zero mean and unit variance before any computations are performed.
#' The final output is in the original scales.
#' 
#' @seealso \code{\link{predict.cv.SplitReg}}, \code{\link{coef.cv.SplitReg}}
#' 
#' @examples 
#' library(MASS)
#' set.seed(1)
#' beta <- c(rep(5, 5), rep(0, 45))
#' Sigma <- matrix(0.5, 50, 50)
#' diag(Sigma) <- 1
#' x <- mvrnorm(50, mu = rep(0, 50), Sigma = Sigma)
#' y <- x %*% beta + rnorm(50)
#' fit <- cv.SplitReg(x, y, num_models=2)
#' coefs <- predict(fit, type="coefficients")
#' 

cv.SplitReg <- function(x, y, num_lambdas_sparsity = 100, num_lambdas_diversity = 100, alpha = 1, num_models = 10,
                       tolerance = 1e-8, max_iter = 1e5, num_folds = 10, num_threads = 1){
  # Some sanity checks on the input
  if (all(!inherits(x, "matrix"), !inherits(x, "data.frame"))) {
    stop("x should belong to one of the following classes: matrix, data.frame")
  } else if (all(!inherits(y, "matrix"), all(!inherits(y, "numeric")))) {
    stop("y should belong to one of the following classes: matrix, numeric")
  } else if (any(anyNA(x), any(is.nan(x)), any(is.infinite(x)))) {
    stop("x should not have missing, infinite or nan values")
  } else if (any(anyNA(y), any(is.nan(y)), any(is.infinite(y)))) {
    stop("y should not have missing, infinite or nan values")
  } else {
    if(inherits(y, "matrix")) {
      if (ncol(y)>1){
        stop("y should be a vector")
      }
      # Force to vector if input was a matrix
      y <- as.numeric(y)
    }
    len_y <- length(y)
    if (len_y != nrow(x)) {
      stop("y and x should have the same number of rows")
    }
  }
  if (!inherits(tolerance, "numeric")) {
    stop("tolerance should be numeric")
  } else if (!all(tolerance < 1, tolerance > 0)) {
    stop("tolerance should be between 0 and 1")
  }
  if (!inherits(alpha, "numeric")) {
    stop("alpha should be numeric")
  } else if (!all(alpha <= 1, alpha > 0)) {
    stop("alpha should be between 0 and 1")
  }
  if (!inherits(max_iter, "numeric")) {
    stop("max_iter should be numeric")
  } else if (any(!max_iter == floor(max_iter), max_iter <= 0)) {
    stop("max_iter should be a positive integer")
  }
  if (!inherits(num_models, "numeric")) {
    stop("num_models should be numeric")
  } else if (any(!num_models == floor(num_models), num_models <= 1)) {
    stop("num_models should be an integer, greater than one")
  }
  if (!inherits(num_lambdas_sparsity, "numeric")) {
    stop("num_lambdas_sparsity should be numeric")
  } else if (any(!num_lambdas_sparsity == floor(num_lambdas_sparsity), num_lambdas_sparsity <= 0)) {
    stop("num_lambdas_sparsity should be a positive integer")
  }
  if (!inherits(num_lambdas_diversity, "numeric")) {
    stop("num_lambdas_diversity should be numeric")
  } else if (any(!num_lambdas_diversity == floor(num_lambdas_diversity), num_lambdas_diversity <= 0)) {
    stop("num_lambdas_diversity should be a positive integer")
  }
  
  # Shuffle the data
  n <- nrow(x)
  random.permutation <- sample(1:n, n)
  x.permutation <- x[random.permutation, ]
  y.permutation <- y[random.permutation]
  
  output <- Main_Ensemble_EN(x.permutation, y.permutation, num_lambdas_sparsity, num_lambdas_diversity, alpha, num_models, 
                             tolerance, max_iter, num_folds, num_threads)
  fn_call <- match.call()
  output <- construct.cv.SplitReg(output, fn_call, x, y)
  return(output)
}


