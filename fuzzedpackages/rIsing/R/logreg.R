#' @title L1 Regularized Logistic Regression
#'
#' @description L1 Regularized logistic regression using OWL-QN L-BFGS-B optimization.
#'
#' @param X The design matrix.
#' @param y Vector of binary observations of length equal to \code{nrow(X)}.
#' @param nlambda (positive integer) The number of parameters in the regularization path (default 50).
#' @param lambda.min.ratio (non-negative double) The ratio of \code{max(lambda) / min(lambda)} (default \code{1e-3}).
#' @param lambda A user-supplied vector of regularization parameters. Under the default option (\code{NULL}), the function computes a regularization path using the input data.
#' @param scale (boolean) Whether to scale \code{X} before running the regression. The output parameters will always be rescaled. Use \code{FALSE} if \code{X} is already scaled.
#' @param type (integer 1 or 2) Type 1 aggregates the input data based on repeated rows in \code{X}. Type 2 (default) uses the data as is, and is generally faster. Use Type 1 if the data contains several repeated rows.
#'
#' @return A list containing the matrix of fitted weights (\code{wmat}), the vector of regularization parameters, sorted in decreasing order (\code{lambda}), and the vector of log-likelihoods corresponding to \code{lambda} (\code{logliks}).
#'
#' @examples
#' # simulate some linear regression data
#' n <- 1e3
#' p <- 100
#' X <- matrix(rnorm(n*p),n,p)
#' wt <- sample(seq(0,9),p+1,replace = TRUE) / 10
#' z <- cbind(1,X) %*% wt + rnorm(n)
#' probs <- 1 / (1 + exp(-z))
#' y <- sapply(probs, function(p) rbinom(1,1,p))
#'
#' m1 <- logreg(X, y)
#' m2 <- logreg(X, y, nlambda = 100, lambda.min.ratio = 1e-4, type = 1)
#'
#' \dontrun{
#' # Performance comparison
#' library(glmnet)
#' library(microbenchmark)
#' nlambda = 50; lambda.min.ratio = 1e-3
#' microbenchmark(
#'   logreg_type1 = logreg(X, y, nlambda = nlambda,
#'                          lambda.min.ratio = lambda.min.ratio, type = 1),
#'   logreg_type2 = logreg(X, y, nlambda = nlambda,
#'                          lambda.min.ratio = lambda.min.ratio, type = 2),
#'   glmnet       = glmnet(X, y, family = "binomial",
#'                          nlambda = nlambda, lambda.min.ratio = lambda.min.ratio),
#'   times = 20L
#' )
#' }
#'
#' @export
logreg <- function(X, y, nlambda = 50, lambda.min.ratio = 1e-3, lambda = NULL, scale = TRUE, type = 2) {
  Xy <- cbind(y,X)
  if (anyNA(Xy)) {
    Xy <- Xy[complete.cases(Xy),]
  }
  X <- Xy[,-1,drop=FALSE]; y <- Xy[,1]
  yu <- unique(y)
  if (length(yu) > 2) stop("y must contain binary data.")
  if (!all(yu %in% c(0,1))) {
    y <- c(0,1)[as.factor(y)]
  }
  if (nlambda < 0) stop("nlambda must be >= 0.")
  if (nlambda == 0) lambda <- 0
  if (lambda.min.ratio < 0) stop("lambda.min.ratio must be >=0.")
  if (!(scale %in% c(0,1))) stop("scale must be a boolean value.")

  if (!(type %in% c(1,2))) stop("Invalid type selected (choose either 1 or 2)")

  if (type == 2) {
    logreg2(X, y, nlambda, lambda.min.ratio, lambda, scale)
  } else {
    logreg1(X, y, nlambda, lambda.min.ratio, lambda, scale)
  }
}

logreg1 <- function(X, y, nlambda = 50, lambda.min.ratio = 1e-3, lambda = NULL, scale = TRUE) {
  regpath <- is.null(lambda)
  lrs <- logreg_setup(X, y, scale, regpath, nlambda, lambda.min.ratio) # logreg setup
  if (regpath) {
    lambda <- lrs$lambda
  } else {
    lambda <- sort(lambda, decreasing = T)
  }

  dt <- as.data.table(cbind(y, lrs$Xs))
  dt_ <- dt[, list(b = .N, y = sum(y == 1)), by = setdiff(names(dt), "y")]
  X <- as.matrix(dt_[,!c("b","y"),with = FALSE])
  y <- dt_[["y"]]
  b <- dt_[["b"]]
  cpp_out <- logreg_cpp(X, y, b, lrs$means, lrs$sds, lambda)

  return ( list(wmat = cpp_out$wmat, logliks = cpp_out$logliks, lambda = lambda) )
}

logreg2 <- function(X, y, nlambda = 50, lambda.min.ratio = 1e-3, lambda = NULL, scale = TRUE) {
  return ( logreg_cpp2(X, y, lambda, nlambda, lambda.min.ratio, scale))
}
