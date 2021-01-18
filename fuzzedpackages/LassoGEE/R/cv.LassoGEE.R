#' @title Cross-validation for LassoGEE.
#' @description  Does k-fold cross-validation for LassoGEE to select tuning parameter
#' value for longitudinal data with working independence structure.
#' @param X A design matrix of dimension  \code{(nm) * p}.
#' @param y A response vector of length  \code{m * n}.
#' @param id A vector for identifying subjects/clusters.
#' @param family A family object: a list of functions and expressions
#' for defining link and variance functions. Families supported here is same as
#'   in \pkg{PGEE} which are binomial, gaussian, gamma and poisson.
#' @param lambda.vec A vector of tuning parameters that will be used in the
#' cross-validation.
#' @param method The algorithms that are available. \code{"CGD"} represents the
#' I-CGD algorithm, and \code{"RWL"} represents re-weighted least square algorithm.
#' @param fold The number of folds used in cross-validation.
#' @param scale.fix A logical variable; if true, the scale parameter is
#' fixed at the value of scale.value. The default value is TRUE.
#' @param scale.value If  \code{scale.fix = TRUE}, this assignes a
#' numeric value to which the scale parameter should be fixed.
#' The default value is 1.
#' @param maxiter The number of iterations that is used in the estimation algorithm.
#' The default value is \code{50}.
#' @param tol The tolerance level that is used in the estimation algorithm.
#' The default value is \code{1e^-3}.
#' @return An object class of cv.LassoGEE.
#' @references Li, Y., Gao, X., and Xu, W. (2020). Statistical consistency for
#' generalized estimating equation with \eqn{L_1} regularization.
#' @seealso LassoGEE
#' @importFrom caret createFolds
#' @export

cv.LassoGEE <- function (X, y, id, family, method = c("CGD", "RWL"),
                         scale.fix, scale.value,
                         fold, lambda.vec, maxiter, tol)
{
  call_cv <- match.call()
  method=match.arg(method)

  if (missing(family))
    family = gaussian(link = "identity")
  if (missing(scale.fix))
    scale.fix <- TRUE
  if (missing(scale.value))
    scale.value <- 1
  if (missing(maxiter))
    maxiter <- 100
  if (missing(tol))
    tol <- 1e-3
  # N <- length(unique(id))
  # nx <- dim(X)[2]
  lam.min <- -1
  cv.min <- Inf
  cv.vect <- NULL
  theFolds =  caret::createFolds(unique(id), k= fold)
  for (j in 1:length(lambda.vec)) {
    lam.temp <- lambda.vec[j]
    cv.value <- 0
    for (k in 1:fold) {

      x.train = X[!id %in% theFolds[[k]],];
      y.train = y[!id %in% theFolds[[k]]];
      id.train <- id[!id %in% theFolds[[k]]]
      # y.train <- y[-index.cv]
      # if (colnames(X)[1] == "(Intercept)")
      #   x.train <- X[-index.cv, -1]
      # else x.train <- X[-index.cv, ]
      # id.train <- id[-index.cv]
      # data.train = data.frame(id = id.train, y = y.train, x.train)
      nCGDfit = LassoGEE(X = x.train, y = y.train, id = id.train,
                         family = family, method = method,
                         lambda = lam.temp, corstr = "independence")
      beta.train <- nCGDfit$betaest

      x.cv= X[ id %in% theFolds[[k]],];
      y.cv = y[ id %in% theFolds[[k]]];
      id.cv <- id[id %in% theFolds[[k]]]

      # y.cv <- y[index.cv]
      # x.cv <- X[index.cv, ]
      # id.cv <- id[index.cv]
      yy = y.cv
      eta = x.cv %*% beta.train
      mu = family$linkinv(eta)
      cv.value <- cv.value + sum((family$dev.resids(yy, mu, wt = 1)))
    }
    cv.vect <- c(cv.vect, cv.value)
    if (cv.value < cv.min) {
      lam.min <- lam.temp
      cv.min <- cv.value
    }
  }
  out <- list()
  attr(out, "class") <- c("cv.LassoGEE")
  out$fold = fold
  out$lam.vect = lambda.vec
  out$cv.vect = cv.vect
  out$lam.opt = lam.min
  out$cv.min = cv.min
  out$call <- call_cv
  out
}
