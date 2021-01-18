#' @title Compute lasso estimator
#'
#' @description Computes lasso, group lasso, scaled lasso, or scaled group lasso solution.
#' The outputs are coefficient-estimate and subgradient. If \code{type = "slasso"}
#' or \code{type = "sgrlasso"}, the output will include estimated standard deviation.
#'
#' @param X predictor matrix.
#' @param Y response vector.
#' @param type type of penalty. Must be specified to be one of the following:
#'  \code{"lasso", "grlasso", "slasso"} or \code{"sgrlasso"}, which correspond to
#'  lasso, group lasso, scaled lasso or scaled group lasso.
#' @param lbd penalty term of lasso. By letting this argument be \code{"cv.1se"} or
#' \code{"cv.min"}, users can have the cross-validated lambda that gives either minimum
#' squared error or that is within 1 std error bound.
#' @param weights weight vector with length equal to the number of groups. Default is
#' \code{weights = rep(1, max(group))}.
#' @param group \code{p} x \code{1} vector of consecutive integers describing the group structure.
#' The number of groups should be the same as max(group). Default is \code{group = 1:p}
#' , where \code{p} is number of covariates.
#' @param verbose logical. Only available for \code{type = "slasso"} or \code{type = "sgrlasso"}.
#' @param ... auxiliary arguments for \code{lbd = "cv.min", lbd = "cv.1se"}.
#' See \code{\link{cv.lasso}} for details.
#' @details
#' Computes lasso, group lasso, scaled lasso, or scaled group lasso solution.
#' Users can specify the value of lbd or choose to run cross-validation to get
#' optimal lambda in term of mean squared error. Coordinate decent algorithm is used
#' to fit scaled lasso and scaled group lasso models.
#'
#' @references
#' Mitra, R. and Zhang, C. H. (2016), "The benefit of group sparsity in group inference with
#' de-biased scaled group lasso," Electronic Journal of Statistics, 10, 1829-1873.
#'
#' Yang, Y. and Zou, H. (2015), “A Fast Unified Algorithm for Computing
#' Group-Lasso Penalized Learning Problems,” Statistics and Computing, 25(6), 1129-1141.
#'
#' @return \item{B0}{coefficient estimator.}
#' @return \item{S0}{subgradient.}
#' @return \item{sigmaHat}{estimated standard deviation.}
#' @return \item{lbd, weights, group}{same as input arguments.}
#' @examples
#' set.seed(123)
#' n <- 50
#' p <- 10
#' X <- matrix(rnorm(n*p), n)
#' Y <- X %*% c(1, 1, rep(0, p-2)) + rnorm(n)
#' #
#' # lasso
#' #
#' lassoFit(X = X, Y = Y, type = "lasso", lbd = .5)
#' #
#' # group lasso
#' #
#' lassoFit(X = X, Y = Y, type = "grlasso", lbd = .5, weights = rep(1,2),
#'            group = rep(1:2, each=5))
#' #
#' # scaled lasso
#' #
#' lassoFit(X = X, Y = Y, type = "slasso", lbd = .5)
#' #
#' # scaled group lasso
#' #
#' lassoFit(X = X, Y = Y, type = "sgrlasso", lbd = .5, weights = rep(1,2),
#'            group = rep(1:2, each=5))
#' @export
lassoFit <- function(X, Y, type, lbd,
                       group=1:ncol(X), weights=rep(1,max(group)), verbose = FALSE, ...)
{
  n <- nrow(X)
  p <- ncol(X)
  X <- as.matrix(X)
  Y <- matrix(Y, , 1)

  if (!type %in% c("lasso", "grlasso", "slasso", "sgrlasso")) {
    stop("type has to be either lasso, grlasso, slasso or sgrlasso.")
  }

  if (!all(group==1:p) && (!type %in% c("grlasso", "sgrlasso"))) {
    stop("Choose type to be either grlasso or sgrlasso if group-structure exists.")
  }

  # if (all(group==1:p) && (!type %in% c("lasso", "slasso"))) {
  #   stop("Choose type to be either lasso or slasso if group-structure does not exist.")
  # }

  if (!lbd  %in% c("cv.1se", "cv.min")) {
    if (!is.numeric(lbd) || lbd <= 0) {stop("invalid lbd input.")}
  }

  if (verbose) {cat("# Cross-validation \n")}
  if (lbd %in% c("cv.1se", "cv.min")) {
    lbdTEMP <- cv.lasso(X = X, Y = Y, group = group, weights = weights,
                        type = type, verbose=verbose, ...)
    if (lbd == "cv.1se") {lbd <- lbdTEMP$lbd.1se} else {
      lbd <- lbdTEMP$lbd.min
    }
  }

  # if (missing(lbd)) {
  #   if (type %in% c("lasso", "grlasso")) {
  #     lbd <- .37
  #   } else {
  #     lbd <- .5
  #   }
  # }

  # if (lbd <= 0) {
  #   stop("lbd has to be positive.")
  # }
  if (length(group) != p) {
    stop("length(group) has to be the same with ncol(X)")
  }
  if (length(weights) != length(unique(group))) {
    stop("length(weights) has to be the same with the number of groups")
  }
  if (any(weights <= 0)) {
    stop("weights should be positive.")
  }
  if (any(!1:max(group) %in% group)) {
    stop("group index has to be a consecutive integer starting from 1.")
  }

  if (length(Y) != n) {
    stop("dimension of X and Y are not conformable.")
  }

  IndWeights <- rep(weights,table(group))
  # scale X with weights
  Xtilde   <- scale(X,FALSE,scale=IndWeights)

  # slassoLoss <- function(X,Y,beta,sig,lbd) {
  #   n <- nrow(X)
  #   crossprod(Y-X%*%beta) / 2 / n / sig + sig / 2 + lbd * sum(abs(beta))
  # }

  if (type %in% c("lasso", "grlasso")) {
    # compute group lasso estimator B0 and S0
    # B0 <- coef(gglasso(x=Xtilde,y=Y,group=group,pf=rep(1,max(group)),lambda=lbd,intercept=FALSE))[-1] / IndWeights
    B0 <- grlassoFit(X = Xtilde, Y = Y, group = group, weights = rep(1, max(group)), lbd = lbd)$coef / IndWeights
    S0 <- (t(Xtilde) %*% Y - t(Xtilde) %*% Xtilde %*%
             (B0 * IndWeights)) / n / lbd
    #A <- which(B0!=0)
    return(list(B0=B0, S0=c(S0), lbd=lbd, weights=weights, group=group))
  } else {
    TEMP <- slassoFit.tilde(Xtilde = Xtilde, Y=Y, lbd=lbd, group=group, weights = weights, verbose = verbose)
    if (sum(TEMP$B0!=0) == (n-1)) {
      warning("Active set is too large. Try to increase the value of lbd.")
    }
    return(list(B0=TEMP$B0, S0=TEMP$S0, sigmaHat=TEMP$hsigma, lbd=lbd, weights=weights, group=group))
  }
}


#' @title Compute K-fold cross-validated mean squared error for lasso
#'
#' @description Computes K-fold cross-validated mean squared error
#' to propose a lambda value for lasso, group lasso, scaled lasso or scaled
#' group lasso.
#'
#' @param X predictor matrix.
#' @param Y response vector.
#' @param group \code{p} x \code{1} vector of consecutive integers describing the group structure.
#' The number of groups should be the same as max(group). Default is \code{group = 1:p}
#' , where \code{p} is number of covariates. See examples for a guideline.
#' @param weights weight vector with length equal to the number of groups. Default is
#' \code{rep(1, max(group))}.
#' @param type type of penalty. Must be specified to be one of the following:
#'  \code{"lasso", "grlasso", "slasso"} or \code{"sgrlasso"}, which correspond to
#'  lasso, group lasso, scaled lasso or scaled group lasso.
#' @param K integer. Number of folds
#' @param minlbd numeric. Minimum value of the lambda sequence.
#' @param maxlbd numeric. Maximum value of the lambda sequence.
#' @param num.lbdseq integer. Length of the lambda sequence.
#' @param parallel logical. If \code{parallel = TRUE}, uses parallelization.
#' Default is \code{parallel = FALSE}.
#' @param ncores integer. The number of cores to use for parallelization.
#' @param plot.it logical. If true, plots the squared error curve.
#' @param verbose logical.
#'
#' @return \item{lbd.min}{a value of lambda which gives a minimum squared error.}
#' @return \item{lbd.1se}{a largest lambda within 1 standard error from \code{lbd.min}.}
#' @return \item{lbd.seq}{lambda sequence.}
#' @return \item{cv}{mean squared error at each lambda value.}
#' @return \item{cvsd}{the standard deviation of cv.}
#'
#' @examples
#' set.seed(123)
#' n <- 30
#' p <- 50
#' group <- rep(1:(p/10),each=10)
#' weights <- rep(1, max(group))
#' X <- matrix(rnorm(n*p),n)
#' truebeta <- c(rep(1,5),rep(0,p-5))
#' Y <- X%*%truebeta + rnorm(n)
#'
#' # To accelerate the computational time, we set K=2 and num.lbdseq=2.
#' # However, in practice, Allowing K=10 and num.lbdseq > 100 is recommended.
#' cv.lasso(X = X, Y = Y, group = group, weights = weights, K = 2,
#' type = "grlasso", num.lbdseq = 2, plot.it = FALSE)
#' cv.lasso(X = X, Y = Y, group = group, weights = weights, K = 2,
#' type = "sgrlasso", num.lbdseq = 2, plot.it = FALSE)
#' @export
cv.lasso <- function( X, Y, group = 1:ncol(X), weights = rep(1,max(group)),
  type, K = 10L, minlbd, maxlbd, num.lbdseq = 100L, parallel = FALSE,
  ncores = 2L, plot.it = FALSE, verbose = FALSE)
{
  n <- nrow(X)
  p <- ncol(X)

  Y <- as.vector(Y)
  K <- as.integer(K)
  num.lbdseq <- as.integer(num.lbdseq)

  if (!parallel) {ncores <- 1}

  if(missing(minlbd)) {
    ifelse(type %in% c("lasso", "grlasso"), minlbd <- 0, minlbd <- .1)
  }
  if(missing(maxlbd)) {
    maxlbd <- if (type == "lasso")
      {
        max(abs(1/weights * t(X) %*% Y))/n
      } else if (type == "grlasso")
      {
        lbdTEMP <- c()
        for (i in 1:max(group)) {
          lbdTEMP[i] <- sqrt(crossprod(t(X[,group==i])%*%Y))/weights[i]
        }
        max(lbdTEMP / n)
      } else {2}
  }
  #--------------------
  # Error Handling
  #--------------------
  if (n < 10) {stop("Sample size is too small to run cross-validation.")}
  if (n < 30) {warning("Insufficient sample size. The result can be unstable.")}

  if (length(Y) != n) {
    stop("dimension of X and Y are not conformable.")
  }
  if (!type %in% c("lasso", "grlasso", "slasso", "sgrlasso")) {
    stop("type has to be either lasso, grlasso, slasso or sgrlasso.")
  }
  if (length(group) != p) {
    stop("group must have a same length with the number of X columns")
  }
  if (!all(group==1:p) && (!type %in% c("grlasso", "sgrlasso"))) {
    stop("Choose type to be either grlasso or sgrlasso if group-structure exists.")
  }
  parallelTemp <- ErrorParallel(parallel,ncores)
  parallel <- parallelTemp$parallel
  ncores <- parallelTemp$ncores
  if (length(weights) != length(unique(group))) {
    stop("weights has to have a same length as the number of groups")
  }
  if (any(weights <= 0)) {
    stop("weights should be positive.")
  }
  if (any(!1:max(group) %in% group)) {
    stop("group index has to be a consecutive integer starting from 1.")
  }
  if (any(c(minlbd,maxlbd) < 0)) {stop("minlbd/maxlbd should be non-negative")}
  if (minlbd >= maxlbd) {stop("minlbd is too large compared to maxlbd.")}
  if (num.lbdseq <= 0) {
    stop("num.lbdseq should be non-negative")
  }
  if (K <= 0) {
    stop("K should be a positive integer.")
  }

  all.folds <- split(sample(1:n), rep(1:K, length = n))

  #index=seq(0,max(t(X)%*%Y),length=100)
  index <- seq(minlbd,maxlbd,length=num.lbdseq)[-1]
  residmat <- matrix(0, length(index), K)

  FF <- function(x,omit) {
    fit <- lassoFit(X=X[-omit,,drop=FALSE],Y=Y[-omit],type=type,
                      lbd=index[x],group=group,weights=weights)$B0
    fit <- X[omit,,drop=FALSE]%*%fit
    return(mean((Y[omit]-fit)^2))
  }

  for (i in seq(K)) {
    omit <- all.folds[[i]]
    residmat[,i] <- do.call(rbind,parallel::mclapply(1:length(index),
      FF,omit = omit, mc.cores = ncores))
    if(verbose) {cat("\n CV Fold", i, "\n\n")}
  }
  #apply(residmat, 2, which.min)
  cv <- apply(residmat, 1, mean)
  cvsd <- apply(residmat, 1, sd)
  index.min.cv <- which.min(cv)

  err.1se <- cvsd[index.min.cv] + cv[index.min.cv]

  lbd.1se <- max(index[cv <= err.1se], na.rm=TRUE)
  lbd.min <- index[index.min.cv]

  if (plot.it) {
    matplot(x=index, y=cbind(cv,cv-cvsd,cv+cvsd), type="l",
            lty=c(1,2,2), col=c(1,2,2),
            xlab="lambda",ylab="mean squared error",
            main="cross validation")
    abline(v=lbd.min,lty=2)
    abline(v=lbd.1se,lty=3)
  }
  return(list(lbd.seq=index, cv=cv, cvsd=cvsd,
              lbd.min=lbd.min, lbd.1se=lbd.1se,
              type=type))
}
