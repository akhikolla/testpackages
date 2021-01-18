#' @title Function to fit penalized GEE by I-CGD algorithm.
#' @description  This function fits a \eqn{L_1} penalized GEE model to longitudinal
#'  data by I-CGD algorithm or re-weighted least square algorithm.
#' @param X A design matrix of dimension  \code{(nm) * p}.
#' @param y A response vector of length  \code{m * n}.
#' @param id A vector for identifying subjects/clusters.
#' @param family A family object  representing one of the built-in families.
#' Families supported here are the same as in \pkg{PGEE}, e.g, \code{binomial},
#' \code{gaussian}, \code{gamma} and \code{poisson}, and the corresponding
#' link functions are supported, e.g, \code{identity}, and \code{probit}.
#' @param lambda A user supplied value for the penalization parameter.
#' @param method The algorithms that are available. \code{"CGD"} represents the
#' I-CGD algorithm, and \code{"RWL"} represents re-weighted least square algorithm.
#' @param corstr A character string that indicates the correlation structure among
#' the repeated measurements of a subject. Structures supported in \code{LassoGEE}
#' are "AR1", "exchangeable", "unstructured", and "independence". The default
#' \code{corstr} type is "independence".
#' @param beta.ini User specified initial values for regression parameters.
#' The default value is \code{NULL}.
#' @param R User specified correlation matrix. The default value is \code{NULL}.
#' @param scale.fix  A logical variable. The default value is \code{TRUE},
#' then the value of the scale parameter is fixed to \code{scale.value}.
#' @param scale.value If  \code{scale.fix = TRUE}, a numeric value will be assigned
#'  to the fixed scale parameter. The default value is 1.
#' @param maxiter The maximum number of iterations used in the algorithm.
#' The default value is 50.
#' @param tol The tolerance level used in the algorithm. The default value is \code{1e-3}.
#' @param Mv If either "stat_M_dep", or "non_stat_M_dep" is specified in corstr,
#' then this assigns a numeric value for Mv. Otherwise, the default value is NULL.
#' @param silent A logical variable; if false, the iteration counts
#' at each iteration of CGD are printed. The default value is TRUE.
#' @param verbose A logical variable; Print the out loop iteration counts. The default value is TRUE.
#' @return A list containing the following components:
#'   \item{betaest}{return final estimation}
#'   \item{beta_all_step}{return estimate in each iteration}
#'   \item{inner.count}{iterative count in each stage}
#'   \item{outer.iter}{iterate number of outer loop}
#' @references Li, Y., Gao, X., and Xu, W. (2020). Statistical consistency for
#' generalized estimating equation with \eqn{L_1} regularization.
#' @export
#' @importFrom PGEE mycor_gee1
#' @importFrom stats binomial gaussian quantile family
#' @import MASS
#' @import Rcpp
#' @useDynLib LassoGEE, .registration = TRUE
#' @import mvtnorm
#' @import SimCorMultRes
#' @seealso cv.LassoGEE
#' @examples
#' # required R package
#' library(mvtnorm)
#' library(SimCorMultRes)
#' #
#' set.seed(123)
#' p <- 200
#' s <- ceiling(p^{1/3})
#' n <- ceiling(10 * s * log(p))
#' m <- 4
#' # covariance matrix of p number of continuous covariates
#' X.sigma <- matrix(0, p, p)
#' {
#'   for (i in 1:p)
#'     X.sigma[i,] <- 0.5^(abs((1:p)-i))
#' }
#'
#' # generate matrix of covariates
#' X <- as.matrix(rmvnorm(n*m, mean = rep(0,p), X.sigma))
#'
#' # true regression parameter associated with the covariate
#' bt <- runif(s, 0.05, 0.5) # = rep(1/s,s)
#' beta.true <- c(bt,rep(0,p-s))
#' # intercept
#' beta_intercepts <- 0
#' # unstructure
#' tt <- runif(m*m,-1,1)
#' Rtmp <- t(matrix(tt, m,m))%*%matrix(tt, m,m)+diag(1,4)
#' R_tr <- diag(diag(Rtmp)^{-1/2})%*%Rtmp%*%diag(diag(Rtmp)^{-1/2})
#' diag(R_tr) = round(diag(R_tr))
#'
#' # library(SimCorMultRes)
#' # simulation of clustered binary responses
#' simulated_binary_dataset <- rbin(clsize = m, intercepts = beta_intercepts,
#'                                  betas = beta.true, xformula = ~X, cor.matrix = R_tr,
#'                                  link = "probit")
#' lambda <- 0.2* s *sqrt(log(p)/n)
#' data = simulated_binary_dataset$simdata
#' y = data$y
#' X = data$X
#' id = data$id
#'
#' ptm <- proc.time()
#' nCGDfit = LassoGEE(X = X, y = y, id = id, family = binomial("probit"),
#'                  lambda = lambda, corstr = "unstructured")
#' proc.time() - ptm
#' betaest <- nCGDfit$betaest
#'
LassoGEE <- function(X, y, id, family = binomial("probit"), lambda,
                     corstr = "independence", method = c("CGD", "RWL"),
                     beta.ini = NULL, R = NULL, scale.fix = TRUE,
                     scale.value = 1, maxiter = 50, tol = 1e-3,
                     silent = TRUE, Mv = NULL, verbose = TRUE)  {

  call <- match.call()
  method=match.arg(method)

  if (is.null(id)) {
    stop("Id variable not found!")
  }
  if (length(id) != length(y))
    stop("Id and y do not have the same length!")
  if (!(is.double(X)))
    X <- as.double(X)
  if (!(is.double(y)))
    y <- as.double(y)
  if (!(is.double(id)))
    id <- as.double(id)


  N<-length(unique(id))
  nx=ncol(X)

  avec <- as.integer(unlist(lapply(split(id, id), "length")))
  maxclsz <- max(avec)
  maxcl <- maxclsz
  nt <- avec
  nobs <- sum(nt)
  xnames <- dimnames(X)[[2]]
  if (is.null(xnames)) {
    xnames <- paste("x", 1:dim(X)[2], sep = "")
    dimnames(X) <- list(NULL, xnames)
  }

  if (!(is.double(N)))
    N <- as.double(N)
  if (!(is.double(maxcl)))
    maxcl <- as.double(maxcl)
  if (!(is.double(nobs)))
    nobs <- as.double(nobs)
  if (missing(lambda))
    stop("A value is not assiged for lambda!")
  # if (missing(pindex))
  #   pindex = NULL
  if (missing(family))
    family = gaussian(link = "identity")
  if (missing(corstr))
    corstr = "independence"
  if (missing(Mv))
    Mv <- NULL
  if (corstr == "stat_M_dep" && is.null(Mv))
    stop("corstr is assumed to be 'stat_M_dep' but Mv is not specified!")
  if (corstr == "non_stat_M_dep" && is.null(Mv))
    stop("corstr is assumed to be 'non_stat_M_dep' but Mv is not specified!")
  if ((corstr != "stat_M_dep" && corstr != "non_stat_M_dep") &&
      !is.null(Mv))
    stop("Mv is specified while corstr is assumed to be neither \n'stat_M_dep' nor 'non_stat_M_dep'!")
  if (corstr == "non_stat_M_dep" && length(unique(nt)) != 1)
    stop("corstr cannot be assumed to be 'non_stat_M_dep' for unbalanced data!")
  if (corstr == "unstructured" && length(unique(nt)) != 1)
    stop("corstr cannot be assumed to be 'unstructured' for unbalanced data!")
  if (missing(R))
    R <- NULL
  if (corstr == "fixed" && is.null(R))
    stop("corstr is assumed to be 'fixed' but R is not specified!")
  if (corstr != "fixed" && !is.null(R))
    stop("R is specified although corstr is not assumed to be 'fixed'!")
  if (!is.null(R)) {
    Rr <- nrow(R)
    if (Rr != ncol(R))
      stop("R is not square!")
    if (Rr < maxclsz) {
      stop("R is not big enough to accommodate some clusters!")
    }
    else if (Rr > maxclsz) {
      stop("R is larger than the maximum cluster!")
    }
  }
  if (missing(scale.fix))
    scale.fix <- TRUE
  scale.fix <- as.integer(scale.fix)
  if (missing(scale.value))
    scale.value = 1
  scale.value <- as.integer(scale.value)
  if (missing(maxiter))
    maxiter <- 100
  maxiter <- as.integer(maxiter)
  if (missing(tol))
    tol = 0.0001
  tol = as.double(tol)
  if (missing(silent))
    silent <- TRUE
  silent <- as.integer(silent)
  if (is.character(family))
    family <- get(family)
  if (is.function(family))
    family <- family()
  links <- c("identity", "log", "logit", "inverse", "probit",
             "cloglog")
  fams <- c("gaussian", "poisson", "binomial", "Gamma", "quasi")
  varfuns <- c("constant", "mu", "mu(1-mu)", "mu^2")
  corstrs <- c("independence", "fixed", "stat_M_dep", "non_stat_M_dep",
               "exchangeable", "AR-1", "unstructured")
  linkv <- as.integer(match(c(family$link), links, -1))
  if (linkv < 1)
    stop("unknown link!")
  famv <- match(family$family, fams, -1)
  if (famv < 1)
    stop("unknown family")
  if (famv <= 4) {
    varfunv <- famv
  } else {
    varfunv <- match(family$varfun, varfuns, -1)
  }

  if (varfunv < 1)
    stop("unknown varfun!")
  corstrv <- as.integer(match(corstr, corstrs, -1))
  if (corstrv < 1)
    stop("unknown corstr!")

  Mv <- as.integer(Mv)
  if (!is.null(beta.ini)) {
    betaest <- matrix(beta.ini, ncol = 1)
    if(nrow(betaest) != nx) {
      stop("Dimension of beta != ncol(X)!")
      }
    #message("user\'s initial regression estimate")
  } else {
    betaest = c(rep(0,nx))
  }

  aindex=cumsum(nt)
  index=c(0,aindex[-length(aindex)])
  diff<-1
  iter<-0
  # maxiter <- 30
  count <- c()
  beta_all_step <- list()

  while(iter < maxiter) {

    R.fi.hat <- PGEE::mycor_gee2( N, nt, y, X, family, beta_new = betaest,
                                  corstr = corstr, Mv = Mv, maxclsz = maxclsz,
                                  R = R, scale.fix = scale.fix,
                                  scale.value = scale.value)
    Rhat <- R.fi.hat$Ehat
    fihat <- R.fi.hat$fi

    eta <- drop(X%*%betaest)
    mu=family$linkinv(eta)
    mu_eta = family$mu.eta(eta)
    vari = family$variance(mu)
    S.H.E.M.val = SHM(X = X, y = y, mu = mu, mu_eta = mu_eta, vari = vari,
                      nt = nt, index = index, Rhat = Rhat,
                      N = N, fihat = fihat)
    S<-S.H.E.M.val$S
    v<-S.H.E.M.val$H
    u<- v%*%betaest + S

    if(method == "CGD") {
      inner.fit <- ac_prox_grad(u = u, v = v, lambda = rep(lambda*N, nx),
                                tol = tol, maxiter = maxiter, silent = silent)
      betaest1 <- inner.fit$beta_k
      diff<-sum(abs(betaest-betaest1))
      betaest<-betaest1

      iter<-iter+1
      count[iter] <- inner.fit$k
      beta_all_step[[iter]] <- inner.fit$beta_inner_step
      if (diff <= tol) {
        Smat <- S.H.E.M.val$Smat
        if(verbose) cat("iter: ",iter, "diff: ",diff,"\n")
        break
      }

    } else {
      betaest1 <- WLreglass(v = v, u = u, lambda = rep(lambda*N, nx), tol = tol)
      diff<-sum(abs(betaest-betaest1))
      betaest<-betaest1

      iter<-iter+1
      if (diff <= tol) {
        Smat <- S.H.E.M.val$Smat
        if(verbose) cat("iter: ",iter, "diff: ",diff,"\n")
        break
      }

    }

    # if (silent == 0)
    #   cat("iter", iter, "betaest", betaest, "diff",
    #       diff, "\n")
    # if (diff <= tol)
    #   break
  }

  fit <- list()
  attr(fit, "class") <- c("LassoGEE")
  fit$title <- paste("LassoGEE by", method, "Algorithm")
  fit$version <- "Version: 1.0"
  links <- c("Identity", "Logarithm", "Logit", "Reciprocal",
             "Probit", "Cloglog")
  varfuns <- c("Gaussian", "Poisson", "Binomial", "Gamma")
  corstrs <- c("Independent", "Fixed", "Stationary M-dependent",
               "Non-Stationary M-dependent", "Exchangeable", "AR-1",
               "Unstructured")
  fit$model <- list()
  fit$model$link <- links[linkv]
  fit$model$varfun <- varfuns[varfunv]
  fit$model$corstr <- corstrs[corstrv]
  if (!is.na(match(c(corstrv), c(3, 4))))
    fit$model$M <- Mv
  fit$call <- call
  fit$nobs <- nobs
  fit$outer.iter <- iter
  fit$betaest <- as.vector(betaest)
  fit$nas <- is.na(fit$betaest)
  if(method == "CGD"){
    fit$beta_all_step <- beta_all_step
    fit$inner.count <- count
  }
  names(fit$betaest) <- xnames
  eta <- as.vector(X %*% fit$betaest)
  fit$linear.predictors <- eta
  mu <- as.vector(family$linkinv(eta))
  fit$fitted.values <- mu
  fit$residuals <- y - mu
  fit$family <- family
  fit$y <- as.vector(y)
  fit$id <- as.vector(id)
  fit$max.id <- maxcl
  fit$working.correlation <- Rhat[1:maxclsz, 1:maxclsz, which(avec ==
                                                                maxclsz)[1]]
  fit$scale <- fihat
  fit$S <- S
  fit$Smat <- Smat
  fit$lambda.value <- lambda
  fit$xnames <- xnames
  fit$error <- diff
  fit

  # return(list(betaest = c(betaest), beta_all_step = beta_all_step,
  #             inner.count = count, outer.iter = iter))
}






# proximal of L1 norm
prox_L1 = function(x, lambda){
  return(sign(x) * pmax(0, abs(x) - lambda))
}

##
ac_prox_grad <- function(u, v, lambda, tol, maxiter, silent) {

  # initilization
  L.max <- max(eigen(v)$values)
  L.min <- min(eigen(v)$values)
  L <- ifelse(L.max > 0, L.max, L.min)

  ## step size
  # nu <- 0.5             ## line search for t paramter

  beta_last <-solve(v)%*%u          ## initilization of x
  z_last <- beta_last
  t_last <- 1

  k <- 0
  # maxiter <- 100           ## number of iterations
  # tol <- 1e-4
  beta_inner_step <- beta_last

  while(k < maxiter) {

    need_project <- z_last-(v%*%z_last)/L + u/L
    beta_new <- prox_L1(need_project, lambda/abs(L))
    distance <- beta_new-beta_last
    k <- k+1
    if (silent == 0) cat(k, fill = TRUE)
    beta_inner_step <- cbind(beta_inner_step, beta_new)

    if (sqrt(sum(distance^2)) <= tol * sqrt(sum(beta_last^2))) {
      break
    }
    t_new <- (1 + sqrt(1 + 4 * t_last * t_last)) / 2
    z_new <- beta_new + (t_last - 1) * distance / t_new
    beta_last <- beta_new
    t_last <- t_new
    z_last <- z_new
  }

  return(list(beta_k = beta_new, beta_inner_step = beta_inner_step, k = k))

}



##
WLreglass<-function(v, u, lambda, tol){

  ###x is the design matrix, y is the response, w is the weight matrix, we do a regression using lasso penalty
  ### using lasso in regression, with v=t(x)w x, u=t(x)w y  using Friendmans's way of iteration

  p<-dim(v)[1]
  diff<-1

  beta<-solve(v)%*%u


  while (diff > tol){
    oldbeta<-beta
    z<-u-v%*%oldbeta


    for (j in 1:p){
      beta[j]<- prox_L1((z[j]+v[j,j]*oldbeta[j]),lambda[j])/(v[j,j])  ###only sentence is changed
      z<-z-(beta[j]-oldbeta[j])*v[,j]
      #print(j)
    }
    diff<-max(abs(beta-oldbeta))
    #print(diff)
  }
  # beta_k = beta
  return(beta)
}
