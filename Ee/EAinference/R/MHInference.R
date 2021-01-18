#' @title Post-inference with lasso estimator
#'
#' @description Provides confidence intervals for the set of active coefficients
#' of lasso using Metropolis-Hastings sampler.
#'
#' @param X predictor matrix.
#' @param Y response vector.
#' @param lbd penalty term of lasso. By letting this argument be \code{"cv.1se"} or
#' \code{"cv.min"}, users can have the cross-validated lambda that gives either minimum
#' squared error or that is within 1 std error bound.
#' @param weights weight vector with length equal to the number of coefficients.
#' Default is \code{rep(1, ncol(X))}.
#' @param tau numeric vector. Standard deviation of proposal distribution
#'  for each beta. Adjust the value to get relevant level of acceptance rate.
#'  Default is \code{rep(1, ncol(X))}.
#' @param sig2.hat variance of error term.
#' @param alpha confidence level for confidence interval.
#' @param nChain the number of chains. For each chain, different plug-in beta will be generated
#' from its confidence region.
#' @param niterPerChain the number of iterations per chain.
#' @param method Type of robust method. Users can choose either \code{"coeff"} or \code{"mu"}.
#' @param parallel logical. If \code{parallel = TRUE}, uses parallelization.
#' Default is \code{parallel = FALSE}.
#' @param ncores integer. The number of cores to use for parallelization.
#' @param returnSamples logical. If \code{returnSamples = TRUE}, print Metropolis-Hastings samples.
#' @param ... auxiliary \code{\link{MHLS}} arguments.
#' @details
#' This function provides post-selection inference for the active coefficients selected by lasso.
#' Uses Metropolis-Hastings sampler with multiple chains to draw from the
#' distribution under a fixed active set and generates \code{(1-alpha)}
#' confidence interval for each active coefficients.
#' Set \code{returnSamples = TRUE} to check the Metropolis-Hastings samples.
#' Check the acceptance rate and adjust \code{tau} accordingly.
#' We recommend to set \code{nChain >= 10} and \code{niterPerChain >= 500}.
#'
#' @return \item{MHsamples}{a list of class MHLS.}
#' @return \item{confidenceInterval}{(1-alpha) confidence interval
#' for each active coefficient.}
#' @examples
#' set.seed(123)
#' n <- 6
#' p <- 10
#' X <- matrix(rnorm(n*p),n)
#' Y <- X %*% rep(1,p) + rnorm(n)
#' sig2 <- 1
#' lbd <- .37
#' weights <- rep(1,p)
#' parallel <- (.Platform$OS.type != "windows")
#' postInference.MHLS(X = X, Y = Y, lbd = lbd, sig2.hat = 1, alpha = .05,
#' nChain = 3, niterPerChain = 20, method = "coeff", parallel = parallel)
#' postInference.MHLS(X = X, Y = Y, lbd = lbd, sig2.hat = 1, alpha = .05,
#' nChain = 3, niterPerChain = 20, method = "coeff", parallel = parallel, returnSamples = TRUE)
#' postInference.MHLS(X = X, Y = Y, lbd = lbd, sig2.hat = 1, alpha = .05,
#' nChain = 3, niterPerChain = 20, method = "mu", parallel = parallel)
#' postInference.MHLS(X = X, Y = Y, lbd = lbd, sig2.hat = 1, alpha = .05,
#' nChain = 3, niterPerChain = 20, method = "mu", parallel = parallel, returnSamples = TRUE)
#' @export
postInference.MHLS <- function(X, Y, lbd, weights = rep(1, ncol(X)),
  tau = rep(1, ncol(X)), sig2.hat, alpha = .05, nChain = 10, method,
  niterPerChain = 500, parallel = FALSE, ncores = 2L, returnSamples=FALSE, ...)
{
  # nChain : the number of MH chains
  # niterPerChain : the number of iteration for each chain
  # B0, S0 : The lasso estimator
  # tau : same as in MHLS function

  LassoEst <- lassoFit(X=X, Y=Y, type = "lasso", lbd=lbd, weights=weights)
  B0 <- LassoEst$B0
  S0 <- LassoEst$S0
  lbd <- LassoEst$lbd
  A <- which(B0!=0)
  if (length(A)==0) {
    stop("Given lbd, active set is empty.")
  }

  Y <- matrix(Y, , 1)
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  beta.refit <- rep(0,p)
  beta.refit[A] <- coef(lm(Y~X[,A]+0))

  nChain <- as.integer(nChain)
  niterPerChain <- as.integer(niterPerChain)
  #--------------------
  # Error Handling
  #--------------------
  parallelTemp <- ErrorParallel(parallel,ncores)
  parallel <- parallelTemp$parallel
  ncores <- parallelTemp$ncores

  if (!method %in% c("coeff", "mu")) {
    stop("Invalide method type.")
  }

  if (length(tau) != p) {
    stop("length(tau) has to be the same with col(X)")
  }

  if (nrow(X) != nrow(Y)) {
    stop("The dimension of X and Y are not conformable.")
  }
  if (sig2.hat <=0 || lbd <= 0) {
    stop("sig2.hat and/or lbd have to be positive.")
  }
  if (length(weights) != ncol(X)) {
    stop("length(weights) has to be the same with col(X).")
  }
  if (any(weights <= 0)) {
    stop("weights should be positive.")
  }
  if (alpha <=0 || alpha >=1) {
    stop("alpha needs to be between 0 and 1.")
  }
  if (any(c(nChain,niterPerChain) <= 0)) {
    stop("nChain & niterPerChain have to be a positive integer.")
  }
  # if (all.equal(coef(gglasso(X, Y, pf = weights, group = 1:p, loss="ls",
  #                             intercept=F, lambda=lbd))[-1],B0) != TRUE ||
  #     all.equal(c(((t(X)/weights)%*%Y - (t(X) /weights) %*% X %*% B0) / n / lbd)
  #                , S0) != TRUE) {
  #   stop("Invalid B0 or S0, use lassoFit to get a valid lasso solution.")
  # }
  # Draw samples of pluginbeta from the 95% confidence
  #  region boundary of restricted lse.
  # If nChain ==1, we just use restricted lse.

  if (missing(sig2.hat)) {
    if (length(A) >= nrow(X)) {
      stop("If size of active set matches with nrow(X), sig2.hat needs to be provided.")
    }
    sig2.hat <- summary((lm(Y~X[,A]+0)))$sigma^2
  }

  if (method == "coeff") {
    Plugin.seq <- Pluginbeta.MHLS(X = X, Y = Y, A = A, nPlugin = nChain,
                                  sigma.hat = sqrt(sig2.hat))
    betaCenter <- beta.refit
  } else {
    Plugin.seq <- PluginMu.MHLS(X = X, Y = Y, lbd = lbd, sigma.hat = sqrt(sig2.hat),
      # ratioSeq = seq(0,1,by=0.01), alpha = 0.05, nChain = nChain, niter = 100,
      alpha = 0.05, nChain = nChain, niter = 100,
      method = "boundary", parallel = parallel, ncores = ncores)
    betaCenter <- rep(0,p)
    betaCenter[A] <- solve(crossprod(X[,A]))%*%t(X[,A])%*% Plugin.seq[nChain+1, ]
  }

  FF <- function(x) {
    MHLS(X = X, PE = Plugin.seq[x,], sig2 = sig2.hat, lbd = lbd,
         weights = weights, niter=niterPerChain,
         burnin = 0, B0 = B0, S0 = S0, tau = tau, PEtype = method, verbose=FALSE, ...)
  }

  if (parallel) {
    TEMP <- parallel::mclapply(1:nChain,FF, mc.cores = ncores)
  } else {
    TEMP <- lapply(1:nChain,FF)
  }
  names(TEMP) <- paste("Chain",1:nChain,sep="")

  MCSAMPLE <- TEMP[[1]]
  if (nChain > 1) {
    for (i in 2:nChain) {
      MCSAMPLE$beta <- rbind(MCSAMPLE$beta, TEMP[[i]]$beta)
      MCSAMPLE$subgrad <- rbind(MCSAMPLE$subgrad, TEMP[[i]]$subgrad)
      MCSAMPLE$acceptHistory <- MCSAMPLE$acceptHistory + TEMP[[i]]$acceptHistory
      MCSAMPLE$niteration <- MCSAMPLE$niteration + TEMP[[i]]$niteration
      MCSAMPLE$burnin <- MCSAMPLE$burnin + TEMP[[i]]$burnin
    }
  }

  # Using MH samples, refit the coeff.
  RefitBeta <- Refit.MHLS(X,weights,lbd,MCSAMPLE)
  if (returnSamples) {
    return(list(MHsamples = TEMP, pluginValue = Plugin.seq, method = method,
            confidenceInterval = CI.MHLS(betaRefitMH = RefitBeta, betaCenter = betaCenter,
                                         betaRefit = beta.refit, alpha = alpha)))
  } else {
    return(CI.MHLS(betaRefitMH = RefitBeta, betaRefit = beta.refit,
                   betaCenter = betaCenter, alpha = alpha))
  }
}
