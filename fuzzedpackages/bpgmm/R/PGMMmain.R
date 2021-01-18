#' bpgmm Model-Based Clustering Using Baysian PGMM Carries out model-based clustering using parsimonious Gaussian mixture models. MCMC are used for parameter estimation. The RJMCMC is used for model selection.
#'
#' @import stats MASS mcmcse pgmm label.switching fabMix
#' @param X the observation matrix with size p * m
#' @param mInit the number of initial clusters
#' @param mVec the range of the number of clusters
#' @param qnew the number of factor for a new cluster
#' @param delta scaler hyperparameters
#' @param ggamma scaler hyperparameters
#' @param burn the number of burn in iterations
#' @param niter the number of iterations
#' @param constraint the pgmm initial constraint, a vector of length three with binary entry. For example, c(1,1,1) means the fully constraint model
#' @param dVec a vector of hyperparameters with length three, shape parameters for alpha1, alpha2 and bbeta respectively
#' @param sVec sVec a vector of hyperparameters with length three, rate parameters for alpha1, alpha2 and bbeta respectively
#' @param Mstep the indicator of whether do model selection on the number of clusters
#' @param Vstep the indicator of whether do model selection on variance structures
#' @param SCind the indicator of whether use split/combine step in Mstep
#' @examples
#'
#'
#' library("fabMix")
#' library("mclust")
#' library("pgmm")
#' library("mvtnorm")
#' library("mcmcse")
#' library("MASS")
#' library("gtools")
#' n <- 500
#' p <- 10
#' q <- 4
#' K <- 10
#' nsim <- 10
#' burn <- 20
#' qnew <- 4
#' Mstep <- 1
#' Vstep <- 1
#' constraint <- c(0, 0, 0)
#' mInit <- 20
#' mVec <- c(1, 20)
#' X <- t(simData(
#'   sameLambda = TRUE,
#'   sameSigma = TRUE,
#'   K.true = K, n = n, q = q, p = p, sINV_values = 1 / ((1:p))
#' )$data)
#' \donttest{
#' pgmmRJMCMC(X,
#'   mInit, mVec, qnew,
#'   niter = nsim, burn = burn,
#'   constraint = constraint, Mstep = Mstep, Vstep = Vstep
#' )
#' }
#'
#' @export

pgmmRJMCMC <- function(X,
                       mInit,
                       mVec,
                       qnew,
                       delta = 2,
                       ggamma = 2,
                       burn = 20,
                       niter = 1000,
                       constraint = C(0, 0, 0),
                       dVec = c(1, 1, 1),
                       sVec = c(1, 1, 1),
                       Mstep = 0,
                       Vstep = 0,
                       SCind = 0) {
  n <- ncol(X)
  p <- nrow(X)

  alpha1 <- rgamma(1, dVec[1], sVec[1])
  alpha2 <- rgamma(1, dVec[2], sVec[2])
  bbeta <- rgamma(1, dVec[3], sVec[3])

  hparam <- new("Hparam", alpha1 = alpha1, alpha2 = alpha2, bbeta = bbeta, delta = delta, ggamma = ggamma)

  hparamInit <- hparam

  muBar <- X[, sample(1:n, 1)]

  ## cluster indicator
  clusInd <- rep(0, mVec[2])
  clusInd[1:mInit] <- 1

  ## qinit
  qVec <- rep(0, mVec[2])
  qVec[1:mInit] <- qnew


  ## priors
  ZOneDim <- kmeans(x = t(X), centers = mInit)$cluster
  thetaYList <- generatePriorThetaY(mInit, n, p, muBar, hparam, qVec, ZOneDim, constraint)

  ## burn in
  for (i in 1:burn) {
    MCMCobj <- stayMCMCupdate(X, thetaYList, ZOneDim, hparam, qVec, qnew, dVec, sVec, constraint, clusInd)
    ZOneDim <- MCMCobj$ZOneDim
    thetaYList <- MCMCobj$thetaYList
    hparam <- MCMCobj$hparam
    hparam@alpha2 <- max(0.01, hparam@alpha2)
  }

  thetaYList <- clearCurrentThetaYlist(thetaYList, clusInd, mVec[2])
  ##
  alpha1Vec <- c()
  alpha2Vec <- c()
  bbetaVec <- c()
  taoList <- list()
  psyList <- list()
  MList <- list()
  lambdaList <- list()
  YList <- list()
  ZmatList <- list()
  constraintList <- list()
  clusIndList <- list()
  ##

  for (h in 1:niter) {
    cat("iter = ", h, "======>\n")

    ## choose m or choose v
    if (Mstep == 1) {
      MCMCobj <- MstepRJMCMCupdate(X, muBar, p, thetaYList, ZOneDim, hparam, hparamInit, qVec, qnew, dVec, sVec, constraint, clusInd, mVec, "BD")
      ZOneDim <- MCMCobj$ZOneDim
      thetaYList <- MCMCobj$thetaYList
      hparam <- MCMCobj$hparam
      qVec <- MCMCobj$qVec
      clusInd <- MCMCobj$clusInd
      ##
      if (SCind == 1) {
        MCMCobj <- MstepRJMCMCupdate(X, muBar, p, thetaYList, ZOneDim, hparam, hparamInit, qVec, qnew, dVec, sVec, constraint, clusInd, mVec, "SC")
        ZOneDim <- MCMCobj$ZOneDim
        thetaYList <- MCMCobj$thetaYList
        hparam <- MCMCobj$hparam
        qVec <- MCMCobj$qVec
        clusInd <- MCMCobj$clusInd
      }
    }

    if (Vstep == 1) {
      MCMCobj <- VstepRJMCMCupdate(X, muBar, p, thetaYList, ZOneDim, hparam, hparamInit, qVec, qnew, ggamma, dVec, sVec, constraint, clusInd)
      ZOneDim <- MCMCobj$ZOneDim
      thetaYList <- MCMCobj$thetaYList
      hparam <- MCMCobj$hparam
      constraint <- MCMCobj$constraint
    }

    # stay step
    MCMCobj <- stayMCMCupdate(X, thetaYList, ZOneDim, hparam, qVec, qnew, dVec, sVec, constraint, clusInd)
    ZOneDim <- MCMCobj$ZOneDim
    thetaYList <- MCMCobj$thetaYList
    hparam <- MCMCobj$hparam
    hparam@alpha2 <- max(0.01, hparam@alpha2)

    ## save
    clusIndList[[h]] <- clusInd
    alpha1Vec[h] <- hparam@alpha1
    alpha2Vec[h] <- hparam@alpha2
    bbetaVec[h] <- hparam@bbeta
    taoList[[h]] <- thetaYList@tao
    psyList[[h]] <- thetaYList@psy
    MList[[h]] <- thetaYList@M
    lambdaList[[h]] <- thetaYList@lambda
    YList[[h]] <- thetaYList@Y
    ZmatList[[h]] <- ZOneDim
    constraintList[[h]] <- constraint
  }

  list(
    taoList = taoList, psyList = psyList, MList = MList, lambdaList = lambdaList,
    YList = YList, ZmatList = ZmatList, constraintList = constraintList,
    alpha1Vec = alpha1Vec, alpha2Vec = alpha2Vec, bbetaVec = bbetaVec,
    clusIndList = clusIndList
  )
}
