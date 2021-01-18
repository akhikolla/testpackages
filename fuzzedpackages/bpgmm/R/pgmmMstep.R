#' MstepRJMCMCupdate
#'
#' @param X X in MstepRJMCMCupdate
#' @param muBar muBar
#' @param p p
#' @param thetaYList thetaYList
#' @param ZOneDim ZOneDim
#' @param hparam hparam
#' @param hparamInit hparamInit
#' @param qVec qVec
#' @param qnew qnew
#' @param dVec dVec
#' @param mVec mVec
#' @param sVec sVec
#' @param constraint constraint
#' @param clusInd clusInd
#' @param Mind Mind
#' @examples
#' set.seed(100)
#' n <- 10
#' p <- 2
#' q <- 1
#' K <- 2
#' m <- 2
#' muBar <- c(0, 0)
#' qVec <- c(1, 1)
#' constraint <- c(0, 0, 0)
#' X <- t(
#'   fabMix::simData(
#'     sameLambda = TRUE,
#'     sameSigma = TRUE,
#'     K.true = K,
#'     n = n,
#'     q = q,
#'     p = p,
#'     sINV_values = 1 / ((1:p))
#'   )$data
#' )
#' hparam <- new(
#'   "Hparam",
#'   alpha1 = 0.567755037123148,
#'   alpha2 = 1.1870201935945,
#'   delta = 2,
#'   ggamma = 2,
#'   bbeta = 3.39466184520673
#' )
#' ZOneDim <- sample(seq_len(m), n, replace = TRUE)
#' thetaYList <-
#'   new(
#'     "ThetaYList",
#'     tao = 0.366618687752634,
#'     psy = list(structure(
#'       c(
#'         4.18375613018654,
#'         0, 0, 5.46215996830771
#'       ),
#'       .Dim = c(2L, 2L)
#'     )),
#'     M = list(structure(
#'       c(
#'         3.27412045866392,
#'         -2.40544145363349
#'       ),
#'       .Dim = 1:2
#'     )),
#'     lambda = list(structure(
#'       c(
#'         2.51015961514781,
#'         -0.0741189919182549
#'       ),
#'       .Dim = 2:1
#'     )),
#'     Y = list(structure(
#'       c(
#'         -0.244239011725104,
#'         -0.26876172736886,
#'         0.193431511203083,
#'         0.41624466812811,
#'         -0.54581548068437,
#'         -0.0479517628308146,
#'         -0.633383997203325,
#'         0.856855296613208,
#'         0.792850576988512,
#'         0.268208848994559
#'       ),
#'       .Dim = c(1L, 10L)
#'     ))
#'   )
#' qnew <- 1
#' dVec <- c(1, 1, 1)
#' sVec <- c(1, 1, 1)
#' constraint <- c(0, 0, 0)
#' clusInd <- rep(1, m)
#' Mind <- "BD"
#' mVec <- c(1, m)
#' \donttest{
#' MstepRJMCMCupdate(
#'   X,
#'   muBar,
#'   p,
#'   thetaYList,
#'   ZOneDim,
#'   hparam,
#'   hparamInit,
#'   qVec,
#'   qnew,
#'   dVec,
#'   sVec,
#'   constraint,
#'   clusInd,
#'   mVec,
#'   Mind
#' )
#' }
MstepRJMCMCupdate <- function(X,
                              muBar,
                              p,
                              thetaYList,
                              ZOneDim,
                              hparam,
                              hparamInit,
                              qVec,
                              qnew,
                              dVec,
                              sVec,
                              constraint,
                              clusInd,
                              mVec,
                              Mind) {


  # print("mStep")
  ## rjmcmc for m

  m <- sum(clusInd)
  n <- ncol(X)
  # empty cluster
  emptyZ <- setdiff(which(clusInd == 1), ZOneDim)

  currentStep <- c()

  ## birth and death
  if (Mind == "BD") {
    if (length(emptyZ) > 1) {
      currentStep <- c("death")
    }
    else if (m == mVec[1] | length(unique(ZOneDim)) == 1) {
      currentStep <- c("birth")
    }
    else if (m == mVec[2]) {
      if (length(emptyZ) > 0) {
        currentStep <- c("death")
      } else {
        currentStep <- c("NA")
      }
    } else {
      if (length(emptyZ) > 0) {
        currentStep <- c("birth", "death")
      } else {
        currentStep <- c("birth")
      }
    }
  }

  if (Mind == "SC") {
    if (m == mVec[1] | length(unique(ZOneDim)) == 1) {
      currentStep <- c("split")
    } else if (m == mVec[2]) {
      currentStep <- c("combine")
    } else {
      currentStep <- c("combine", "split")
    }
  }


  indChar <- sample(currentStep, 1)
  # cat("======",indChar, "======>\n")

  if (indChar == "birth") {
    birthInd <- which(clusInd == 0)[1]

    newClusterThetaY <- generateNewClusPara(m, n, p, hparam, muBar, qnew, thetaYList, constraint, clusInd)

    ## adjust tao in the combine
    combinedThetaList <- combineClusterPara(thetaYList, newClusterThetaY, birthInd)
    jacobianAdj <- m * log(1 - newClusterThetaY@tao)

    ## adjust q
    qVecNew <- qVec
    qVecNew[birthInd] <- qnew


    ## adjust clusInd
    clusIndNew <- clusInd
    clusIndNew[birthInd] <- 1


    oldDensity <- likelihood(thetaYList, ZOneDim, qVec, muBar, X) + evaluatePrior(m, p, muBar, hparam, thetaYList, ZOneDim, qVec, constraint, clusInd)
    newDensity <- likelihood(combinedThetaList, ZOneDim, qVecNew, muBar, X) + evaluatePrior(m + 1, p, muBar, hparam, combinedThetaList, ZOneDim, qVecNew, constraint, clusIndNew)

    newClusterThetaYEval <- evaluateNewClusPara(m, p, hparam, muBar, newClusterThetaY, constraint)

    emptyZ <- setdiff(which(clusIndNew == 1), ZOneDim)

    numeratorAlpha <- newDensity + jacobianAdj
    denominatorAlpha <- oldDensity + newClusterThetaYEval + log(length(emptyZ))

    probAlpha <- calculateRatio(numeratorAlpha, denominatorAlpha)
    acceptP <- min(1, probAlpha)

    res <- rbinom(1, size = 1, prob = acceptP)
    # cat("birth prob = ", probAlpha, "====>\n")


    if (res == 1) {
      ## birth
      # print("birth succeed=====>")
      # cat("birth prob = ", probAlpha, "====>\n")

      thetaYList <- combinedThetaList
      qVec <- qVecNew
      clusInd <- clusIndNew
      thetaYList <- clearCurrentThetaYlist(thetaYList, clusInd, mVec[2])
    } else {
      ## stay
      # print("birth reject =====>")
    }
  }
  else if (indChar == "death") {
    ## no update

    emptyZ <- setdiff(which(clusInd == 1), ZOneDim)
    chooseZprob <- 0
    if (length(emptyZ) == 1) {
      deathInd <- emptyZ
    } else {
      deathInd <- sample(x = emptyZ, size = 1)
      chooseZprob <- log(1 / length(emptyZ))
    }

    ## adj qVec and clusInd
    qVecNew <- qVec
    clusIndNew <- clusInd
    qVecNew[deathInd] <- 0
    clusIndNew[deathInd] <- 0


    dropCluster <- getIndThetaY(thetaYList, deathInd)
    dropClusterDval <- evaluateNewClusPara(m - 1, p, hparam, muBar, dropCluster, constraint)

    ## adjust tao

    leftThetaY <- getRemovedIndThetaY(thetaYList, deathInd)

    jacobianAdj <- (m - 1) * log(1 - dropCluster@tao)


    oldDensity <- evaluatePrior(m, p, muBar, hparam, thetaYList, ZOneDim, qVec, constraint, clusInd)
    newDensity <- evaluatePrior(m - 1, p, muBar, hparam, leftThetaY, ZOneDim, qVecNew, constraint, clusIndNew)


    numeratorAlpha <- newDensity + dropClusterDval + chooseZprob
    denominatorAlpha <- oldDensity + jacobianAdj

    # cat("diff  = ", numeratorAlpha - denominatorAlpha, "=====> \n" )

    probAlpha <- calculateRatio(numeratorAlpha, denominatorAlpha)
    acceptP <- min(1, probAlpha)

    res <- rbinom(1, size = 1, prob = acceptP)

    if (res == 1) {
      ## contract
      # print("death succeed=====>")
      # cat("death prob = ", probAlpha, "====>\n")

      # stop("=====>")

      thetaYList <- leftThetaY
      qVec <- qVecNew
      clusInd <- clusIndNew
      thetaYList <- clearCurrentThetaYlist(thetaYList, clusInd, mVec[2])
    } else {
      # print("death stay =====>")
    }
  }

  else if (indChar == "combine") {
    ## choose two  positions to combine
    combineClusInd <- sample(x = sort(unique(ZOneDim)), size = 2, prob = 1 / table(ZOneDim))

    logProbCombineZ <- log(1 / choose(length(unique(ZOneDim)), 2))

    clusIndNew <- clusInd
    qVecNew <- qVec
    clusIndNew[combineClusInd] <- 0
    qVecNew[combineClusInd] <- 0

    ## choose one position to fill combined clus
    combinedClusInd <- which(clusIndNew == 0)[1]
    clusIndNew[combinedClusInd] <- 1
    qVecNew[combinedClusInd] <- qnew

    ## propose values for new clusters
    combinedThetaYList <- proposeCombinedClusters(X, thetaYList, hparam, combineClusInd, combinedClusInd, qVecNew, constraint)
    combinedThetaYList <- clearCurrentThetaYlist(combinedThetaYList, clusIndNew, mVec[2])
    combineEval <- evaluateCombinedClusters(X, combinedThetaYList, hparam, combinedClusInd, qVecNew, constraint)
    ## evaluate values for old cluster
    splitEval <- evaluateSplitedClusters(X, combinedThetaYList, thetaYList, hparam, combinedClusInd, combineClusInd, qVec, constraint)


    ## calculate Jacobian

    JacEval <- calculateJacobian(p, combinedThetaYList, thetaYList, combinedClusInd, combineClusInd)


    ## combine Z

    ZOneDimNew <- combineZOneDim(ZOneDim, combineClusInd, combinedClusInd)

    ## eval split

    evalZProb <- splitEvalZOneDim(ZOneDimNew, ZOneDim, thetaYList, combinedClusInd, combineClusInd)

    # update Y
    combinedThetaYList@Y <- updateY(X, combinedThetaYList, ZOneDimNew, clusIndNew, qVecNew)



    combinedThetaYListTEST <- combinedThetaYList
    ZOneDimNewTEST <- ZOneDimNew
    for (i in 1:5) {
      MCMCobj <- stayMCMCupdate(X, combinedThetaYListTEST, ZOneDimNewTEST, hparam, qVecNew, qnew, dVec, sVec, constraint, clusIndNew)
      ZOneDimNewTEST <- MCMCobj$ZOneDim
      combinedThetaYListTEST <- MCMCobj$thetaYList
      # hparam = MCMCobj$hparam
    }


    oldDensity <- likelihood(thetaYList, ZOneDim, qVec, muBar, X) + evaluatePrior(m, p, muBar, hparam, thetaYList, ZOneDim, qVec, constraint, clusInd)
    newDensity <- likelihood(combinedThetaYListTEST, ZOneDimNewTEST, qVecNew, muBar, X) + evaluatePrior(m - 1, p, muBar, hparam, combinedThetaYListTEST, ZOneDimNewTEST, qVecNew, constraint, clusIndNew)

    numer <- newDensity + splitEval + JacEval
    denom <- oldDensity + combineEval + logProbCombineZ

    # cat("combine diff = ", newDensity - oldDensity, "====>\n")
    # cat("combine diff = ", numer - denom, "====>\n")

    probAlpha <- calculateRatio(numer, denom)
    acceptP <- min(1, probAlpha)
    res <- rbinom(1, size = 1, prob = acceptP)

    if (res == 1) {
      # print("combine success=====>")

      thetaYList <- combinedThetaYListTEST
      ZOneDim <- ZOneDimNewTEST
      clusInd <- clusIndNew
      qVec <- qVecNew
      thetaYList <- clearCurrentThetaYlist(thetaYList, clusInd, mVec[2])
    } else {
      # print("combine fail=====>")
    }
  }

  else if (indChar == "split") {
    ## choose one existing position to split
    if (length(unique(ZOneDim)) == 1) {
      splitClusInd <- unique(ZOneDim)
    } else {
      splitClusInd <- sample(sort(unique(ZOneDim)), size = 1, prob = table(ZOneDim))
    }
    logProbSplitZ <- log(1 / length(unique(ZOneDim)))

    clusIndNew <- clusInd
    qVecNew <- qVec
    clusIndNew[splitClusInd] <- 0
    qVecNew[splitClusInd] <- 0
    ## choose two positions to fill the two clusters
    splitedClusInd <- which(clusIndNew == 0)[1:2]
    clusIndNew[splitedClusInd] <- 1
    qVecNew[splitedClusInd] <- qnew
    ## propose values for new clusters

    splitedThetaYList <- proposeSplitedClusters(X, thetaYList, hparam, splitClusInd, splitedClusInd, qVecNew, constraint)
    splitEval <- evaluateSplitedClusters(X, thetaYList, splitedThetaYList, hparam, splitClusInd, splitedClusInd, qVecNew, constraint)

    ## calculate Jacobian

    JacEval <- calculateJacobian(p, thetaYList, splitedThetaYList, splitClusInd, splitedClusInd)

    ## evaluate values for old cluster
    combineEval <- evaluateCombinedClusters(X, thetaYList, hparam, splitClusInd, qVec, constraint)

    ## split Z

    splitZObj <- splitZOneDim(ZOneDim, splitedThetaYList, splitClusInd, splitedClusInd)

    ZOneDimNew <- splitZObj$ZOneDim
    evalZProb <- splitZObj$evalProb

    # update Y
    splitedThetaYList@Y <- updateY(X, splitedThetaYList, ZOneDimNew, clusIndNew, qVecNew)

    splitedThetaYListTEST <- splitedThetaYList
    ZOneDimNewTEST <- ZOneDimNew
    for (i in 1:5) {
      MCMCobj <- stayMCMCupdate(X, splitedThetaYListTEST, ZOneDimNewTEST, hparam, qVecNew, qnew, dVec, sVec, constraint, clusIndNew)
      ZOneDimNewTEST <- MCMCobj$ZOneDim
      splitedThetaYListTEST <- MCMCobj$thetaYList
    }

    ## calculate acceptancy propability

    oldDensity <- likelihood(thetaYList, ZOneDim, qVec, muBar, X) + evaluatePrior(m, p, muBar, hparam, thetaYList, ZOneDim, qVec, constraint, clusInd)
    newDensity <- likelihood(splitedThetaYListTEST, ZOneDimNewTEST, qVecNew, muBar, X) + evaluatePrior(m + 1, p, muBar, hparam, splitedThetaYListTEST, ZOneDimNewTEST, qVecNew, constraint, clusIndNew)

    numer <- newDensity + combineEval
    denom <- oldDensity + splitEval + logProbSplitZ + JacEval

    # cat("diff  = ", numer - denom, "====>\n")
    probAlpha <- calculateRatio(numer, denom)
    acceptP <- min(1, probAlpha)
    res <- rbinom(1, size = 1, prob = acceptP)

    if (res == 1) {
      # print("split success=====>")
      # cat("split prob = ", probAlpha, "====>\n")


      thetaYList <- splitedThetaYListTEST
      ZOneDim <- ZOneDimNewTEST
      clusInd <- clusIndNew
      qVec <- qVecNew
      thetaYList <- clearCurrentThetaYlist(thetaYList, clusInd, mVec[2])
    } else {
      # print("split fail=====>")
    }
  }

  return(list(thetaYList = thetaYList, ZOneDim = ZOneDim, hparam = hparam, qVec = qVec, clusInd = clusInd))
}

utils::globalVariables(c(
  "mVec", "generateNewClusPara", "evaluateNewClusPara",
  "proposeCombinedClusters", "evaluateCombinedClusters",
  "evaluateSplitedClusters", "calculateJacobian",
  "combineZOneDim", "splitEvalZOneDim", "updateY",
  "proposeSplitedClusters", "splitZOneDim", "res"
))
