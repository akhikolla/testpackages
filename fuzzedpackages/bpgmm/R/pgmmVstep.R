#' VstepRJMCMCupdate
#'
#' @param X X
#' @param muBar muBar
#' @param p p
#' @param thetaYList thetaYList
#' @param ZOneDim ZOneDim
#' @param hparam hparam
#' @param hparamInit hparamInit
#' @param qVec qVec
#' @param qnew qnew
#' @param ggamma ggamma
#' @param dVec dVec
#' @param sVec sVec
#' @param constraint constraint
#' @param clusInd clusInd
#'
#' @examples
#' set.seed(100)
#' n <- 10
#' p <- 2
#' q <- 1
#' K <- 2
#' m <- 1
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
#' \donttest{
#' VstepRJMCMCupdate(
#'   X,
#'   muBar,
#'   p,
#'   thetaYList,
#'   ZOneDim,
#'   hparam,
#'   hparamInit,
#'   qVec,
#'   qnew,
#'   ggamma,
#'   dVec,
#'   sVec,
#'   constraint,
#'   clusInd
#' )
#' }
VstepRJMCMCupdate <- function(X,
                              muBar,
                              p,
                              thetaYList,
                              ZOneDim,
                              hparam,
                              hparamInit,
                              qVec,
                              qnew,
                              ggamma,
                              dVec,
                              sVec,
                              constraint,
                              clusInd) {



  ## rjmcmc for V
  # print("vStep")
  n <- ncol(X)
  m <- sum(clusInd)
  # steps = c("stay", "lambda", "psi1", "psi2")
  steps <- c("lambda", "psi1", "psi2")
  currentStep <- sample(size = 1, x = steps)
  # currentStep = "psi2"
  # currentStep = "stay"


  ## trans to non empty obj
  NEobj <- toNEthetaYlist(thetaYList, ZOneDim, qVec, clusInd)
  thetaYList <- NEobj$thetaYList
  ZOneDim <- NEobj$ZOneDim
  qVec <- NEobj$qVec
  # if(currentStep == "stay"){
  #   print("stay")
  #   ZOneDim = update_PostZ(X, m, n, thetaYList)
  #   # ZOneDim = updatePostZ(m, n, thetaYList)
  #   thetaYList = updatePostThetaY(m,  n,hparam, thetaYList, ZOneDim,  qVec,constraint, X, ggamma)
  #   hparam = update_Hyperparameter(m, p, qnew, hparam,thetaYList, dVec, sVec)
  #
  # }
  if (currentStep == "lambda") {

    # print("lambda step")
    proposeConstraint <- constraint
    proposeConstraint[1] <- (constraint[1] + 1) %% 2

    CxyList <- Calculate_Cxy(m, n, hparam, thetaYList, ZOneDim, qVec, X)

    newLambda <- CalculateProposalLambda(hparam, thetaYList, CxyList, proposeConstraint, m, p, qVec)

    newthetaYList <- thetaYList
    newthetaYList@lambda <- newLambda

    newCxyList <- Calculate_Cxy(m, n, hparam, newthetaYList, ZOneDim, qVec, X)

    oldDensity <- likelihood(thetaYList, ZOneDim, qVec, muBar, X) + evaluatePrior(m, p, muBar, hparam, thetaYList, ZOneDim, qVec, constraint, sort(clusInd, decreasing = T))
    oldLambdaEval <- EvaluateProposalLambda(hparam, newthetaYList, newCxyList, constraint, thetaYList@lambda, m, qVec, p)
    newLambdaEval <- EvaluateProposalLambda(hparam, thetaYList, CxyList, proposeConstraint, newthetaYList@lambda, m, qVec, p)




    ## Gibbs
    newhparam <- hparam
    for (i in 1:10) {
      newZOneDim <- update_PostZ(X, m, n, newthetaYList)
      # newZOneDim = updatePostZ(m, n, newthetaYList)

      newthetaYList <- updatePostThetaY(m = m, n = n, p, newhparam, newthetaYList, ZOneDim = newZOneDim, qVec = qVec, constraint = proposeConstraint, X, ggamma)
      newhparam <- update_Hyperparameter(m, p, qnew, newhparam, newthetaYList, dVec, sVec)
    }

    newDensity <- likelihood(newthetaYList, newZOneDim, qVec, muBar, X) + evaluatePrior(m, p, muBar, newhparam, newthetaYList, newZOneDim, qVec, proposeConstraint, sort(clusInd, decreasing = T))

    numer <- newDensity + oldLambdaEval
    denom <- oldDensity + newLambdaEval
    # diffVec[h] = numer - denom
    # cat("accepting prob = exp ", numer - denom, " \n")

    probAlpha <- calculateRatio(numer, denom)
    acceptP <- min(1, probAlpha)
    res <- rbinom(1, size = 1, prob = acceptP)

    if (res == 1) {
      # print("lambda success=====>")
      # cat("lambda prob = ", probAlpha, "====>\n")

      constraint <- proposeConstraint
      thetaYList <- newthetaYList
      ZOneDim <- newZOneDim
      hparam <- newhparam
    } else {
      # print("lambda fail=====>")
    }
  }
  else if (currentStep == "psi1" | currentStep == "psi2") {

    # print(currentStep)
    proposeConstraint <- constraint
    if (currentStep == "psi1") {
      proposeConstraint[2] <- (constraint[2] + 1) %% 2
    } else {
      proposeConstraint[3] <- (constraint[3] + 1) %% 2
    }

    CxyList <- Calculate_Cxy(m, n, hparam, thetaYList, ZOneDim, qVec, X)
    newPsy <- CalculateProposalPsy(hparam, thetaYList, CxyList, proposeConstraint, m, p, qVec)

    newthetaYList <- thetaYList
    newthetaYList@psy <- newPsy
    newCxyList <- Calculate_Cxy(m, n, hparam, newthetaYList, ZOneDim, qVec, X)


    # oldDensity = likelihood(thetaYList, ZOneDim,qVec,muBar, X) + evaluatePrior(m, p, muBar, hparam, thetaYList, ZOneDim, qVec, constraint,clusInd)
    # newDensity = likelihood(newthetaYList, ZOneDim,qVec,muBar, X) + evaluatePrior(m, p, muBar,hparam, newthetaYList, ZOneDim, qVec, proposeConstraint,clusInd)

    oldDensity <- likelihood(thetaYList, ZOneDim, qVec, muBar, X) + evaluatePrior(m, p, muBar, hparam, thetaYList, ZOneDim, qVec, constraint, sort(clusInd, decreasing = T))
    newDensity <- likelihood(newthetaYList, ZOneDim, qVec, muBar, X) + evaluatePrior(m, p, muBar, hparam, newthetaYList, ZOneDim, qVec, proposeConstraint, sort(clusInd, decreasing = T))


    oldPsyEval <- EvaluateProposalPsy(hparam, newthetaYList, newCxyList, constraint, thetaYList@psy, m, p, qVec)
    newPsyEval <- EvaluateProposalPsy(hparam, thetaYList, CxyList, proposeConstraint, newthetaYList@psy, m, p, qVec)

    numer <- newDensity + oldPsyEval
    denom <- oldDensity + newPsyEval
    # print(numer - denom)
    # diffVec[h] = numer - denom

    probAlpha <- calculateRatio(numer, denom)
    acceptP <- min(1, probAlpha)
    res <- rbinom(1, size = 1, prob = acceptP)

    if (res == 1) {
      # print("psi success=====>")
      # cat("psi prob = ", probAlpha, "====>\n")
      constraint <- proposeConstraint
      thetaYList <- newthetaYList
      # psy = newPsy
    } else {

      # print("psi fail=====>")
    }
  }

  ## trans to empty obj
  Eobj <- toEthetaYlist(thetaYList, ZOneDim, qnew, clusInd)
  thetaYList <- Eobj$thetaYList
  ZOneDim <- Eobj$ZOneDim
  qVec <- Eobj$qVec



  return(list(thetaYList = thetaYList, ZOneDim = ZOneDim, hparam = hparam, constraint = constraint))
}
