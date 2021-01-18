
#' proposeCombinedClusters
#'
#' @param X X
#' @param thetaYList thetaYList
#' @param hparam hparam
#' @param combineClusInd combineClusInd
#' @param combinedClusInd combinedClusInd
#' @param qVec qVec
#' @param constraint constraint
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
#' combineClusInd <- c(1, 1)
#' combinedClusInd <- 1
#' proposeCombinedClusters(X, thetaYList, hparam, combineClusInd, combinedClusInd, qVec, constraint)
proposeCombinedClusters <- function(X, thetaYList, hparam, combineClusInd, combinedClusInd, qVec, constraint) {
  p <- nrow(X)
  n <- ncol(X)

  ggamma <- hparam@ggamma
  delta <- hparam@delta
  bbeta <- hparam@bbeta
  alpha1 <- hparam@alpha1
  alpha2 <- hparam@alpha2


  resThetaYList <- thetaYList

  ## propose w
  w1 <- thetaYList@tao[combineClusInd[1]]
  w2 <- thetaYList@tao[combineClusInd[2]]
  w <- w1 + w2

  if (w > 1) {
    w <- 1
  }

  resThetaYList@tao[combinedClusInd] <- w


  ## propose mu

  mu1 <- thetaYList@M[[combineClusInd[1]]]
  mu2 <- thetaYList@M[[combineClusInd[2]]]

  mu <- w1 / w * mu1 + w2 / w * mu2

  # mu =   (mu1 + mu2)/2

  resThetaYList@M[[combinedClusInd]] <- mu

  ## propose psy

  ## propose lambda and psy
  if (constraint[2] == 1) {
    psy <- thetaYList@psy[combineClusInd[1]]
  } else {
    psy <- generatePriorPsi(p, 1, delta, bbeta, constraint)
  }

  resThetaYList@psy[[combinedClusInd]] <- psy[[1]]

  if (constraint[1] == 1) {
    lambda <- thetaYList@lambda[combineClusInd[1]]
  } else {
    lambda <- generatePriorLambda(p, 1, alpha2, qVec, psy, constraint)
  }

  resThetaYList@lambda[combinedClusInd] <- lambda


  ## proposed Y
  for (k in combinedClusInd) {
    resThetaYList@Y[[k]] <- matrix(NA, qVec[k], n)
    for (i in 1:n) {
      resThetaYList@Y[[k]][, i] <- mvtnorm::rmvnorm(1, mean = rep(0, qVec[k]), sigma = diag(qVec[k]))
    }
  }



  return(resThetaYList)
}

#' evaluateCombinedClusters
#'
#' @param X X
#' @param combinedThetaYList combinedThetaYList
#' @param hparam hparam
#' @param combinedClusInd combinedClusInd
#' @param qVec qVec
#' @param constraint constraint
#'
#' @return
#' @export
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
#' combinedThetaYList <-
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
#' combinedClusInd <- 1
#' evaluateCombinedClusters(X, combinedThetaYList, hparam, combinedClusInd, qVec, constraint)
evaluateCombinedClusters <- function(X, combinedThetaYList, hparam, combinedClusInd, qVec, constraint) {
  p <- nrow(X)
  n <- ncol(X)

  ggamma <- hparam@ggamma
  delta <- hparam@delta
  bbeta <- hparam@bbeta
  alpha1 <- hparam@alpha1
  alpha2 <- hparam@alpha2


  ## eval  psy
  psy <- combinedThetaYList@psy[combinedClusInd]
  evalpsy <- 0
  if (constraint[2] == 0) {
    evalpsy <- evaluatePriorPsi(psy, p, 1, delta, bbeta, constraint, 1)
  }

  ## eval lambda
  lambda <- combinedThetaYList@lambda[combinedClusInd]
  evallambda <- 0
  if (constraint[1] == 0) {
    evallambda <- evaluatePriorLambda(p, 1, alpha2, qVec, psy, lambda, constraint, 1)
  }


  return(sum(evalpsy, evallambda))
}




#' combineZOneDim
#'
#' @param ZOneDim ZOneDim
#' @param combineClusInd combineClusInd
#' @param combinedClusInd combinedClusInd
#'
#' @examples
#' ZOneDim <- c(1, 1, 1, 2, 3, 4)
#' combineClusInd <- c(1, 2)
#' combinedClusInd <- 1
#' combineZOneDim(ZOneDim, combineClusInd, combinedClusInd)
combineZOneDim <- function(ZOneDim, combineClusInd, combinedClusInd) {
  ZOneDimNew <- ZOneDim
  ZOneDimNew[ZOneDimNew == combineClusInd[1]] <- combinedClusInd
  ZOneDimNew[ZOneDimNew == combineClusInd[2]] <- combinedClusInd
  return(ZOneDimNew)
}
