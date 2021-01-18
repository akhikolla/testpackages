#' generate a new extra cluster Theta and Y
#'
#' @param m m
#' @param n n
#' @param p p
#' @param hparam hparam
#' @param muBar muBar
#' @param qnew qnew
#' @param thetaYList thetaYList
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
#' clusInd <- rep(1, m)
#' generateNewClusPara(m, n, p, hparam, muBar, qnew, thetaYList, constraint, clusInd)
generateNewClusPara <- function(m, n, p, hparam, muBar, qnew, thetaYList, constraint, clusInd) {
  ggamma <- hparam@ggamma
  delta <- hparam@delta
  bbeta <- hparam@bbeta
  alpha1 <- hparam@alpha1
  alpha2 <- hparam@alpha2


  lambdaPsyList <- getLambdaPsiForNewClus(m, p, qnew, thetaYList, hparam, constraint, clusInd)

  lambda <- lambdaPsyList$lambda
  psy <- lambdaPsyList$psy

  ## proposal tao
  tao <- rbeta(1, shape1 = 1, shape2 = m)
  # tao = rbeta(1, shape1 = 1, shape2 = 100)
  # taoDval = log(dbeta(tao, shape1 = 1, shape2 = m))

  ## proposal M
  Mvar <- 1 / (alpha1) * psy ## prior
  M <- t(c(mvtnorm::rmvnorm(1, muBar, Mvar)))
  # MDval = mvtnorm::dmvnorm(M[[1]], muBar, Mvar, log = T)


  ## proposal Y
  # YvecDval = c()
  Y <- matrix(NA, qnew, n)
  for (i in 1:n) {
    Y[, i] <- mvtnorm::rmvnorm(1, mean = rep(0, qnew), sigma = diag(qnew))
    # YvecDval[i] = mvtnorm::dmvnorm(Y[[1]][,i], mean = rep(0,qNew), sigma = diag(qNew), log = T)
  }
  # YDval = sum(YvecDval)
  # ThetaYGenDval = sum(taoDval, psyDval, MDval, lambdaDval, YDval)
  # ThetaYGenDval = sum(psyDval, MDval, lambdaDval, YDval)


  new("ThetaYList",
    tao = tao,
    psy = list(psy),
    M = list(M),
    lambda = list(lambda),
    Y = list(Y)
  )
}


#' getLambdaPsiForNewClus
#'
#' @param m m
#' @param p p
#' @param qnew qnew
#' @param thetaYList thetaYList
#' @param hparam hparam
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
#' clusInd <- rep(1, m)
#' getLambdaPsiForNewClus(
#'   m,
#'   p,
#'   qnew,
#'   thetaYList,
#'   hparam,
#'   constraint,
#'   clusInd
#' )
getLambdaPsiForNewClus <- function(m,
                                   p,
                                   qnew,
                                   thetaYList,
                                   hparam,
                                   constraint,
                                   clusInd) {
  ggamma <- hparam@ggamma
  delta <- hparam@delta
  bbeta <- hparam@bbeta
  alpha1 <- hparam@alpha1
  alpha2 <- hparam@alpha2

  ## psy
  psy <- c()
  if (constraint[2] == T) {
    psy <- thetaYList@psy[[which(clusInd == 1)[1]]]
  } else if (constraint[2] == F & constraint[3] == T) {
    psyValue <- 1 / rgamma(1, shape = hparam@delta, rate = hparam@bbeta)
    psy <- diag(rep(psyValue, p))
  } else if (constraint[2] == F & constraint[3] == F) {
    psy <- diag(1 / rgamma(p, shape = hparam@delta, rate = hparam@bbeta))
  }


  ## lambda
  lambda <- c()
  if (constraint[1] == T) {
    lambda <- thetaYList@lambda[[which(clusInd == 1)[1]]]
  } else if (constraint[1] == F) {
    lambda <- matrix(0, p, qnew)
    for (j in 1:qnew) {
      lambda[, j] <- mvtnorm::rmvnorm(1, rep(0, p), 1 / alpha2 * psy)
    }
  }
  return(list(psy = psy, lambda = lambda))
}


## evaluate a new extra cluster Theta and Y
#'
#' @param m m
#' @param p p
#' @param hparam hparam
#' @param muBar muBar
#' @param newClusterThetaY newClusterThetaY
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
#' clusInd <- rep(1, m)
#' \donttest{
#' evaluateNewClusPara(m, p, hparam, muBar, thetaYList, constraint)
#' }
#'
evaluateNewClusPara <- function(m, p, hparam, muBar, newClusterThetaY, constraint) {
  ggamma <- hparam@ggamma
  delta <- hparam@delta
  bbeta <- hparam@bbeta
  alpha1 <- hparam@alpha1
  alpha2 <- hparam@alpha2

  ## eval lambda psy
  lambdaPsyEval <- evalLambdaPsiForNewClus(m, p, qnew, thetaYList, newClusterThetaY, hparam, constraint)

  ## eval tao
  tao <- newClusterThetaY@tao
  taoDval <- log(dbeta(tao, shape1 = 1, shape2 = m))

  # taoDval = log(dbeta(tao, shape1 = 1, shape2 = 100))

  # hist(rbeta(10000, shape1 = 1, shape2 = 100))
  ## test
  # taoDval = 0
  ##

  ## eval M
  M <- newClusterThetaY@M[[1]]
  MDval <- mvtnorm::dmvnorm(M, muBar, 1 / alpha1 * newClusterThetaY@psy[[1]], log = T)

  return(sum(taoDval, MDval, lambdaPsyEval$lambdaEval, lambdaPsyEval$psyEval))
}




#' evalLambdaPsiForNewClus
#'
#' @param m m
#' @param p p
#' @param qnew qnew
#' @param thetaYList thetaYList
#' @param newClusterThetaY newClusterThetaY
#' @param hparam hparam
#' @param constraint constraint
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
#' newClusterThetaY <- thetaYList <-
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
#' clusInd <- rep(1, m)
#' evalLambdaPsiForNewClus(
#'   m,
#'   p,
#'   qnew,
#'   thetaYList,
#'   newClusterThetaY,
#'   hparam,
#'   constraint
#' )
evalLambdaPsiForNewClus <- function(m,
                                    p,
                                    qnew,
                                    thetaYList,
                                    newClusterThetaY,
                                    hparam,
                                    constraint) {
  ggamma <- hparam@ggamma
  delta <- hparam@delta
  bbeta <- hparam@bbeta
  alpha1 <- hparam@alpha1
  alpha2 <- hparam@alpha2

  ## psy
  psy <- newClusterThetaY@psy[[1]]
  psyEval <- 0
  if (constraint[2] == T) {
    ## psyVal = 0 in this case
  } else if (constraint[2] == F & constraint[3] == T) {
    psyEval <- dgamma((diag(1 / psy)[1]), shape = hparam@delta, rate = hparam@bbeta, log = T)
  } else if (constraint[2] == F & constraint[3] == F) {
    psyEval <- sum(dgamma(1 / diag(psy), shape = hparam@delta, rate = hparam@bbeta, log = T))
  }


  ## lambda
  lambda <- newClusterThetaY@lambda[[1]]
  lambdaEval <- 0
  if (constraint[1] == T) {
    ## lambdaEval = 0
  } else if (constraint[1] == F) {
    for (j in 1:qnew) {
      lambdaEval[j] <- mvtnorm::dmvnorm(lambda[, j], rep(0, p), 1 / alpha2 * psy, log = T)
    }
    lambdaEval <- sum(lambdaEval)
  }
  return(list(psyEval = psyEval, lambdaEval = lambdaEval))
}
