
#' proposeSplitedClusters
#'
#' @param X X
#' @param thetaYList thetaYList
#' @param hparam hparam
#' @param splitClusInd splitClusInd
#' @param splitedClusInd splitedClusInd
#' @param qVec qVec
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
#' splitClusInd <- 1
#' splitedClusInd <- c(1, 2)
#' proposeSplitedClusters(X, thetaYList, hparam, splitClusInd, splitedClusInd, qVec, constraint)
proposeSplitedClusters <- function(X, thetaYList, hparam, splitClusInd, splitedClusInd, qVec, constraint) {
  p <- nrow(X)
  n <- ncol(X)

  ggamma <- hparam@ggamma
  delta <- hparam@delta
  bbeta <- hparam@bbeta
  alpha1 <- hparam@alpha1
  alpha2 <- hparam@alpha2


  resThetaYList <- thetaYList

  ## beta(2,2) is symatric between (0,1)
  a1 <- rbeta(1, 2, 2)
  # a2 = rbeta(p, 2, 2)
  a2 <- rgamma(p, 1, 2)


  # print(a1)
  # print("===========propose>")
  # print(a2)
  # print("===========propose>")
  ## propose weight
  w1 <- thetaYList@tao[splitClusInd] * a1
  w2 <- thetaYList@tao[splitClusInd] * (1 - a1)

  resThetaYList@tao[splitedClusInd] <- c(w1, w2)

  ## propose mean
  splitMu <- thetaYList@M[[splitClusInd]]
  ## calculate split clus Variance
  splitClusVar <- thetaYList@psy[[splitClusInd]] + thetaYList@lambda[[splitClusInd]] %*% t(thetaYList@lambda[[splitClusInd]])
  splitClusSDvec <- sqrt(diag(splitClusVar))
  # print("=======")
  # print(splitClusSDvec)
  mu1 <- c()
  mu2 <- c()

  for (i in 1:p) {
    signInd <- sample(c(-1, 1), 1)
    mu1[i] <- splitMu[i] + signInd * a2[i] * splitClusSDvec[i] * sqrt(w2 / w1)
    mu2[i] <- splitMu[i] + -1 * signInd * a2[i] * splitClusSDvec[i] * sqrt(w1 / w2)
  }

  # print(mu1)
  # print(mu2)
  # print(splitMu)
  # print("=======>")

  resThetaYList@M[[splitedClusInd[1]]] <- t(mu1)
  resThetaYList@M[[splitedClusInd[2]]] <- t(mu2)

  ## propose lambda and psy
  if (constraint[2] == 1) {
    psy1 <- thetaYList@psy[[splitClusInd]]
    psy2 <- thetaYList@psy[[splitClusInd]]
    psy <- list(psy1, psy2)
  } else {
    psy <- generatePriorPsi(p, 2, delta, bbeta, constraint)
    psy1 <- psy[[1]]
    psy2 <- psy[[2]]
  }

  resThetaYList@psy[[splitedClusInd[1]]] <- psy1
  resThetaYList@psy[[splitedClusInd[2]]] <- psy2


  if (constraint[1] == 1) {
    lambda1 <- thetaYList@lambda[[splitClusInd]]
    lambda2 <- thetaYList@lambda[[splitClusInd]]
  } else {
    lambda <- generatePriorLambda(p, 2, alpha2, qVec[splitedClusInd], psy, constraint)
    lambda1 <- lambda[[1]]
    lambda2 <- lambda[[2]]
  }

  resThetaYList@lambda[[splitedClusInd[1]]] <- lambda1
  resThetaYList@lambda[[splitedClusInd[2]]] <- lambda2

  ## propose Y
  for (k in splitedClusInd) {
    resThetaYList@Y[[k]] <- matrix(NA, qVec[k], n)
    for (i in 1:n) {
      resThetaYList@Y[[k]][, i] <- mvtnorm::rmvnorm(1, mean = rep(0, qVec[k]), sigma = diag(qVec[k]))
    }
  }
  return(resThetaYList)
}



#' evaluateSplitedClusters
#'
#' @param X X
#' @param thetaYList thetaYList
#' @param splitedThetaYList splitedThetaYList
#' @param hparam hparam
#' @param splitClusInd splitClusInd
#' @param splitedClusInd splitedClusInd
#' @param qVec qVec
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
#' splitedThetaYList <- thetaYList <-
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
#' splitClusInd <- 1
#' splitedClusInd <- c(1, 1)
#' evaluateSplitedClusters(X, thetaYList, splitedThetaYList, hparam, splitClusInd, splitedClusInd, qVec, constraint)
evaluateSplitedClusters <- function(X, thetaYList, splitedThetaYList, hparam, splitClusInd, splitedClusInd, qVec, constraint) {
  p <- nrow(X)
  n <- ncol(X)

  ggamma <- hparam@ggamma
  delta <- hparam@delta
  bbeta <- hparam@bbeta
  alpha1 <- hparam@alpha1
  alpha2 <- hparam@alpha2

  ## evaluate a1
  w1 <- splitedThetaYList@tao[splitedClusInd[1]]
  w2 <- splitedThetaYList@tao[splitedClusInd[2]]
  w <- thetaYList@tao[splitClusInd]


  EVALa1 <- dbeta(w1 / w, 2, 2, log = T)

  ## evaluate a2

  splitMu <- thetaYList@M[[splitClusInd]]
  ## calculate split clus Variance

  splitClusVar <- thetaYList@psy[[splitClusInd]] + thetaYList@lambda[[splitClusInd]] %*% t(thetaYList@lambda[[splitClusInd]])
  splitClusSDvec <- sqrt(diag(splitClusVar))
  # print("=======")
  # print(splitClusSDvec)

  mu1 <- splitedThetaYList@M[[splitedClusInd[1]]]
  mu2 <- splitedThetaYList@M[[splitedClusInd[2]]]

  # print(mu1)
  # print(mu2)
  # print(splitMu)
  # print("=======>")
  a2 <- c()
  EVALa2 <- c()
  for (i in 1:p) {
    a2[i] <- abs((mu1[i] - splitMu[i]) / (splitClusSDvec[i] * sqrt(w2 / w1)))
    # a2[i] = abs( (mu2[i] - splitMu[i]) /(splitClusSDvec[i] * sqrt(w1/w2)) )

    # EVALa2[i] = dbeta(a2[i],2,2, log = T)
    # EVALa2[i] = dbeta(a2[i],1,1, log = T)

    EVALa2[i] <- dgamma(a2[i], 1, 2, log = T)
  }
  # print(a2)
  # print(EVALa2)
  # x = rbeta(10000, 2,2)
  # y = x/(1-x)* rgamma(10000,1,2)
  #
  # hist(sqrt(y))
  #
  # hist(rgamma(10000, 2,4))
  # hist(rbeta(10000, 2,2))

  # print(w1/w)
  # print("===========evaluate>")
  # print(a2)
  # print("===========evaluate>")
  ## eval  psy
  psy1 <- splitedThetaYList@psy[[splitedClusInd[1]]]
  psy2 <- splitedThetaYList@psy[[splitedClusInd[2]]]
  psy <- list(psy1, psy2)
  evalpsy <- 0
  if (constraint[2] == 0) {
    evalpsy <- evaluatePriorPsi(psy, p, 2, delta, bbeta, constraint, c(1, 1))
  }

  ## eval lambda
  lambda1 <- splitedThetaYList@lambda[[splitedClusInd[1]]]
  lambda2 <- splitedThetaYList@lambda[[splitedClusInd[2]]]
  lambda <- list(lambda1, lambda2)
  evallambda <- 0
  if (constraint[1] == 0) {
    evallambda <- evaluatePriorLambda(p, 2, alpha2, qVec[splitedClusInd], psy, lambda, constraint, c(1, 1))
  }

  # signEval = log(0.5 ^ p)
  # print(c(EVALa1, EVALa2, evalpsy, evallambda))
  return(sum(EVALa1, EVALa2, evalpsy, evallambda))
}


#' splitZOneDim
#'
#' @param ZOneDim ZOneDim
#' @param splitedThetaYList splitedThetaYList
#' @param splitClusInd splitClusInd
#' @param splitedClusInd splitedClusInd
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
#' splitedThetaYList <-
#'   new("ThetaYList",
#'     tao = c(0.90162050961987, 0.0983794903801295),
#'     psy = list(
#'       structure(c(3.68472841602225, 0, 0, 8.34691978354054), .Dim = c(2L, 2L)),
#'       structure(c(0.785011896130842, 0, 0, 1.19022383323437), .Dim = c(2L, 2L))
#'     ),
#'     M = list(structure(c(
#'       2.96424305287004,
#'       1.08454861414306
#'     ), .Dim = 1:2), structure(c(
#'       -0.232625450433964,
#'       0.984505960868685
#'     ), .Dim = 1:2)), lambda = list(structure(c(
#'       -0.964026624054337,
#'       0.89378616732449
#'     ), .Dim = 2:1), structure(c(
#'       0.533334148228635,
#'       -1.80033696090263
#'     ), .Dim = 2:1)), Y = list(structure(c(
#'       -0.15346475266988,
#'       1.6584112693271, 0.409294936277862, -1.46628591247549, -0.532753243163142,
#'       -0.332143130316749, 0.307558110800446, -0.525374243612587, 0.527667526535661,
#'       0.748193650431916
#'     ), .Dim = c(1L, 10L)), structure(c(
#'       0.571325118638535,
#'       0.542462985882966, 0.559971315637159, -1.73905343105432, -0.583549598471542,
#'       1.71264245945391, -0.327119395945831, 1.02464651767821, -1.11462280255215,
#'       0.81095592501554
#'     ), .Dim = c(1L, 10L)))
#'   )
#' splitClusInd <- 1
#' splitedClusInd <- c(1, 2)
#' splitZOneDim(ZOneDim, splitedThetaYList, splitClusInd, splitedClusInd)
splitZOneDim <- function(ZOneDim, splitedThetaYList, splitClusInd, splitedClusInd) {
  resZOneDim <- ZOneDim
  evalProb <- c()

  splitObsInd <- which(ZOneDim == splitClusInd)
  for (j in 1:length(splitObsInd)) {
    # print(j)
    i <- splitObsInd[j]

    p1 <- log(splitedThetaYList@tao[splitedClusInd[1]]) + mvtnorm::dmvnorm(X[, i],
      mean = splitedThetaYList@M[[splitedClusInd[1]]],
      sigma = splitedThetaYList@psy[[splitedClusInd[1]]] + splitedThetaYList@lambda[[splitedClusInd[1]]] %*% t(splitedThetaYList@lambda[[splitedClusInd[1]]]),
      log = T
    )
    # p1 =   mvtnorm::dmvnorm(X[,i], mean = splitedThetaYList@M[[splitedClusInd[1]]],
    #                                                               sigma = splitedThetaYList@psy[[splitedClusInd[1]]] + splitedThetaYList@lambda[[splitedClusInd[1]]]%*%t(splitedThetaYList@lambda[[splitedClusInd[1]]]),
    #                                                               log = T)

    p2 <- log(splitedThetaYList@tao[splitedClusInd[2]]) + mvtnorm::dmvnorm(X[, i],
      mean = splitedThetaYList@M[[splitedClusInd[2]]],
      sigma = splitedThetaYList@psy[[splitedClusInd[2]]] + splitedThetaYList@lambda[[splitedClusInd[2]]] %*% t(splitedThetaYList@lambda[[splitedClusInd[2]]]),
      log = T
    )
    # p2 =  mvtnorm::dmvnorm(X[,i], mean = splitedThetaYList@M[[splitedClusInd[2]]],
    #                                                               sigma = splitedThetaYList@psy[[splitedClusInd[2]]] + splitedThetaYList@lambda[[splitedClusInd[2]]]%*%t(splitedThetaYList@lambda[[splitedClusInd[2]]]),
    #                                                               log = T)
    # cat(" for sub", i,"===\n")
    # cat("p1 = ", p1, "===\n")
    # cat("p2 = ", p2, "===\n")
    # print("=======>")
    if (is.infinite(p1)) {
      print("===infi")
      p1 <- -1e-10
    }
    if (is.infinite(p2)) {
      print("===infi")
      p2 <- -1e-10
    }
    # print(c(p1,p2))
    prob <- calculateRatio(p1, c(p1, p2))
    # cat("prob = ", prob, "====>\n")
    probVec <- c(prob, 1 - prob)

    Zind <- sample(x = 1:2, size = 1, prob = probVec)
    evalProb[j] <- probVec[Zind]
    resZOneDim[i] <- splitedClusInd[Zind]
  }
  return(list(ZOneDim = resZOneDim, evalProb = sum(log(evalProb))))
}



#' splitEvalZOneDim
#'
#' @param ZOneDimSplit ZOneDimSplit
#' @param ZOneDimSplited ZOneDimSplited
#' @param splitedThetaYList splitedThetaYList
#' @param splitClusInd splitClusInd
#' @param splitedClusInd splitedClusInd
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
#' ZOneDimSplit <- sample(seq_len(m), n, replace = TRUE)
#' ZOneDimSplit[ZOneDimSplit == 2] <- 1
#' constraint <- c(0, 0, 0)
#' splitedThetaYList <- new("ThetaYList", tao = c(0.90162050961987, 0.0983794903801295), psy = list(structure(c(3.68472841602225, 0, 0, 8.34691978354054), .Dim = c(2L, 2L)), structure(c(0.785011896130842, 0, 0, 1.19022383323437), .Dim = c(2L, 2L))), M = list(structure(c(
#'   2.96424305287004,
#'   1.08454861414306
#' ), .Dim = 1:2), structure(c(
#'   -0.232625450433964,
#'   0.984505960868685
#' ), .Dim = 1:2)), lambda = list(structure(c(
#'   -0.964026624054337,
#'   0.89378616732449
#' ), .Dim = 2:1), structure(c(
#'   0.533334148228635,
#'   -1.80033696090263
#' ), .Dim = 2:1)), Y = list(structure(c(
#'   -0.15346475266988,
#'   1.6584112693271, 0.409294936277862, -1.46628591247549, -0.532753243163142,
#'   -0.332143130316749, 0.307558110800446, -0.525374243612587, 0.527667526535661,
#'   0.748193650431916
#' ), .Dim = c(1L, 10L)), structure(c(
#'   0.571325118638535,
#'   0.542462985882966, 0.559971315637159, -1.73905343105432, -0.583549598471542,
#'   1.71264245945391, -0.327119395945831, 1.02464651767821, -1.11462280255215,
#'   0.81095592501554
#' ), .Dim = c(1L, 10L))))
#' splitClusInd <- 1
#' splitedClusInd <- c(1, 2)
#' ZOneDimSplited <- sample(seq_len(m), n, replace = TRUE)
#' splitEvalZOneDim(ZOneDimSplit, ZOneDimSplited, splitedThetaYList, splitClusInd, splitedClusInd)
splitEvalZOneDim <- function(ZOneDimSplit, ZOneDimSplited, splitedThetaYList, splitClusInd, splitedClusInd) {
  evalProb <- c()

  splitObsInd <- which(ZOneDimSplit == splitClusInd)
  for (j in 1:length(splitObsInd)) {
    i <- splitObsInd[j]

    p1 <- log(splitedThetaYList@tao[splitedClusInd[1]]) + mvtnorm::dmvnorm(X[, i],
      mean = splitedThetaYList@M[[splitedClusInd[1]]],
      sigma = splitedThetaYList@psy[[splitedClusInd[1]]] + splitedThetaYList@lambda[[splitedClusInd[1]]] %*% t(splitedThetaYList@lambda[[splitedClusInd[1]]]),
      log = T
    )


    p2 <- log(splitedThetaYList@tao[splitedClusInd[2]]) + mvtnorm::dmvnorm(X[, i],
      mean = splitedThetaYList@M[[splitedClusInd[2]]],
      sigma = splitedThetaYList@psy[[splitedClusInd[2]]] + splitedThetaYList@lambda[[splitedClusInd[2]]] %*% t(splitedThetaYList@lambda[[splitedClusInd[2]]]),
      log = T
    )

    prob <- calculateRatio(p1, c(p1, p2))
    # cat("prob = ", prob, "====>\n")
    probVec <- c(prob, 1 - prob)
    # print(prob)
    Zind <- which(splitedClusInd == ZOneDimSplited[i])
    # Zind = sample(x = 1:2,size = 1, prob = probVec)
    evalProb[j] <- probVec[Zind]
  }
  return(sum(log(evalProb)))
}


#' calculateJacobian
#'
#' @param p p
#' @param thetaYList thetaYList
#' @param splitedThetaYList splitedThetaYList
#' @param splitClusInd splitClusInd
#' @param splitedClusInd splitedClusInd
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
#' splitedThetaYList <- new("ThetaYList", tao = c(0.90162050961987, 0.0983794903801295), psy = list(structure(c(3.68472841602225, 0, 0, 8.34691978354054), .Dim = c(2L, 2L)), structure(c(0.785011896130842, 0, 0, 1.19022383323437), .Dim = c(2L, 2L))), M = list(structure(c(
#'   2.96424305287004,
#'   1.08454861414306
#' ), .Dim = 1:2), structure(c(
#'   -0.232625450433964,
#'   0.984505960868685
#' ), .Dim = 1:2)), lambda = list(structure(c(
#'   -0.964026624054337,
#'   0.89378616732449
#' ), .Dim = 2:1), structure(c(
#'   0.533334148228635,
#'   -1.80033696090263
#' ), .Dim = 2:1)), Y = list(structure(c(
#'   -0.15346475266988,
#'   1.6584112693271, 0.409294936277862, -1.46628591247549, -0.532753243163142,
#'   -0.332143130316749, 0.307558110800446, -0.525374243612587, 0.527667526535661,
#'   0.748193650431916
#' ), .Dim = c(1L, 10L)), structure(c(
#'   0.571325118638535,
#'   0.542462985882966, 0.559971315637159, -1.73905343105432, -0.583549598471542,
#'   1.71264245945391, -0.327119395945831, 1.02464651767821, -1.11462280255215,
#'   0.81095592501554
#' ), .Dim = c(1L, 10L))))
#' splitClusInd <- 1
#' splitedClusInd <- c(1, 2)
#' calculateJacobian(p, thetaYList, splitedThetaYList, splitClusInd, splitedClusInd)
calculateJacobian <- function(p, thetaYList, splitedThetaYList, splitClusInd, splitedClusInd) {
  w1 <- splitedThetaYList@tao[splitedClusInd[1]]
  w2 <- splitedThetaYList@tao[splitedClusInd[2]]
  w <- thetaYList@tao[splitClusInd]

  splitClusVar <- thetaYList@psy[[splitClusInd]] + thetaYList@lambda[[splitClusInd]] %*% t(thetaYList@lambda[[splitClusInd]])
  splitClusSDvec <- sqrt(diag(splitClusVar))

  mu1 <- splitedThetaYList@M[[splitedClusInd[1]]]
  mu2 <- splitedThetaYList@M[[splitedClusInd[2]]]
  splitMu <- thetaYList@M[[splitClusInd]]

  a2 <- c()
  for (i in 1:p) {
    a2[i] <- abs((mu1[i] - splitMu[i]) / (splitClusSDvec[i] * sqrt(w2 / w1)))
  }

  res <- log(w) + sum(log(abs((mu1 - mu2) / a2)))
  return(res)
}
