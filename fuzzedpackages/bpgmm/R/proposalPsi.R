#' CalculateProposalPsy
#'
#' @param hparam hparam
#' @param thetaYList thetaYList
#' @param CxyList CxyList
#' @param constraint constraint
#' @param m the number of clusters
#' @param qVec the vector of the number of factors in each clusters
#' @param p the number of features
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
#' constraint <- c(0, 0, 0)
#' CxyList <-
#'   list(
#'     A = list(structure(
#'       c(0.567755037123148, 0, 0, 1.1870201935945),
#'       .Dim = c(2L, 2L)
#'     )),
#'     nVec = structure(10, .Dim = c(1L, 1L)),
#'     Cxxk = list(structure(
#'       c(
#'         739.129405647622,
#'         671.040583460732,
#'         671.040583460732,
#'         618.754338945564
#'       ),
#'       .Dim = c(2L, 2L)
#'     )),
#'     Cxyk = list(structure(
#'       c(-18.5170828875512, -16.5748393456787),
#'       .Dim = 2:1
#'     )),
#'     Cyyk = list(structure(2.4786991560888, .Dim = c(
#'       1L,
#'       1L
#'     ))),
#'     Cytytk = list(structure(
#'       c(
#'         10, 0.787438922114998, 0.787438922114998,
#'         2.4786991560888
#'       ),
#'       .Dim = c(2L, 2L)
#'     )),
#'     Cxtytk = list(structure(
#'       c(
#'         -57.5402230447872,
#'         -54.6677145995824,
#'         -18.5170828875512,
#'         -16.5748393456787
#'       ),
#'       .Dim = c(
#'         2L,
#'         2L
#'       )
#'     )),
#'     CxL1k = list(structure(
#'       c(-59.5168204264758, -54.6093504204781),
#'       .Dim = 2:1
#'     )),
#'     Cxmyk = list(structure(
#'       c(
#'         -21.0952527723962,
#'         -14.6807011202188
#'       ),
#'       .Dim = 2:1
#'     )),
#'     sumCxmyk = structure(c(
#'       -21.0952527723962,
#'       -14.6807011202188
#'     ), .Dim = 2:1),
#'     sumCyyk = structure(3.6657193496833, .Dim = c(
#'       1L,
#'       1L
#'     ))
#'   )
#' #'
#' \donttest{
#' CalculateProposalPsy(hparam, thetaYList, CxyList, constraint, m, p, qVec)
#' }
CalculateProposalPsy <- function(hparam, thetaYList, CxyList, constraint, m, p, qVec) {
  alpha1 <- hparam@alpha1
  alpha2 <- hparam@alpha2
  bbeta <- hparam@bbeta
  delta <- hparam@delta
  M <- thetaYList@M
  psy <- thetaYList@psy
  lambda <- thetaYList@lambda


  Cxxk <- CxyList$Cxxk
  Cxyk <- CxyList$Cxyk
  Cyyk <- CxyList$Cyyk
  Cytytk <- CxyList$Cytytk
  Cxtytk <- CxyList$Cxtytk
  CxL1k <- CxyList$CxL1k
  Cxmyk <- CxyList$Cxmyk
  sumCxmyk <- CxyList$sumCxmyk
  sumCyyk <- CxyList$sumCyyk
  A <- CxyList$A
  nVec <- CxyList$nVec

  # post tilda lambda_k = {mu_k, lambda_k}, first column is mu_k
  tildaLambda <- list()
  # print(M)
  # print(lambda)
  for (k in 1:m) {
    tildaLambda[[k]] <- cbind(t(M[[k]]), lambda[[k]])
  }

  if (constraint[1] == T & constraint[2] == T & constraint[3] == T) {
    ## model 1

    ## post psy
    psy <- list()

    shapePara <- 0
    ratePara <- 0
    for (k in 1:m) {
      shapePara <- shapePara + p / 2 * (nVec[k] + qVec[k] / m + (2 * delta - 2) / (m * p) + 1)
      ratePara <- ratePara + 1 / 2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
        + tildaLambda[[k]] %*% (Cytytk[[k]] + A[[k]] / m) %*% t(tildaLambda[[k]])
        + (2 * bbeta) / (m * p) * diag(rep(1, p)))
    }
    shapePara <- shapePara + 1
    ratePara <- sum(ratePara)

    invpsy <- rgamma(1, shape = shapePara, rate = ratePara)
    for (k in 1:m) {
      psy[[k]] <- diag(rep(1 / invpsy, p))
    }
    ## model 1 end
  } else if (constraint[1] == T & constraint[2] == T & constraint[3] == F) {
    ## model 2

    ## post psy
    psy <- list()

    shapePara <- 0
    ratePara <- 0
    for (k in 1:m) {
      shapePara <- shapePara + 1 / 2 * (nVec[k] + qVec[k] / m + (2 * delta - 2) / m + 1)
      ratePara <- ratePara + 1 / 2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
        + tildaLambda[[k]] %*% (Cytytk[[k]] + A[[k]] / m) %*% t(tildaLambda[[k]])
        + 2 * bbeta / m * diag(rep(1, p)))
    }
    shapePara <- shapePara + 1

    invpsy <- c()
    for (j in 1:p) {
      invpsy[j] <- rgamma(1, shape = shapePara, rate = ratePara[j])
    }
    for (k in 1:m) {
      psy[[k]] <- diag(1 / invpsy)
    }
    ## end model 2
  } else if (constraint[1] == T & constraint[2] == F & constraint[3] == T) {
    ## model 3

    ## post psy
    psy <- list()

    shapePara <- 0
    ratePara <- 0
    for (k in 1:m) {
      shapePara <- p / 2 * (nVec[k] + qVec[k] / m + (2 * delta - 2) / p + 1) + 1
      ratePara <- 1 / 2 * sum(diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
        + tildaLambda[[k]] %*% (Cytytk[[k]] + A[[k]]) %*% t(tildaLambda[[k]])
        + 2 * bbeta / p * diag(rep(1, p))))

      invpsy <- rgamma(1, shape = shapePara, rate = ratePara)
      psy[[k]] <- diag(rep(1 / invpsy, p))
    }


    ## end model 3
  } else if (constraint[1] == T & constraint[2] == F & constraint[3] == F) {
    ## model 4
    ## post psy
    psy <- list()

    shapePara <- 0
    ratePara <- 0
    for (k in 1:m) {
      shapePara <- 1 / 2 * (nVec[k] + qVec[k] / m + 2 * delta - 1) + 1
      ratePara <- 1 / 2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
        + tildaLambda[[k]] %*% (Cytytk[[k]] + A[[k]] / m) %*% t(tildaLambda[[k]])
        + 2 * bbeta * diag(rep(1, p)))


      invpsy <- c()
      for (j in 1:p) {
        invpsy[j] <- rgamma(1, shape = shapePara, rate = ratePara[j])
      }
      # invpsy = rgamma(p, shape = shapePara, rate = ratePara)
      psy[[k]] <- diag(1 / invpsy)
    }

    ## end model 4
  } else if (constraint[1] == F & constraint[2] == T & constraint[3] == T) {
    ## model 5
    ## post psy
    psy <- list()

    shapePara <- 0
    ratePara <- 0
    for (k in 1:m) {
      shapePara <- shapePara + p / 2 * (nVec[k] + qVec[k] + (2 * delta - 2) / (m * p) + 1)
      ratePara <- ratePara + 1 / 2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
        + tildaLambda[[k]] %*% (Cytytk[[k]] + A[[k]]) %*% t(tildaLambda[[k]])
        + 2 * bbeta / (m * p) * diag(rep(1, p)))
    }
    shapePara <- shapePara + 1
    ratePara <- sum(ratePara)

    invpsy <- rgamma(1, shape = shapePara, rate = ratePara)

    for (k in 1:m) {
      psy[[k]] <- diag(rep(1 / invpsy), p)
    }
    ## end model 5
  } else if (constraint[1] == F & constraint[2] == T & constraint[3] == F) {
    ## model 6
    ## post psy
    psy <- list()

    shapePara <- 0
    ratePara <- 0
    for (k in 1:m) {
      shapePara <- shapePara + 1 / 2 * (nVec[k] + qVec[k] + (2 * delta - 2) / m + 1)
      ratePara <- ratePara + 1 / 2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
        + tildaLambda[[k]] %*% (Cytytk[[k]] + A[[k]]) %*% t(tildaLambda[[k]])
        + 2 * bbeta / m * diag(rep(1, p)))
    }
    shapePara <- shapePara + 1

    invpsy <- c()
    for (j in 1:p) {
      invpsy[j] <- rgamma(1, shape = shapePara, rate = ratePara[j])
    }
    for (k in 1:m) {
      psy[[k]] <- diag(1 / invpsy)
    }
    ## end model 6
  } else if (constraint[1] == F & constraint[2] == F & constraint[3] == T) {
    ## model 7

    ## post psy
    psy <- list()

    shapePara <- 0
    ratePara <- 0
    for (k in 1:m) {
      shapePara <- p / 2 * (nVec[k] + qVec[k] + (2 * delta - 2) / p + 1) + 1
      ratePara <- 1 / 2 * sum(diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
        + tildaLambda[[k]] %*% (Cytytk[[k]] + A[[k]]) %*% t(tildaLambda[[k]])
        + 2 * bbeta / p * diag(rep(1, p))))

      invpsy <- rgamma(1, shape = shapePara, rate = ratePara)
      psy[[k]] <- diag(rep(1 / invpsy, p))
    }
    ## end model 7
  } else if (constraint[1] == F & constraint[2] == F & constraint[3] == F) {
    ## model 8
    ## post psy
    psy <- list()

    shapePara <- 0
    ratePara <- 0
    for (k in 1:m) {
      shapePara <- 1 / 2 * (nVec[k] + qVec[k] + 2 * delta - 1) + 1
      ratePara <- 1 / 2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
        + tildaLambda[[k]] %*% (Cytytk[[k]] + A[[k]]) %*% t(tildaLambda[[k]])
        + 2 * bbeta * diag(rep(1, p)))

      invpsy <- rgamma(p, shape = shapePara, rate = ratePara)
      psy[[k]] <- diag(1 / invpsy)
    }

    ## end model 8
  }
  return(psy)
}


#' EvaluateProposalPsy
#'
#' @param hparam hparam
#' @param thetaYList thetaYList
#' @param CxyList CxyList
#' @param constraint constraint
#' @param newpsy newpsy
#' @param m the number of clusters
#' @param qVec the vector of the number of factors in each clusters
#' @param p the number of features
#' @param delta hyperparameters
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
#' constraint <- c(0, 0, 0)
#' CxyList <-
#'   list(
#'     A = list(structure(
#'       c(0.567755037123148, 0, 0, 1.1870201935945),
#'       .Dim = c(2L, 2L)
#'     )),
#'     nVec = structure(10, .Dim = c(1L, 1L)),
#'     Cxxk = list(structure(
#'       c(
#'         739.129405647622,
#'         671.040583460732,
#'         671.040583460732,
#'         618.754338945564
#'       ),
#'       .Dim = c(2L, 2L)
#'     )),
#'     Cxyk = list(structure(
#'       c(-18.5170828875512, -16.5748393456787),
#'       .Dim = 2:1
#'     )),
#'     Cyyk = list(structure(2.4786991560888, .Dim = c(
#'       1L,
#'       1L
#'     ))),
#'     Cytytk = list(structure(
#'       c(
#'         10, 0.787438922114998, 0.787438922114998,
#'         2.4786991560888
#'       ),
#'       .Dim = c(2L, 2L)
#'     )),
#'     Cxtytk = list(structure(
#'       c(
#'         -57.5402230447872,
#'         -54.6677145995824,
#'         -18.5170828875512,
#'         -16.5748393456787
#'       ),
#'       .Dim = c(
#'         2L,
#'         2L
#'       )
#'     )),
#'     CxL1k = list(structure(
#'       c(-59.5168204264758, -54.6093504204781),
#'       .Dim = 2:1
#'     )),
#'     Cxmyk = list(structure(
#'       c(
#'         -21.0952527723962,
#'         -14.6807011202188
#'       ),
#'       .Dim = 2:1
#'     )),
#'     sumCxmyk = structure(c(
#'       -21.0952527723962,
#'       -14.6807011202188
#'     ), .Dim = 2:1),
#'     sumCyyk = structure(3.6657193496833, .Dim = c(
#'       1L,
#'       1L
#'     ))
#'   )
#' \donttest{
#' EvaluateProposalPsy(hparam, thetaYList, CxyList, constraint, thetaYList@psy, m, p, qVec, delta)
#' }
EvaluateProposalPsy <- function(hparam, thetaYList, CxyList, constraint, newpsy, m, p, qVec, delta) {
  alpha1 <- hparam@alpha1
  alpha2 <- hparam@alpha2
  bbeta <- hparam@bbeta
  delta <- hparam@delta
  M <- thetaYList@M
  lambda <- thetaYList@lambda

  ##
  Cxxk <- CxyList$Cxxk
  Cxyk <- CxyList$Cxyk
  Cyyk <- CxyList$Cyyk
  Cytytk <- CxyList$Cytytk
  Cxtytk <- CxyList$Cxtytk
  CxL1k <- CxyList$CxL1k
  Cxmyk <- CxyList$Cxmyk
  sumCxmyk <- CxyList$sumCxmyk
  sumCyyk <- CxyList$sumCyyk
  A <- CxyList$A
  nVec <- CxyList$nVec

  ##
  psyEval <- c()
  psy <- newpsy
  tildaLambda <- list()
  for (k in 1:m) {
    tildaLambda[[k]] <- cbind(t(M[[k]]), lambda[[k]])
  }
  ##
  if (constraint[1] == T & constraint[2] == T & constraint[3] == T) {
    ## model 1

    ## post psy
    shapePara <- 0
    ratePara <- 0
    for (k in 1:m) {
      shapePara <- shapePara + p / 2 * (nVec[k] + qVec[k] / m + (2 * delta - 2) / (m * p) + 1)
      ratePara <- ratePara + 1 / 2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
        + tildaLambda[[k]] %*% (Cytytk[[k]] + A[[k]] / m) %*% t(tildaLambda[[k]])
        + (2 * bbeta) / (m * p) * diag(rep(1, p)))
    }
    shapePara <- shapePara + 1
    ratePara <- sum(ratePara)

    psyEval <- dgamma(1 / psy[[1]][1, 1], shape = shapePara, rate = ratePara, log = T)
    ## model 1 end
  } else if (constraint[1] == T & constraint[2] == T & constraint[3] == F) {
    ## model 2

    shapePara <- 0
    ratePara <- 0
    for (k in 1:m) {
      shapePara <- shapePara + 1 / 2 * (nVec[k] + qVec[k] / m + (2 * delta - 2) / m + 1)
      ratePara <- ratePara + 1 / 2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
        + tildaLambda[[k]] %*% (Cytytk[[k]] + A[[k]] / m) %*% t(tildaLambda[[k]])
        + 2 * bbeta / m * diag(rep(1, p)))
    }
    shapePara <- shapePara + 1

    invpsy <- 1 / diag(psy[[1]])
    for (j in 1:p) {
      psyEval[j] <- dgamma(invpsy[j], shape = shapePara, rate = ratePara[j], log = T)
    }
    ## end model 2
  } else if (constraint[1] == T & constraint[2] == F & constraint[3] == T) {
    ## post psy

    shapePara <- 0
    ratePara <- 0
    for (k in 1:m) {
      shapePara <- p / 2 * (nVec[k] + qVec[k] / m + (2 * delta - 2) / p + 1) + 1
      ratePara <- 1 / 2 * sum(diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
        + tildaLambda[[k]] %*% (Cytytk[[k]] + A[[k]]) %*% t(tildaLambda[[k]])
        + 2 * bbeta / p * diag(rep(1, p))))

      psyEval[k] <- dgamma(1 / psy[[k]][1, 1], shape = shapePara, rate = ratePara, log = T)
    }


    ## end model 3
  } else if (constraint[1] == T & constraint[2] == F & constraint[3] == F) {
    ## model 4
    ## post psy

    shapePara <- 0
    ratePara <- 0
    psyEval <- matrix(NA, nrow = m, ncol = p)
    for (k in 1:m) {
      shapePara <- 1 / 2 * (nVec[k] + qVec[k] / m + 2 * delta - 1) + 1
      ratePara <- 1 / 2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
        + tildaLambda[[k]] %*% (Cytytk[[k]] + A[[k]] / m) %*% t(tildaLambda[[k]])
        + 2 * bbeta * diag(rep(1, p)))


      invpsy <- 1 / diag(psy[[k]])
      for (j in 1:p) {
        psyEval[k, j] <- dgamma(invpsy[j], shape = shapePara, rate = ratePara[j], log = T)
      }
    }

    ## end model 4
  } else if (constraint[1] == F & constraint[2] == T & constraint[3] == T) {
    ## model 5

    ## post psy
    shapePara <- 0
    ratePara <- 0
    for (k in 1:m) {
      shapePara <- shapePara + p / 2 * (nVec[k] + qVec[k] + (2 * delta - 2) / (m * p) + 1)
      ratePara <- ratePara + 1 / 2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
        + tildaLambda[[k]] %*% (Cytytk[[k]] + A[[k]]) %*% t(tildaLambda[[k]])
        + 2 * bbeta / (m * p) * diag(rep(1, p)))
    }
    shapePara <- shapePara + 1
    ratePara <- sum(ratePara)

    psyEval <- dgamma(1 / psy[[1]][1, 1], shape = shapePara, rate = ratePara, log = T)

    ## end model 5
  } else if (constraint[1] == F & constraint[2] == T & constraint[3] == F) {
    ## model 6


    ## post psy

    shapePara <- 0
    ratePara <- 0
    for (k in 1:m) {
      shapePara <- shapePara + 1 / 2 * (nVec[k] + qVec[k] + (2 * delta - 2) / m + 1)
      ratePara <- ratePara + 1 / 2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
        + tildaLambda[[k]] %*% (Cytytk[[k]] + A[[k]]) %*% t(tildaLambda[[k]])
        + 2 * bbeta / m * diag(rep(1, p)))
    }
    shapePara <- shapePara + 1

    invpsy <- 1 / diag(psy[[1]])
    for (j in 1:p) {
      psyEval[j] <- dgamma(invpsy[j], shape = shapePara, rate = ratePara[j], log = T)
    }

    ## end model 6
  } else if (constraint[1] == F & constraint[2] == F & constraint[3] == T) {
    ## model 7


    ## post psy

    shapePara <- 0
    ratePara <- 0
    for (k in 1:m) {
      shapePara <- p / 2 * (nVec[k] + qVec[k] + (2 * delta - 2) / p + 1) + 1
      ratePara <- 1 / 2 * sum(diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
        + tildaLambda[[k]] %*% (Cytytk[[k]] + A[[k]]) %*% t(tildaLambda[[k]])
        + 2 * bbeta / p * diag(rep(1, p))))

      psyEval[k] <- dgamma(1 / psy[[k]][1, 1], shape = shapePara, rate = ratePara, log = T)
    }
    ## end model 7
  } else if (constraint[1] == F & constraint[2] == F & constraint[3] == F) {
    ## model 8

    ## post psy
    psyEval <- matrix(NA, nrow = m, ncol = p)
    shapePara <- 0
    ratePara <- 0
    for (k in 1:m) {
      shapePara <- 1 / 2 * (nVec[k] + qVec[k] + 2 * delta - 1) + 1
      ratePara <- 1 / 2 * diag(Cxxk[[k]] - 2 * Cxtytk[[k]] %*% t(tildaLambda[[k]])
        + tildaLambda[[k]] %*% (Cytytk[[k]] + A[[k]]) %*% t(tildaLambda[[k]])
        + 2 * bbeta * diag(rep(1, p)))


      invpsy <- 1 / diag(psy[[k]])
      for (j in 1:p) {
        psyEval[k, j] <- dgamma(invpsy[j], shape = shapePara, rate = ratePara[j], log = T)
      }
    }
    ## end model 8
  }
  return(sum(psyEval))
}
