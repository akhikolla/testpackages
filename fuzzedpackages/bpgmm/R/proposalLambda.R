#' CalculateProposalLambda
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
#' m <- 1
#' muBar <- c(0, 0)
#'
#' hparam <- new(
#'   "Hparam",
#'   alpha1 = 0.567755037123148,
#'   alpha2 = 1.1870201935945,
#'   delta = 2,
#'   ggamma = 2,
#'   bbeta = 3.39466184520673
#' )
#' qVec <- c(1, 1)
#' constraint <- c(0, 0, 0)
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
#' CalculateProposalLambda(hparam, thetaYList, CxyList, constraint, m, p, qVec)
#' }
#'
CalculateProposalLambda <- function(hparam, thetaYList, CxyList, constraint, m, p, qVec) {
  alpha1 <- hparam@alpha1
  alpha2 <- hparam@alpha2

  M <- thetaYList@M
  psy <- thetaYList@psy
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
  lambda <- list()
  if (constraint[1] == T & constraint[2] == T & constraint[3] == T) {
    ## model 1

    sumCxmyk <- 0
    sumCyyk <- 0
    for (k in 1:m) {
      sumCxmyk <- sumCxmyk + Cxmyk[[k]]
      sumCyyk <- sumCyyk + Cyyk[[k]] + alpha2 / m * diag(qVec[k])
    }

    for (k in 1:m) {
      if (k == 1) {
        lambda[[k]] <- mvtnorm::rmvnorm(1,
          mean = c(sumCxmyk %*% solve(sumCyyk)),
          sigma = kronecker(solve(sumCyyk), psy[[k]])
        )
        lambda[[k]] <- matrix(lambda[[k]], p, qVec[k])
      } else {
        lambda[[k]] <- lambda[[1]]
      }
    }
    ## model 1 end
  } else if (constraint[1] == T & constraint[2] == T & constraint[3] == F) {
    ## model 2

    sumCxmyk <- 0
    sumCyyk <- 0
    for (k in 1:m) {
      sumCxmyk <- sumCxmyk + Cxmyk[[k]]
      sumCyyk <- sumCyyk + Cyyk[[k]] + alpha2 / m * diag(qVec[k])
    }
    for (k in 1:m) {
      if (k == 1) {
        lambda[[k]] <- mvtnorm::rmvnorm(1,
          mean = c(sumCxmyk %*% solve(sumCyyk)),
          sigma = kronecker(solve(sumCyyk), psy[[k]])
        )
        lambda[[k]] <- matrix(lambda[[k]], p, qVec[k])
      } else {
        lambda[[k]] <- lambda[[1]]
      }
    }

    ## end model 2
  } else if (constraint[1] == T & constraint[2] == F & constraint[3] == T) {
    ## model 3


    sumPhiCxy <- 0
    sumPhiCyy <- 0

    for (k in 1:m) {
      sumPhiCxy <- sumPhiCxy + 1 / psy[[k]][1, 1] * Cxmyk[[k]]
      sumPhiCyy <- sumPhiCyy + 1 / psy[[k]][1, 1] * (Cyyk[[k]] + alpha2 / m * diag(qVec[k]))
    }

    for (k in 1:m) {
      if (k == 1) {
        lambda[[k]] <- mvtnorm::rmvnorm(1,
          mean = c(sumPhiCxy %*% solve(sumPhiCyy)),
          sigma = kronecker(solve(sumPhiCyy), diag(p))
        )

        lambda[[k]] <- matrix(lambda[[k]], p, qVec[k])
      } else {
        lambda[[k]] <- lambda[[1]]
      }
    }

    ## end model 3
  } else if (constraint[1] == T & constraint[2] == F & constraint[3] == F) {
    ## model 4
    sumVar <- 0
    B <- 0
    for (k in 1:m) {
      sumVar <- sumVar + kronecker(
        Cyyk[[k]] + alpha2 / m * diag(qVec[k]),
        solve(psy[[k]])
      )
      B <- B + solve(psy[[k]]) %*% Cxmyk[[k]]
    }
    lambdaVar <- solve(sumVar)
    lambdaMean <- t(c(B)) %*% lambdaVar
    for (k in 1:m) {
      if (k == 1) {
        lambda[[k]] <- mvtnorm::rmvnorm(1,
          mean = lambdaMean,
          sigma = lambdaVar
        )

        lambda[[k]] <- matrix(lambda[[k]], p, qVec[k])
      } else {
        lambda[[k]] <- lambda[[1]]
      }
    }
    ## end model 4
  } else if (constraint[1] == F & constraint[2] == T & constraint[3] == T) {
    ## model 5
    for (k in 1:m) {
      lambda[[k]] <- mvtnorm::rmvnorm(1,
        mean = c(Cxmyk[[k]] %*% solve(Cyyk[[k]] + alpha2 * diag(qVec[k]))),
        sigma = kronecker(solve(Cyyk[[k]] + alpha2 * diag(qVec[k])), psy[[k]])
      )
      lambda[[k]] <- matrix(lambda[[k]], p, qVec[k])
    }

    ## end model 5
  } else if (constraint[1] == F & constraint[2] == T & constraint[3] == F) {
    ## model 6
    for (k in 1:m) {
      lambda[[k]] <- mvtnorm::rmvnorm(1,
        mean = c(Cxmyk[[k]] %*% solve(Cyyk[[k]] + alpha2 * diag(qVec[k]))),
        sigma = kronecker(solve(sumCyyk), psy[[k]])
      )
      lambda[[k]] <- matrix(lambda[[k]], p, qVec[k])
    }

    ## end model 6
  } else if (constraint[1] == F & constraint[2] == F & constraint[3] == T) {
    ## model 7
    for (k in 1:m) {
      lambda[[k]] <- mvtnorm::rmvnorm(1,
        mean = c(Cxmyk[[k]] %*% solve(Cyyk[[k]] + alpha2 * diag(qVec[k]))),
        sigma = kronecker(solve(sumCyyk), psy[[k]])
      )
      lambda[[k]] <- matrix(lambda[[k]], p, qVec[k])
    }
    ## end model 7
  } else if (constraint[1] == F & constraint[2] == F & constraint[3] == F) {
    ## model 8
    for (k in 1:m) {
      lambda[[k]] <- mvtnorm::rmvnorm(1,
        mean = c(Cxmyk[[k]] %*% solve(Cyyk[[k]] + alpha2 * diag(qVec[k]))),
        sigma = kronecker(solve(sumCyyk), psy[[k]])
      )
      lambda[[k]] <- matrix(lambda[[k]], p, qVec[k])
    }
    ## end model 8
  }
  return(lambda)
}

#' EvaluateProposalLambda
#'
#' @param hparam hparam
#' @param thetaYList thetaYList
#' @param CxyList CxyList
#' @param constraint constraint
#' @param newlambda newlambda
#' @param m the number of clusters
#' @param qVec the vector of the number of factors in each clusters
#' @param p the number of features
#'
#' @export
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
#' EvaluateProposalLambda(hparam, thetaYList, CxyList, constraint, thetaYList@lambda, m, qVec, p)
#' }
EvaluateProposalLambda <- function(hparam, thetaYList, CxyList, constraint, newlambda, m, qVec, p) {
  alpha1 <- hparam@alpha1
  alpha2 <- hparam@alpha2

  M <- thetaYList@M
  psy <- thetaYList@psy
  lambda <- newlambda
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
  lambdaEval <- c()
  if (constraint[1] == T & constraint[2] == T & constraint[3] == T) {
    ## model 1

    sumCxmyk <- 0
    sumCyyk <- 0
    for (k in 1:m) {
      sumCxmyk <- sumCxmyk + Cxmyk[[k]]
      sumCyyk <- sumCyyk + Cyyk[[k]] + alpha2 / m * diag(qVec[k])
    }

    for (k in 1:m) {
      if (k == 1) {
        lambdaEval[k] <- mvtnorm::dmvnorm(
          x = c(lambda[[k]]), mean = c(sumCxmyk %*% solve(sumCyyk)),
          sigma = kronecker(solve(sumCyyk), psy[[k]]), log = T
        )
      } else {
        lambdaEval[k] <- 0
      }
    }
    ## model 1 end
  } else if (constraint[1] == T & constraint[2] == T & constraint[3] == F) {
    ## model 2

    sumCxmyk <- 0
    sumCyyk <- 0
    for (k in 1:m) {
      sumCxmyk <- sumCxmyk + Cxmyk[[k]]
      sumCyyk <- sumCyyk + Cyyk[[k]] + alpha2 / m * diag(qVec[k])
    }

    for (k in 1:m) {
      if (k == 1) {
        lambdaEval[k] <- mvtnorm::dmvnorm(
          x = c(lambda[[k]]), mean = c(sumCxmyk %*% solve(sumCyyk)),
          sigma = kronecker(solve(sumCyyk), psy[[k]]),
          log = T
        )
      } else {
        lambdaEval[k] <- 0
      }
    }

    ## end model 2
  } else if (constraint[1] == T & constraint[2] == F & constraint[3] == T) {
    ## model 3

    sumPhiCxy <- 0
    sumPhiCyy <- 0

    for (k in 1:m) {
      sumPhiCxy <- sumPhiCxy + 1 / psy[[k]][1, 1] * Cxmyk[[k]]
      sumPhiCyy <- sumPhiCyy + 1 / psy[[k]][1, 1] * (Cyyk[[k]] + alpha2 / m * diag(qVec[k]))
    }

    for (k in 1:m) {
      if (k == 1) {
        lambdaEval[k] <- mvtnorm::dmvnorm(
          x = c(lambda[[k]]), mean = c(sumPhiCxy %*% solve(sumPhiCyy)),
          sigma = kronecker(solve(sumPhiCyy), diag(p)), log = T
        )
      } else {
        lambdaEval[k] <- 0
      }
    }
    ## end model 3
  } else if (constraint[1] == T & constraint[2] == F & constraint[3] == F) {
    ## model 4
    sumVar <- 0
    B <- 0
    for (k in 1:m) {
      sumVar <- sumVar + kronecker(
        Cyyk[[k]] + alpha2 / m * diag(qVec[k]),
        solve(psy[[k]])
      )
      B <- B + solve(psy[[k]]) %*% Cxmyk[[k]]
    }
    lambdaVar <- solve(sumVar)
    lambdaMean <- t(c(B)) %*% lambdaVar
    for (k in 1:m) {
      if (k == 1) {
        lambdaEval[k] <- mvtnorm::dmvnorm(
          x = c(lambda[[k]]), mean = lambdaMean,
          sigma = lambdaVar, log = T
        )
      } else {
        lambdaEval[k] <- 0
      }
    }
    ## end model 4
  } else if (constraint[1] == F & constraint[2] == T & constraint[3] == T) {
    ## model 5
    for (k in 1:m) {
      lambdaEval[k] <- mvtnorm::dmvnorm(
        x = c(lambda[[k]]), mean = c(Cxmyk[[k]] %*% solve(Cyyk[[k]] + alpha2 * diag(qVec[k]))),
        sigma = kronecker(solve(sumCyyk), psy[[k]]), log = T
      )
    }

    ## end model 5
  } else if (constraint[1] == F & constraint[2] == T & constraint[3] == F) {
    ## model 6
    for (k in 1:m) {
      lambdaEval[k] <- mvtnorm::dmvnorm(
        x = c(lambda[[k]]), mean = c(Cxmyk[[k]] %*% solve(Cyyk[[k]] + alpha2 * diag(qVec[k]))),
        sigma = kronecker(solve(sumCyyk), psy[[k]]), log = T
      )
    }

    ## end model 6
  } else if (constraint[1] == F & constraint[2] == F & constraint[3] == T) {
    ## model 7
    for (k in 1:m) {
      lambdaEval[k] <- mvtnorm::dmvnorm(
        x = c(lambda[[k]]), mean = c(Cxmyk[[k]] %*% solve(Cyyk[[k]] + alpha2 * diag(qVec[k]))),
        sigma = kronecker(solve(sumCyyk), psy[[k]]), log = T
      )
    }
    ## end model 7
  } else if (constraint[1] == F & constraint[2] == F & constraint[3] == F) {
    ## model 8
    for (k in 1:m) {
      lambdaEval[k] <- mvtnorm::dmvnorm(c(lambda[[k]]),
        mean = c(Cxmyk[[k]] %*% solve(Cyyk[[k]] + alpha2 * diag(qVec[k]))),
        sigma = kronecker(solve(sumCyyk), psy[[k]]), log = T
      )
    }
    ## end model 8
  }
  return(sum(lambdaEval))
}
