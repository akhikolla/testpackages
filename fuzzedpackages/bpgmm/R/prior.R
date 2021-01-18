#' PriorThetaY list
#' @description generate prior value for parameter Theta and Y.
#' @importFrom gtools rdirichlet
#' @importFrom mvtnorm rmvnorm
#' @import stats
#' @param m the number of cluster
#' @param n sample size
#' @param p number of covariates
#' @param muBar parameter
#' @param hparam hyperparameters
#' @param qVec the vector of the number of factors in each clusters
#' @param ZOneDim ZOneDim
#' @param constraint constraint
#' @examples
#' m <- 20
#' n <- 500
#' p <- 10
#' muBar <- c(
#'   -33.1342706763595, -35.2699639183419, 48.276928009445, 16.2370659601218,
#'   19.0023163870536, -23.4900965314972, 37.1081269873873, 4.74944562930846,
#'   4.6918997353449, -4.55088073255655
#' )
#' hparam <- new("Hparam",
#'   alpha1 = 0.567755037123148, alpha2 = 1.1870201935945,
#'   delta = 2, ggamma = 2, bbeta = 3.39466184520673
#' )
#' qVec <- c(4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4)
#' ZOneDim <- sample(seq_len(m), n, replace = TRUE)
#' constraint <- c(0, 0, 0)
#' #'
#' \donttest{
#' generatePriorThetaY(m, n, p, muBar, hparam, qVec, ZOneDim, constraint)
#' }
#'
#' @export
generatePriorThetaY <- function(m,
                                n,
                                p,
                                muBar,
                                hparam,
                                qVec,
                                ZOneDim,
                                constraint) {
  ggamma <- hparam@ggamma
  delta <- hparam@delta
  bbeta <- hparam@bbeta
  alpha1 <- hparam@alpha1
  alpha2 <- hparam@alpha2

  # prior tao
  if (m == 1) {
    tao <- rbeta(n = 1, shape1 = 1, shape2 = m)
  } else {
    tao <- gtools::rdirichlet(n = 1, alpha = rep(ggamma, m))
  }

  # prior psy
  psy <- generatePriorPsi(p, m, delta, bbeta, constraint)

  # prior M
  M <- list()
  for (i in 1:m) {
    M[[i]] <- (mvtnorm::rmvnorm(n = 1, mean = muBar, sigma = 1 / alpha1 * psy[[i]]))
  }

  # prior lambda
  lambda <- generatePriorLambda(p, m, alpha2, qVec, psy, constraint)


  # Zmat = getZmat(ZOneDim, m, n)
  Zmat <- get_Z_mat(ZOneDim, m, n)

  # post Y
  Y <- list()
  for (k in 1:m) {
    Y[[k]] <- matrix(NA, qVec[k], n)
    for (i in 1:n) {
      Y[[k]][, i] <- mvtnorm::rmvnorm(1, mean = rep(0, qVec[k]), sigma = diag(qVec[k]))
    }
  }

  new("ThetaYList",
    tao = tao,
    psy = psy,
    M = M,
    lambda = lambda,
    Y = Y
  )
}

#' evaluate Prior
#' @description evaluate prior value for parameter Theta and Y.
#' @importFrom gtools ddirichlet
#' @importFrom mvtnorm dmvnorm
#' @import stats
#' @param m m
#' @param p p
#' @param muBar mu_bar
#' @param hparam hyper parameter class
#' @param thetaYList theta Y list
#' @param ZOneDim one dim of z
#' @param qVec q vector
#' @param constraint type of constraint
#' @param clusInd cluster indicator vector
#' @examples
#' m <- 20
#' n <- 500
#' p <- 10
#' muBar <- c(
#'   -33.1342706763595, -35.2699639183419, 48.276928009445, 16.2370659601218,
#'   19.0023163870536, -23.4900965314972, 37.1081269873873, 4.74944562930846,
#'   4.6918997353449, -4.55088073255655
#' )
#' hparam <- new("Hparam",
#'   alpha1 = 0.567755037123148, alpha2 = 1.1870201935945,
#'   delta = 2, ggamma = 2, bbeta = 3.39466184520673
#' )
#' qVec <- c(4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4)
#' ZOneDim <- sample(seq_len(m), n, replace = TRUE)
#' constraint <- c(0, 0, 0)
#' thetaYList <- generatePriorThetaY(m, n, p, muBar, hparam, qVec, ZOneDim, constraint)
#' clusInd <- rep(1, m)
#' #'
#' \donttest{
#' evaluatePrior(
#'   m,
#'   p,
#'   muBar,
#'   hparam,
#'   thetaYList,
#'   ZOneDim,
#'   qVec,
#'   constraint,
#'   clusInd
#' )
#' }
evaluatePrior <- function(m,
                          p,
                          muBar,
                          hparam,
                          thetaYList,
                          ZOneDim,
                          qVec,
                          constraint,
                          clusInd) {
  ggamma <- hparam@ggamma
  alpha1 <- hparam@alpha1
  alpha2 <- hparam@alpha2
  bbeta <- hparam@bbeta
  delta <- hparam@delta

  loopm <- which(clusInd == 1)

  ## 2.2: Y
  # Yval = 0
  # for(k in loopm){
  #   Yval = Yval + sum(dnorm(thetaYList@Y[[k]], log = T))
  # }


  ## 2.3: Z
  # adjustTao = thetaYList@tao/sum(thetaYList@tao, na.rm = T)

  Zval <- log(sum(thetaYList@tao[ZOneDim], na.rm = T))
  # Zval = log(sum(adjustTao[ZOneDim], na.rm = T))
  ## 2.6: tao
  # taoVal = log(gtools::ddirichlet(x = adjustTao[!is.na(adjustTao)], alpha = rep(ggamma, m)))
  # print(thetaYList@tao[clusInd == 1])
  # print(rep(ggamma, m))

  taoVal <- log(gtools::ddirichlet(x = thetaYList@tao[clusInd == 1], alpha = rep(ggamma, m)))

  ## test
  # Zval = 0;taoVal = 0
  ##


  ## 2.7: M
  Mval <- 0
  for (k in loopm) {
    # cat("k =", k , "===> \n")
    Mval <- Mval + mvtnorm::dmvnorm(c(thetaYList@M[[k]]), mean = c(muBar), sigma = 1 / alpha1 * thetaYList@psy[[k]], log = T)
  }

  ## 2.8: lambda
  lambdaVal <- evaluatePriorLambda(
    p, m, alpha2, qVec, thetaYList@psy, thetaYList@lambda, constraint,
    clusInd
  )

  ## 2.9: psy
  psyVal <- evaluatePriorPsi(thetaYList@psy, p, m, delta, bbeta, constraint, clusInd)


  # print(Zval)
  # print(taoVal)
  # print(Mval)
  # print(lambdaVal)
  # print(psyVal)
  # totalVal = sum(Yval + Zval + taoVal + Mval + lambdaVal + psyVal)


  # print(c( Zval,taoVal,Mval,lambdaVal,psyVal))
  totalVal <- sum(Zval + taoVal + Mval + lambdaVal + psyVal)
  return(totalVal)
}
#' generatePriorPsi
#'
#' @description generate prior value for parameter Psi
#' @import stats
#' @param p the number of features
#' @param m the number of clusters
#' @param delta hyperparameters
#' @param bbeta hyperparameters
#' @param constraint the pgmm constraint, a vector of length three with binary entry. For example, c(1,1,1) means the fully constraint model
#' @export
#'
#' @examples
#' p <- 10
#' m <- 20
#' delta <- 2
#' bbeta <- 2
#' constraint <- c(0, 0, 0)
#'
#' \donttest{
#' generatePriorPsi(
#'   p,
#'   m,
#'   delta,
#'   bbeta,
#'   constraint
#' )
#' }
generatePriorPsi <- function(p,
                             m,
                             delta,
                             bbeta,
                             constraint) {
  psy <- list()

  if (constraint[2] == T & constraint[3] == T) {
    for (i in 1:m) {
      if (i == 1) {
        psyValue <- 1 / rgamma(1, shape = delta, rate = bbeta)
        psy[[i]] <- diag(rep(psyValue, p))
      } else {
        psy[[i]] <- psy[[1]]
      }
    }
  } else if (constraint[2] == T & constraint[3] == F) {
    for (i in 1:m) {
      if (i == 1) {
        psy[[i]] <- diag(1 / rgamma(p, shape = delta, rate = bbeta))
      } else {
        psy[[i]] <- psy[[1]]
      }
    }
  } else if (constraint[2] == F & constraint[3] == T) {
    for (i in 1:m) {
      psyValue <- 1 / rgamma(1, shape = delta, rate = bbeta)
      psy[[i]] <- diag(rep(psyValue, p))
    }
  } else {
    for (i in 1:m) {
      psy[[i]] <- diag(1 / rgamma(p, shape = delta, rate = bbeta))
    }
  }
  return(psy)
}

#' evaluatePriorPsi
#'
#' @description evaluate prior value for parameter Psi
#' @import stats
#' @param psy parameter
#' @param p the number of features
#' @param m the number of clusters
#' @param delta parameter
#' @param bbeta parameter
#' @param constraint parameter
#' @param clusInd cluster indicator vector
#' @examples
#' p <- 10
#' m <- 20
#' delta <- 2
#' bbeta <- 2
#' constraint <- c(0, 0, 0)
#' psy <- generatePriorPsi(
#'   p,
#'   m,
#'   delta,
#'   bbeta,
#'   constraint
#' )
#' clusInd <- rep(1, m)
#' #'
#' \donttest{
#' evaluatePriorPsi(
#'   psy,
#'   p,
#'   m,
#'   delta,
#'   bbeta,
#'   constraint,
#'   clusInd
#' )
#' }
evaluatePriorPsi <- function(psy,
                             p,
                             m,
                             delta,
                             bbeta,
                             constraint,
                             clusInd) {
  loopm <- which(clusInd == 1)
  psyeval <- 0
  if (constraint[2] == T & constraint[3] == T) {
    for (i in loopm) {
      if (i == 1) {
        psyValue <- 1 / psy[[i]][1, 1]
        psyeval <- psyeval + dgamma(psyValue, shape = delta, rate = bbeta, log = T)
      }
    }
  } else if (constraint[2] == T & constraint[3] == F) {
    for (i in loopm) {
      if (i == 1) {
        psyValue <- 1 / diag(psy[[i]])
        psyeval <- psyeval + sum(dgamma(psyValue, shape = delta, rate = bbeta, log = T))
      }
    }
  } else if (constraint[2] == F & constraint[3] == T) {
    for (i in loopm) {
      psyValue <- 1 / psy[[i]][1, 1]
      psyeval <- psyeval + dgamma(psyValue, shape = delta, rate = bbeta, log = T)
    }
  } else {
    for (i in loopm) {
      psyValue <- 1 / diag(psy[[i]])
      psyeval <- psyeval + sum(dgamma(psyValue, shape = delta, rate = bbeta, log = T))
    }
  }
  return(psyeval)
}


#' generatePriorLambda
#'
#' @description evaluate prior value for parameter Lambda
#' @importFrom mvtnorm rmvnorm
#' @param p the number of features
#' @param m the number of clusters
#' @param alpha2 hyper parameter
#' @param qVec parameter
#' @param psy parameter
#' @param constraint parameter
#' @export
#'
#' @examples
#' p <- 10
#' m <- 20
#' alpha2 <- 1.18
#' qVec <- rep(4, m)
#' delta <- 2
#' bbeta <- 2
#' constraint <- c(0, 0, 0)
#' psy <- generatePriorPsi(
#'   p,
#'   m,
#'   delta,
#'   bbeta,
#'   constraint
#' )
#' #'
#' \donttest{
#' generatePriorLambda(
#'   p,
#'   m,
#'   alpha2,
#'   qVec,
#'   psy,
#'   constraint
#' )
#' }
generatePriorLambda <- function(p,
                                m,
                                alpha2,
                                qVec,
                                psy,
                                constraint) {
  lambda <- list()
  if (constraint[1] == T & constraint[2] == T) {
    for (k in 1:m) {
      if (k == 1) {
        qk <- qVec[k]
        lambda[[k]] <- matrix(0, p, qk)
        for (j in 1:qk) {
          lambda[[k]][, j] <- mvtnorm::rmvnorm(1, rep(0, p), 1 / alpha2 * psy[[k]])
        }
      } else {
        lambda[[k]] <- lambda[[1]]
      }
    }
  } else if (constraint[1] == T & constraint[2] == F) {
    psyAve <- matrix(0, nrow = p, ncol = p)
    for (k in 1:m) {
      psyAve <- psyAve + solve(psy[[k]])
    }
    psyAve <- solve(1 / m * psyAve)

    for (k in 1:m) {
      if (k == 1) {
        qk <- qVec[k]
        lambda[[k]] <- matrix(0, p, qk)
        for (j in 1:qk) {
          lambda[[k]][, j] <- mvtnorm::rmvnorm(1, rep(0, p), 1 / alpha2 * psyAve)
        }
      } else {
        lambda[[k]] <- lambda[[1]]
      }
    }
  } else {
    for (k in 1:m) {
      qk <- qVec[k]
      lambda[[k]] <- matrix(0, p, qk)
      for (j in 1:qk) {
        lambda[[k]][, j] <- mvtnorm::rmvnorm(1, rep(0, p), 1 / alpha2 * psy[[k]])
      }
    }
  }
  return(lambda)
}

#' evaluatePriorLambda
#'
#' @description evaluate prior value for parameter Lambda
#' @importFrom mvtnorm dmvnorm
#' @param p the number of features
#' @param m the number of clusters
#' @param alpha2 hyper parameter
#' @param qVec the vector of the number of factors in each clusters
#' @param psy parameter
#' @param lambda parameter
#' @param constraint the pgmm constraint, a vector of length three with binary entry. For example, c(1,1,1) means the fully constraint model
#' @param clusInd cluster indicator vector
#' @examples
#' p <- 10
#' m <- 20
#' alpha2 <- 1.18
#' qVec <- rep(4, m)
#' delta <- 2
#' bbeta <- 2
#' constraint <- c(0, 0, 0)
#' psy <- generatePriorPsi(
#'   p,
#'   m,
#'   delta,
#'   bbeta,
#'   constraint
#' )
#' lambda <- generatePriorLambda(
#'   p,
#'   m,
#'   alpha2,
#'   qVec,
#'   psy,
#'   constraint
#' )
#' clusInd <- rep(1, m)
#' #'
#' \donttest{
#' evaluatePriorLambda(
#'   p,
#'   m,
#'   alpha2,
#'   qVec,
#'   psy,
#'   lambda,
#'   constraint,
#'   clusInd
#' )
#' }
evaluatePriorLambda <- function(p,
                                m,
                                alpha2,
                                qVec,
                                psy,
                                lambda,
                                constraint,
                                clusInd) {
  loopm <- which(clusInd == 1)
  evallambda <- 0
  if (constraint[1] == T & constraint[2] == T) {
    for (k in loopm) {
      if (k == 1) {
        # qk = qVec[k]
        qk <- ncol(lambda[[k]])
        for (j in 1:qk) {
          evallambda <- evallambda + mvtnorm::dmvnorm(lambda[[k]][, j], rep(0, p), 1 / alpha2 * psy[[k]], log = T)
        }
      }
    }
  } else if (constraint[1] == T & constraint[2] == F) {
    psyAve <- matrix(0, nrow = p, ncol = p)
    for (k in loopm) {
      psyAve <- psyAve + solve(psy[[k]])
    }
    psyAve <- solve(1 / m * psyAve)

    for (k in loopm) {
      if (k == 1) {
        # qk = qVec[k]
        qk <- ncol(lambda[[k]])
        for (j in 1:qk) {
          evallambda <- evallambda + mvtnorm::dmvnorm(lambda[[k]][, j], rep(0, p), 1 / alpha2 * psyAve, log = T)
        }
      }
    }
  } else {
    for (k in loopm) {
      # qk = qVec[k]
      qk <- ncol(lambda[[k]])
      for (j in 1:qk) {
        # print(k)
        # print("=====")
        # print(lambda[[k]][, j])
        # print(psy[[k]])
        evallambda <- evallambda + mvtnorm::dmvnorm(lambda[[k]][, j], rep(0, p), 1 / alpha2 * psy[[k]], log = T)
      }
    }
  }
  return(evallambda)
}
