#' @title Compute importance weights for  lasso, group lasso, scaled lasso or
#' scaled group lasso estimator under high-dimensional setting
#'
#' @description \code{hdIS} computes importance weights using samples
#' drawn by \code{\link{PBsampler}}. See the examples
#' below for details.
#'
#' @param PBsample bootstrap samples of class \code{PB} from \code{\link{PBsampler}}.
#' @param PETarget,sig2Target,lbdTarget parameters of target distribution.
#' (point estimate of beta or \code{E(y)}, estimated variance of error and lambda)
#' @param TsA.method method to construct \code{T(eta(s),A)} matrix. See Zhou and Min(2017)
#' for details.
#' @param log logical. If \code{log = TRUE}, importance weight is computed in log scale.
#' @param parallel logical. If \code{parallel = TRUE}, uses parallelization.
#' Default is \code{parallel = FALSE}.
#' @param ncores integer. The number of cores to use for parallelization.
#'
#' @details computes importance weights which is defined as (target density)/(proposal density),
#'  when the samples are drawn from the proposal
#'  distribution with the function \code{\link{PBsampler}} while the parameters of
#'  the target distribution are (PETarget, sig2Target, lbdTarget). \cr
#'  Say that we are interested in computing the expectation of a function of a random variable, \code{h(X)}.
#'  Let \code{f(x)} be the true or target distribution and \code{g(x)} be the proposal distribution.
#'  We can approximate the expectation, \code{E[h(X)]}, by a weighted average of samples, \code{x_i}, drawn from
#'  the proposal distribution as follows, \code{E[h(X)] = mean( h(x_i) * f(x_i)/h(x_i) )}.
#'
#' @references
#' Zhou, Q. (2014), "Monte Carlo simulation for Lasso-type problems by estimator augmentation,"
#' Journal of the American Statistical Association, 109, 1495-1516.
#'
#' Zhou, Q. and Min, S. (2017), "Estimator augmentation with applications in
#' high-dimensional group inference," Electronic Journal of Statistics, 11(2), 3039-3080.
#'
#' @return importance weights of the proposed samples.
#'
#' @examples
#' set.seed(1234)
#' n <- 10
#' p <- 30
#' Niter <-  10
#' Group <- rep(1:(p/10), each = 10)
#' Weights <- rep(1, p/10)
#' x <- matrix(rnorm(n*p), n)
#'
#' # Target distribution parameter
#' PETarget <- rep(0, p)
#' sig2Target <- .5
#' lbdTarget <- .37
#'
#' #
#' # Using non-mixture distribution
#' # ------------------------------
#' ## Proposal distribution parameter
#' PEProp1 <- rep(1, p)
#' sig2Prop1 <- .5
#' lbdProp1 <- 1

#' PB <- PBsampler(X = x, PE_1 = PEProp1, sig2_1 = sig2Prop1,
#'  lbd_1 = lbdProp1, weights = Weights, group = Group, niter = Niter,
#'  type="grlasso", PEtype = "coeff")
#'
#' hdIS(PB, PETarget = PETarget, sig2Target = sig2Target, lbdTarget = lbdTarget,
#'  log = TRUE)
#'
#' #
#' # Using mixture distribution
#' # ------------------------------
#' # Target distribution parameters (coeff, sig2, lbd) = (rep(0,p), .5, .37)
#' # Proposal distribution parameters
#' #  (coeff, sig2, lbd) = (rep(0,p), .5, .37) & (rep(1,p), 1, .5)
#' #
#' #
#' PEProp1 <- rep(0,p); PEProp2 <- rep(1,p)
#' sig2Prop1 <- .5; sig2Prop2 <- 1
#' lbdProp1 <- .37; lbdProp2 <- .5
#'
#' PBMixture <- PBsampler(X = x, PE_1 = PEProp1,
#'  sig2_1 = sig2Prop1, lbd_1 = lbdProp1, PE_2 = PEProp2,
#'  sig2_2 = sig2Prop2, lbd_2 = lbdProp2, weights = Weights, group = Group,
#'  niter = Niter, type = "grlasso", PEtype = "coeff")
#' hdIS(PBMixture, PETarget = PETarget, sig2Target = sig2Target, lbdTarget = lbdTarget,
#'  log = TRUE)
#' @export
hdIS <- function(PBsample, PETarget, sig2Target, lbdTarget,
            TsA.method = "default", log = TRUE, parallel = FALSE, ncores = 2L)
{

  if (class(PBsample) != "PB") {
    stop("Use EAlasso::PBsampler to generate Bootstrap samples")
  }

  if (any(missing(PETarget), missing(sig2Target), missing(lbdTarget))) {
    stop("provide all the parameters for the target distribution")
  }

  parallelTemp <- ErrorParallel(parallel,ncores)
  parallel <- parallelTemp[[1]]
  ncores <- parallelTemp[[2]]

  X <- PBsample$X
  n <- nrow(X)
  p <- ncol(X)

  if (n >= p) {
    stop("High dimensional setting is required, i.e. nrow(X) < ncol(X) is required.")
  }

  Mixture <- PBsample$mixture
  Btype <- PBsample$Btype
  type <- PBsample$type
  PEtype <- PBsample$PEtype
  method <- PBsample$method
  group <- PBsample$group
  weights <- PBsample$weights
  B <- PBsample$beta
  S <- PBsample$subgrad

  if (Mixture) {
    PEProp1 <- PBsample$PE[1,]
    PEProp2 <- PBsample$PE[2,]
    sig2Prop1 <- PBsample$sig2[1]
    sig2Prop2 <- PBsample$sig2[2]
    lbdProp1 <- PBsample$lbd[1]
    lbdProp2 <- PBsample$lbd[2]
  } else {
    PEProp1 <- PBsample$PE
    sig2Prop1 <- PBsample$sig2
    lbdProp1 <- lbdProp2 <- PBsample$lbd
  }

  if (Btype == "wild") {
    Y <- PBsample$Y
    if (PEtype == "coeff") {
      resTarget <- Y - X%*%PETarget
      resProp1 <- Y - X%*%PEProp1
      if (Mixture) {
        resProp2 <- Y - X%*%PEProp2
      }
    } else {
      resTarget <- Y - PETarget
      resProp1 <- Y - PEProp1
      if (Mixture) {
        resProp2 <- Y - PEProp2
      }
    }
    resTarget <- resTarget - mean(resTarget)
    resProp1 <- resProp1 - mean(resProp1)
    if (Mixture) {
      resProp2 <- resProp2 - mean(resProp2)
    }
  }

  if (type %in% c("slasso", "sgrlasso")) {
    hSigma <- PBsample$hsigma
  }
  niter <- nrow(B)


  if (type %in% c("lasso", "grlasso")) {
    if (all(group==1:p)) {
      #=========================================================================
      #-------------------
      # Lasso
      #-------------------
      #=========================================================================
      ## precalculation
      C <- t(X) %*% X / n
      egC <- eigen(C)
      V <- egC$vectors
      R <- 1:n
      N <- (n+1):p
      InvVarR     <- 1 / (egC$values[R] * sig2Target / n) #inverse of (sig2Target*Lambda_i/n)
      InvVarRprop1 <- 1 / (egC$values[R] * sig2Prop1 / n) #inverse of (sig2Prop1*Lambda_i/n)
      if (Mixture) {
        InvVarRprop2 <- 1 / (egC$values[R] * sig2Prop2 / n) #inverse of (sig2Prop1*Lambda_i/n)
      }
      VR <-matrix(V[, R], p, n)
      VRC <- VRCprop <- t(VR)%*%C
      W <- diag(weights)
      LBD <- LBDprop <- diag(egC$values[R])
      VRW <- VRWprop <- t(VR)%*%W
      if (PEtype == "coeff") {
        VRCB     <- t(VR) %*% C %*% PETarget
        VRCBprop1 <- t(VR) %*% C %*% PEProp1
        if (Mixture) {
          VRCBprop2 <- t(VR) %*% C %*% PEProp2
        }
      } else {
        VRCB     <- t(VR) %*% t(X) %*% PETarget / n
        VRCBprop1 <- t(VR) %*% t(X) %*% PEProp1 / n
        if (Mixture) {
          VRCBprop2 <- t(VR) %*% t(X) %*% PEProp2 / n
        }
      }

      if (Btype == "wild") {
        XVR <- X%*%VR
        VARTarget <- t(XVR) %*% diag(resTarget * sig2Target) %*% XVR / n^2
        VARProp1 <- t(XVR) %*% diag(resProp1 * sig2Prop1) %*% XVR / n^2
        if (Mixture) {
          VARProp2 <- t(XVR) %*% diag(resProp2 * sig2Prop2) %*% XVR / n^2
        }
      }

      Weight <- numeric(niter)
      FF <- function(x) {
        if (Btype == "gaussian") {
          log.f0 <- -0.5 * sum((VRC %*% B[x, ] + lbdTarget * VRW %*% S[x, ] -
                                  VRCB)^2 * InvVarR) + log(lbdTarget) *
            (n - sum(B[x,] != 0)) - 0.5 * ( n * log(sig2Target/n) )
          log.f1 <- -0.5 * sum((VRCprop %*% B[x, ] + lbdProp1 *
                                  VRWprop %*% S[x, ] - VRCBprop1)^2 * InvVarRprop1) + log(lbdProp1) *
            (n - sum(B[x,] != 0)) - 0.5 * ( n * log(sig2Prop1/n) )
          if (Mixture) {
            log.f2 <- -0.5 * sum((VRCprop %*% B[x, ] + lbdProp2 *
                                    VRWprop %*% S[x, ] - VRCBprop2)^2 * InvVarRprop2) + log(lbdProp2) *
              (n - sum(B[x,] != 0)) - 0.5 * ( n * log(sig2Prop2/n) )
          }
          # (n - sum(B[x,] != 0)) * log(lbdTarget / lbdProp1) -
          #   0.5 * sum((VRC %*% B[x, ] + lbdTarget * VRW %*% S[x, ] -
          #   VRCB)^2 * InvVarR) + 0.5 * sum((VRCprop %*% B[x, ] +
          #   lbdProp1 * VRWprop %*% S[x, ] - VRCBprop)^2 * InvVarRprop)
        } else {
          ##############################
          # wild multiplier bootstrap
          ##############################
          log.f0 <- dmvnorm(c(VRC %*% B[x, ] + lbdTarget * VRW %*% S[x, ] - VRCB),
                            mean=rep(0,n),sigma=VARTarget,log=TRUE) +
            log(lbdTarget) * (n - sum(B[x,] != 0))
          log.f1 <- dmvnorm(c(VRCprop %*% B[x, ] + lbdProp1 *
                                VRWprop %*% S[x, ] - VRCBprop1),
                            mean=rep(0,n),sigma=VARProp1,log=TRUE) +
            log(lbdProp1) * (n - sum(B[x,] != 0))
          if (Mixture) {
            log.f2 <- dmvnorm(c(VRCprop %*% B[x, ] + lbdProp2 *
                                  VRWprop %*% S[x, ] - VRCBprop2),
                              mean=rep(0,n),sigma=VARProp2,log=TRUE) +
              log(lbdProp2) * (n - sum(B[x,] != 0))
          }
        }
        if (!Mixture) {
          if (log) {Weight <- log.f0 - log.f1} else {
            Weight <- exp(log.f0 - log.f1)
          }
        } else {
          if (log) {
            Weight <- - log(exp(-log(2) + log.f1 - log.f0) +
                              exp(-log(2) + log.f2 - log.f0))
          } else {
            Weight <- 1/ (exp(-log(2) + log.f1 - log.f0) +
                            exp(-log(2) + log.f2 - log.f0))
          }
        }
        return(Weight)
      }
    } else {
      #=========================================================================
      #-------------------
      # Group Lasso
      #-------------------
      #=========================================================================
      if (!TsA.method %in% c("default", "qr")) {
        stop("TsA.method should be either \"default\" or \"qr\"")
      }

      #TsA.select <- switch(TsA.method, null = TsA.null, qr = TsA.qr, default = TsA)
      TsA.select <- switch(TsA.method, qr = TsA.qr, default = TsA)
      Psi <- t(X)%*%X / n
      ginv.tX <- solve(tcrossprod(X)) %*% X
      #precaculation
      W <- rep(weights, table(group))

      if (!all(lbdTarget == c(lbdProp1, lbdProp2)) && TsA.method != "null") {
        Q <- MASS::Null(t(X)/W)#round(Null(t(X)/W),5)
      }
      t.XWinv <- t(X)/W
      Weight <- numeric(niter)

      if (Btype == "gaussian") {
        f0sd <- sqrt(sig2Target/n)
        f1sd <- sqrt(sig2Prop1/n)
        if (Mixture) {
          f2sd <- sqrt(sig2Prop2/n)
        }
      } else {
        f0sd <- sqrt(sig2Target/n) * abs(resTarget)
        f1sd <- sqrt(sig2Prop1/n) * abs(resProp1)
        if (Mixture) {
          f2sd <- sqrt(sig2Prop2/n) * abs(resProp2)
        }
      }

      FF <- function(x) {
        Beta <- B[x,]
        Subgrad <- S[x,]
        if (n < p) {
          if (PEtype == "coeff") {
            H.tilde.target <- sqrt(n) * ginv.tX %*% (Psi %*% (Beta - PETarget) + lbdTarget * W * Subgrad) #H.tilde
            H.tilde.prop1 <- sqrt(n) * ginv.tX %*% (Psi %*% (Beta - PEProp1) + lbdProp1 * W * Subgrad) #H.tilde proposed1
            if (Mixture) {
              H.tilde.prop2 <- sqrt(n) * ginv.tX %*% (Psi %*% (Beta - PEProp2) + lbdProp2 * W * Subgrad) #H.tilde proposed2
            }
          } else {
            H.tilde.target <- sqrt(n) * ginv.tX %*% (Psi %*% Beta + lbdTarget * W * Subgrad) - PETarget / sqrt(n) #H.tilde
            H.tilde.prop1 <- sqrt(n) * ginv.tX %*% (Psi %*% Beta + lbdProp1 * W * Subgrad) - PEProp1 / sqrt(n)  #H.tilde proposed1
            if (Mixture) {
              H.tilde.prop2 <- sqrt(n) * ginv.tX %*% (Psi %*% Beta + lbdProp2 * W * Subgrad) - PEProp2 / sqrt(n)   #H.tilde proposed2
            }
          }

          r <- group.norm2(Beta, group)
          A <- unique(group[Beta != 0])

          if (!all(lbdTarget == c(lbdProp1, lbdProp2))) {
            #           if (TsA.method == "null") {
            #             TSA <- TsA.select(t.XWinv, Subgrad, group, A, n, p)
            #           } else {
            #             TSA <- TsA.select(Q, Subgrad, group, A, n, p)
            #           }
            TSA <- TsA.select(Q, Subgrad, group, A, n, p)

            log.f1 <- sum(dnorm(H.tilde.prop1, 0, f1sd, log = T)) +
              (logJacobiPartial(X, Subgrad, r, Psi, group, A, lbdProp1, W, TSA) )
            if (Mixture) {
              log.f2 <- sum(dnorm(H.tilde.prop2, 0, f2sd, log = T)) +
                (logJacobiPartial(X, Subgrad, r, Psi, group, A, lbdProp2, W, TSA) )
            }
            log.f0 <- sum(dnorm(H.tilde.target, 0, f0sd, log = T)) +
              (logJacobiPartial(X, Subgrad, r, Psi, group, A, lbdTarget, W, TSA) )

            if (!Mixture) {
              if (log) {Weight <- log.f0 - log.f1} else {
                Weight <- exp(log.f0 - log.f1)
              }
            } else {
              if (log) {
                Weight <- - log(exp(-log(2) + log.f1 - log.f0) + exp(-log(2) + log.f2 - log.f0))
              } else {
                Weight <- 1 / (exp(-log(2) + log.f1 - log.f0) + exp(-log(2) + log.f2 - log.f0))
              }
            }
          } else {
            log.f1 <- sum(dnorm(H.tilde.prop1, 0, f1sd, log = TRUE))
            if (Mixture) {log.f2 <- sum(dnorm(H.tilde.prop2, 0, f2sd, log = TRUE))}
            log.f0 <- sum(dnorm(H.tilde.target, 0, f0sd, log = TRUE))

            if (!Mixture) {
              if (log) {Weight <- log.f0 - log.f1} else {
                Weight <- exp(log.f0 - log.f1)
              }
            } else {
              if (log) {
                Weight <- - log(exp(-log(2) + log.f1 - log.f0) + exp(-log(2) + log.f2 - log.f0))
              } else {
                Weight <- 1/ (exp(-log(2) + log.f1 - log.f0) + exp(-log(2) + log.f2 - log.f0))
              }
            }
          }
        }
        # else {
        #   if (PEtype == "coeff") {
        #     H.target <- Psi %*% (Beta - PETarget) + lbdTarget * W * Subgrad #H.tilde
        #     H.prop1 <-  Psi %*% (Beta - PEProp1) + lbdProp1 * W * Subgrad #H.tilde proposed1
        #     if (Mixture) H.prop2 <-  Psi %*% (Beta - PEProp2) + lbdProp2 * W * Subgrad #H.tilde proposed2
        #   } else {
        #     H.target <- Psi %*% Beta + lbdTarget * W * Subgrad - t(X) %*% PETarget / n #H.tilde
        #     H.prop1 <-  Psi %*% Beta + lbdProp1 * W * Subgrad - t(X) %*% PEProp1 / n  #H.tilde proposed1
        #     if (Mixture) H.prop2 <-  Psi %*% Beta + lbdProp2 * W * Subgrad - t(X) %*% PEProp2 / n  #H.tilde proposed2
        #   }
        #
        #   r <- group.norm2(Beta, group)
        #   A <- unique(group[Beta != 0])
        #
        #   if (!all(lbdTarget == c(lbdProp1, lbdProp2))) {
        #
        #     if (TsA.method == "null") {
        #       TSA <- TsA.select(t.XWinv, Subgrad, group, A, n, p)
        #     } else {
        #       TSA <- TsA.select(Q, Subgrad, group, A, n, p)
        #     }
        #
        #     log.f1 <- dmvnorm(drop(H.prop1), , sig2Prop1/n * Psi, log = T) +
        #       (logJacobiPartial(X, Subgrad, r, Psi, group, A, lbdProp1, weights, TSA) )
        #     if (Mixture) {
        #       log.f2 <- dmvnorm(drop(H.prop2), , sig2Prop2/n * Psi, log = T) +
        #         (logJacobiPartial(X, Subgrad, r, Psi, group, A, lbdProp2, weights, TSA) )
        #     }
        #     log.f0 <- dmvnorm(drop(H.target), , sig2Target/n * Psi, log = T) +
        #       (logJacobiPartial(X, Subgrad, r, Psi, group, A, lbdTarget, weights, TSA) )
        #
        #     if (!Mixture) {
        #       if (log) {Weight <- log.f0 - log.f1} else {
        #         Weight <- exp(log.f0 - log.f1)
        #       }
        #     } else {
        #       if (log) {
        #         Weight <- - log(exp(-log(2) + log.f1 - log.f0) + exp(-log(2) + log.f2 - log.f0))
        #       } else {
        #         Weight <- 1/ (exp(-log(2) + log.f1 - log.f0) + exp(-log(2) + log.f2 - log.f0))
        #       }
        #     }
        #   } else {
        #     log.f1 <- sum(dmvnorm(drop(H.prop1), , sig2Prop1/n * Psi, log = TRUE))
        #     if (Mixture) {
        #       log.f2 <- sum(dmvnorm(drop(H.prop2), , sig2Prop2/n * Psi, log = TRUE))
        #     }
        #     log.f0 <- sum(dmvnorm(drop(H.target), , sig2Target/n * Psi, log = TRUE))
        #
        #     if (!Mixture) {
        #       if (log) {Weight <- log.f0 - log.f1} else {
        #         Weight <- exp(log.f0 - log.f1)
        #       }
        #     } else {
        #       if (log) {
        #         Weight <- - log(exp(-log(2) + log.f1 - log.f0) + exp(-log(2) + log.f2 - log.f0))
        #       } else {
        #         Weight <- 1/ (exp(-log(2) + log.f1 - log.f0) + exp(-log(2) + log.f2 - log.f0))
        #       }
        #     }
        #   }
        # }
        return(Weight)
      }
    }
  } else {
    #===========================================================================
    #-------------------
    # scaled lasso, scaled group lasso
    #-------------------
    #===========================================================================
    SVD <- svd(X)
    Psi <- t(X)%*%X / n
    ginv.tX <- solve(tcrossprod(X)) %*% X
    #precaculation
    W <- rep(weights, table(group))

    if (!all(lbdTarget == c(lbdProp1, lbdProp2))) {
      Q <- MASS::Null(t(X)/W)#round(Null(t(X)/W),5)
    }
    t.XWinv <- t(X)/W
    Weight <- numeric(niter)
    SVD.temp <- SVD$v %*% diag(1/SVD$d^2)%*%t(SVD$v)
    if (Btype == "gaussian") {
      f0sd <- sqrt(sig2Target/n)
      f1sd <- sqrt(sig2Prop1/n)
      if (Mixture) {
        f2sd <- sqrt(sig2Prop2/n)
      }
    } else {
      f0sd <- sqrt(sig2Target/n) * resTarget
      f1sd <- sqrt(sig2Prop1/n) * resProp1
      if (Mixture) {
        f2sd <- sqrt(sig2Prop2/n) * resProp2
      }
    }
    FF <- function(x) {
      Beta <- B[x,]
      Subgrad <- S[x,]
      hatSigma <- hSigma[x]
      if (PEtype == "coeff") {
        H.tilde.target <- sqrt(n) * ginv.tX %*% (Psi %*% (Beta - PETarget) +
                                                   lbdTarget * hatSigma * W * Subgrad) #H.tilde
        H.tilde.prop1 <- sqrt(n) * ginv.tX %*% (Psi %*% (Beta - PEProp1) +
                                                  lbdProp1 * hatSigma * W * Subgrad) #H.tilde proposed1
        if (Mixture) {
          H.tilde.prop2 <- sqrt(n) * ginv.tX %*% (Psi %*% (Beta - PEProp2) +
                                                    lbdProp2 * hatSigma * W * Subgrad) #H.tilde proposed2
        }
      } else {
        H.tilde.target <- sqrt(n) * ginv.tX %*% (Psi %*% Beta +
                                                   lbdTarget * hatSigma * W * Subgrad) - PETarget / sqrt(n) #H.tilde
        H.tilde.prop1 <- sqrt(n) * ginv.tX %*% (Psi %*% Beta +
                                                  lbdProp1 * hatSigma * W * Subgrad) - PEProp1 / sqrt(n)  #H.tilde proposed1
        if (Mixture) {
          H.tilde.prop2 <- sqrt(n) * ginv.tX %*% (Psi %*% Beta +
                                                    lbdProp2 * hatSigma * W * Subgrad) - PEProp2 / sqrt(n)   #H.tilde proposed2
        }
      }

      r <- group.norm2(Beta, group)
      A <- unique(group[Beta != 0])

      if (!all(lbdTarget == c(lbdProp1, lbdProp2))) { # Jacobian terms stay
        TSA <- TsA.slasso(SVD.temp = SVD.temp, Q = Q, s = Subgrad, W = W, group = group,
                          A = A, n = n, p = p)
        log.f1 <- sum(dnorm(H.tilde.prop1, 0, f1sd, log = TRUE)) +
          (logJacobiPartial.slasso(X = X, s = Subgrad, r = r, Psi = Psi, group = group, A = A, lam = lbdProp1, hsigma = hatSigma, W = W, TSA = TSA))
        if (Mixture) {
          log.f2 <- sum(dnorm(H.tilde.prop2, 0, f2sd, log = TRUE)) +
            (logJacobiPartial.slasso(X = X, s = Subgrad, r = r, Psi = Psi, group = group, A = A, lam = lbdProp2, hsigma = hatSigma, W = W, TSA = TSA))
        }
        log.f0 <- sum(dnorm(H.tilde.target, 0, f0sd, log = TRUE)) +
          (logJacobiPartial.slasso(X = X, s = Subgrad, r = r, Psi = Psi, group = group, A = A, lam = lbdTarget, hsigma = hatSigma, W = W, TSA = TSA))
      } else { # Jacobian terms cancel out
        log.f1 <- sum(dnorm(H.tilde.prop1, 0, f1sd, log = TRUE))
        if (Mixture) {log.f2 <- sum(dnorm(H.tilde.prop2, 0, f2sd, log = TRUE))}
        log.f0 <- sum(dnorm(H.tilde.target, 0, f0sd, log = TRUE))
      }
      if (!Mixture) {
        if (log) {Weight <- log.f0 - log.f1} else {
          Weight <- exp(log.f0 - log.f1)
        }
      } else {
        if (log) {
          Weight <- - log(exp(-log(2) + log.f1 - log.f0) + exp(-log(2) + log.f2 - log.f0))
        } else {
          Weight <- 1 / (exp(-log(2) + log.f1 - log.f0) + exp(-log(2) + log.f2 - log.f0))
        }
      }
      return(Weight)
    }
  }

  if (!parallel) {
    for (t in 1:niter) {
      Weight[t] <- FF(t)
    }
  } else {
    Weight <- parallel::mclapply(1:niter, FF, mc.cores = ncores)
    Weight <- do.call(c,Weight)
  }
  return(Weight)
}
