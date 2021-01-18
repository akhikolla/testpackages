#' @title Metropolis-Hastings lasso sampler under a fixed active set.
#'
#' @description Metropolis-Hastings sampler to simulate from the sampling
#' distribution of lasso given a fixed active set.
#'
#' @param X predictor matrix.
#' @param PE,sig2,lbd parameters of target distribution.
#' (point estimate of beta or \code{E(y)} depends on \code{PEtype}, variance estimate of error and lambda).
#' @param weights weight vector with length \code{p}(the number of covariates).
#' Default is \code{weights = rep(1, p)}.
#' @param B0 numeric vector with length \code{p}.
#' Initial value of lasso estimator.
#' @param S0 numeric vector with length \code{p}.
#' Initial value of subgradients.
#' If not given, this will be generated in a default way.
#' @param A numeric vector. Active coefficient index.
#' Every active coefficient index in \code{B0} must be included.
#' Default is \code{A = which(B0 != 0)}.
#' @param tau numeric vector with length \code{p}.
#' Standard deviation of proposal distribution for each coefficient.
#' @param niter integer. The number of iterations. Default is \code{niter = 2000}
#' @param burnin integer. The length of burin-in periods. Default is \code{burnin = 0}
#' @param updateS.itv integer. Update subgradients every \code{updateS.itv} iterations.
#' Set this value larger than \code{niter} if one wants to skip updating subgradients.
#' @param PEtype Type of \code{PE} which is needed to characterize the target distribution.
#' Users can choose either \code{"coeff"} or \code{"mu"}.
#' @param verbose logical. If true, print out the progress step.
#' @param ... complementary arguments.
#' \itemize{
#'  \item{\code{FlipSA :}}{ optional parameter.
#' This has to be a subset of active set, A. If the index is not listed in FlipSA,
#' the sign of coefficients which correspond to the listed index will remain fixed.
#' The default is \code{FlipSA=A}}
#'  \item{\code{SFindex :}}{ optional parameter. subgradient index for the free coordinate.}
#'  \item{\code{randomSFindex :}}{ logical. If \code{true}, resample \code{SFindex} every
#' \code{updateSF.itv} iterations.}
#'  \item{\code{updateSF.itv :}}{ integer. In every \code{updateSF.itv} iterations,
#'  randomize \code{SFindex}.}
#' }
#'
#' @details Given appropriate initial value, provides Metropolis-Hastings samples
#'  under the fixed active set. \cr
#' From the initial values, \code{B0} and {S0}, \code{\link{MHLS}} draws \code{beta} and \code{subgrad} samples.
#' In every iteration, given \code{t}-th iteration values, \code{t}-th \code{beta} and \code{t}-th \code{subgrad},
#' a new set of proposed beta and subgradient is sampled. We either accept the proposed sample
#' and use that as \code{(t+1)}-th iteration values or reuse \code{t}-th iteration values. \cr
#' See Zhou(2014) for more details.
#'
#' @return \code{\link{MHLS}} returns an object of class \code{"MHLS"}.
#' The functions \code{\link{summary.MHLS}} and \code{\link{plot.MHLS}}
#' provide a brief summary and generate plots.
#' \item{beta}{lasso samples.}
#' \item{subgrad}{subgradient samples.}
#' \item{acceptHistory}{numbers of acceptance and proposal.}
#' \item{niter, burnin, PE, type}{same as function arguments.}
#'
#' @references
#' Zhou, Q. (2014), "Monte Carlo simulation for Lasso-type problems by estimator augmentation,"
#' Journal of the American Statistical Association, 109, 1495-1516.
#'
#' @examples
#' #-------------------------
#' # Low dim
#' #-------------------------
#' set.seed(123)
#' n <- 10
#' p <- 5
#' X <- matrix(rnorm(n * p), n)
#' Y <- X %*% rep(1, p) + rnorm(n)
#' sigma2 <- 1
#' lbd <- .37
#' weights <- rep(1, p)
#' LassoResult <- lassoFit(X = X, Y = Y, lbd = lbd, type = "lasso", weights = weights)
#' B0 <- LassoResult$B0
#' S0 <- LassoResult$S0
#' MHLS(X = X, PE = rep(0, p), sig2 = 1, lbd = 1,
#'      weights = weights, B0 = B0, S0 = S0, niter = 50, burnin = 0,
#'      PEtype = "coeff")
#' MHLS(X = X, PE = rep(0, n), sig2 = 1, lbd = 1,
#'      weights = weights, B0 = B0, S0 = S0, niter = 50, burnin = 0,
#'      PEtype = "mu")
#'
#' #-------------------------
#' # High dim
#' #-------------------------
#' set.seed(123)
#' n <- 5
#' p <- 10
#' X <- matrix(rnorm(n*p),n)
#' Y <- X %*% rep(1,p) + rnorm(n)
#' weights <- rep(1,p)
#' LassoResult <- lassoFit(X = X,Y = Y,lbd = lbd, type = "lasso", weights = weights)
#' B0 <- LassoResult$B0
#' S0 <- LassoResult$S0
#' MHLS(X = X, PE = rep(0, p), sig2 = 1, lbd = 1,
#'      weights = weights, B0 = B0, S0 = S0, niter = 50, burnin = 0,
#'      PEtype = "coeff")
#' MHLS(X = X, PE = rep(0, n), sig2 = 1, lbd = 1,
#'      weights = weights, B0 = B0, S0 = S0, niter = 50, burnin = 0,
#'      PEtype = "mu")
#' @export
MHLS <-  function(X, PE, sig2, lbd,
                   weights = rep(1, ncol(X)), B0, S0, A = which(B0 != 0),
                   tau = rep(1, ncol(X)), niter = 2000, burnin = 0,
                   PEtype = "coeff", updateS.itv = 1, verbose = FALSE, ...)
{
  MHLSmain(X = X, PE = PE, sig2 = sig2, lbd = lbd,
      weights = weights, B0 = B0, S0 = S0, A = A,
      tau = tau, niter = niter, burnin = burnin, PEtype = PEtype,
      updateS.itv = updateS.itv, verbose = verbose, ...)
}

MHLSmain <- function (X, PE, sig2, lbd,
  weights, B0, S0, A, tau, niter, burnin, PEtype, updateS.itv, verbose, ...)
{
  #------------------
  # Error handling
  #------------------
  n <- nrow(X)
  p <- ncol(X)
  if (PEtype == "coeff" && length(PE) != p) {
    stop("length(PE) must be the same with ncol(X), if PEtype = \"coeff\"")
  }
  if (PEtype == "mu" && length(PE) != n) {
    stop("length(PE) must be the same with nrow(X), if PEtype = \"mu\"")
  }
  if (length(B0) != p || (!missing(S0) && length(S0) != p)) {
    stop("length(B0) and/or length(S0) has to be the same with ncol(X)")
  }
  if (n < p && length(A) > n) {
    stop("Invalid active set index, A. Cannot be larger than min(nrow(X),ncol(X)).")
  }
  if (length(A) != length(unique(A))) {
    stop("Invalid active set index, A.")
  }
  if (!PEtype %in% c("coeff", "mu")) {
    stop("Invalide PEtype input.")
  }
  if (length(weights) != p) {
    stop("length(weights) has to be the same with the number of coefficients")
  }
  if (any(weights <= 0)) {
    stop("weights should be positive.")
  }
  if (sig2 <=0 || lbd <= 0) {
    stop("sig2 and/or lbd have to be positive.")
  }
  # if (!all(group == 1:ncol(X)) && missing(S0)) {
  #   stop("Missing S0. Use LassoMHLS for a good initial value.")
  # }
  if (any(missing(PE), missing(sig2), missing(lbd))) {
    stop("provide all the parameters for the distribution")
  }
  if (burnin >= niter) {
    stop("burnin has to be greater than niter")
  }
  if (niter <= 1) {
    stop("niter should be a integer greater than 1.")
  }

  est <- MHLSswp(X = X, PE = PE, sig2 = sig2,
                 lbd = lbd, weights = weights, B0 = B0, S0 = S0, A = A,
                 tau = tau, niter = niter, burnin = burnin, PEtype = PEtype,
                 updateS.itv = updateS.itv, verbose = verbose, ...)
  class(est) <- "MHLS"
  est$niter <- niter
  est$burnin <- burnin
  est$PE <- PE
  est$PEtype <- PEtype
  return(est)
}

MHLSswp <- function(X, PE, sig2, lbd, weights,
  B0, S0, A, tau, niter,
  burnin, PEtype, FlipSA = A, SFindex,
  randomSFindex = TRUE, updateSF.itv = max(round(niter/20), 50), updateS.itv,
  verbose, ...)
{
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  NN <- floor((1:10) * niter / 10) # verbose index

  A <- unique(A) # active set
  Ac <- setdiff(1:p,A) # inactive set
  nA <- length(A) # size of the active set
  nI <- length(Ac) # size of the inactive set
  C <- crossprod(X) / n #Gram matrix

  if (!all(which(B0 != 0) %in% A)) {
    stop("Invalid active set index, A. The active set, A, has to include every index of nonzero B0.")
  }
  if (any(!A %in% 1:p)) {
    stop("Invalid active set index, A. The active set, A, has to be a subset of 1:ncol(X)")
  }
  if (!missing(S0) && !all(round(S0[which(B0 != 0)], 3) == sign(B0[B0 != 0]))) {
    stop("Invalid S0. Leave S0 blank, if S0 is unknown.")
  }
  if (length(tau) != ncol(X)) {
    stop("tau must have a same length with the active set, A.")
  }
  if (n >= p) {   # Low-dim MH
    #precalculation
    #for(j in 1:p){X[,j]=X[,j]-mean(X[,j])}
    #If X to be centered, we need to re-compute B0 and S0 using centered X.
    Cinv <- solve(C) #Inverse Gram matrix
    logdiagC <- log(diag(C))
    Vinv <- n / (2 * sig2) * Cinv # non exponential part of pdf of U
    if (PEtype == "coeff") {
      CB <- C %*% PE    # PE : True beta
    } else {
      CB <- crossprod(X, PE) / n    # PE : True beta
    }

    lbdwgt <- lbd * weights
    loglbd <- log(lbdwgt)

    #initialization
    B <- matrix(0, niter, p) #Store Coefficient
    S <- matrix(0, niter, p) #Store Sub-grad
    B[1,] <- B0 #B0 initial value of beta

    if (missing(S0)) {
      S[1, A] <- sign(B[1, A])
      S[1, Ac] <- runif(length(Ac), -1, 1)
    } else {
      S[1, ] <- S0
    }

    Ucur <- C %*% B[1, ] + lbdwgt * S[1, ] - CB #Current U
    ldetDRatio <- 0
    if (nA >= 1) {negCAAinv <- -solve(C[A, A])} else {negCAAinv <- NULL}

    nAccept <- numeric(2)
    nProp <- (niter - burnin) * c(nA, nI)

    if (length(setdiff(FlipSA, A))!=0)
      stop("FlipSA has to be a subset of active set, A.")
    A2 <- setdiff(A,FlipSA)
    # if (any(B0[A2]==0))
    #   stop("To fix the sign of beta_j, use non-zero B0_j.")

    if (length(A2) != 0) {
      LUbounds <- matrix(0, p, 2);
      LUbounds[B0 < 0, 1] <- -Inf;
      LUbounds[B0 > 0, 2] <- Inf;
    }


    for(t in 2:niter)
    {
      if(nA >= 1){
        if (length(FlipSA)!=0) {
          for (j in FlipSA) {
            b_prop <- rnorm(1, mean = B[t - 1, j], sd = tau[j])
            s_prop <- sign(b_prop)
            DiffU <- (b_prop - B[t - 1, j]) * C[, j]
            DiffU[j] <- DiffU[j] + lbdwgt[j] * (s_prop - S[t - 1, j])
            Uprop <- Ucur + DiffU
            logMH <- -t(Ucur + Uprop) %*% Vinv %*% DiffU
            u <- runif(1)
            if(log(u) < logMH)
            {
              B[t, j] <- b_prop
              S[t, j] <- s_prop
              Ucur <- Uprop
              if (t > burnin) {nAccept[1] <- nAccept[1] + 1}
              #nAccept[2]=nAccept[2]+1
            }else{
              B[t, j] <- B[t - 1, j]
              S[t, j] <- S[t - 1, j]
            }
          }
        }

        if (length(A2)!=0) {
          for (j in A2) {
            b_prop <- rtnorm(1, mean = B[t - 1, j], sd = tau[j],
                             lower = LUbounds[j, 1],
                             upper = LUbounds[j, 2])
            Ccur <- pnorm(0,mean=B[t-1, j],sd=tau[j],lower.tail=(B[t-1,j]<0),log.p=FALSE);
            Cnew <- pnorm(0,mean=b_prop,sd=tau[j],lower.tail=(b_prop<0),log.p=FALSE);
            lqratio=log(Ccur/Cnew);


            DiffU <- (b_prop - B[t - 1, j]) * C[, j]
            Uprop <- Ucur + DiffU
            logMH <- -t(Ucur + Uprop) %*% Vinv %*% DiffU + lqratio
            u <- runif(1)
            if(log(u) < logMH)
            {
              B[t, j] <- b_prop
              S[t, j] <- s_prop
              Ucur <- Uprop
              if (t > burnin) {nAccept[1] <- nAccept[1] + 1}
              #nAccept[2]=nAccept[2]+1
            }else{
              B[t, j] <- B[t - 1, j]
              S[t, j] <- S[t - 1, j]
            }
          }
        }

      }
      if(nI >= 1 && (t %% updateS.itv == 0)){
        for(j in Ac)
        {
          s_prop <- runif(1, -1, 1)
          diffu <- lbdwgt[j] * (s_prop - S[t - 1, j])
          Uprop <- Ucur
          Uprop[j] <- Ucur[j] + diffu
          logMH <- -t(Ucur + Uprop) %*% Vinv[, j] * diffu
          u <- runif(1)
          if(log(u) < logMH)
          {
            S[t,j] <- s_prop
            Ucur <- Uprop
            if (t > burnin) {nAccept[2] <- nAccept[2] + 1}
            #nAccept[3]=nAccept[3]+1
          } else {
            S[t, j] <- S[t - 1, j]
          }
        }
      } else {
        S[t, Ac] <- S[t - 1, Ac]
      }
      if (verbose && sum(t == NN)==1) {
        aa <- which(NN==t)
        cat(paste("Updating : ", aa * 10,"%" ,sep = ""), "\n")
      }
    }
    #nAccept=nAccept/c((niter-1)*selectsize,nProp)
    #nAccept=nAccept/nProp
  }
  if (n < p) {
    #precalculation---------------------------
    #for (j in 1:p) {X[,j]=X[,j]-mean(X[,j])}
    #If X to be centered, we need to recompute B0 and S0 using centered X.
    C <- t(X) %*% X / n
    egC <- eigen(C)
    V <- egC$vectors
    R <- 1:n
    N <- (n + 1):p
    InvVarR <- 1 / (egC$values[R] * sig2 / n) #inverse of (sig2*Lambda_i/n)
    VR <- matrix(V[,R], p, n)
    VRC <- t(VR) %*% C
    W <- diag(weights)
    LBD <- diag(egC$values[R])
    lbdVRW <- lbd * t(VR) %*% W
    if (PEtype == "coeff") {
      VRCB <- t(VR) %*% C %*% PE
    } else {
      VRCB <- t(VR) %*% crossprod(X, PE) / n
    }


    # if (is.missing(B0)) {
    #   if (signBA == NULL)
    #     stop("If initial value of 'B0' is not given, 'signBA' has to be given.")
    #   B0 = rep(0,p)
    #   B0[A] = abs(rnorm(nA)) * signBA
    # }

    if (length(setdiff(FlipSA, A))!=0)
      stop("FlipSA has to be a subset of active set, A.")
    A2 <- setdiff(A,FlipSA)

    V_AN <- V[A,N,drop=FALSE]
    V_IN <- V[Ac,N,drop=FALSE]
    V_AR <- V[A,R,drop=FALSE]
    V_IR <- V[Ac,R,drop=FALSE]

    BB <- t(V_IN) %*% W[Ac,Ac]
    tVAN.WA <- t(V_AN)%*%W[A,A]

    if (!missing(S0) && !all.equal(t(V[,N])%*%S0,matrix(0,length(N),1))) {
      warning("Invalid S0. Regenerate S0 with a default way")
      S0 <- NULL
    }

    if (missing(S0) || is.null(S0)) {
      E <- BB
      H <- rep(-1,2*nI)
      G <- rbind(diag(rep(1,nI)),diag(rep(-1,nI)))
      F1 <- -tVAN.WA%*%sign(B0[A])
      S0.prop <- limSolve::lsei(G=G,H=H,E=E,F=F1)
      if (S0.prop$IsError) {
        stop("There exist no solution for the given 'B0'.")
      } else {
        S0 <- rep(0,p)
        S0[A] <- sign(B0[A])
        S0[-A] <- S0.prop$X
      }
    }

    if (n-nA != 0) {
      # SFindex : S_F, free coordinate, index
      if (missing(SFindex)) {SFindex <- 1:(n-nA)}

      if (length(SFindex) != n-nA) {
        warning("Length of SFindex has to be same as 'n-length(A)'. Automatically set 'SFindex <- 1:(n-nA)'.")
        SFindex <- 1:(n-nA)
      }

      B_F <- BB[, SFindex, drop=FALSE]
      B_D <- BB[,-SFindex, drop=FALSE]
      invB_D <- solve(B_D)
      BDF <- invB_D%*%B_F
    } else if (n-nA == 0) {
      if (missing(SFindex)) {SFindex <- NULL}
      if (!is.null(SFindex)) {
        warning("If 'n-nA == 0', SFindex has to be set to NULL. Automatically set 'SFindex <- NULL'.")
        SFindex <- NULL
      }
      BDF <- invB_D <- solve(BB)  # NEED TO CHECK
      B_F <- NULL
    }

    #initialization--------------------------
    B <- matrix(0,niter,p)
    S <- matrix(0,niter,p)
    B[1,] <- B0
    S[1,] <- S0

    Rcur <- VRC[ , A, drop=FALSE] %*% t(B[1, A, drop=FALSE]) + lbdVRW %*% S[1, ] - VRCB
    logfRcur <- -0.5 * sum(Rcur^2 * InvVarR)
    Tcur <- CalTmat(p,n,V,LBD,W,lbd,R,N,A,Ac)
    logdensitycur <- logfRcur+Tcur$logdetT

    #record acceptance rates
    nAccept <- numeric(2)
    #nProp=numeric(2)
    nProp <- c(nA*(niter-max(1,burnin)),(n-nA)*(niter-max(1,burnin)))
    #Change sign count
    nSignChange <- numeric(3)
    for(t in 2:niter )
    {
      #P1: update b_A
      if(nA>=1)
      {
        S[t,] <- S[t-1,]
        if (length(FlipSA)!=0) {
          MoveBA <- UpdateBA(B[t-1,],S[t-1,],tau,A,Ac,Rcur,logfRcur,VRC,lbdVRW,InvVarR,
                          tVAN.WA,invB_D,B_F,FlipSA,SFindex)
          B[t,] <- MoveBA$B
          S[t,] <- MoveBA$S
          nSignChange <- nSignChange+MoveBA$nChangeSign
          Rcur <- MoveBA$Rvec
          logfRcur <- MoveBA$logf
          nAccept[1] <- nAccept[1]+MoveBA$nAccept
          #nProp[1]=nProp[1]+nA
        }

        if (length(A2)!=0) {
          MoveBA <- UpdateBA.fixedSA(B[t-1,],tau,A2,Rcur,logfRcur,VRC,InvVarR)
          B[t,A2] <- MoveBA$B[A2]
          #S[t,]=S[t-1,]
          Rcur <- MoveBA$Rvec
          logfRcur <- MoveBA$logf
          nAccept[1] <- nAccept[1] + MoveBA$nAccept
        }
      }else{ B[t,] <- B[t-1,]; S[t,] <- S[t-1,]}

      # P2: update S_I
      if(nA<n && (t %% updateS.itv == 0))
      {
        MoveSI <- UpdateSI(S[t,],A,Ac,Rcur,n,p,logfRcur,lbdVRW,InvVarR,
                        tVAN.WA,invB_D,BDF,B_F,SFindex,...)
        S[t,] <- MoveSI$S
        Rcur <- MoveSI$Rvec
        logfRcur <- MoveSI$logf
        nAccept[2] <- nAccept[2]+MoveSI$nAccept
        #nProp[2]=nProp[2]+(n-nA)
      }#else{S[t,]=S[t-1,];S[t,A]=sign(B[t,A])}

      if (!is.null(SFindex) && randomSFindex && (t%%updateSF.itv==0) ) {
        SFindex <- sort(sample(1:(p-nA),n-nA))
        if (verbose) {cat("New SFindex : [",paste(SFindex,collapse=", "),"]\n")}
        B_F <- BB[,SFindex,drop=FALSE]
        B_D <- BB[,-SFindex,drop=FALSE]
        invB_D <- solve(B_D)
        BDF <- invB_D%*%B_F
      }

      if (verbose && sum(t == NN)==1) {
        aa <- which(NN==t)
        cat(paste("Updating : ", aa*10  ,"%",sep=""),"\n")
      }
    }
  }
  return(list(beta = B[if (burnin != 0){-c(1:burnin)} else {1:niter}, ],
              subgrad = S[if (burnin != 0){-c(1:burnin)} else {1:niter}, ],
              acceptHistory = rbind(nAccept, nProp)))
  #
  # return(list(beta = B[if (burnin != 0){-c(1:burnin)} else {1:niter}, ],
  #             subgrad = S[if (burnin != 0){-c(1:burnin)} else {1:niter}, ],
  #             PE = PE,
  #             signchange=nSignChange,
  #             acceptHistory = rbind(nAccept, nProp)))
}

# MHLSgroup <- function(X, PE, sig2, lbd,
#  weights, group, B0, S0, A, tau, niter, burnin, PEtype = "coeff", updateS.itv, verbose)
# {
#   if ( all.equal(group.norm2(S0, group)[A], rep(1, length(A)), tolerance = 1e-04) != TRUE ) {
#     stop("Invalid S0. Use LassoMHLS for a good initial value.")
#   }
#   n <- nrow(X)
#   p <- ncol(X)
#
#   K <- 10
#   W <- rep(weights,table(group))
#   Psi <- 1/n * crossprod(X)
#   if (n > p) {
#     inv.Var <- n/sig2 * solve(Psi)
#   }
#   nA <- length(A)
#
#   r.seq <- matrix(, niter, max(group))
#   S.seq <- matrix(, niter ,p)
#   nAccept <- numeric(2)
#   nProp <- c(nA*(niter-max(1,burnin)), max(group)*(niter-max(1,burnin)))
#
#   rcur <- group.norm2(B0,group)
#   r.seq[1, ] <- rcur
#   S.seq[1, ] <- Scur <- S0
#
#   if (PEtype == "coeff") {
#     Hcur <- drop(Psi %*% drop(B0 - PE) + lbd * W * drop(S0))
#   } else {
#     Hcur <- drop(Psi %*% drop(B0) - t(X) %*% PE / n + lbd * W * drop(S0))
#   }
#
#   if (n >= p) {
#     for (i in 2:niter) {
#       r.new <- ld.Update.r(rcur,Scur,A,Hcur,X,PE,Psi,W,lbd,group,inv.Var,tau,PEtype,n,p)
#       r.seq[i,] <- rcur <- r.new$r
#       Hcur <- r.new$Hcur
#       if (i > burnin) {nAccept[1] <- nAccept[1] + r.new$nrUpdate}
#
#       if (i %% updateS.itv == 0) {
#         S.new <- ld.Update.S (rcur,Scur,A,Hcur,X,PE,Psi,W,lbd,group,inv.Var,PEtype,n,p)
#         S.seq[i,] <- Scur <- S.new$S
#         Hcur <- S.new$Hcur
#       } else {
#         S.seq[i,] <- Scur
#       }
#       if (i > burnin) {nAccept[2] <- nAccept[2] + S.new$nSUpdate}
#
#       if (verbose && (i %% round(niter/10) == 0)) {
#         cat("MCMC step,", K, "% Finished\n")
#         K <- K+10
#       }
#     }
#   # } else {
#   #   for (i in 2:niter) {
#   #     r.new <- hd.Update.r(rcur,Scur,A,Hcur,X,PE,Psi,W,lbd,group,inv.Var,1)
#   #     r.seq[i,] <- rcur <- r.new$r
#   #     Hcur <- r.new$Hcur
#   #     nAccept[1] <- nAccept[1] + r.new$nrUpdate
#   #
#   #     S.new <- hd.Update.S (rcur,Scur,A,Hcur,X,PE,Psi,W,lbd,group,inv.Var,p)
#   #     S.seq[i,] <- Scur <- S.new$S
#   #     Hcur <- S.new$Hcur
#   #     nAccept[2] <- nAccept[2] + S.new$nSUpdate
#   #
#   #     if (verbose && (i %% round(niter/10) == 0)) {
#   #       cat("MCMC step,", K, "% Finished\n")
#   #       K <- K+10
#   #     }
#   #   }
#   }
#   return(list(
#     group.l2.norm = r.seq[if (burnin != 0){-c(1:burnin)} else {1:niter}, ],
#     subgrad = S.seq[if (burnin != 0){-c(1:burnin)} else {1:niter}, ],
#     acceptHistory = rbind(nAccept, nProp)))
# }

#' @method print MHLS
#' @title Print Metropolis-Hastings sampler outputs
#'
#' @description Print a brief summary of the MH sampler outputs.
#'
#' @param x an object of class "MHLS", which is a result of \code{\link{MHLS}}.
#' @param ...	additional print arguments.
#' @details
#' \code{\link{print.MHLS}} prints out last 10 iterations and a brief summary
#' of the simulation; number of iterations, number of burn-in periods, PE, PEtype and
#' acceptance rate.
#'
#' @return Above results are silently returned.
#'
#' @examples
#' set.seed(123)
#' n <- 10
#' p <- 5
#' X <- matrix(rnorm(n * p), n)
#' Y <- X %*% rep(1, p) + rnorm(n)
#' sigma2 <- 1
#' lbd <- .37
#' weights <- rep(1, p)
#' LassoResult <- lassoFit(X = X, Y = Y, lbd = lbd, type="lasso", weights = weights)
#' B0 <- LassoResult$B0
#' S0 <- LassoResult$S0
#' Result <- MHLS(X = X, PE = rep(0, p), sig2 = sigma2, lbd = lbd, group = 1:p,
#'      weights = weights, B0 = B0, S0 = S0, niter = 50, burnin = 0,
#'      type = "coeff")
#' print(Result)
#' @export
print.MHLS <- function (x, ...) {
  cat ("===========================\n")
  cat ("Number of iteration: ", x$niter,"\n\n")
  cat ("Burn-in period: ", x$burnin,"\n\n")
  cat ("Plug-in PE: \n")
  print(x$PE)
  cat ("PEtype: \n")
  print(x$PEtype)

  # if (inherits(x,"Group")) {
  #   Group.matrix <- matrix(0, length(unique(x$group)), p)
  #   Group.matrix[cbind(x$group,1:p)] <- 1
  #   beta <- (x$group.l2.norm %*% Group.matrix) * x$subgrad
  # }
  cat ("\nLast 10 steps of beta's:\n")

  # if (inherits(x,"GroupLasso")) {
  #   if (x$niter-x$burnin <= 9) {
  #     print(x$group.l2.norm)
  #   } else {
  #     print(x$group.l2.norm[(x$niter-x$burnin-9):(x$niter-x$burnin),])
  #   }
  # }

  if (x$niter-x$burnin <= 9) {
    # if (inherits(x,"Group")) {
    #   #print(x$group.l2.norm)
    #   print(beta)
    # } else {
    print(x$beta)
    # }
  } else {
    # if (inherits(x,"Group")) {
    #   print(beta[(x$niter-x$burnin-9):(x$niter-x$burnin),])
    #   #print(x$group.l2.norm[(x$niter-x$burnin-9):(x$niter-x$burnin),])
    # } else {
    print(x$beta[(x$niter-x$burnin-9):(x$niter-x$burnin),])
    # }
  }

  cat ("\nlast 10 steps of subgradients:\n")
  if (x$niter-x$burnin <= 9) {
    print(x$subgrad)
  } else {
    print(x$subgrad[(x$niter-x$burnin-9):(x$niter-x$burnin),])
  }

  cat ("\nAcceptance rate:\n")
  cat("-----------------------------\n")
  cat("\t \t \t beta \t subgrad\n")
  # if (inherits(x,"GroupLasso")) {
  #   cat("\t \t l_2 group norm\t subgrad\n")
  # } else {
  #   cat("\t \t \t beta \t subgrad\n")
  # }
  cat("# Accepted\t : \t", paste(x$acceptHistory[1,],"\t"),"\n")
  cat("# Moved\t\t : \t", paste(x$acceptHistory[2,],"\t"),"\n")
  cat("Acceptance rate\t : \t", paste(round(x$acceptHistory[1,]/x$acceptHistory[2,],3),"\t"),"\n")
  # cat ("\nSignChange rate:\n")
  # cat("-----------------------------\n")
  # cat("# Accepted\t : \t", paste(x$signchange[1],"\t"),"\n")
  # cat("# Moved\t\t : \t", paste(x$signchange[2],"\t"),"\n")
  # cat("# Cdt Accept \t : \t", paste(x$signchange[3],"\t"),"\n")
  # cat("Acceptance rate\t : \t", paste(round(x$signchange[1]/x$signchange[2],3),"\t"),"\n")
}

#' @method summary MHLS
#' @title Summarizing Metropolis-Hastings sampler outputs
#'
#' @description Summary method for class "MHLS".
#'
#' @param object an object of class "MHLS", which is a result of \code{\link{MHLS}}.
#' @param ... additional arguments affecting the summary produced.
#'
#' @details
#' This function provides a summary of each sampled beta and subgradient.
#' @return mean, median, standard deviation, 2.5\% quantile and 97.5\% quantile
#' for each beta and its subgradient.
#' @examples
#' #' set.seed(123)
#' n <- 10
#' p <- 5
#' X <- matrix(rnorm(n * p), n)
#' Y <- X %*% rep(1, p) + rnorm(n)
#' sigma2 <- 1
#' lbd <- .37
#' weights <- rep(1, p)
#' LassoResult <- lassoFit(X = X, Y = Y, lbd = lbd, type = "lasso", weights = weights)
#' B0 <- LassoResult$B0
#' S0 <- LassoResult$S0
#' summary(MHLS(X = X, PE = rep(0, p), sig2 = sigma2, lbd = lbd,
#'      weights = weights, B0 = B0, S0 = S0, niter = 50, burnin = 0,
#'      type = "coeff"))
#' @export
summary.MHLS <- function (object, ...) {
  betaSummary <- t(apply(object$beta,2,SummBeta))
  #signsummary <- t(apply(object$beta,2,SummSign))
  subgradSummary <- t(apply(object$subgrad,2,SummBeta))
  result <- list(beta=betaSummary,subgradient=subgradSummary)
  class(result) <- "summary.MHLS"
  return(result)
}

#' @method plot MHLS
#' @title Plot Metropolis-Hastings sampler outputs
#'
#' @description Provides six plots for each covariate index;
#'  histogram, path plot and acf plot for beta and for its subgradient.
#'
#' @param x an object of class "MHLS", which is an output of \code{\link{MHLS}}.
#' @param index an index of covariates to plot.
#' @param skipS logical. If \code{skipS = TRUE}, plots beta only.
#' @param ...	additional arguments passed to or from other methods.
#' @details
#' \code{\link{plot.MHLS}} provides summary plots of beta and subgradient.
#'  The first column provides histogram of beta and subgradient, while the second
#'  and the third columns provide path and acf plots, respectively.
#'  If \code{skipS = TRUE}, this function provides summary plots for beta only.
#' @examples
#' #' set.seed(123)
#' n <- 10
#' p <- 5
#' X <- matrix(rnorm(n * p), n)
#' Y <- X %*% rep(1, p) + rnorm(n)
#' sigma2 <- 1
#' lbd <- .37
#' weights <- rep(1, p)
#' LassoResult <- lassoFit(X = X, Y = Y, lbd = lbd, type="lasso", weights = weights)
#' B0 <- LassoResult$B0
#' S0 <- LassoResult$S0
#' plot(MHLS(X = X, PE = rep(0, p), sig2 = 1, lbd = 1, group = 1:p,
#'      weights = weights, B0 = B0, S0 = S0, niter = 50, burnin = 0,
#'      type = "coeff"))
#' @export
plot.MHLS <- function(x, index = 1:ncol(x$beta), skipS = FALSE, ...) {
  #	n=nrow(x$beta)
  if (any(!index %in% 1:ncol(x$beta))) {
    stop("Invalid index.")
  }

  niter <- x$niter
  burnin <- x$burnin

  if (!skipS) {par(mfrow = c(2,3))} else {par(mfrow = c(1,3))}

  if (!skipS)	{
    for (i in index) {
      hist(x$beta[,i],breaks=20,prob=T,xlab=paste("Beta_",i,sep=""),ylab="Density",main="")
      #ts.plot(x$beta[,i],xlab="Iterations",ylab="Samples")
      plot((burnin+1):niter,x$beta[,i],xlab="Iterations",ylab="Path",type="l")
      if ( sum(abs(diff(x$beta[,i]))) == 0 ) { plot( 0,type="n",axes=F,xlab="",ylab="")
        text(1,0,"Auto correlation plot \n not available",cex=1)} else {
          acf(x$beta[,i],xlab="Lag",main="")
        }
      hist(x$subgrad[,i],breaks=seq(-1-1/10,1,by=1/10)+1/20,prob=T,xlim=c(-1-1/20,1+1/20),xlab=paste("Subgradient_",i,sep=""),ylab="Density",main="")
      #ts.plot(x$subgrad[,i],xlab=Iterations,ylab="Samples")
      plot((burnin+1):niter,x$subgrad[,i],xlab="Iterations",ylab="Samples",type="l")
      if ( sum(abs(diff(x$subgrad[,i]))) == 0 ) { plot( 0,type="n",axes=F,xlab="",ylab="")
        text(1,0,"Auto correlation plot \n not available",cex=1)} else {
          acf(x$subgrad[,i],xlab="Lag",main="")
        }
      readline("Hit <Return> to see the next plot: ")
    }
  } else {
    for (i in index) {
      hist(x$beta[,i],breaks=20,prob=T,xlab=paste("Beta_",i,sep=""),ylab="Density",main="")
      #ts.plot(x$beta[,i],xlab="Iterations",ylab="Samples")
      plot((burnin+1):niter,x$beta[,i],xlab="Iterations",ylab="Path",type="l")
      if ( sum(abs(diff(x$beta[,i]))) == 0 ) { plot( 0,type="n",axes=F,xlab="",ylab="")
        text(1,0,"Auto correlation plot \n not available",cex=1)} else {
          acf(x$beta[,i],xlab="Lag",main="")
        }
      readline("Hit <Return> to see the next plot: ")
    }
  }
}



