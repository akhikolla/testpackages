`%do%` <- foreach::`%do%`
`%dopar%` <- foreach::`%dopar%`
utils::globalVariables("i")

#' @importFrom stats binomial dnorm glm quantile sd
NULL

# SCORNET.R: Contains SCORNET survival estimator function.
# Author: Yuri Ahuja
# Last Updated: 9/15/2020
#
# Semi-Supervised Calibration of Risk with Noisy Event Times (SCORNET) is a consistent, non-parametric
# survival curve estimator that boosts efficiency over existing non-parametric estimators
# by (1) utilizing unlabeled patients in a semi-supervised fashion, and (2) leveraging
# information-dense engineered EHR features to maximize unlabeled set imputation precision
# See Ahuja et al. (2020) BioArxiv for details

expit <- function(x){
  1/(1+exp(-x))
}

logit <- function(x){
  log(x/(1-x))
}

# Default kernel function: standard normal PDF
Knorm <- function(t0,t,b=1){
  dnorm(abs(t-t0),sd=b)
}

# Kernel-Smoothed Cox/Breslow estimator of C|Z0
estimate.h <- function(C,Z=NULL,b=NULL,nCores=1){
  N <- length(C)
  if (is.null(Z)){Z <- matrix(1,nrow=length(C),ncol=1)}
  if (is.null(b)){b <- N^(-1/4) * min(sd(C), (quantile(C,0.75)-quantile(C,0.25))/1.34)}
  
  beta_C.Z <- survival::coxph(survival::Surv(C)~Z)$coefficients; beta_C.Z[is.na(beta_C.Z)] <- 0
  
  ord <- order(C); ordInv <- order(ord)
  Z <- as.matrix(Z[ord,]); C <- C[ord]
  denom <- sapply(1:N,function(i){sum(exp(as.matrix(Z[i:N,]) %*% beta_C.Z))})
  hC0 <- 1 / denom
  hC0k <- c(kernelSmoothen(hC0,C,b))
  
  if (nCores == 1){
    HC0 <- foreach::foreach(i=1:N, .combine=c) %do% {
      pracma::trapz(C[1:i],hC0k[1:i])
    }
  }
  else{
    HC0 <- foreach::foreach(i=1:N, .combine=c, .export='trapz') %dopar% {
      pracma::trapz(C[1:i],hC0k[1:i])
    }
  }
  SC0 <- exp(-HC0)
  h <- exp(Z %*% beta_C.Z)
  
  as.vector(h * hC0k * SC0^h)[ordInv]
}

#' SCORNET Estimator
#' @param Delta Labeled set current status labels (I(T<C))
#' @param C Labeled set censoring times
#' @param t0.all Times at which to estimate survival 
#' @param C.UL Unlabeled set censoring times
#' @param filter Labeled set binary filter indicators
#' @param filter.UL Unlabeled set filter indicators
#' @param Z0 Labeled set baseline feature matrix
#' @param Z0.UL Unlabeled set baseline feature matrix
#' @param Zehr Labeled set EHR-derived feature matrix
#' @param Zehr.UL Unlabeled set EHR-derived feature matrix
#' @param K Kernel function (defaults to standard normal) 
#' @param b bandwidth (default set heuristically)
#' @param h_hat N^1/4-consistent pdf estimator of C|Z0 (defaults to Kernel-Smoothed Cox/Breslow estimator)
#' @param nCores Number of cores to use for parallelization (defaults to 1)
#' @return S_hat: Survival function estimates at times t0.all; StdErrs: Asymptotically consistent standard error estimates corresponding to S_hat
#' @export
scornet <- function(Delta, C, t0.all, C.UL = NULL, filter = NULL, filter.UL = NULL, Z0 = NULL, Z0.UL = NULL,
                    Zehr = NULL, Zehr.UL = NULL, K = Knorm, b = NULL, h_hat = NULL, nCores = 1) {
  Ctot <- c(C,C.UL)
  N <- length(C)
  Ntot <- length(Ctot)
  if (is.null(filter)){filter <- rep(TRUE,N)}
  if (is.null(filter.UL)){filter.UL <- rep(TRUE,length(C.UL))}
  filtertot <- c(filter,filter.UL)
  if (!is.null(Z0)){Z0 <- as.matrix(Z0)} else{Z0 <- matrix(1,N,1)}
  if (is.null(Z0.UL)){Z0.UL <- rep(1,length(C.UL))}
  Z0tot <- rbind(Z0,as.matrix(Z0.UL))
  if (!is.null(Zehr)){Zehr <- as.matrix(Zehr)}
  Zehrtot <- Zehr; if (!is.null(Zehr.UL)){Zehrtot <- rbind(Zehrtot,as.matrix(Zehr.UL))}
  Cfp <- Ctot[filtertot]
  Z0fp <- Z0tot[filtertot,]
  Zehrfp <- Zehrtot[filtertot,]
  Nfp <- sum(filtertot)
  
  if (nCores > 1){
    logfile <- "SCORNET.log"
    writeLines(c(""), file(logfile,'w'))
    clust <- parallel::makeCluster(nCores, outfile=logfile)
    doParallel::registerDoParallel(clust)
  }
  
  if (is.null(b)){
    b <- N^(-1/3) * min(sd(Cfp), (quantile(Cfp,0.75)-quantile(Cfp,0.25))/1.34)
  }
  nu <- Ntot^(-1/4) * min(sd(Ctot), (quantile(Ctot,0.75)-quantile(Ctot,0.25))/1.34)
  
  
  # STEP 1A: Estimate density h(C|Z0)
  
  if (is.null(h_hat)){
    h_hat <- estimate.h(Ctot,Z0tot,nu,nCores)
  }
  
  
  # STEP 1B: Estimate probability P(T<=t|Z,Z0,C)
  
  Kmat <- outer(Ctot,t0.all,function(x,y){K(x,y,b)})
  beta_T.Z <- sapply(1:length(t0.all),function(i){
    glm(Delta~cbind(Z0,Zehr), family=quasibinomial, weights=filter*Kmat[1:N,i])$coef
  })
  beta_T.Z[is.na(beta_T.Z)] <- 0
  
  
  # STEP 2: Estimate marginal survival function S_T(t)
  
  Fi_hat <- expit(cbind(1,Z0tot,Zehrtot) %*% beta_T.Z)
  Fi_hat[!filtertot,] <- 0
  S_hat <- 1 - (colSums(Fi_hat*Kmat/h_hat) / colSums(Kmat/h_hat))
  
  
  # Estimate standard errors for S_T(t)
  
  tau_sq <- 1 / (2*sqrt(pi))
  StdErrs <- sapply(1:length(t0.all),function(i){
    tryCatch({
      W <- cbind(1,Z0fp,Zehrfp)
      gTW <- expit(W %*% beta_T.Z[,i])
      gTWprime <- gTW * (1-gTW)
      if (sum(gTWprime)==0){NA}
      else{
        At <- matrix(0,ncol(W),ncol(W))
        fps <- which(filtertot)
        for (j in 1:Nfp){
          At <- At + gTWprime[j] * Kmat[fps[j],i] * (W[j,] %*% t(W[j,]))
        }
        At <- At / sum(Kmat[fps,i])
        ABA <- solve(At) / (b * sum(Kmat[1:N,i]) / tau_sq)
        
        gTWprime <- Fi_hat[filtertot,i] * (1-Fi_hat[filtertot,i])
        P <- t(W) %*% (Kmat[filtertot,i]*gTWprime/h_hat[filtertot]) / sum(Kmat[,i]/h_hat)
        ABAP <- c(t(P) %*% ABA %*% P)
        
        sqrt(ABAP)
      }      
    }, error=function(e){NA})
  })
  
  return(list('S_hat'=S_hat, 'StdErrs'=StdErrs))
}
