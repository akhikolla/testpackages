#' @useDynLib EAinference
#' @importFrom Rcpp evalCpp
NULL

#' @importFrom msm rtnorm
#' @importFrom mvtnorm dmvnorm
#' @import graphics
#' @import stats
#-------------------------------------------
# Utility functions for lasso
#-------------------------------------------
logFlipTR <- function(x,t)
{
  logr <- log(2) + dnorm(x, mean=0, sd=t,log=TRUE);
  # theta(b_j;0,tau_j^2)/(1/2)
  return(logr);
}

SummBeta <- function(x) {
  c( mean=mean(x), median = median(x), s.d = sd(x), quantile(x,c(.025,.975)) )
}

SummSign <- function(x) {
  n=length(x)
  return ( c(Positive.proportion=sum(x>0)/n , Zero.proportion=sum(x==0)/n, Negative.proportion=sum(x<0)/n  ) )
}

UpdateBA <- function(Bcur,Scur,tau,A,I,Rcur,logfRcur,VRC,lbdVRW,InvVarR,
                  tVAN.WA,invB_D,B_F,FlipSA,IndexSF,nA=length(A),nI=length(I))
{
  Bnew=Bcur;
  Bnew[A]=rnorm(nA,Bcur[A],tau[A])

  diffB=Bnew-Bcur;
  nAcceptBA=0;
  nChangeSign=c(0,0,0);
  for(j in FlipSA)
  {
    if (sign(Bnew[j])!=sign(Bcur[j])) {
      nChangeSign[2]=nChangeSign[2]+1
      Snew=Scur;
      Snew[j]=sign(Bnew[j])
      if (is.null(IndexSF)) {
        Snew[I] = -invB_D%*%tVAN.WA%*%Snew[A]
      } else {
        Snew[I][-IndexSF]=invB_D %*% (-B_F %*% Snew[I][IndexSF] - tVAN.WA%*%Snew[A])
      }

      if (if (is.null(IndexSF)) {all(abs(Snew[I]) <= 1)} else {all(abs(Snew[I][-IndexSF]) <= 1)}) {
        nChangeSign[3]=nChangeSign[3]+1
        diffR=diffB[j]*VRC[,j] + lbdVRW%*%(Snew-Scur);
        Rnew=Rcur+diffR;
        logfRnew=-0.5*sum(Rnew^2*InvVarR);
        logTargetRatio=logfRnew-logfRcur;

        logMH=logTargetRatio;
        u=runif(1);
        if(log(u)<logMH)
        {
          # loc = which(apply(S.temp$Sign,1,function(x) all(x==Scur[A])))
          # if(length(loc)!=0) {
          # 	S.temp$Sign[loc,]=sign(Bcur[A])
          # 	S.temp$Subgrad[loc,]=Scur[I]
          # } else {
          # 	S.temp$Sign=rbind(S.temp$Sign,sign(Bcur[A]))
          # 	S.temp$Subgrad=rbind(S.temp$Subgrad,Scur[I])
          # }
          Rcur=Rnew;
          logfRcur=logfRnew;
          Bcur[j]=Bnew[j];
          Scur=Snew
          nAcceptBA=nAcceptBA+1;
          nChangeSign[1]=nChangeSign[1]+1;
        }
      }
    } else {
      diffR=diffB[j]*VRC[,j];
      Rnew=Rcur+diffR;
      logfRnew=-0.5*sum(Rnew^2*InvVarR);
      logTargetRatio=logfRnew-logfRcur;

      logMH=logTargetRatio;
      u=runif(1);
      if(log(u)<logMH)
      {
        Rcur=Rnew;
        logfRcur=logfRnew;
        Bcur[j]=Bnew[j];
        nAcceptBA=nAcceptBA+1;
      }
    }
  }
  return(list(B=Bcur,S=Scur,Rvec=Rcur,logf=logfRcur,nAccept=nAcceptBA,nChangeSign=nChangeSign));
}

UpdateBA.fixedSA <- function(Bcur,tau,A,Rcur,logfRcur,VRC,InvVarR) # Fix sign(beta[A])
{
  nA=length(A);
  LUbounds=matrix(0,nA,2);
  LUbounds[Bcur[A]>0,2]=Inf;
  LUbounds[Bcur[A]<0,1]=-Inf;
  Bnew=Bcur;
  Bnew[A]=rtnorm(nA,Bcur[A],tau[A],lower=LUbounds[,1],upper=LUbounds[,2]);

  Ccur=pnorm(0,mean=Bcur[A],sd=tau[A],lower.tail=TRUE,log.p=FALSE);
  Ccur[Bcur[A]>0]=1-Ccur[Bcur[A]>0];
  Cnew=pnorm(0,mean=Bnew[A],sd=tau[A],lower.tail=TRUE,log.p=FALSE);
  Cnew[Bcur[A]>0]=1-Cnew[Bcur[A]>0];
  lqratio=log(Ccur/Cnew);

  diffB=Bnew-Bcur;
  i=1;
  nAcceptBA=0;
  for(j in A)
  {
    diffR=diffB[j]*VRC[,j];
    Rnew=Rcur+diffR;
    logfRnew=-0.5*sum(Rnew^2*InvVarR);
    logTargetRatio=logfRnew-logfRcur;

    logMH=logTargetRatio+lqratio[i];
    u=runif(1);
    if(log(u)<logMH)
    {
      Rcur=Rnew;
      logfRcur=logfRnew;
      Bcur[j]=Bnew[j];
      nAcceptBA=nAcceptBA+1;
    }
    i=i+1;
  }
  return(list(B=Bcur,Rvec=Rcur,logf=logfRcur,nAccept=nAcceptBA));
}

# Subfunction for UpdateBA : In high-dim, when sign(BA) changes, update SI.
# ProposeSI=function(PropSA, CurrSI, S.temp, tVAN.WA, nI, E, G, H, iter=10){
# 	loc = which(apply(S.temp$Sign,1,function(x) all(x==PropSA)))
# 	if(length(loc)!=0) return(S.temp$Subgrad[loc,])
# 	# There already is a valid SI
#
# 	F1 <- -tVAN.WA%*%PropSA
#
# 	x0 <- suppressWarnings(lsei(A=diag(nI),B=c(CurrSI),G=G,H=H,E=E,F=F1))
# 	if (x0$IsError == TRUE) return(rep(9,nI)) # Solution does not exist
# 	#return(x0$X[1:nI])
# 	#return(suppressWarnings(xsample(A=diag(nI),B=c(CurrSI),G=G,H=H,E=E,F=F1,iter=iter,x0=x0$X))$X[100,1:nI])
# 	return(suppressWarnings(xsample(G=G,H=H,E=E,F=F1,iter=iter,x0=x0$X))$X[iter,1:nI])
# 	# Draw valid SI
# }

UpdateSI <- function(Scur,A,I,Rcur,n,p,logfRcur,lbdVRW,InvVarR,
                  tVAN.WA,invB_D,BDF, B_F, IndexSF, SIscale=1, nA=length(A))
{
  if (SIscale <1) stop("SIscale has to be greater than 1")

  nAcceptSI=0;
  U=-tVAN.WA%*%Scur[A]
  #V_AN=V[A,N,drop=FALSE];
  #V_IN=V[I,N,drop=FALSE];
  #V_AR=V[A,R,drop=FALSE];
  #V_IR=V[I,R,drop=FALSE];
  #t(V_AN)%*%W[A,A]%*%Snew[A]+t(V_IN)%*%W[I,I]%*%Snew[I]

  a=invB_D%*%U + rep(-1,p-n)
  b=a + rep(2,p-n)

  for ( i in 1:(n-nA)) {

    Scur_F=Scur[I][IndexSF]
    Scur_D=setdiff(Scur[I],Scur_F)

    bd=cbind(a-BDF[,-i,drop=FALSE]%*%Scur_F[-i], b-BDF[,-i,drop=FALSE]%*%Scur_F[-i])
    #		BDF[,i]*Scur_F[i]
    bd=bd/BDF[,i]

    bd=t(apply(bd,1,sort))

    lowbd=max(c(bd[,1],-1))
    highbd=min(c(bd[,2],1))

    if(lowbd>highbd) {
      if ( lowbd-highbd > 1e-10) {
        stop("Lowerbound is greater than Uppderbound?")
      } else {
        lowbd = highbd
      }
    }

    lnth=(highbd-lowbd)/SIscale

    Snew_F=Scur_F
    Snew_F[i]=runif(1,(lowbd+highbd)/2-lnth/2,(lowbd+highbd)/2+lnth/2)
    Snew_D = invB_D%*%(U-B_F%*%Snew_F)

    Snew_I=Scur[I]
    Snew_I[IndexSF]=Snew_F
    Snew_I[-IndexSF]=Snew_D

    diffSI=Snew_I-Scur[I]

    diffR=lbdVRW[,I,drop=FALSE]%*%diffSI;
    Rnew=Rcur+diffR;
    logfRnew=-0.5*sum(Rnew^2*InvVarR);
    logTargetRatio=logfRnew-logfRcur;

    #logPropRatio = log(lnth) -log(2)
    # Proposal density should be the same in both directions.

    #logMH=logTargetRatio + logPropRatio;
    logMH=logTargetRatio;

    u=runif(1);
    if(log(u)<logMH)
    {
      Rcur=Rnew;
      logfRcur=logfRnew;
      Scur[I]=Snew_I;
      nAcceptSI=nAcceptSI+1;
    }
  }
  return(list(S=Scur,Rvec=Rcur,logf=logfRcur,nAccept=nAcceptSI));
}


CalTmat <- function(p,n,V,LBD,W,lbd,R,N,A,I)
{
  V_IN=V[I,N,drop=FALSE];
  V_AR=V[A,R,drop=FALSE];
  V_IR=V[I,R,drop=FALSE];
  if(length(A)<n)
  {
    BNul=as.matrix(svd(t(V_IN)%*%as.matrix(W[I,I]),nv=length(I))$v[,(p-n+1):length(I)]);
    Tmat=cbind(LBD%*%t(V_AR),lbd*t(V_IR)%*%as.matrix(W[I,I])%*%BNul);
    #LBD%*%t(V_AR) == t(VR)%*%C[,A]
  }else{
    BNul=NULL;
    Tmat=LBD%*%t(V_AR);
  }
  logdetT=determinant(as.matrix(Tmat),logarithm=TRUE)$modulus[1];
  return(list(Basisset=BNul,T=Tmat,logdetT=logdetT));
}
#-------------------------------------------
# Utility functions for Group Lasso
#-------------------------------------------
group.norm <- function(x, group, al = 2) {
  tapply(x,group, function(x) sum(abs(x)^al)^(1/al))
}

group.norm2 <- function(x, group) {
  result <- c()
  for (i in 1:max(group)) {
    result[i] <- sqrt(crossprod(x[group==i]))
  }
  return(result)
}

# T(s,A) : p x (n-|A|) matrix s.t. ds = T(s,A)ds_F
TsA <- function(Q, s, group, A, n, p) {
  if (n < p && missing(Q)) {
    stop("High dimensional setting needs Q")
  }
  nA <- length(A)
  rankX <- min(c(n,p))
  Subgradient.group.matix <- matrix(0, length(unique(group)), p)

  for (i in 1:length(unique(group))) {
    Subgradient.group.matix[i, which(group == i)] = s[group == i]
  }
  Subgradient.group.matix <- Subgradient.group.matix[A, ,drop=FALSE]

  if (n >= p) {
    B <- Subgradient.group.matix
  } else {
    B <- rbind(t(Q), Subgradient.group.matix)
  }

  if (nA != 0 ) {
    P <- matrix(0, p, p)
    Permute.order <- 1:p

    if (n < p) {nrowB <- p - n + nA} else
    {nrowB <- nA}
    for (i in 1:nrowB) {
      if (B[i, rankX - nA + i] == 0) {
        W1 <- which(B[i,] !=0)[1]
        #Permute.order[c(W1,n-length(A)+i)] = Permute.order[c(n-length(A)+i,W1)]
        Permute.order[c(which(Permute.order == W1),rankX-nA+i)] =
          Permute.order[c(rankX-nA+i,which(Permute.order == W1))]
        #print(Permute.order)
      }
    }
    for ( i in 1:p) {
      P[Permute.order[i], i] <- 1
    }
    B <- (B%*%P)
  } else {
    P <- diag(p)
  }

  if (nA == 0 && n >= p) return(P);

  if (rankX-nA >= 1) {
    BF <- B[, 1:(rankX - nA),drop=FALSE]
    BD <- B[, -c(1:(rankX - nA)),drop=FALSE]
    #Result <- P %*% rbind(diag(n-length(A)), -solve(BD)%*%BF)
    Result <- P %*% rbind(diag(rankX-nA), -solve(BD,BF))
  } else {
    Result <- -solve(B)
  }
  return(Result)
}
#------------------------------
TsA.qr <- function(Q, s, group, A, n, p) { # T(s,A)
  if (n == length(A)) {stop("|A| should be smaller than n")}
  if (n < p && missing(Q)) {
    stop("High dimensional setting needs Q")
  }

  nA <- length(A)

  Subgradient.group.matix <- matrix(0, length(unique(group)), p)
  for (i in 1:length(unique(group))) {
    Subgradient.group.matix[i, which(group == i)] = s[group == i]
  }
  Subgradient.group.matix <- Subgradient.group.matix[A, , drop=FALSE]
  if (n >= p) {
    B <- Subgradient.group.matix
  } else {
    B <- rbind(t(Q), Subgradient.group.matix)
  }

  if (n < p) { nrowB <- p - n + nA ; rankX <- n} else
  {nrowB <- nA; rankX <- p}

  QR.B <- qr(B)
  Pivot <- sort.list(qr(B)$pivot)
  B.Q <- qr.Q(QR.B)
  #B.R <- qr.R(QR.B)[,Pivot]
  B.R <- qr.R(QR.B)
  tP <- matrix(0, p, p)

  for (i in 1:p) {
    tP[Pivot[i], i] = 1
  }

  P <- t(tP)

  if (rankX-nA >= 1) {
    RF <- B.R[, 1:(rankX - nA), drop=FALSE]
    RD <- B.R[, -c(1:(rankX - nA)), drop=FALSE]
    Result <- P %*% rbind(diag(rankX-nA), -solve(RD)%*%RF)
  }
  return(Result)
}

# TsA.null <- function(t.XWinv, s, group, A, n, p) { # T(s,A)
#   # even if length(A) == 0, everything will work just fine !!
#   # when length(A) == n, we only compute F2 function.
#   # Updated for Low-dim case
#   if (missing(t.XWinv) && n < p) {
#     stop("When n < p, t.XWinv is needed")
#   }
#
#   nA <- length(A)
#   if (n < p) { # High-dim
#     if (nA !=0) {
#       Subgradient.group.matix <- matrix(0, nA, p)
#       for (i in 1:nA) {
#         Subgradient.group.matix[i, which(group == A[i])] <- s[group == A[i]]
#       }
#       t.XWinv %*% Null(t(Subgradient.group.matix %*% t.XWinv))
#     } else {
#       t.XWinv
#     }
#   } else { # Low-dim
#     if (nA !=0) {
#       Subgradient.group.matix <- matrix(0, nA, p)
#       for (i in 1:nA) {
#         Subgradient.group.matix[i, which(group == A[i])] <- s[group == A[i]]
#       }
#       Null(t(Subgradient.group.matix))
#     } else {
#       diag(p)   # if |A| = 0, T should be pXp identity matrix.
#     }
#   }
# }

# F1 = r \circ \psi , eq(3.6), p x p matrix
F1 <- function(r, Psi, group) {
  Result <- Psi
  for (i in 1:length(unique(group))) {
    Result[, group == i] <- Result[, group == i] * r[i]
  }
  return(Result)
}

# F2 = \psi \circ \eta , eq(3.7), p x J matrix, where J is the number of groups
F2 <- function(s, Psi, group) {
  Result <- matrix(, nrow(Psi), length(unique(group)))
  for (i in 1:length(unique(group))) {
    Result[, i] = crossprod(t(Psi[, group == i]), s[group == i])
  }
  return(Result)
}

logJacobiPartial <- function(X, s, r, Psi, group, A, lam, W, TSA) { # log(abs(det(X %*% [F2[,A] | (F1 + lam * W) %*% TsA])))
  n <- nrow(X)
  p <- ncol(X)
  table.group <- table(group)

  # W <- c()
  # for (i in 1:length(table.group)) {
  #   W <- c(W, rep(weights[i], table.group[i]))
  # }

  if (n < p) { # High-dim
    if (n == length(A)) {
      log.Det <- determinant(X %*% F2(s, Psi, group)[,A])
    } else {
      log.Det <- determinant(X %*% cbind(F2(s, Psi, group)[, A], (F1(r, Psi, group) + lam * diag(W)) %*% TSA))
    }
    return(log.Det[[1]][1]);
  } else { # Low-dim
    if (p == length(A)) {
      log.Det <- determinant(F2(s,Psi,group))
    } else {
      log.Det <- determinant(cbind(F2(s, Psi, group)[, A], (F1(r, Psi, group) + lam * diag(W)) %*% TSA))
    }
    return(log.Det[[1]][1]);
  }
}

ld.Update.r <- function(rcur,Scur,A,Hcur,X,pointEstimate,Psi,W,lbd,group,inv.Var,tau,PEtype,n,p) {
  rprop <- rcur;
  nrUpdate <- 0;
  Bcur <- Bprop <- Scur * rep(rcur,table(group));
  TSA.cur <- TSA.prop <- TsA(,Scur,group,A,n,p);
  for (i in A) {
      #rprop[i] <- truncnorm::rtruncnorm(1, 0, , rcur[i], sqrt(tau[i] * ifelse(rcur[i]!=0,rcur[i],1)))
    rprop[i] <- rtnorm(n = 1, mean = rcur[i], sd = sqrt(tau[i] * ifelse(rcur[i]!=0,rcur[i],1)), lower = 0, upper = Inf)
    Bprop[group==i] <- rprop[i] * Scur[group==i]

    if (PEtype == "coeff") {
      Hprop <- drop(Psi %*% drop(Bprop - pointEstimate) + lbd * W * drop(Scur))
    } else {
      Hprop <- drop(Psi %*% drop(Bprop) - t(X) %*% pointEstimate / n + lbd * W * drop(Scur))
    }

    Hdiff <- Hcur - Hprop

    lNormalRatio <- drop(t(Hdiff)%*% inv.Var %*% (Hprop + Hdiff/2))
    #dmvnorm(Hprop,,sig2 / n * Psi,log=T) - dmvnorm(Hcur,,sig2 / n * Psi,log=T)
    lJacobianRatio <- logJacobiPartial(X,Scur,rprop,Psi,group,A,lbd,W,TSA.prop) -
      logJacobiPartial(X,Scur,rcur,Psi,group,A,lbd,W,TSA.cur)
    lProposalRatio <- pnorm(0,rcur[i],sqrt(tau[i] * rcur[i]), lower.tail=FALSE, log.p=TRUE) -
      pnorm(0,rprop[i],sqrt(tau[i] * rprop[i]), lower.tail=FALSE, log.p=TRUE)
    lAcceptanceRatio <-  lNormalRatio + lJacobianRatio + lProposalRatio
    if (lAcceptanceRatio <= log(runif(1))) { # Reject
      rprop[i] <- rcur[i];
      Bprop[group==i] <- Bcur[group==i]
    } else { # Accept
      nrUpdate <- nrUpdate + 1;
      Hcur <- Hprop;
    }
  }
  return(list(r = rprop, Hcur = Hcur, nrUpdate = nrUpdate))
}

ld.Update.S <- function(rcur,Scur,A,Hcur,X,pointEstimate,Psi,W,lbd,group,inv.Var,PEtype,n,p) {
  Sprop <- Scur;
  nSUpdate <- 0;
  #p <- ncol(X)
  for (i in 1:max(group)) {
    if (i %in% A) {Sprop[group == i] <- rUnitBall.surface(sum(group == i))} else {
      Sprop[group ==i] <- rUnitBall(sum(group==i))
    }

    if (PEtype == "coeff") {
      Hprop <- drop(Psi %*% drop(Sprop * rep(rcur,table(group)) - pointEstimate) + lbd * W * drop(Sprop))
    } else {
      Hprop <- drop(Psi %*% drop(Sprop * rep(rcur,table(group))) - t(X) %*% pointEstimate / n + lbd * W * drop(Sprop))
    }
    Hdiff <- Hcur - Hprop

    lNormalRatio <- drop(t(Hdiff)%*% inv.Var %*% (Hprop + Hdiff/2))
    #dmvnorm(Hprop,,sig2 / n * Psi,log=T) - dmvnorm(Hcur,,sig2 / n * Psi,log=T)
    lJacobianRatio <- logJacobiPartial(X,Sprop,rcur,Psi,group,A,lbd,W,TsA(,Sprop,group,A,n,p)) -
      logJacobiPartial(X,Scur,rcur,Psi,group,A,lbd,W,TsA(,Scur,group,A,n,p))
    lAcceptanceRatio <-  lNormalRatio + lJacobianRatio
    if (lAcceptanceRatio <= log(runif(1))) { # Reject
      Sprop[group == i] <- Scur[group == i];
    } else { # Accept
      nSUpdate <- nSUpdate + 1;
      Hcur <- Hprop;
    }
  }
  return(list(S = Sprop, Hcur = Hcur, nSUpdate = nSUpdate))
}
rUnitBall.surface <- function(p, n = 1) {
  # p : dimension
  # n : the number of samples you want to draw
  Result <- matrix(0,n,p)
  for (i in 1:n) {
    x <- rnorm(p)
    Result[i,] <- x / c(sqrt(crossprod(x)))
  }
  return(Result)
}
rUnitBall <- function(p, n = 1) {
  # p : dimension
  # n : the number of samples you want to draw
  Result <- matrix(0,n,p)
  # for (i in 1:n) {
  #   x <- rnorm(p,,1/sqrt(2));
  #   y <- rexp(1)
  #   Result[i,] <- x / c(sqrt(y+crossprod(x)))
  # }
  for (i in 1:n) {
    x <- rnorm(p);
    y <- rexp(1)
    Result[i,] <- x * runif(1)^(1/p) / c(sqrt(crossprod(x)))
  }
  return(Result)
}
# hd.Update.r <- function(rcur,Scur,A,Hcur,X,coeff,Psi,W,lbd,group,inv.Var,tau) {}
# hd.Update.S <- function(rcur,Scur,A,Hcur,X,coeff,Psi,W,lbd,group,inv.Var,p) {}
#-------------------------------------------
# Utility functions for MHLS summary
#-------------------------------------------
SummBeta <- function ( x ) {
  c( mean=mean(x) , median = median(x) , s.d = sd(x) , quantile(x,c(.025,.975)) )
}

SummSign <- function ( x ) {
  n=length(x)
  return ( c(Positive.proportion=sum(x>0)/n , Zero.proportion=sum(x==0)/n, Negative.proportion=sum(x<0)/n  ) )
}
#-------------------------------------------
# Utility functions for scaled lasso / scaled group lasso
#-------------------------------------------
TsA.slasso <- function(SVD.temp, Q, s, W, group, A, n, p) {
  if (n < p && missing(Q)) {
    stop("High dimensional setting needs Q")
  }
  nA <- length(A)
  Subgradient.group.matix <- matrix(0, nA, p)

  IndWeights <- W
  if (nA != 0) {
    for (i in 1:nA) {
      Subgradient.group.matix[i, which(group == A[i])] = s[group == A[i]]
    }
  }
  #Subgradient.group.matix <- Subgradient.group.matix[A, ,drop=FALSE]

  #all.equal(V%*%diag(1/D^2)%*%t(V) , t(X) %*% solve(tcrossprod(X)) %*% solve(tcrossprod(X)) %*% X)
  #all.equal(U%*%diag(D)%*%t(V) , X)

  B <- rbind(t(Q), Subgradient.group.matix,
             t(s * IndWeights) %*% SVD.temp %*% diag(IndWeights))

  if (nA != 0 ) {
    P <- matrix(0, p, p)
    Permute.order <- 1:p

    # if (n < p) {nrowB <- p - n + nA} else
    # {nrowB <- nA}
    for (i in 1:nrow(B)) {
      if (B[i, n - nA - 1 + i] == 0) {
        W1 <- which(B[i,] !=0)[1]
        #Permute.order[c(W1,n-length(A)+i)] = Permute.order[c(n-length(A)+i,W1)]
        Permute.order[c(which(Permute.order == W1),n-nA-1+i)] =
          Permute.order[c(n-nA-1+i,which(Permute.order == W1))]
        #print(Permute.order)
      }
    }
    for ( i in 1:p) {
      P[Permute.order[i], i] <- 1
    }
    B <- (B%*%P)
  } else {
    P <- diag(p)
  }

  if (n-nA-1 >= 1) { # n-nA-1 : # of free coordinate
    BF <- B[, 1:(n - nA - 1),drop=FALSE]
    BD <- B[, -c(1:(n - nA - 1)),drop=FALSE]
    #Result <- P %*% rbind(diag(n-length(A)), -solve(BD)%*%BF)
    Result <- P %*% rbind(diag(n-nA-1), -solve(BD,BF))
  } else {
    Result <- -solve(B)
  }
  return(Result)
}

logJacobiPartial.slasso <- function(X, s, r, Psi, group, A, lam, hsigma, W, TSA) { # log(abs(det(X %*% [F2[,A] | (F1 + lam * W) %*% TsA])))
  # This function is only for high-dimensional cases.
  n <- nrow(X)
  p <- ncol(X)
  table.group <- table(group)
  if (n > p) stop("High dimensional setting is required.")
  if (n == (length(A)-1)) {
    log.Det <- determinant(X %*% cbind(F2(s, Psi, group)[,A], diag(lam * hsigma * W) %*% s))
  } else {
    log.Det <- determinant(X %*% cbind(F2(s, Psi, group)[, A], (F1(r, Psi, group)
                                                                + diag(lam * hsigma * W)) %*% TSA,diag(lam * W) %*% s))
  }
  return(log.Det[[1]][1]);
}

slassoLoss <- function(X,Y,beta,sig,lbd) {
  n <- nrow(X)
  crossprod(Y-X%*%beta) / 2 / n / sig + sig / 2 + lbd * sum(abs(beta))
}

# Scaled lasso / group lasso function, use scaled X matrix.
slassoFit.tilde <- function(Xtilde, Y, lbd, group, weights, Gamma, verbose=FALSE){
  n <- nrow(Xtilde)
  p <- ncol(Xtilde)

  if (verbose) {
    cat("niter \t Loss \t hat.sigma \n")
    cat("-------------------------------\n")
  }
  sig <- signew <- .1
  K <- 1 ; niter <- 0
  if (missing(Gamma)) {
    Gamma <- groupMaxEigen(X = Xtilde, group = group)
  }
  B0 <- rep(1, p)
  while(K == 1 & niter < 500){
    sig <- signew;
    lam <- lbd * sig
    B0 <- grlassoFit(X = Xtilde, Y = Y, group = group, weights = rep(1, max(group)), Gamma = Gamma, lbd = lam, initBeta = B0)$coef
    # B0 <- coef(gglasso(Xtilde,Y,loss="ls",group=group,pf=rep(1,max(group)),lambda=lam,intercept = FALSE))[-1]
    signew <- sqrt(crossprod(Y-Xtilde %*% B0) / n)

    niter <- niter + 1
    if (abs(signew - sig) < 1e-04) {K <- 0}
    if (verbose) {
      cat(niter, "\t", sprintf("%.3f", slassoLoss(Xtilde,Y,B0,signew,lbd)),"\t",
          sprintf("%.3f", signew), "\n")
    }
  }
  lam <- lbd * signew
  B0 <- grlassoFit(X = Xtilde, Y = Y, group = group, weights = rep(1, max(group)),  Gamma = Gamma, lbd = lam, initBeta = B0)$coef
  # B0 <- coef(gglasso(Xtilde,Y,loss="ls",group=group,pf=rep(1,max(group)),lambda=lam,intercept = FALSE))[-1]
  hsigma <- c(signew)
  S0 <- t(Xtilde) %*% (Y - Xtilde %*% B0) / n / lbd / hsigma
  B0 <- B0 / rep(weights,table(group))
  return(list(B0=B0, S0=c(S0), hsigma=hsigma,lbd=lbd))
}


#---------------------
# Error handling
#---------------------
ErrorParallel <- function(parallel, ncores) {
  if(.Platform$OS.type == "windows" && (parallel == TRUE | ncores > 1)){
    ncores <- 1L
    parallel <- FALSE
    warning("Under Windows platform, parallel computing cannot be executed.")
  }

  if (parallel && ncores == 1) {
    ncores <- 2L
    warning("If parallel=TRUE, ncores needs to be greater than 1. Automatically
            set ncores to 2.")
  }
  if (parallel && ncores > parallel::detectCores()){
    ncores <- parallel::detectCores()
    warning("ncores is greater than the maximum number of available processores.
            Set it to the maximum possible value.")
  }
  return(list(parallel=parallel,ncores=ncores))
}


#-------------------------------------------
# Utility functions for MHInference
#-------------------------------------------
#Propose mu hat from 95% region
PluginMu.MHLS <- function(X, Y, lbd, alpha
                          , sigma.hat
                          , nChain, niter, method = "boundary"
                          , parallel, ncores) {
  # method can be either "boundary" or "unif"
  Sample <- switch(method, boundary = rUnitBall.surface, unif = rUnitBall)

  n <- length(Y)
  B0 <- lassoFit(X = X, Y = Y, type = "lasso", lbd = lbd)$B0
  TEMP <- projStein(X, Y, B0, sig2 = sigma.hat^2
                    , alpha=alpha, niter=niter, fixSeed = FALSE)
  Result <- matrix(0,nChain + 1, n)
  Result[1:nChain, ] <- tcrossprod(rep(1,nChain), TEMP$mu_s) + Sample(p = n, n = nChain) * TEMP$r_s * sqrt(n) +
    tcrossprod(rep(1,nChain), TEMP$mu_w) + Sample(p = n, n = nChain) * TEMP$r_w  * sqrt(n)
  Result[nChain+1,] <- TEMP$mu_s + TEMP$mu_w

  # returns nChain+1 by n matrix
  return(Result)
}

# PluginMu.MHLS <- function(X, Y, lbd, ratioSeq = seq(0,1,by=0.01), alpha = 0.05
#                           , nChain, niter = 100, method = "boundary"
#                           , parallel = FALSE, ncores = 2) {
#   # method can be either "boundary" or "unif"
#   Sample <- switch(method, boundary = rUnitBall.surface, unif = rUnitBall)
#
#   n <- length(Y)
#   B0 <- lassoFit(X = X, Y = Y, type = "lasso", lbd = lbd)$B0
#
#   FF <- function(x) {
#     muhat <- numeric(2 * n)
#     TEMP <- projStein(X, Y, B0, lbd = lbd, ratios = ratioSeq,
#                       alpha=alpha, niter=niter, fixSeed = FALSE)
#     muhat[1:n] <- TEMP$mu_s + Sample(n) * TEMP$r_s + TEMP$mu_w + Sample(n) * TEMP$r_w
#     muhat[1:n + n] <- TEMP$mu_s + TEMP$mu_w
#     return(muhat)
#   }
#   TEMP <- t(parallel::mcmapply(FF,1:nChain,mc.cores = ncores))
#
#   # returns nChain by 2n matrix
#   return(rbind(TEMP[,1:n], colMeans(TEMP[,1:n+n])))
# }
####################################################

# Shrink Y towards zero by Stein estimate
# Y: n*1 matrix
# df: degress of freedom
# return: list(mu: estimate of true mu, L: sure (stein unbiased risk estimate))
ShrinkToZeroStein <- function(Y, df, sig2)
{
  B <- sig2 * df / sum(Y^2)
  return(list(mu = (1 - B) * Y, L = sig2 * max(1 - B, 0)))
}

####################################################

# calculate c(alpha)
# mu: n*1 matrix
# df: degress of freedom
# quantile:
# proj: identity matrix by default or a n*n matrix
# niter
# return: the quantile of c(alpha) in Professor's note

####################################################
quantile_stein <- function(mu, df, Quantile, sig2, proj = NULL, niter = 100)
{
  #n = dim(mu)[1]
  n <- length(mu)
  if (is.null(proj)) {proj <- diag(n)}
  z <- numeric(niter)
  for (i in 1:niter)
  {
    Y_prime <- mu + proj %*% matrix(rnorm(n,sd = sqrt(sig2)), ncol = 1)
    res <- ShrinkToZeroStein(Y = Y_prime, df = df, sig2 = sig2)
    z[i] <- res$L - sum((res$mu - mu)^2) / df
  }
  c <- quantile(abs(z) * sqrt(df) / sig2, probs = Quantile)
  return(c)
}
####################################################

# We split X and Y into 2 pieces. (X1, X2) (Y1, Y2).
# Use (X1,Y1) to get the active set A
# Then draw inference using (X2, Y2).

# Given a sequence of ratios, X, Y, lbd, calculate r_s r_w, mu_s mu_w which minimize volumne
# X1, X2: n*p matrx
# Y1, Y2: n*1 matrix
# lbd: tunning parameter, lambda, in lasso
# ratios: a sequence of ratios of tau against lambda (tau the threshold value of \hat{beta})
# alpha:
# niter:
# return: list(r_s, r_w, mu_s, mu_w)
projStein <- function(X, Y, B0, sig2, alpha,
                      niter, fixSeed)
{
  n <- length(Y)
  # B0 <- lassoFit(X = X, Y = Y, type = "lasso", lbd = lbd)$B0

  # estimate c by A=0
  stein <- ShrinkToZeroStein(Y = Y, df = n, sig2 = sig2)
  if (fixSeed) {set.seed(123)}
  c <- quantile_stein(mu = stein$mu, df = n, Quantile = 1-alpha/2, sig2 = sig2, niter = niter)
    # Why 1-alpha/2 instead of 1-alpha?, Bonferroni

  # # information along ratios
  # r_s_ratio <- r_w_ratio <- logVol_ratio <- numeric(length(ratios))
  #
  # for (i in 1:length(ratios))
  # {
  #   A <- which(abs(B0) > ratios[i] * lbd)
  #   if (length(A) != 0)
  #   {
  #     r_s_ratio[i] <- sqrt(qchisq(1-alpha/2, df = length(A)) / n)
  #     P_A <- X[,A] %*% solve(t(X[,A]) %*% X[,A]) %*% t(X[,A])
  #     P_perp <- diag(n) - P_A
  #   } else
  #   {
  #     r_s_ratio[i] <- 0
  #     P_perp <- diag(n)
  #   }
  #   nMinusA <- n - length(A)
  #
  #   stein <- ShrinkToZeroStein(P_perp %*% Y, nMinusA)
  #   r_w_ratio[i] <- sqrt(nMinusA / n * (c / sqrt(nMinusA) + stein$L))
  #
  #   #update r_s_ratio r_w_ratio if r_s_ratio exists and calculate volume
  #   if (length(A) != 0)
  #   {
  #     r_s_ratio[i] <- r_s_ratio[i] * sqrt(n / length(A))
  #     r_w_ratio[i] <- r_w_ratio[i] * sqrt(n / (n - length(A)))
  #     logVol_ratio[i] <- length(A) * log(r_s_ratio[i]) + nMinusA * log(r_w_ratio[i])
  #   } else
  #   {
  #     logVol_ratio[i] <- nMinusA * log(r_w_ratio[i])
  #   }
  # }
  # # information of minimum volume
  # i <- which.min(logVol_ratio)
  A <- which(B0 != 0)
  # compute mu_s
  if (length(A) != 0)
  {
    P_A <- X[,A] %*% solve(t(X[,A]) %*% X[,A]) %*% t(X[,A])
    mu_s <- P_A %*% Y
    P_perp <- diag(n) - P_A
  } else
  {
    mu_s <- 0
    P_perp <- diag(n)
  }
  # compute mu_w
  stein <- ShrinkToZeroStein(Y = P_perp %*% Y, df = n - length(A), sig2 = sig2)
  mu_w <- stein$mu


    if (length(A) != 0)
    {
      r_s <- sqrt(qchisq(1-alpha/2, df = length(A)) / n) * sqrt(sig2)
    } else
    {
      r_s <- 0
    }
    r_w <- sqrt(n - length(A) / n * (c  * sig2 / sqrt(n - length(A)) + stein$L))

    #update r_s_ratio r_w_ratio if r_s_ratio exists and calculate volume
    if (length(A) != 0)
    {
      r_s <- r_s * sqrt(n / length(A))
      r_w <- r_w * sqrt(n / (n - length(A)))
    }


  result <- list(r_s = r_s, r_w = r_w, mu_s = mu_s, mu_w = mu_w)
  return(result)
}
# #Propose mu hat from 95% region
# PluginMu.MHLS <- function(X, Y, lbd, ratioSeq = seq(0,1,by=0.01), alpha = 0.05
#                           , nChain, niter = 100, method = "boundary"
#                           , parallel = FALSE, ncores = 2) {
#   # method can be either "boundary" or "unif"
#   Sample <- switch(method, boundary = rUnitBall.surface, unif = rUnitBall)
#
#   n <- length(Y)
#   n1 <- round(n/2)
#   n2 <- n - n1
#
#   FF <- function(x) {
#     muhat <- numeric(2 * n)
#     Index <- sort(sample(n, n1))
#     #print(Index)
#     X1 <- X[Index,]
#     X2 <- X[-Index,]
#     Y1 <- Y[Index]
#     Y2 <- Y[-Index]
#     TEMP <- projStein(X1, X2, Y1, Y2, lbd = lbd, ratios = ratioSeq,
#                       alpha=alpha, niter=niter, fixSeed = FALSE)
#     muhat[(1:n)[-Index]] <- TEMP$mu_s + Sample(n2) * TEMP$r_s + TEMP$mu_w + Sample(n2) * TEMP$r_w
#     muhat[(1:n + n)[-Index]] <- TEMP$mu_s + TEMP$mu_w
#     TEMP <- projStein(X2, X1, Y2, Y1, lbd = lbd, ratios = ratioSeq,
#                       alpha=alpha, niter=niter, fixSeed = FALSE)
#     muhat[Index] <- TEMP$mu_s + Sample(n1) * TEMP$r_s + TEMP$mu_w + Sample(n1) * TEMP$r_w
#     muhat[(1:n + n)[Index]] <- TEMP$mu_s + TEMP$mu_w
#     return(muhat)
#   }
#   TEMP <- t(parallel::mcmapply(FF,1:nChain,mc.cores = ncores))
#
#   # returns nChain by p matrix
#   return(rbind(TEMP[,1:n], colMeans(TEMP[,1:n+n])))
#
#   # for (i in 1:nChain) {
#   #   Index <- sort(sample(n, n1))
#   #   #print(Index)
#   #   X1 <- X[Index,]
#   #   X2 <- X[-Index,]
#   #   Y1 <- Y[Index]
#   #   Y2 <- Y[-Index]
#   #   TEMP <- projStein(X1, X2, Y1, Y2, lbd = lbd, ratios = ratioSeq,
#   #                     alpha=alpha, niter=niter, fixSeed = FALSE)
#   #   muhat[i,-Index] <- TEMP$mu_s + rUnitBall(n1) * TEMP$r_s + TEMP$mu_w + rUnitBall(n1) * TEMP$r_w
#   #   TEMP <- projStein(X2, X1, Y2, Y1, lbd = lbd, ratios = ratioSeq,
#   #                     alpha=alpha, niter=niter, fixSeed = FALSE)
#   #   muhat[i,Index] <- TEMP$mu_s + rUnitBall(n2) * TEMP$r_s + TEMP$mu_w + rUnitBall(n2) * TEMP$r_w
#   # }
# }
#
# #library(glmnet)
# ####################################################
#
# # Shrink Y towards zero by Stein estimate
# # Y: n*1 matrix
# # df: degress of freedom
# # return: list(mu: estimate of true mu, L: sure (stein unbiased risk estimate))
#
# ShrinkToZeroStein <- function(Y, df)
# {
#   B <- df / sum(Y^2)
#   return(list(mu = (1 - B) * Y, L = max(1 - B, 0)))
# }
#
# ####################################################
#
# # calculate c(alpha)
# # mu: n*1 matrix
# # df: degress of freedom
# # quantile:
# # proj: identity matrix by default or a n*n matrix
# # niter
# # return: the quantile of c(alpha) in Professor's note
#
# ####################################################
# quantile_stein <- function(mu, df, quantile, proj = NULL, niter = 100)
# {
#   #n = dim(mu)[1]
#   n <- length(mu)
#   if (is.null(proj)) {proj <- diag(n)}
#   z <- numeric(niter)
#   for (i in 1:niter)
#   {
#     Y_prime <- mu + proj %*% matrix(rnorm(n), ncol = 1)
#     res <- ShrinkToZeroStein(Y_prime, df)
#     z[i] <- res$L - sum((res$mu - mu)^2) / df
#   }
#   c <- quantile(abs(z) * sqrt(df), probs = quantile)
#   return(c)
# }
# ####################################################
#
# # We split X and Y into 2 pieces. (X1, X2) (Y1, Y2).
# # Use (X1,Y1) to get the active set A
# # Then draw inference using (X2, Y2).
#
# # Given a sequence of ratios, X, Y, lbd, calculate r_s r_w, mu_s mu_w which minimize volumne
# # X1, X2: n*p matrx
# # Y1, Y2: n*1 matrix
# # lbd: tunning parameter, lambda, in lasso
# # ratios: a sequence of ratios of tau against lambda (tau the threshold value of \hat{beta})
# # alpha:
# # niter:
# # return: list(r_s, r_w, mu_s, mu_w)
# projStein <- function(X1, X2, Y1, Y2, lbd, ratios, alpha = 0.05,
#                       niter = 100, fixSeed = TRUE)
# {
#   # apply lasso regression of Y onto X
#   n1 <- nrow(X1)
#   n2 <- nrow(X2)
#   #lassomodel <- glmnet::glmnet(X1, Y1, family = 'gaussian', intercept = F, lambda = lbd)
#
#   #lassomodel <- glmnet::glmnet(X1, Y1, family = 'gaussian', intercept = F, lambda = lbd,
#   #                              standardize = FALSE)$beta
#   B0 <- lassoFit(X = X1, Y = Y1, type = "lasso", lbd = lbd)$B0
#
#   # estimate c by A=0
#   stein <- ShrinkToZeroStein(Y2, n2)
#   if (fixSeed) {set.seed(123)}
#   c <- quantile_stein(stein$mu, n2, 1-alpha/2, niter = niter) # Why 1-alpha/2 instead of 1-alpha? Bonferroni
#
#   # information along ratios
#   r_s_ratio <- r_w_ratio <- logVol_ratio <- numeric(length(ratios))
#
#   for (i in 1:length(ratios))
#   {
#     A <- which(abs(B0) > ratios[i] * lbd)
#     if (length(A) != 0)
#     {
#       r_s_ratio[i] <- sqrt(qchisq(1-alpha/2, df = length(A)) / n2)
#       P_A <- X2[,A] %*% solve(t(X2[,A]) %*% X2[,A]) %*% t(X2[,A])
#       P_perp <- diag(n2) - P_A
#     } else
#     {
#       r_s_ratio[i] <- 0
#       P_perp <- diag(n2)
#     }
#     nMinusA <- n2 - length(A)
#
#     stein <- ShrinkToZeroStein(P_perp %*% Y2, nMinusA)
#     r_w_ratio[i] <- sqrt(nMinusA / n2 * (c / sqrt(nMinusA) + stein$L))
#
#     #update r_s_ratio r_w_ratio if r_s_ratio exists and calculate volume
#     if (length(A) != 0)
#     {
#       r_s_ratio[i] <- r_s_ratio[i] * sqrt(n2 / length(A))
#       r_w_ratio[i] <- r_w_ratio[i] * sqrt(n2 / (n2 - length(A)))
#       logVol_ratio[i] <- length(A) * log(r_s_ratio[i]) + nMinusA * log(r_w_ratio[i])
#     } else
#     {
#       logVol_ratio[i] <- nMinusA * log(r_w_ratio[i])
#     }
#   }
#   # information of minimum volume
#   i <- which.min(logVol_ratio)
#   A <- which(abs(B0) > ratios[i] * lbd)
#   # compute mu_s
#   if (length(A) != 0)
#   {
#     P_A <- X2[,A] %*% solve(t(X2[,A]) %*% X2[,A]) %*% t(X2[,A])
#     mu_s <- P_A %*% Y2
#     P_perp <- diag(n2) - P_A
#   } else
#   {
#     mu_s <- 0
#     P_perp <- diag(n2)
#   }
#   # compute mu_w
#   stein <- ShrinkToZeroStein(P_perp %*% Y2, n2 - length(A))
#   mu_w <- stein$mu
#
#   result <- list(r_s = r_s_ratio[i], r_w = r_w_ratio[i], mu_s = mu_s, mu_w = mu_w)
#   return(result)
# }


# Refit the beta estimator to remove the bias
Refit.MHLS <- function(X,weights,lbd,MHsample) {
  # W : diag(weights)
  # MHsample : MH samples from MHLS function
  n <- nrow(X)
  # Active set
  A <- which(MHsample$beta[1,]!=0)
  # Recover y using KKT condition
  Y.MH <- solve(X%*%t(X))%*%X %*% (crossprod(X) %*% t(MHsample$beta) / n +
                                     lbd * weights * t(MHsample$subgrad)) * n
  # Refit y to the restricted set of X
  beta.refit <- solve(t(X[,A])%*%X[,A])%*%t(X[,A]) %*% Y.MH
  return(beta.refit)
}

# Generate 1-alpha Confidence Interval based on the deviation
CI.MHLS <- function(betaRefitMH, betaCenter, betaRefit, alpha=.05) {
  # betaCenter : ginv(X[,A]) %*% X %*% beta0 or ginv(X[,A]) %*% hatMu
  #               , vector of length |A|
  # betaRefitMH : refitted beta via Refit.MHLS, a niter x |A| matrix.
  # betaRefit : refitted beta from original data, a vector with length p.
  # alpha : significant level.
  A <- which(betaRefit!=0)
  Quantile <- apply(betaRefitMH - betaCenter[A], 1, function(x)
  {quantile(x,prob=c(alpha/2, 1 - alpha/2))})
  Result <- rbind(LowerBound = -Quantile[2,] + betaRefit[A],
                  UpperBound = -Quantile[1,] + betaRefit[A])
  colnames(Result) <- paste("beta", A, sep="")
  return(Result)
}

# Generate pluginbeta's from 95% confidence region
Pluginbeta.MHLS <- function(X,Y,A,nPlugin,sigma.hat) {
  # nPlugin : number of pluginbeta's want to generate
  # sigma.hat : estimator of sigma , \epsilon ~ N(0, sigma^2)
  #             If missing, use default way to generate it
  if (length(A)==0) stop("The lasso solution has an empty active set.")

  if (missing(sigma.hat)) {
    sigma.hat <- summary((lm(Y~X[,A]+0)))$sigma
  }
  betaRefit <- coef(lm(Y~X[,A]+0))

  if (nPlugin == 1) {
    return(matrix(betaRefit,1))
  } else {
    xy <- matrix(rnorm(length(A)*(nPlugin)), nPlugin)
    lambda <- 1 / sqrt(rowSums(xy^2))
    xy <- xy * lambda * sqrt(qchisq(0.95, df=length(A)))
    coeff.seq <- matrix(0,nPlugin,ncol(X))
    coeff.seq[,A]  <- rbind(t(t(xy%*%chol(solve(crossprod(X[,A])))) *
                                           sigma.hat + betaRefit))
    return(coeff.seq=coeff.seq)
  }
}

# Compute group-wize max eigen-value
groupMaxEigen <- function(X, group) {
  Psi <- crossprod(X) / nrow(X)
  J <- max(group)
  Gamma <- numeric(J)
  for (i in 1:J) {
    Gamma[i] <- max(eigen(Psi[group==i, group==i], only.values= TRUE)$values)
  }
  Gamma <- Gamma * (1 + 1e-10)
  return(Gamma)
}

