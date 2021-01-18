#CorrHetSNSol <- function(Nres,SNres,CvCase,Conf,Xscld,Xmean,Xsd,XsdOutPrd,grouping,Mxt,lglikdif,limlnk2,OptCntrl,bordertol=1e-2,maxsk=0.99527)
CorrHetSNSol <- function(Nres,SNres,CvCase,Conf,Xscld,Xmean,Xsd,XsdOutPrd,grouping,Mxt,lglikdif,limlnk2,OptCntrl,getvcov=TRUE,bordertol=1e-2,maxsk=0.99527)
{

  GetvcovHetD <- function(Res)
  {
    HoMxtparnam <- paste("mu_",Xnames,"_",grplvls[1],sep="")
    for (g in 2:k) HoMxtparnam <- c(HoMxtparnam,paste("mu_",Xnames,"_",grplvls[g],sep=""))
    for (i in 1:p) HoMxtparnam <- c(HoMxtparnam,paste("Sigma_",Xnames[i],"_",Xnames[i:p],sep=""))
    HoMxtparnam <- c(HoMxtparnam,paste("gamma1_",Xnames,sep=""))
    npar <- SKnpar(Conf,p,p/2,Ngrps=k)
    grpModMat <- model.matrix(~ grouping)
    Gmat <- model.matrix(~ grplvls)
    InFData <- try( sn.infoMv( dp=list(beta=as.matrix(rbind.data.frame(Res$ksi[1,],Res$beta2k)),Omega=Res$Omega,alpha=Res$alpha), y=Xscld, x=grpModMat ) )
    if ( is.null(InFData) || class(InFData)[1] == "try-error" || is.null(InFData$asyvar.cp) )  {
      return( list(mleCPvcov=NULL,muEse=NULL,SigmaEse=NULL,gamma1Ese=NULL,status="Invalid") )
    }
    if (Conf==1)  {
      betavcov <- InFData$asyvar.cp
      nSpar <- p*(p+1)/2
    }  else  {
      nparC1 <- (k+1)*p + p*(p+1)/2  
      parind <- c(1:(k*p),(k*p)+SigCind(Conf,p/2),(nparC1-p+1):nparC1)
      nSpar <- length(parind) - (k+1)*p
      betavcov <- Safepdsolve(InFData$info.cp[parind,parind],maxlnk2=limlnk2,scale=TRUE)
    }
    if ( !is.null(betavcov) )
    {
      mleCPvcov <- matrix(nrow=npar,ncol=npar)		
      muind <- 1:(p*k)
      nmugind <- NULL
      for (g in 1:k) nmugind <- c(nmugind,(0:(p-1))*k+g)
      Sind <- p*k + 1:nSpar
      gamma1ind <- p*k + nSpar + 1:p
      Sgamma1ind <- c(Sind,gamma1ind)
      M <- kronecker(diag(p),Gmat)
      mleCPvcov[muind,muind] <- (M %*% betavcov[muind,muind] %*% t(M))[nmugind,nmugind]
      mleCPvcov[muind,Sgamma1ind] <- (M %*% betavcov[muind,Sgamma1ind])[nmugind,]
      mleCPvcov[Sgamma1ind,muind] <- t(mleCPvcov[muind,Sgamma1ind])
      mleCPvcov[Sgamma1ind,Sgamma1ind] <- betavcov[Sgamma1ind,Sgamma1ind]
      if (Conf==1)  {
        rownames(mleCPvcov) <- colnames(mleCPvcov) <- HoMxtparnam
      } else  {
        rownames(mleCPvcov) <- colnames(mleCPvcov) <- HoMxtparnam[parind]
      }
    } else { 
      mleCPvcov <- NULL
    }

    if (is.null(mleCPvcov))  return( list(mleCPvcov=NULL,muEse=NULL,SigmaEse=NULL,gamma1Ese=NULL) )
    mleCPvcov <- mleCPvcov * SNVCovscaling(Conf,p,Xsd,k=k) # scale back vcov matrix   
    CPStderr <- sqrt(diag(mleCPvcov))
    muEse <- matrix(CPStderr[1:(k*p)],nrow=k,ncol=p,byrow=TRUE)
    gammaind <- (npar-p+1):npar
    gamma1Ese <- CPStderr[gammaind]
    SigmaEse <- matrix(nrow=p,ncol=p) 
    cnt <- k*p
    for (j1 in 1:p) for (j2 in j1:p)
    {   
      if (FreePar(q,j1,j2,Conf))
      {
        cnt <- cnt+1
        SigmaEse[j1,j2] <- SigmaEse[j2,j1] <- CPStderr[cnt]
      }
    }
    names(gamma1Ese) <- rownames(SigmaEse) <- colnames(SigmaEse) <- Xnames      
    rownames(muEse) <- grplvls
    colnames(muEse) <- Xnames

    list(mleCPvcov=mleCPvcov,muEse=muEse,SigmaEse=SigmaEse,gamma1Ese=gamma1Ese,status="Regular")
  }

  GetvcovSingD <- function(Res)
  {
    SngDparnam <- paste("mu_",Xnames,sep="")
    for (i in 1:p)  SngDparnam <- c(SngDparnam,paste("Sigma_",Xnames[i],"_",Xnames[i:p],sep=""))
    SngDparnam <- c(SngDparnam,paste("gamma1_",Xnames,sep=""))
    npar <- SKnpar(Conf,p,p/2)
    InFData <- try( sn.infoMv( dp=list(xi=Res$ksi,Omega=Res$Omega,alpha=Res$alpha), y=Xscld, x=matrix(1,nrow=n,ncol=1) ) )
    if ( is.null(InFData) || class(InFData)[1] == "try-error" || is.null(InFData$asyvar.cp) )  {
      return( list(mleCPvcov=NULL,muEse=NULL,SigmaEse=NULL,gamma1Ese=NULL,status="Invalid") )
    }
    if (Conf==1)  {
      mleCPvcov <- InFData$asyvar.cp
      rownames(mleCPvcov) <- colnames(mleCPvcov) <- SngDparnam
    }  else  {
      nparC1 <- 2*p + p*(p+1)/2  
      parind <- c(1:p,p+SigCind(Conf,p/2),(nparC1-p+1):nparC1)
      mleCPvcov <- Safepdsolve(InFData$info.cp[parind,parind],maxlnk2=limlnk2,scale=TRUE)
      if ( !is.null(mleCPvcov) ) { rownames(mleCPvcov) <- colnames(mleCPvcov) <- SngDparnam[parind] }
    }
    if (is.null(mleCPvcov))  return( list(mleCPvcov=NULL,muEse=NULL,SigmaEse=NULL,gamma1Ese=NULL) )
    mleCPvcov <- mleCPvcov * SNVCovscaling(Conf,p,Xsd,k=k) # scale back vcov matrix   
    CPStderr <- sqrt(diag(mleCPvcov))
    muEse <- CPStderr[1:p]
    gammaind <- (npar-p+1):npar
    gamma1Ese <- CPStderr[gammaind]
    SigmaEse <- matrix(nrow=p,ncol=p) 
    cnt <- k*p
    for (j1 in 1:p) for (j2 in j1:p)
    {   
      if (FreePar(q,j1,j2,Conf))
      {
        cnt <- cnt+1
        SigmaEse[j1,j2] <- SigmaEse[j2,j1] <- CPStderr[cnt]
      }
    }
    names(gamma1Ese) <- rownames(SigmaEse) <- colnames(SigmaEse) <- names(muEse) <- Xnames      

    list(mleCPvcov=mleCPvcov,muEse=muEse,SigmaEse=SigmaEse,gamma1Ese=gamma1Ese,status="Regular")
  }


##########################################################################
### End of auxiliary functions
### Start of main CorrHetSNSol function
##########################################################################

  p <- ncol(Xscld)
  q <- p/2 
  n <- nrow(Xscld)
  grplvls <- levels(grouping)
  k <- length(grplvls)
  n1scvct <- rep(1,n)
  Xnames <- names(Xmean)

  zmu <- scale(Nres@mleNmuE,center=Xmean,scale=Xsd)
  displcmnt <- -Xmean
  sclfactor <- 1./Xsd
  NewSNres <- list()

  if (Mxt=="Hom" || Mxt=="Loc")  {
    zmug1 <- zmu[1,]
    if (k>2) {
      beta2k <- scale(zmu[-1,],center=zmug1,scale=FALSE)
    } else {
      beta2k <- zmu[-1,] - zmug1
    }
    zSigma <- Nres@CovConfCases[[CvCase]]$mleSigE/XsdOutPrd 
    if (Conf==1) {
      DP <- cnvCPtoDP(p,zmu[1,],zSigma,rep(0.,p),limlnk2=limlnk2) # check if this is correct (it not -- we need the betas !!!
      newpar <- c(DP$ksi,beta2k,DP$alpha/DP$omega)
      SNStdDtRes <- SNCnf1MaxLik(Xscld,initpar=newpar,grouping=grouping,limlnk2=limlnk2,OptCntrl=OptCntrl)
    } else {
      SigmaSrpar <- GetCovPar(zSigma,Conf,test=FALSE) 
      newpar <- c(zmug1,beta2k,SigmaSrpar,rep(0.,p))
      SNStdDtRes <- SNCMaxLik(Xscld,Config=Conf,initpar=newpar,grouping=grouping,limlnk2=limlnk2,OptCntrl=OptCntrl)    
    }
    NewSNres$muE <- scale( scale(SNStdDtRes$mu,center=FALSE,scale=sclfactor),center=displcmnt,scale=FALSE )
    NewSNres$ksiE <- scale( scale(SNStdDtRes$ksi,center=FALSE,scale=sclfactor),center=displcmnt,scale=FALSE )
    attr(NewSNres$muE,"scaled:center") <- attr(NewSNres$ksiE,"scaled:center") <-
      attr(NewSNres$muE,"scaled:scale") <- attr(NewSNres$ksiE,"scaled:scale") <- NULL
    NewSNres$SigmaE <- XsdOutPrd*SNStdDtRes$Sigma
    NewSNres$gamma1E <- SNStdDtRes$gamma1
    NewSNres$OmegaE <- XsdOutPrd*SNStdDtRes$Omega
    NewSNres$alphaE <- SNStdDtRes$alpha
    NewSNres$logLik <- SNStdDtRes$lnLik + lglikdif

    colnames(NewSNres$muE) <-  colnames(NewSNres$ksiE) <- Xnames
    rownames(NewSNres$muE) <-  rownames(NewSNres$ksiE) <- grplvls
    if ( !is.null(NewSNres$gamma1E) && !is.null(NewSNres$alphaE) && !is.null(NewSNres$SigmaE) && !is.null(NewSNres$OmegaE) )
    {
      names(NewSNres$gamma1E) <- names(NewSNres$alphaE) <-
      dimnames(NewSNres$SigmaE)[[1]] <- dimnames(NewSNres$SigmaE)[[2]] <- 
      dimnames(NewSNres$OmegaE)[[1]] <- dimnames(NewSNres$OmegaE)[[2]] <- Xnames
    }  

  } else if (Mxt=="Het" || Mxt=="Gen")  {
    NewSNres$muE <- matrix(nrow=k,ncol=p)
    NewSNres$gamma1E <- matrix(nrow=k,ncol=p)
    NewSNres$ksiE <- matrix(nrow=k,ncol=p)
    NewSNres$alphaE <- matrix(nrow=k,ncol=p)
    nparbyg <- SKnpar(Conf,p,q)
    anams <- list(Xnames,Xnames,grplvls)
    mnams <- list(grplvls,Xnames)
    NewSNres$SigmaE <- array(dim=c(p,p,k),dimnames=anams)
    NewSNres$OmegaE <- array(dim=c(p,p,k),dimnames=anams)
    NewSNres$logLik <- 0. 

    if (getvcov) {
      NewSNres$status <- "Regular"
      NewSNres$mleCPvcov <- array(dim=c(nparbyg,nparbyg,k),dimnames=list(NULL,NULL,grplvls))
      NewSNres$muEse <- matrix(nrow=k,ncol=p,dimnames=mnams)
      NewSNres$SigmaEse <- array(dim=c(p,p,k),dimnames=anams)
      NewSNres$gamma1Ese <- matrix(nrow=k,ncol=p,dimnames=mnams)
     } else {
       NewSNres$status <- "OnHold"
       NewSNres$mleCPvcov <- NewSNres$muEse <- NewSNres$SigmaEse <- NewSNres$gamma1Ese <- NULL
     }

    for (g in 1:k) {
      Xscldg <- Xscld[grouping==grplvls[g],] 
      zSigma <- Nres@CovConfCases[[CvCase]]$mleSigE[,,g]/XsdOutPrd 
      if (Conf==1) {
        DP <- cnvCPtoDP(p,zmu[g,],zSigma,rep(0.,p),limlnk2=limlnk2) # check if this is correct !!!
        newpar <- c(DP$ksi,DP$alpha/DP$omega)
        SNStdDtRes <- SNCnf1MaxLik(Xscldg,initpar=newpar,grouping=NULL,limlnk2=limlnk2,OptCntrl=OptCntrl)
      } else {
        SigmaSrpar <- GetCovPar(zSigma,Conf,test=FALSE) 
        newpar <- c(zmu[g,],SigmaSrpar,rep(0.,p))
        SNStdDtRes <- SNCMaxLik(Xscldg,Config=Conf,initpar=newpar,grouping=NULL,limlnk2=limlnk2,OptCntrl=OptCntrl)
      }
      NewSNres$muE[g,] <- SNStdDtRes$mu/sclfactor-displcmnt
      NewSNres$ksiE[g,] <- SNStdDtRes$ksi/sclfactor-displcmnt
      NewSNres$SigmaE[,,g] <- XsdOutPrd*SNStdDtRes$Sigma
      NewSNres$gamma1E[g,] <- SNStdDtRes$gamma1
      NewSNres$OmegaE[,,g] <- XsdOutPrd*SNStdDtRes$Omega
      if (!is.null(SNStdDtRes$alpha)) NewSNres$alphaE[g,] <- SNStdDtRes$alpha
      NewSNres$logLik <- NewSNres$logLik + SNStdDtRes$lnLik

      if (getvcov) {
        if ( (!is.null(SNStdDtRes$c2) && SNStdDtRes$c2 > bordertol) || (maxsk-max(abs(SNStdDtRes$gamma1)) < bordertol) )
        {
          vcovl <- GetvcovSingD(SNStdDtRes)
          if (!is.null(vcovl$muEse))  NewSNres$muEse[g,] <- vcovl$muEse 
          if (!is.null(vcovl$gamma1Ese))NewSNres$gamma1Ese[g,] <- vcovl$gamma1Ese 
          if (!is.null(vcovl$SigmaEse)) NewSNres$SigmaEse[,,g] <- vcovl$SigmaEse
          if (!is.null(vcovl$mleCPvcov))  {
            NewSNres$mleCPvcov[,,g] <- vcovl$mleCPvcov
            if (g==1) dimnames(NewSNres$mleCPvcov)[[1]]  <- dimnames(NewSNres$mleCPvcov)[[2]] <- rownames(NewSNres$mleCPvcov)
          }
        } else {
          NewSNres$status <- "Onborder"
        }
     }

    }
    NewSNres$logLik <- NewSNres$logLik + lglikdif
    colnames(NewSNres$muE) <-  colnames(NewSNres$ksiE) <- Xnames
    rownames(NewSNres$muE) <-  rownames(NewSNres$ksiE) <- grplvls
    if ( !is.null(NewSNres$gamma1E) && !is.null(NewSNres$alphaE) && !is.null(NewSNres$SigmaE) && !is.null(NewSNres$OmegaE) )
    {
      names(NewSNres$gamma1E) <- names(NewSNres$alphaE) <-
        dimnames(NewSNres$SigmaE)[[1]] <- dimnames(NewSNres$SigmaE)[[2]] <- 
        dimnames(NewSNres$OmegaE)[[1]] <- dimnames(NewSNres$OmegaE)[[2]] <- Xnames
      dimnames(NewSNres$SigmaE)[[3]] <- dimnames(NewSNres$OmegaE)[[3]] <- grplvls
    }  
  }


  NewSNres  #  return(NewSNres)
}
