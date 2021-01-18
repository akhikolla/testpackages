#CorrSNSol <- function(Nres,SNres,CvCase,Conf,Xscld,Xmean,Xsd,XsdOutPrd,lglikdif,limlnk2,OptCntrl,bordertol=1e-2,maxsk=0.99527)
CorrSNSol <- function(Nres,SNres,CvCase,Conf,Xscld,Xmean,Xsd,XsdOutPrd,lglikdif,limlnk2,OptCntrl,getvcov=TRUE,bordertol=1e-2,maxsk=0.99527)
{
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
    mleCPvcov <- mleCPvcov * SNVCovscaling(Conf,p,Xsd,k=1) # scale back vcov matrix   
    CPStderr <- sqrt(diag(mleCPvcov))
    muEse <- CPStderr[1:p]
    gammaind <- (npar-p+1):npar
    gamma1Ese <- CPStderr[gammaind]
    SigmaEse <- matrix(nrow=p,ncol=p) 
    cnt <- p
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
### Start of main CorrSNSol function
##########################################################################

  p <- ncol(Xscld)
  q <- p/2 
  n <- nrow(Xscld)
  n1scvct <- rep(1,n)
  Xnames <- names(Xmean)

  zmu <- (Nres@mleNmuE-Xmean)/Xsd
  if (Conf==1) {
    zSigma <- Nres@CovConfCases[[CvCase]]$mleSigE/XsdOutPrd 
    DP <- cnvCPtoDP(p,zmu,zSigma,rep(0.,p),limlnk2=limlnk2)
    newpar <- c(DP$ksi,DP$alpha/DP$omega)
#    newpar <- c(zmu,rep(0.,p))  rep(0.,p) # This should replace the three previous lines with the same results, but faster -- check it !!!
    SNStdDtRes <- SNCnf1MaxLik(Xscld,initpar=newpar,grouping=NULL,limlnk2=limlnk2,OptCntrl=OptCntrl)
  } else {
    zSigma <- Nres@CovConfCases[[CvCase]]$mleSigE/XsdOutPrd
    SigmaSrpar <- GetCovPar(zSigma,Conf,test=FALSE) 
    newpar <- c(zmu,SigmaSrpar,rep(0.,p))
    SNStdDtRes <- SNCMaxLik(Xscld,Config=Conf,initpar=newpar,grouping=NULL,limlnk2=limlnk2,OptCntrl=OptCntrl)
  }

  NewSNres <- list()
  NewSNres$muE <- Xmean + Xsd*SNStdDtRes$mu
  NewSNres$ksiE <- Xmean + Xsd*SNStdDtRes$ksi   
  NewSNres$SigmaE <- XsdOutPrd*SNStdDtRes$Sigma
  NewSNres$gamma1E <- SNStdDtRes$gamma1
  NewSNres$OmegaE <- XsdOutPrd*SNStdDtRes$Omega
  NewSNres$alphaE <- SNStdDtRes$alpha
  NewSNres$logLik <- SNStdDtRes$lnLik + lglikdif

  names(NewSNres$muE) <-  names(NewSNres$ksiE) <- Xnames
  if ( !is.null(NewSNres$gamma1E) && !is.null(NewSNres$alphaE) && !is.null(NewSNres$SigmaE) && !is.null(NewSNres$OmegaE) )
  {
    names(NewSNres$gamma1E) <- names(NewSNres$alphaE) <-
    dimnames(NewSNres$SigmaE)[[1]] <- dimnames(NewSNres$SigmaE)[[2]] <- 
    dimnames(NewSNres$OmegaE)[[1]] <- dimnames(NewSNres$OmegaE)[[2]] <- Xnames
  }  

  if (getvcov) {
    if ( (SNStdDtRes$c2 > bordertol) || (maxsk-max(abs(SNStdDtRes$gamma1)) < bordertol) )
    {
      vcovl <- GetvcovSingD(SNStdDtRes)
      NewSNres$status <- vcovl$status
      NewSNres$mleCPvcov <- vcovl$mleCPvcov
      NewSNres$muEse <- vcovl$muEse
      NewSNres$SigmaEse <- vcovl$SigmaEse
      NewSNres$gamma1Ese <- vcovl$gamma1Ese
    } else {
      NewSNres$status <- "Onborder"
    }
  } else {
    NewSNres$status <- "OnHold"
    NewSNres$mleCPvcov <- NewSNres$muEse <- NewSNres$SigmaEse <- NewSNres$gamma1Ese <- NULL
  }

  NewSNres  #  return(NewSNres)
}
