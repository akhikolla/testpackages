error <- function(onerror,msg)     # Things to do: replace this by an exception-based mechanism !!
{
  if (onerror=="stop") {
    stop(msg)
  }  else if (onerror=="warning")  { 
    warning(msg)
    return(NULL)
  }  else if (onerror=="silentNull")  {
    return(NULL)
  }
} 

IdtNmle <- function(Idt, grouping=NULL, Type=c("SingDst","HomMxt"), CVtol=1.0e-5, limlnk2=log(1e8),
  OptCntrl=list(), onerror=c("stop","warning","silentNull"), CovCaseArg, Config, SelCrit)
{                                   
  onerror <- match.arg(onerror)
  Type <- match.arg(Type)
  q <- Idt@NIVar
  p <- 2*q
  n <- Idt@NObs
  if (Type=="SingDst")  { 
    k <- 1
  }  else if (Type=="HomMxt")  {
    if (is.null(grouping))  { error(onerror,"Argument grouping is missing from MANOVA method\n") }
    if (!is.factor(grouping))  { error(onerror,"'grouping' is not a factor") }
    nk <- as.numeric(table(grouping))
    if (n != sum(nk))  {
      error(onerror,"Number of observations in IData object and grouping factor do not agree with each other\n")
    }
    k <- length(nk) 
    if (k==1)  {
      error(onerror,"The data belongs to one single group. A partition into at least two different groups is required\n")
    }
  }

  if (CovCaseArg)  {
    nCovCases <- 4 
    CovCaseMap <- c(1,NA,2,3,4)
    modnames <- paste("NModCovC",1:4,sep="")
  } else {
    nCovCases <- 5 
    CovCaseMap <- 1:5
    modnames <- paste("NC",1:5,sep="")
  }
  logLiks <- AICs <- BICs <- rep(NA_real_,nCovCases)
  CovConfCases <- vector("list",nCovCases)
  names(logLiks) <- names(AICs) <- names(BICs) <- names(CovConfCases) <- modnames
  X <- cbind.data.frame(Idt@MidP,Idt@LogR)

  for (model in Config)  {  
    if (model!=2) { 
      CovConfCases[[CovCaseMap[model]]] <- list(mleSigE=NULL,mleSigEse=NULL,mlevcov=NULL,logLik=NULL,AIC=NULL,BIC=NULL)
    }  else   {
      CovConfCases[[2]] <- list(mleSigE=NULL,mleSigEse=NULL,mlevcov=NULL,logLik=NULL,AIC=NULL,BIC=NULL,optres=NULL)
    }
  }
  if (Type=="SingDst")
  {
    mleNmuE <- colMeans(X)
    mleNsdE <- sapply(X,sd)
    mleNmuEse <- mleNsdE/sqrt(n)
    if (is.element(1,Config))
    {
      CovConfCases[[1]]$mleSigE <- var(X) * (n-1)/n
      XVar <- diag(CovConfCases[[1]]$mleSigE)
      CnstV <- which(mleNsdE/abs(mleNmuE)<CVtol)
      if (length(CnstV)==1)  { 
        error(onerror,paste("Variable",names(CnstV),"appears to be constant \n"))
      }  else if (length(CnstV)>0)  {
        error( onerror,paste("Variables",paste(names(CnstV),collapse=" "),"appear to be constant\n") )
      }
      CovConfCases[[1]]$mlevcov <- NmleVCov(n,p,CovConfCases[[1]]$mleSigE)
      CovConfCases[[1]]$mleSigEse  <- sqrt((CovConfCases[[1]]$mleSigE^2+outer(XVar,XVar))/n)
      dimnames(CovConfCases[[1]]$mleSigEse) <- dimnames(CovConfCases[[1]]$mleSigE)
    }  else {
      SigNEC1 <- var(X) * (n-1)/n
      XVar <- diag(SigNEC1)
      CnstV <- which(mleNsdE/abs(mleNmuE)<CVtol)
      if (length(CnstV)==1)  {
        error(onerror,paste("Variable",names(CnstV),"appears to be constant \n"))
      }  else if (length(CnstV)>0)  {
        error(onerror,paste("Variables",paste(names(CnstV),collapse=" "),"appear to be constant\n"))
      }
      vcovNEC1 <- NmleVCov(n,p,SigNEC1)
      SigNEseC1 <- sqrt((SigNEC1^2+outer(XVar,XVar))/n)
      dimnames(SigNEseC1) <- dimnames(SigNEC1)
    }
  }
  else if (Type=="HomMxt")
  {
    mleNmuE <- apply(X,2,function(x) tapply(x,grouping,mean))
    mleNmuEse <- apply(X,2,function(x) tapply(x,grouping,sd))/sqrt(nk)
    rownames(mleNmuE) <- rownames(mleNmuEse) <- levels(grouping)
    if (is.element(1,Config))
    {
      CovConfCases[[1]]$mleSigE <- var(X[grouping==levels(grouping)[1],]) * (nk[1]-1)/n
      for (g in 2:k)
        CovConfCases[[1]]$mleSigE <- CovConfCases[[1]]$mleSigE + var(X[grouping==levels(grouping)[g],]) * (nk[g]-1)/n
      XVar <- diag(CovConfCases[[1]]$mleSigE)
      CnstV <- which(XVar/mleNmuE^2<CVtol^2)
      if (length(CnstV)==1)  { 
        error(onerror,paste("Variable",names(CnstV),"appears to be constant within groups\n"))
      }  else if (length(CnstV)>0)  { 
        error(onerror,paste("Variables",paste(names(CnstV),collapse=" "),"appear to be constant within groups\n"))
      }
      CovConfCases[[1]]$mlevcov <- NmleVCov(nk,p,CovConfCases[[1]]$mleSigE,k,grpnames=levels(grouping))
      CovConfCases[[1]]$mleSigEse  <- sqrt((CovConfCases[[1]]$mleSigE^2+outer(XVar,XVar))/n)
      dimnames(CovConfCases[[1]]$mleSigEse) <- dimnames(CovConfCases[[1]]$mleSigE)
    }  else  {
      SigNEC1 <- var(X[grouping==levels(grouping)[1],]) * (nk[1]-1)/n
      for (g in 2:k) SigNEC1 <- SigNEC1 + var(X[grouping==levels(grouping)[g],]) * (nk[g]-1)/n
      XVar <- diag(SigNEC1)
      CnstV <- which(XVar/mleNmuE^2<CVtol^2)
      if (length(CnstV)==1)  { 
        error(onerror,paste("Variable",names(CnstV),"appears to be constant within groups\n"))
      }  else if (length(CnstV)>0)  { 
        error(onerror,paste("Variables",paste(names(CnstV),collapse=" "),"appear to be constant within groups\n"))
      }
      vcovNEC1 <- NmleVCov(nk,p,SigNEC1,k,grpnames=levels(grouping))
      SigNEseC1  <- sqrt((SigNEC1^2+outer(XVar,XVar))/n)
      dimnames(SigNEseC1) <- dimnames(SigNEC1)
    }
  }
  if (any(Config!=2))
  {
    Cnf35 <- Config[Config>2]
    if (length(Cnf35)!=0)
    {  
      for (Conf in Cnf35)
      {
        if (is.element(1,Config))
        {
          CovConfCases[[CovCaseMap[Conf]]]$mleSigE <- mleNSigC35(Conf,CovConfCases[[1]]$mleSigE,q,p,0.)
          CovConfCases[[CovCaseMap[Conf]]]$mleSigEse <- mleNSigC35(Conf,CovConfCases[[1]]$mleSigEse,q,p,NA)
          CovConfCases[[CovCaseMap[Conf]]]$mlevcov <- mleNvcovC35(Conf,CovConfCases[[1]]$mlevcov,p,k)
        } else {
          CovConfCases[[CovCaseMap[Conf]]]$mleSigE <- mleNSigC35(Conf,SigNEC1,q,p,0.)
          CovConfCases[[CovCaseMap[Conf]]]$mleSigEse <- mleNSigC35(Conf,SigNEseC1,q,p,NA)
          CovConfCases[[CovCaseMap[Conf]]]$mlevcov <- mleNvcovC35(Conf,vcovNEC1,p,k)
        }
      }
    }
    if (Type=="SingDst")  {
      Xdev <- scale(X,scale=FALSE)
    }  else if (Type=="HomMxt")  {
      Xdev <- scale(X[grouping==levels(grouping)[1],],center=mleNmuE[1,],scale=FALSE)
      for (g in 2:k)
        Xdev <- rbind.data.frame(Xdev,scale(X[grouping==levels(grouping)[g],],center=mleNmuE[g,],scale=FALSE))
    }
    Conf134 <- Config[Config !=2 & Config!=5]
    if (length(Conf134)!=0)  {
      for (Conf in Conf134)
      {
	CvCase <- CovCaseMap[Conf]	
        logdet <- pdwt.solve(CovConfCases[[CvCase]]$mleSigE,silent=TRUE,onlylogdet=TRUE)
        if (is.null(logdet))  {
          logLiks[CvCase] <- CovConfCases[[CvCase]]$logLik <- -Inf
        }  else  {
          logLiks[CvCase] <- CovConfCases[[CvCase]]$logLik <- -n*(p*(log(2*pi)+1)+logdet)/2
        }		
	AICs[CvCase] <- CovConfCases[[CvCase]]$AIC <- -2*CovConfCases[[CvCase]]$logLik + 2*npar(Conf,p,q,Ngrps=k)
        BICs[CvCase] <- CovConfCases[[CvCase]]$BIC <- -2*CovConfCases[[CvCase]]$logLik + log(n)*npar(Conf,p,q,Ngrps=k)
      }
    }
    if (is.element(5,Config))
    {
      CvCase <- CovCaseMap[5]	
      logLiks[CvCase] <- CovConfCases[[CvCase]]$logLik <- -n*(p*(log(2*pi)+1)+sum(log(XVar)))/2 
      AICs[CvCase] <- CovConfCases[[CvCase]]$AIC <- -2*CovConfCases[[CvCase]]$logLik + 2*npar(5,p,q,Ngrps=k)
      BICs[CvCase] <- CovConfCases[[CvCase]]$BIC <- -2*CovConfCases[[CvCase]]$logLik + log(n)*npar(5,p,q,Ngrps=k)
    }
  } 
  if (is.element(2,Config))
  {
    Xsd <- sqrt(XVar)
    lglikdif <- -n*sum(log(Xsd)) 
    if (Type=="SingDst")  {
      Xscld <- scale(X,scale=Xsd)
    }  else if (Type=="HomMxt")  {
      Xscld <- scale(X[grouping==levels(grouping)[1],],center=mleNmuE[1,],scale=Xsd)
      for (g in 2:k)
        Xscld <- rbind.data.frame(Xscld,scale(X[grouping==levels(grouping)[g],],center=mleNmuE[g,],scale=Xsd))
    }
    C2res <- Cnf2MaxLik(Xscld,OptCntrl=OptCntrl)
    if ( is.element(4,Config) && C2res$lnLik < logLiks[CovCaseMap[4]] )  {
      C2res <- Cnf2MaxLik(Xscld,initpar=initparconf2(CovConfCases[[CovCaseMap[4]]]$mleSigE,n,q))
    }
    if ( is.element(3,Config) && C2res$lnLik < logLiks[CovCaseMap[3]] )  {
      C2res <- Cnf2MaxLik(Xscld,initpar=initparconf2(CovConfCases[[CovCaseMap[3]]]$mleSigE,n,q))
    }    
    CovConfCases[[2]]$mleSigE <- C2GetCov(C2res$SigmaSr,outer(Xsd,Xsd),q)
    CovConfCases[[2]]$mleSigEse <- C2GetCovStderr(CovConfCases[[2]]$mleSigE,X,q,limlnk2=limlnk2,ue=colMeans(X))
    dimnames(CovConfCases[[2]]$mleSigEse) <- dimnames(CovConfCases[[2]]$mleSigE)
    logLiks[2] <- CovConfCases[[2]]$logLik <- C2res$lnLik + lglikdif 
    AICs[2] <- CovConfCases[[2]]$AIC <- -2*CovConfCases[[2]]$logLik + 2*npar(2,p,q,Ngrps=k)
    BICs[2] <- CovConfCases[[2]]$BIC <- -2*CovConfCases[[2]]$logLik + log(n)*npar(2,p,q,Ngrps=k)
    CovConfCases[[2]]$optres <- C2res$optres
  }
  if (SelCrit=="AIC")  {
    bestmod <- which.min(AICs)
  }  else if (SelCrit=="BIC")  {
    bestmod <- which.min(BICs)
  }

  if (Type=="SingDst")  {   
    new("IdtSngNDE",ModelNames=modnames,ModelType=rep("Normal",nCovCases),ModelConfig=1:nCovCases,
      NIVar=q,mleNmuE=mleNmuE,mleNmuEse=mleNmuEse,CovConfCases=CovConfCases,SelCrit=SelCrit,
      logLiks=logLiks,AICs=AICs,BICs=BICs,BestModel=bestmod,SngD=TRUE)
  }  else  {
    new("IdtMxNDE",ModelNames=modnames,ModelType=rep("Normal",nCovCases),ModelConfig=1:nCovCases,
      NIVar=q,grouping=grouping,Hmcdt=TRUE,mleNmuE=mleNmuE,mleNmuEse=mleNmuEse,CovConfCases=CovConfCases,
      SelCrit=SelCrit,logLiks=logLiks,AICs=AICs,BICs=BICs,BestModel=bestmod,SngD=FALSE,Ngrps=k)
  }
}

IdtHetMxtNmle <- function( Idt,grouping, CVtol=1.0e-5, OptCntrl=OptCntrl,
	onerror=c("stop","warning","silentNull"), CovCaseArg, Config, SelCrit )
{
  onerror <- match.arg(onerror)
  n <- Idt@NObs
  q <- Idt@NIVar
  p <- 2*q
  nk <- as.numeric(table(grouping))
  k <- length(nk) 
  if (k==1)  { 
    error(onerror,"The data belongs to one single group. A partition into at least two different groups is required\n")
  }
  mleNmuE <- matrix(nrow=k,ncol=p)
  rownames(mleNmuE) <- levels(grouping)
  colnames(mleNmuE) <- c(names(Idt@MidP),names(Idt@LogR))
  mleNmuEse <- matrix(nrow=k,ncol=p)
  rownames(mleNmuEse) <- rownames(mleNmuE)
  colnames(mleNmuEse) <- colnames(mleNmuE)
  Xnams <- colnames(mleNmuE)
  anams <- list(Xnams,Xnams,levels(grouping))
  if (CovCaseArg)  {
    nCovCases <- 4 
    CovCaseMap <- c(1,NA,2,3,4)
    modnames <- paste("NModCovC",1:4,sep="")
  }  else  {
    nCovCases <- 5 
    CovCaseMap <- 1:5
    modnames <- paste("NC",1:5,sep="")
  }
  logLiks <- AICs <- BICs <- rep(NA_real_,nCovCases)
  CovConfCases <- vector("list",nCovCases)
  names(logLiks) <- names(AICs) <- names(BICs) <- names(CovConfCases) <- modnames

  for (model in Config)
  { 
    nparbyg <- npar(model,p,q)
    if (model==2)
    {
      CovConfCases[[2]] <- list( mleSigE=array(dim=c(p,p,k),dimnames=anams),
        mleSigEse=array(dim=c(p,p,k),dimnames=anams),
        mlevcov=array(dim=c(nparbyg,nparbyg,k),dimnames=list(NULL,NULL,levels(grouping))),
        logLik=0.,AIC=NULL,BIC=NULL,optres=vector("list",k) )
    }  else  {
      CovConfCases[[CovCaseMap[model]]] <- list( mleSigE=array(dim=c(p,p,k),dimnames=anams),
        mleSigEse=array(dim=c(p,p,k),dimnames=anams),
        mlevcov=array(dim=c(nparbyg,nparbyg,k),dimnames=list(NULL,NULL,levels(grouping))),
        logLik=0.,AIC=NULL,BIC=NULL )
    }
  }
  if (is.element(2,Config)) { names(CovConfCases[[2]]$optres) <- levels(grouping) }
  for (g in 1:k)
  {
    Idtg <- Idt[grouping==levels(grouping)[g],]
    IdtgDF <- cbind.data.frame(Idtg@MidP,Idtg@LogR)
    Xbar <- colMeans(IdtgDF)
    Xstdev <- sapply(IdtgDF,sd)
    CnstV <- which(Xstdev/abs(Xbar)<CVtol)
    if (length(CnstV)==1)  { 
      error( onerror,paste("Variable",names(CnstV),"appears to be constant in group",levels(grouping)[g],"\n") )
    }  else if (length(CnstV)>1)  {
      error( onerror,paste("Variables",paste(names(CnstV),collapse=" "),"appear to be constant in group ",levels(grouping)[g],"\n") )
    }
    pres <- IdtNmle(Idtg,grouping,Type="SingDst",CVtol=CVtol,OptCntrl=OptCntrl,CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit)
    mleNmuE[g,] <- pres@mleNmuE 
    mleNmuEse[g,] <- pres@mleNmuEse
    for (model in Config)
    { 
      CvCase <- CovCaseMap[model]	
      CovConfCases[[CvCase]]$mleSigE[,,g] <- pres@CovConfCases[[CvCase]]$mleSigE
      CovConfCases[[CvCase]]$mleSigEse[,,g] <- pres@CovConfCases[[CvCase]]$mleSigEse
      CovConfCases[[CvCase]]$mlevcov[,,g] <- pres@CovConfCases[[CvCase]]$mlevcov
      if (g==1) {
        dimnames(CovConfCases[[CvCase]]$mlevcov)[[1]]  <- dimnames(CovConfCases[[CvCase]]$mlevcov)[[2]] <- 
          rownames(pres@CovConfCases[[CvCase]]$mlevcov)
      }
      CovConfCases[[CvCase]]$logLik <- CovConfCases[[CvCase]]$logLik + pres@CovConfCases[[CvCase]]$logLik
      if (model==2)  { CovConfCases[[2]]$optres[[g]] <- pres@CovConfCases[[2]]$optres }
    }            
  }
  for (model in Config)
  {
    CvCase <- CovCaseMap[model]	
    logLiks[CvCase] <- CovConfCases[[CvCase]]$logLik 
    nmodelfreepar <- npar(model,p,q,Ngrps=k,Mxt="Het")
    AICs[CvCase] <- CovConfCases[[CvCase]]$AIC <- -2*CovConfCases[[CvCase]]$logLik + 2 * nmodelfreepar
    BICs[CvCase] <- CovConfCases[[CvCase]]$BIC <- -2*CovConfCases[[CvCase]]$logLik + log(n) * nmodelfreepar
  }
  if (SelCrit=="AIC")  {
    bestmod = which.min(AICs)
  }  else if (SelCrit=="BIC")  {
    bestmod = which.min(BICs)
  }
 
  new("IdtMxNDE",ModelNames=modnames,ModelType=rep("Normal",nCovCases),ModelConfig=1:nCovCases,grouping=grouping,
    Hmcdt=FALSE,mleNmuE=mleNmuE,mleNmuEse=mleNmuEse,CovConfCases=CovConfCases,SelCrit=SelCrit,NIVar=q,
    logLiks=logLiks,AICs=AICs,BICs=BICs,BestModel=bestmod,SngD=FALSE,Ngrps=k)
}

mleNSigC35 <- function(Conf,mat,q,p,defval)
{
  Newmat <- matrix(defval,p,p)
  dimnames(Newmat) <- dimnames(mat)
  if (Conf==3)  {
    for (i in 1:q) Newmat[c(i,q+i),c(i,q+i)] <- mat[c(i,q+i),c(i,q+i)]
  }  
  if (Conf==4)
  {
    Newmat[1:q,1:q] <- mat[1:q,1:q]
    Newmat[(q+1):(2*q),(q+1):(2*q)] <- mat[(q+1):(2*q),(q+1):(2*q)]
  }
  if (Conf==5) { diag(Newmat) <- diag(mat) }
  Newmat  #  return(Newmat)
}

mleNvcovC35 <- function(Conf,mat,p,k)
{
  kp <- k*p
  mind <- c(1:kp,kp+vcovCind(Conf,p))
  mat[mind,mind]  #  return(mat[mind,mind])
}

ILogLikNC <- function(Xdev,SigmaInv,const)   const -0.5 * Xdev%*%SigmaInv%*%Xdev
ILogLikNC1 <- function(Xdev,SigmaSrInv,const)   const -0.5 * sum((SigmaSrInv%*%Xdev)^2)
ILogLikDNC <- function(Xdev,IVar,const)  const -0.5 * sum(Xdev^2 * IVar)

SigCind <- function(Conf,q)
{
  if (!is.element(Conf,1:5)) {
    stop("Wrong value for Conf element\n")
  }
  if (Conf==1) { return( 1:(q*(2*q+1)) ) }
  mind <- NULL
  cind <- 0
  for (col in 1:q)  {
    nmidpv <- q-col+1  
    if (Conf==2) mind <- c(mind,cind+1:nmidpv,cind+nmidpv+col)
    if (Conf==3) mind <- c(mind,cind+1,cind+nmidpv+col)
    if (Conf==4) mind <- c(mind,cind+1:nmidpv)
    if (Conf==5) mind <- c(mind,cind+1)
    cind <- cind + nmidpv + q
  }
  for (col in 1:q)  {
    nlogrv <- q-col+1  
    if (Conf==3 || Conf==5) mind <- c(mind,cind+1)
    if (Conf==2 || Conf==4) mind <- c(mind,cind+1:nlogrv)
    cind <- cind + nlogrv
  }
  mind  #  return(mind)
}

npar <- function(Conf,p,q,Ngrps=1,Mxt=c("Hom","Het"))
{
  Mxt <- match.arg(Mxt)
  if (Mxt=="Hom")
  {
    if (Conf==1)  { return(Ngrps*p + p*(p+1)/2) }
    if (Conf==2)  { return(Ngrps*p + q*(q+1) + p) }
    if (Conf==3)  { return(Ngrps*p + p + q) }
    if (Conf==4)  { return(Ngrps*p + q*(q+1)) }
    if (Conf==5)  { return(Ngrps*p + p) }
  }  else if (Mxt=="Het")  {
    if (Conf==1)  { return(Ngrps*(p + p*(p+1)/2)) }
    if (Conf==2)  { return(Ngrps*(p + q*(q+1) + p)) }
    if (Conf==3)  { return(Ngrps*(p + p + q)) }
    if (Conf==4)  { return(Ngrps*(p + q*(q+1))) }
    if (Conf==5)  { return(Ngrps*(p + p)) }
  }
}
