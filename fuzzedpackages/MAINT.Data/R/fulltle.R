setMethod("fulltle",
  signature(Idt = "IData"),
  function(Idt, CovCase=1:4, SelCrit=c("BIC","AIC"), alpha=0.75, use.correction=TRUE, getalpha="TwoStep", 
    rawMD2Dist=c("ChiSq","HardRockeAsF","HardRockeAdjF"), MD2Dist=c("ChiSq","CerioliBetaF"),
    eta=0.025, multiCmpCor=c("never","always","iterstep"), outlin=c("MidPandLogR","MidP","LogR"), reweighted=TRUE, k2max=1e6,
    force=FALSE, ...)
  {

    if (!requireNamespace("robustbase",quietly=TRUE)) 
      stop("fulltle needs the robustbase package to work. Please install it.\n")

    SelCrit <- match.arg(SelCrit)
    rawMD2Dist <- match.arg(rawMD2Dist)
    MD2Dist <- match.arg(MD2Dist)
    multiCmpCor <- match.arg(multiCmpCor)
    outlin <- match.arg(outlin)
    limlnk2 <- log(k2max)
    q <- Idt@NIVar
    if (q==1) CovCase <- q1CovCase(CovCase) 
    n <- Idt@NObs    

    if (getalpha=="TwoStep")
    {
#      X <- cbind(Idt@MidP,Idt@LogR)
      X <- cbind.data.frame(Idt@MidP,Idt@LogR)
      if (MD2Dist=="ChiSq") {
        fstsol <- fulltle1(Idt,CovCase,SelCrit,alpha,use.correction,rawMD2Dist,eta,multiCmpCor,outlin,reweighted,limlnk2,force)
        nOtls <- MDOtlDet(X,coef(fstsol)$mu,coef(fstsol)$Sigma,eta=eta,RefDist="ChiSq",multiCmpCor=multiCmpCor,otp="onlycnt")
      }  else if (MD2Dist=="CerioliBetaF")  {
        fstsol <- fulltle1(Idt,CovCase,SelCrit,alpha,use.correction,rawMD2Dist,eta,multiCmpCor,outlin,reweighted,limlnk2,force)
        nOtls <- MDOtlDet(X,coef(fstsol)$mu,coef(fstsol)$Sigma,eta=eta,RefDist="CerioliBetaF",Rewind=fstsol@RewghtdSet,multiCmpCor=multiCmpCor,otp="onlycnt")
      }
      if (nOtls==0) {
        Config <- getConfig(...)
        if (is.null(Config))  
        {
          Config <- ifelse(CovCase==1,1,CovCase+1)
          CovCaseArg <- TRUE	
        } else {  
          CovCaseArg <- FALSE
        }	
        finalsol <- IdtNmle(Idt,CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit)
        warning(paste("fulltle returned the classical maximum likelihood estimates because the data does not appear to include any outlier.\n",
          "If you want to force a trimmed likelihood estimator run fasttle with the argument getalpha=FALSE.\n"))
          rawSet <- RewghtdSet <- 1: Idt@NObs
        names(rawSet) <- names(RewghtdSet) <- Idt@ObsNames
        if (outlin=="MidPandLogR") {
          RobMD2 <- GetMD2(X,coef(finalsol)$mu,coef(finalsol)$Sigma)
        } else if (outlin=="MidP") {
          RobMD2 <- GetMD2(X[,1:q],coef(finalsol)$mu[1:q],coef(finalsol)$Sigma[1:q,1:q])
        } else if (outlin=="LogR") {
          RobMD2 <- GetMD2(X[,(q+1):(2*q)],coef(finalsol)$mu[(q+1):(2*q)],coef(finalsol)$Sigma[(q+1):(2*q),(q+1):(2*q)])
        } 
        for (case in 1:length(finalsol@CovConfCases)) {
          if (!is.null(finalsol@CovConfCases[[case]])) {
            names(finalsol@CovConfCases[[case]])[1] <- "RobSigE"
            finalsol@CovConfCases[[case]][2] <- finalsol@CovConfCases[[case]][3] <- NULL
          }
        }
        return(new("IdtSngNDRE",ModelNames=finalsol@ModelNames,ModelType=finalsol@ModelType,ModelConfig=finalsol@ModelConfig,
                   NIVar=finalsol@NIVar,SelCrit=finalsol@SelCrit,logLiks=finalsol@logLiks,BICs=finalsol@BICs,AICs=finalsol@AICs,
                   BestModel=finalsol@BestModel,RobNmuE=finalsol@mleNmuE,CovConfCases=finalsol@CovConfCases,SngD=TRUE,
                   rawSet=rawSet,RewghtdSet=RewghtdSet,RobMD2=RobMD2,
                   cnp2=c(1.,1.),raw.cov=coef(fstsol)$Sigma,raw.cnp2=c(1.,1.))
        )
      }  
      newalpha <- 1. - nOtls/n
      return( fulltle1(Idt,CovCase,SelCrit,newalpha,use.correction,rawMD2Dist,eta,multiCmpCor,outlin,reweighted,limlnk2,force) )
    }  else {
      return( fulltle1(Idt,CovCase,SelCrit,alpha,use.correction,rawMD2Dist,eta,multiCmpCor,outlin,reweighted,limlnk2,force) )
    } 
  }
)

fulltle1 <-  function(Idt,CovCase,SelCrit,alpha,use.correction,rawMD2Dist,eta,multiCmpCor,outlin,reweighted,limlnk2,force,...)
{
  Config <- getConfig(...)
  if (is.null(Config))  
  {
    Config <- ifelse(CovCase==1,1,CovCase+1)
    CovCaseArg <- TRUE	
    nCovCases <- 4 
    CovCaseMap <- c(1,NA,2,3,4)
    modnames <- paste("NModCovC",1:4,sep="")
  } else {  
    CovCaseArg <- FALSE
    nCovCases <- 5 
    CovCaseMap <- 1:5
    modnames <- paste("NC",1:5,sep="")
  }	

  n <- Idt@NObs    
  if (outlin=="MidPandLogR")  {
#    X <- as.matrix(cbind(Idt@MidP,Idt@LogR))
    X <- as.matrix(cbind.data.frame(Idt@MidP,Idt@LogR))
    p <- 2*Idt@NIVar
  }
  else {
    p <- Idt@NIVar
    if (outlin=="MidP") {
      X <- as.matrix(Idt@MidP)
      Vind <- 1:p
    }
    if (outlin=="LogR") {
      X <- as.matrix(Idt@LogR)
      Vind <- (p+1):(2*p)
    }
  }
  rownames(X) <- Idt@ObsNames
  k <- robustbase::h.alpha.n(alpha,n,p)
  if (k >=n || k <= p) stop("fulltle stopped because the choice of the trimming parameter implies trimmed samples so small\n", 
   "that the resulting covariance matrix estimates would be singular.\n")
  if (!force) {
    maxnCk <- 10000000
    nCk <- choose(n,k) 
    if (nCk> maxnCk) stop(paste("fulltle might take too long since",nCk,"different subsets",
                      "need to be evaluated.\nTo proceed anyway set the 'force'",
                      "argument to TRUE, otherwise try the fasttle method instead.\n"))
  }                    

  c0 <- -n*(p*(log(2*pi)+1))/2
  bestCrt <- Inf
  if (SelCrit=="BIC") penC <- log(n)
  else if (SelCrit=="AIC") penC <- 2
  if (outlin=="MidPandLogR")
  {
    for (Cnf in Config) {
      if (Cnf==2) {
        Cftmpsol <- Rfulltle(Idt,k,2,SelCrit,TRUE,...)
      }  else {
        Cftmpsol <- .Call("Cfulltle", X, n, p, k, Cnf, limlnk2, c0, PACKAGE = "MAINT.Data" )
      }
     CmpCrt <- -2*Cftmpsol$LogLik + penC*npar(Cnf,p,p/2)
     Cftmpsol$Set <- Cftmpsol$Set + 1 # Conversion of C 0-ind convention to R 1-ind convention
      if (CmpCrt < bestCrt) {
        bestCrt <- CmpCrt
        bestSet <- Cftmpsol$Set
      }
    }
  } else {  
    Cftmpsol <- .Call("Cfulltle", X, n, p, k, 1, limlnk2, c0, PACKAGE = "MAINT.Data" )
    CmpCrt <- -2*Cftmpsol$LogLik + penC*npar(1,2*p,p)
    Cftmpsol$Set <- Cftmpsol$Set + 1 # Conversion of C 0-ind convention to R 1-ind convention
    if (CmpCrt < bestCrt) {
      bestCrt <- CmpCrt
      bestSet <- Cftmpsol$Set
    }
  }
  bestsol <- IdtNmle(Idt[bestSet,],CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit)
  for (Cnf in Config) bestsol@CovConfCases[[CovCaseMap[Cnf]]]$mleSigEse <- bestsol@CovConfCases[[CovCaseMap[Cnf]]]$mlevcov <- NULL   
  
  dhn1 <- robustbase::.MCDcons(p,k/n)
  if (use.correction) {
    dhn <- dhn1*MCDcnp2(p,n,alpha)
  } else {
    dhn <- rep(dhn1,4)
  }

  rawcov <- bestsol@CovConfCases[[bestsol@BestModel]]$mleSigE
  if (outlin=="MidPandLogR")
  {
    Xdev <- scale(X,center=bestsol@mleNmuE,scale=FALSE)
    Sigma <- rawcov <- dhn[bestsol@BestModel] * rawcov       
  } else { 
    Xdev <- scale(X,center=bestsol@mleNmuE[Vind],scale=FALSE)
    Sigma <- rawcov[Vind,Vind] <- dhn[bestsol@BestModel] * rawcov[Vind,Vind]
  }
  SigmaI <- pdwt.solve(Sigma)
  MD2 <- apply(Xdev,1,function(x) x%*%SigmaI%*%x)

  if (!reweighted)  {
    if (outlin=="MidPandLogR") {
      for (Cnf in Config) {
        CvCase <- CovCaseMap[Cnf]	
        bestsol@CovConfCases[[CvCase]]$RobSigE <- dhn[CvCase] * bestsol@CovConfCases[[CvCase]]$mleSigE
        bestsol@CovConfCases[[CvCase]]$mleSigE <- NULL
      }
    } else {
      for (Cnf in Config) {
        CvCase <- CovCaseMap[Cnf]	
        bestsol@CovConfCases[[CvCase]]$RobSigE <- bestsol@CovConfCases[[CvCase]]$mleSigE
        bestsol@CovConfCases[[CvCase]]$RobSigE[Vind,Vind] <- dhn[CvCase] * bestsol@CovConfCases[[CvCase]]$mleSigE[Vind,Vind]
        bestsol@CovConfCases[[CvCase]]$mleSigE <- NULL
      }
    }  
    RewghtdSet <- bestSet
    raw.cnp2 <- cnp2 <- c(dhn1,dhn[CvCase]/dhn1)
    names(cnp2) <- names(raw.cnp2) <- NULL 
  } else {
    for (Cnf in Config) bestsol@CovConfCases[[CovCaseMap[Cnf]]]$mleSigE <- NULL
    oneminuseta <- 1-eta
    if (multiCmpCor=="always" || multiCmpCor=="iterstep") {
      oneminusalpha <- oneminuseta^(1/n)
    } else {
      oneminusalpha <- oneminuseta
    }
    h <- length(bestSet)
    findotl <- TRUE
    iter <- 1
    while (findotl && iter<3)  { 
      if (rawMD2Dist=="ChiSq")  { 
        MD2trshld <- qchisq(oneminusalpha,p)
      } else if (rawMD2Dist=="HardRockeAdjF")  {
        MD2trshld <- qHardRoqF(oneminusalpha,n,p,h)
      } else if (rawMD2Dist=="HardRockeAsF")  {
        MD2trshld <- qHardRoqF(oneminusalpha,n,p,h,adj=FALSE)
      }
      RewghtdSet <- which(MD2<=MD2trshld)
      if (multiCmpCor=="iterstep")
      {
        if (length(RewghtdSet)==n) {
          findotl <- FALSE
        } else {
          oneminusalpha <- oneminuseta
          iter <- iter+1
        }
      } else {
        findotl <- FALSE
      }
    }
    bestsol <- IdtNmle(Idt[RewghtdSet,],CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit)
    for (Cnf in Config) {
      CvCase <- CovCaseMap[Cnf]	
      names(bestsol@CovConfCases[[CvCase]])[which(names(bestsol@CovConfCases[[CvCase]])=="mleSigE")] <- "RobSigE" 
      bestsol@CovConfCases[[CvCase]]$mleSigEse <- bestsol@CovConfCases[[CvCase]]$mlevcov <- NULL   
    }

    m <- length(RewghtdSet)
    rdmhn1 <- robustbase::.MCDcons(p,m/n)
    if (outlin=="MidPandLogR")  {
      if (use.correction) {
        rdmhn <- rdmhn1 * MCDRewcnp2(p,n,m/n,alpha)
      } else { 
        rdmhn <- rep(rdmhn1,4)
      }
      for (Cnf in Config) {
        CvCase <- CovCaseMap[Cnf]	
        bestsol@CovConfCases[[CvCase]]$RobSigE <- rdmhn[CvCase] * bestsol@CovConfCases[[CvCase]]$RobSigE
      }
    } else {
      if (use.correction) {
        rdmhn <- rdmhn1 * MCDRewcnp2(Idt@NIVar,n,m/n,alpha)
      } else { 
        rdmhn <- rep(rdmhn1,4)
      }
      for (Cnf in Config) {
        CvCase <- CovCaseMap[Cnf]	
        bestsol@CovConfCases[[CvCase]]$RobSigE[Vind,Vind] <- rdmhn[CvCase] * bestsol@CovConfCases[[CvCase]]$RobSigE[Vind,Vind]
      }
    }
    cnp2 <- c(rdmhn1,rdmhn[bestsol@BestModel]/rdmhn1)
    raw.cnp2 <- c(dhn1,dhn[bestsol@BestModel]/dhn1)
    names(cnp2) <- names(raw.cnp2) <- NULL 
  }

  rawSet <- sort(bestSet)

  finalsol <- new("IdtSngNDRE",ModelNames=bestsol@ModelNames,ModelType=bestsol@ModelType,ModelConfig=bestsol@ModelConfig,
    NIVar=bestsol@NIVar,SelCrit=bestsol@SelCrit,logLiks=bestsol@logLiks,BICs=bestsol@BICs,AICs=bestsol@AICs,
    BestModel=bestsol@BestModel,RobNmuE=bestsol@mleNmuE,CovConfCases=bestsol@CovConfCases,SngD=TRUE,
    rawSet=rawSet,RewghtdSet=RewghtdSet,RobMD2=MD2,cnp2=cnp2,raw.cov=rawcov,raw.cnp2=raw.cnp2,PerfSt=list()
  )

  finalsol  # return(finalsol)
}

