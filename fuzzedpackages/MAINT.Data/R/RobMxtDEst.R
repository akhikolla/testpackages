setMethod("RobMxtDEst",
  signature(Idt = "IData"),
  function(Idt,grouping,Mxt=c("Hom","Het"),CovEstMet=c("Pooled","Globdev"),
    CovCase=1:4,SelCrit=c("BIC","AIC"),Robcontrol=RobEstControl(), l1medpar=NULL, ...)
  {
    if (!requireNamespace("robustbase",quietly=TRUE)) 
      stop("fasttle needs the robustbase package to work. Please install it\n")

    Mxt <- match.arg(Mxt)
    CovEstMet <- match.arg(CovEstMet)
    SelCrit <- match.arg(SelCrit)

    getalpha <- Robcontrol@getalpha
    alpha <- Robcontrol@alpha
    reweighted <- Robcontrol@reweighted
    getkdblstar <- Robcontrol@getkdblstar
    nsamp <- Robcontrol@nsamp
    ncsteps <- Robcontrol@ncsteps
    trace <- Robcontrol@trace
    use.correction <- Robcontrol@use.correction
    outlin <- Robcontrol@outlin
    trialmethod <- Robcontrol@trialmethod
    m <- Robcontrol@m
    otpType <- Robcontrol@otpType

    n <- Idt@NObs  
    p <- 2*Idt@NIVar  
    if (p==2) CovCase <- q1CovCase(CovCase) 
    if (n <=p) {
      stop("The number of observations is too small (not larger than the total number of variables) and would lead to singular covariance estimates.\n")
    }
    q <- p/2
    if (length(grouping)!=n)  {
      stop("The size of the grouping factor does not agree with the number of observations in the data set supplied")
    }

    grouping <- factor(grouping,exclude=NULL)
    grplvls <- levels(grouping)
    ng <- as.integer(table(grouping))
    if (Mxt=="Het" && any(ng <=p)) {
      stop("The number of observations is too small (for some group(s) not larger than the total number of variables) and would lead to singular covariance estimates.\n")
    }
    ngrps <- length(ng)
    grpRobE <- vector("list",ngrps)
#    Xnams <- names(cbind(Idt@MidP,Idt@LogR))
    Xnams <- c(names(Idt@MidP),names(Idt@LogR))
    anams <- list(Xnams,Xnams,grplvls)
    RobNmuE <- matrix(nrow=ngrps,ncol=p,dimnames=list(grplvls,Xnams))

    nCovC <- length(CovCase)
    CovConfC <- vector("list",nCovC)
    logLiks <- numeric(nCovC) 
    AICs <- numeric(nCovC) 
    BICs <- numeric(nCovC) 
    names(logLiks) <- names(AICs) <- names(BICs) <- names(CovConfC) <- modnames <- paste("NModCovC",CovCase,sep="")

    rawSet <- NULL
    RewghtdSet <- NULL
    RobMD2 <- NULL
    cnp2 <- matrix(nrow=ngrps,ncol=2,dimnames=list(grplvls,NULL))
    raw.cnp2 <- matrix(nrow=ngrps,ncol=2,dimnames=list(grplvls,NULL))
    raw.cov <- array(dim=c(p,p,ngrps),dimnames=anams)
    if (otpType!="SetMD2EstandPrfSt") {
      PerfSt <- list()
    } else {
      PerfSt <- list(RepSteps=vector("list",nCovC),RepLogLik=vector("list",nCovC),StpLogLik=vector("list",nCovC))
      names(PerfSt$RepSteps) <- names(PerfSt$RepLogLik) <- names(PerfSt$StpLogLik) <- modnames 
      Repnames <- paste("Rep",1:nsamp,sep="")
      Stepnames <- paste("Stp",1:ncsteps,sep="")
      for (CovC in 1:nCovC) {
        PerfSt$RepSteps[[CovC]] <- matrix(nrow=nsamp,ncol=ngrps,dimnames=list(Repnames,grplvls))
        PerfSt$RepLogLik[[CovC]] <- matrix(nrow=nsamp,ncol=ngrps,dimnames=list(Repnames,grplvls))
        PerfSt$StpLogLik[[CovC]] <- array(dim=c(nsamp,ncsteps,ngrps),dimnames=list(Repnames,Stepnames,grplvls))
      }  
    }

    if (Mxt=="Hom" && CovEstMet=="Globdev")
    {
#      X <- cbind(Idt@MidP,Idt@LogR)
      X <- cbind.data.frame(Idt@MidP,Idt@LogR)
      Xdev <- matrix(nrow=n,ncol=p) 
      Xgl1med <- matrix(nrow=ngrps,ncol=p) 
      if (!is.null(l1medpar)) {
        MaxStep <- ifelse(is.null(l1medpar$MaxStep),200,l1medpar$MaxStep)
        ItTol <- ifelse(is.null(l1medpar$ItTol),10^-8,l1medpar$ItTol)
        trace <- ifelse(is.null(l1medpar$trace),0,l1medpar$trace)
      }
      for (g in 1:ngrps)
      {
        gind <- which(grouping==grplvls[g])
        Xg <- X[gind,]
        if (is.null(l1medpar)) {
          Xgl1med[g,] <- l1median(Xg)
        }  else {
          if (is.null(l1medpar$m.init)) {
            m.init <- robustbase::colMedians(as.matrix(Xg))
          } else {
            m.init <- l1medpar$m.init
          }
          Xgl1med[g,] <- l1median(Xg,MaxStep,ItTol,trace,m.init)
        }
        Xdev[gind,] <- scale(Xg,center=Xgl1med[g,],scale=FALSE)
      }
    }
    nSteps <- ifelse(getalpha=="TwoStep",2,1)
    for (Steps in 1:nSteps)
    {
      for (CovC in 1:nCovC)
      {
#       maxnSteps <- 0 
        if (Mxt=="Hom")
        {
          CovConfC[[CovC]] <- list(
            RobSigE=matrix(0.,nrow=p,ncol=p,dimnames=list(Xnams,Xnams)),logLik=NULL,AIC=NULL,BIC=NULL
          )
          trmdn <- round(alpha*n)
          if (CovEstMet=="Pooled")  {
            X <- NULL
            for (g in 1:ngrps)  {      # Things to do: add a fulltle option !!!
              Idtg <- Idt[grouping==grplvls[g],]
              grpRobE[[g]] <- fasttle(Idtg,CovC,SelCrit,
                 getalpha="NO",control=Robcontrol, ...)
              RobNmuE[g,] <- grpRobE[[g]]@RobNmuE
              CovConfC[[CovC]]$RobSigE <- 
                CovConfC[[CovC]]$RobSigE + (ng[g]/n) * grpRobE[[g]]@CovConfCases[[CovC]]$RobSigE 
            }
            regset <- ifelse(reweighted,grpRobE[[g]]@RewghtdSet,grpRobE[[g]]@rawSet)
#            X <- rbind(X,cbind(Idtg@MidP[regset,],Idtg@LogR[regset,]))
            X <- rbind.data.frame(X,cbind.data.frame(Idtg@MidP[regset,],Idtg@LogR[regset,]))
          }
          else if (CovEstMet=="Globdev") {
            if  (getkdblstar=="Twopplusone") { 
              kdblstar <- 2*Idt@NIVar+1
            }  else {
              if (!is.finite(getkdblstar)) {
                stop("Wrong value for Robcontrol parameter getkdblstar\n")
              }
              kdblstar <- getkdblstar 
            }
            RobE <- fasttle1(Xdev,CovC,SelCrit,alpha,nsamp,ncsteps,trace,use.correction,kdblstar,outlin,trialmethod,m,
              reweighted,otpType,Idt@VarNames,...)
            for (g in 1:ngrps)  RobNmuE[g,] <- Xgl1med[g,] + RobE@RobNmuE
            CovConfC[[CovC]]$RobSigE <- RobE@CovConfCases[[CovC]]$RobSigE 
          }
          Xdev <- scale(X[grouping==levels(grouping)[1],],center=RobNmuE[1,],scale=FALSE)
          for (g in 2:ngrps)
#            Xdev <- rbind(Xdev,scale(X[grouping==levels(grouping)[g],],center=RobNmuE[g,],scale=FALSE))
            Xdev <- rbind.data.frame(Xdev,scale(X[grouping==levels(grouping)[g],],center=RobNmuE[g,],scale=FALSE))
          logdet <- pdwt.solve(CovConfC[[CovC]]$RobSigE,silent=TRUE,onlylogdet=TRUE)
          if (is.null(logdet))  {
            logLiks[CovC] <- CovConfC[[CovC]]$logLik <- -Inf
          }  else  {
            logLiks[CovC] <- CovConfC[[CovC]]$logLik <- -trmdn*(p*(log(2*pi)+1)+logdet)/2
          }		
        }  else if (Mxt=="Het") {
          CovConfC[[CovC]] <- list( RobSigE=array(dim=c(p,p,ngrps),dimnames=anams),logLik=NULL,AIC=NULL,BIC=NULL )
          for (g in 1:ngrps) {
            grpRobE[[g]] <- fasttle(Idt[grouping==grplvls[g],],CovC,SelCrit,
              getalpha="NO",control=Robcontrol, ...)
            RobNmuE[g,] <- grpRobE[[g]]@RobNmuE
            CovConfC[[CovC]]$RobSigE[,,g] <- grpRobE[[g]]@CovConfCases[[CovC]]$RobSigE
          } 
          logLiks[CovC] <- grpRobE[[1]]@logLiks[CovC]
          for (g in 2:ngrps) logLiks[CovC] <- logLiks[CovC] + grpRobE[[g]]@logLiks[CovC]
        }
        CovConfC[[CovC]]$logLik <- logLiks[CovC]
        nmodelfreepar <- npar(CovC,p,q,Ngrps=ngrps,Mxt=Mxt)
        CovConfC[[CovC]]$AIC <- AICs[CovC] <- -2*logLiks[CovC] + 2*nmodelfreepar
        trmdn <- round(alpha*n)
        CovConfC[[CovC]]$BIC <- BICs[CovC] <- -2*logLiks[CovC] + log(trmdn)*nmodelfreepar

        if (Steps==nSteps) for (g in 1:ngrps)
        {
          rawSet <- c(rawSet,grpRobE[[g]]@rawSet)
          RewghtdSet <- c(RewghtdSet,grpRobE[[g]]@RewghtdSet)
          RobMD2 <- c(RobMD2,grpRobE[[g]]@RobMD2)
          cnp2[g,] <- grpRobE[[g]]@cnp2
          raw.cnp2[g,] <- grpRobE[[g]]@raw.cnp2
          raw.cov[,,g] <- grpRobE[[g]]@raw.cov
          if (otpType=="SetMD2EstandPrfSt") {
            PerfSt$RepSteps[[CovC]][,g] <- grpRobE[[g]]@PerfSt$RepSteps[[CovC]]
            PerfSt$RepLogLik[[CovC]][,g] <- grpRobE[[g]]@PerfSt$RepLogLik[[CovC]]
            maxnSteps <- ncol(grpRobE[[g]]@PerfSt$StpLogLik[[CovC]])
            PerfSt$StpLogLik[[CovC]][,1:maxnSteps,g] <- grpRobE[[g]]@PerfSt$StpLogLik[[CovC]]
#            nSteps <- ncol(grpRobE[[g]]@PerfSt$StpLogLik[[CovC]])
#            tmpStpLogLik[,1:nSteps,g] <- grpRobE[[g]]@PerfSt$StpLogLik[[CovC]]
#            maxnSteps <- max(maxnSteps,nSteps)
          }
        }
#        if (Steps==nSteps && otpType=="SetMD2EstandPrfSt") {
#          PerfSt$StpLogLik[[CovC]] <- array(tmpStpLogLik[,1:maxnSteps,],dim=c(nsamp,maxnSteps,ngrps),dimnames=list(Repnames,Stepnames[1:maxnSteps],grplvls)
#        }
      }
      if (SelCrit=="AIC")  {
        bestmod = which.min(AICs)
      }  else if (SelCrit=="BIC")  {
        bestmod = which.min(BICs)
      }

      if(getalpha=="TwoStep" && Steps==1)
      {
#        X <- data.frame(cbind(Idt@MidP,Idt@LogR))
        X <- data.frame(cbind.data.frame(Idt@MidP,Idt@LogR))
        nOtls <- 0.
        for (g in 1:ngrps)
        {
          if (Mxt=="Hom")  {
            nOtls <- nOtls + MDOtlDet(X[grouping==grplvls[g],],RobNmuE[g,],CovConfC[[bestmod]]$RobSigE,0.025,otp="onlycnt") 
          } else if (Mxt=="Het") {
            nOtls <- nOtls + MDOtlDet(X[grouping==grplvls[g],],RobNmuE[g,],CovConfC[[bestmod]]$RobSigE[,,g],0.025,otp="onlycnt") 
          }
        }
        alpha <- 1.-nOtls/n
      }

    }

    new( "IdtMxNDRE", ModelNames=modnames,ModelType=rep("Normal",nCovC),ModelConfig=1:nCovC,
      grouping=grouping,Hmcdt=(Mxt=="Hom"),RobNmuE=RobNmuE,CovConfCases=CovConfC,
      SelCrit=SelCrit,NIVar=q,logLiks=logLiks,AICs=AICs,BICs=BICs,BestModel=bestmod,SngD=FALSE,Ngrps=ngrps,
      rawSet=rawSet,RewghtdSet=RewghtdSet,RobMD2=RobMD2,cnp2=cnp2,raw.cov=raw.cov,raw.cnp2=raw.cnp2,PerfSt=PerfSt )
  }
)





