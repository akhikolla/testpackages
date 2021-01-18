getConfig <- function(Config=NULL,...) Config 

setMethod("mle",
  signature(Idt = "IData"),
  function(Idt,Model=c("Normal","SKNormal","NrmandSKN"),CovCase=1:4,
  SelCrit=c("BIC","AIC"), k2max=1e6, OptCntrl=list(),...)
  {
    Model <- match.arg(Model)
    SelCrit <- match.arg(SelCrit)
    q <- Idt@NIVar
    p <- 2*q
    n <- Idt@NObs
    limlnk2 <- log(k2max)

    Config <- getConfig(...)
    if (is.null(Config))  
    {
      Config <- ifelse(CovCase==1,1,CovCase+1)
      CovCaseArg <- TRUE
      CovCaseMap <- c(1,NA,2,3,4)
    } else {
      CovCaseArg <- FALSE
      CovCaseMap <- 1:5
    }	
    if (Idt@NIVar==1) {
      CovCase <- q1CovCase(CovCase) 
      Config <- q1Config(Config)
    }  
    if (Model!="SKNormal") { 
      Nres <- IdtNmle(Idt,OptCntrl=OptCntrl,CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit)
    }
    if (Model!="Normal") { 
#      SNres <- IdtSNmle(Idt,OptCntrl=OptCntrl,limlnk2=limlnk2,CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit)
      SNres <- IdtSNmle(Idt,OptCntrl=OptCntrl,limlnk2=limlnk2,CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit,...)
    }
    if (Model=="Normal") { return(Nres) }
    if (Model=="SKNormal") { return(SNres) }
    if (Model=="NrmandSKN")
    {

      Xscld <- NULL 
      for (Conf in Config) {
        CvCase <- CovCaseMap[Conf]
        if (SNres@logLiks[CvCase] < Nres@logLiks[CvCase]) {
          if (is.null(Xscld)) {
            Xscld <- scale(cbind.data.frame(Idt@MidP,Idt@LogR))
            Xmean <- attr(Xscld,"scaled:center")  
            Xsd <- attr(Xscld,"scaled:scale")  
            XsdOutPrd <- outer(Xsd,Xsd)
            lglikdif <- -n*sum(log(Xsd))		# difference between the log-likelihoods for the original and standardized data  
          } 
#          SNres@CovConfCases[[CvCase]] <- CorrSNSol(Nres,SNres,CvCase,Conf,Xscld,Xmean,Xsd,XsdOutPrd,lglikdif,limlnk2,OptCntrl)
          SNres@CovConfCases[[CvCase]] <- CorrSNSol(Nres,SNres,CvCase,Conf,Xscld,Xmean,Xsd,XsdOutPrd,lglikdif,limlnk2,OptCntrl,...)
          SNres@logLiks[CvCase] <- LogLik <- SNres@CovConfCases[[CvCase]]$logLik
          nfreepar <- SKnpar(Conf,p,q)
          SNres@AICs[CvCase] <- -2*LogLik + 2*nfreepar
          SNres@BICs[CvCase] <- -2*LogLik + log(n)*nfreepar
        }
      }

      AICs <- c(Nres@AICs,SNres@AICs)
      BICs <- c(Nres@BICs,SNres@BICs)    			
      if (SelCrit=="AIC") {
        bestmod <- which.min(AICs)
      } else { 
        if (SelCrit=="BIC") { bestmod <- which.min(BICs) }
      }
      return( new("IdtSngNandSNDE",NMod=Nres,SNMod=SNres,ModelNames=c(Nres@ModelNames,SNres@ModelNames),
        ModelType=c(Nres@ModelType,SNres@ModelType),ModelConfig=c(Nres@ModelConfig,SNres@ModelConfig),
        SelCrit=SelCrit,NIVar=Idt@NIVar,logLiks=c(Nres@logLiks,SNres@logLiks),AICs=AICs,BICs=BICs,
        BestModel=bestmod,SngD=TRUE)  
      )
    }
  }
)

setMethod("MANOVA",
  signature(Idt = "IData"),
  function(Idt, grouping, Model=c("Normal","SKNormal","NrmandSKN"), CovCase=1:4,
    SelCrit=c("BIC","AIC"), Mxt=c("Hom","Het","Loc","Gen"), CVtol=1.0e-5, k2max=1e6, 
    OptCntrl=list(), onerror=c("stop","warning","silentNull"), ...)
  {
    onerror <- match.arg(onerror)
    limlnk2 <- log(k2max)
    if (!is.factor(grouping)) { return(error(onerror,"'grouping' is not a factor\n")) }
    if ( Idt@NObs != length(grouping)) { 
      return(error(onerror,"The numbers of data and partition observations are different\n"))
    }
    Mxt <- match.arg(Mxt)
    Model <- match.arg(Model)
    SelCrit <- match.arg(SelCrit)

#    H0res <- mle(Idt,Model=Model,CovCase=CovCase,SelCrit=SelCrit,OptCntrl=OptCntrl)
    H0res <- mle(Idt,Model=Model,CovCase=CovCase,SelCrit=SelCrit,OptCntrl=OptCntrl,getvcov=FALSE,...)

    Config <- getConfig(...)
    if (is.null(Config))  
    {
      Config <- ifelse(CovCase==1,1,CovCase+1)
      CovCaseArg <- TRUE	
      CovCaseMap <- c(1,NA,2,3,4)
    } else {  
      CovCaseArg <- FALSE
      CovCaseMap <- 1:5
    }	
    if (Idt@NIVar==1) {
      CovCase <- q1CovCase(CovCase) 
      Config <- q1Config(Config)
    }  

    grouping <- factor(grouping,exclude=NULL)
    nk <- as.numeric(table(grouping))
    q <- Idt@NIVar 
    p <- 2*q
    n <- Idt@NObs 
    k <- length(nk) 
    if (k==1) {
      return( error(onerror,
        "The data belongs to one single group. A partition into at least two different groups is required\n")
      )
    }
    if (length(H0res@BestModel)==0) {
      return( error(onerror,"Procedure MANOVA failed to find a valid null model mle estimate\n") )
    }

    if (Model=="Normal" && (Mxt=="Hom" || Mxt=="Loc") ) {  
      H1res <- IdtNmle(Idt,grouping,Type="HomMxt",CVtol=CVtol,
        OptCntrl=OptCntrl,CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit)
    }  else if (Model=="SKNormal" && (Mxt=="Hom" || Mxt=="Loc") )  {  
      H1res <- IdtSNmle(Idt,grouping,Type="HomMxt",CVtol=CVtol,limlnk2=limlnk2,
#        OptCntrl=OptCntrl,CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit)
        OptCntrl=OptCntrl,CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit,getvcov=FALSE,...)
    }  else if (Model=="Normal" && (Mxt=="Het" || Mxt=="Gen") ) { 
      H1res <- IdtHetMxtNmle(Idt,grouping,CVtol=CVtol,
        OptCntrl=OptCntrl,CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit)
    }  else if (Model=="SKNormal" && (Mxt=="Het" || Mxt=="Gen") ) { 
      H1res <- IdtFDMxtSNmle(Idt,grouping,CVtol=CVtol,limlnk2=limlnk2,
#        OptCntrl=OptCntrl,CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit)
        OptCntrl=OptCntrl,CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit,getvcov=FALSE,...)

    }  else if (Model=="NrmandSKN")  {
      if (Mxt=="Hom" || Mxt=="Loc")  {
        SNH1res <- IdtSNmle(Idt,grouping,Type="HomMxt",CVtol=CVtol,limlnk2=limlnk2,
#          OptCntrl=OptCntrl,CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit)
          OptCntrl=OptCntrl,CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit,getvcov=FALSE,...)
        NH1res <- IdtNmle(Idt,grouping,Type="HomMxt",CVtol=CVtol,
          OptCntrl=OptCntrl,CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit)
        SNH1res <- IdtSNmle(Idt,grouping,Type="HomMxt",CVtol=CVtol,limlnk2=limlnk2,
#          OptCntrl=OptCntrl,CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit)
          OptCntrl=OptCntrl,CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit,getvcov=FALSE,...)
      } else if (Mxt=="Het" || Mxt=="Gen")  {
        NH1res <- IdtHetMxtNmle(Idt,grouping,CVtol=CVtol,
          OptCntrl=OptCntrl,CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit)
        SNH1res <- IdtFDMxtSNmle(Idt,grouping,CVtol=CVtol,limlnk2=limlnk2,
#          OptCntrl=OptCntrl,CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit)
          OptCntrl=OptCntrl,CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit,getvcov=FALSE,...)
      } 
      Xscld <- NULL 
      for (Conf in Config) {
        CvCase <- CovCaseMap[Conf]
        if (SNH1res@logLiks[CvCase] < NH1res@logLiks[CvCase]) {
          if (is.null(Xscld)) {
            Xscld <- scale(cbind.data.frame(Idt@MidP,Idt@LogR))
            Xmean <- attr(Xscld,"scaled:center")  
            Xsd <- attr(Xscld,"scaled:scale")  
            XsdOutPrd <- outer(Xsd,Xsd)
            lglikdif <- -n*sum(log(Xsd))		# difference between the log-likelihoods for the original and standardized data  
          } 
#          SNH1res@CovConfCases[[CvCase]] <- CorrHetSNSol(NH1res,SNH1res,CvCase,Conf,Xscld,Xmean,Xsd,XsdOutPrd,grouping,Mxt,lglikdif,limlnk2,OptCntrl)
          SNH1res@CovConfCases[[CvCase]] <- 
            CorrHetSNSol(NH1res,SNH1res,CvCase,Conf,Xscld,Xmean,Xsd,XsdOutPrd,grouping,Mxt,lglikdif,limlnk2,OptCntrl,getvcov=FALSE,...)
          SNH1res@logLiks[CvCase] <- LogLik <- SNH1res@CovConfCases[[CvCase]]$logLik
          nfreepar <- SKnpar(Conf,p,q)
          SNH1res@AICs[CvCase] <- -2*LogLik + 2*nfreepar
          SNH1res@BICs[CvCase] <- -2*LogLik + log(n)*nfreepar
        }
      }
      AICs <- c(NH1res@AICs,SNH1res@AICs)
      BICs <- c(NH1res@BICs,SNH1res@BICs)    			
      if (SelCrit=="AIC") { 
        bestmod <- which.min(AICs)
      }  else if (SelCrit=="BIC") {
        bestmod <- which.min(BICs)
      }
      H1res <- new("IdtMxNandSNDE",grouping=grouping,NMod=NH1res,SNMod=SNH1res,
        ModelNames=c(NH1res@ModelNames,SNH1res@ModelNames),ModelType=c(NH1res@ModelType,SNH1res@ModelType),
        ModelConfig=c(NH1res@ModelConfig,SNH1res@ModelConfig),SelCrit=SelCrit,NIVar=Idt@NIVar,
	logLiks=c(NH1res@logLiks,SNH1res@logLiks),AICs=AICs,BICs=BICs,BestModel=bestmod,SngD=FALSE,Ngrps=k)  
    }

    if (length(H1res@BestModel)==0)  { 
      return(error(onerror,"Procedure MANOVA failed to find a valid alternative model mle estimate\n"))
    }
    H0Ll <- H0res@logLiks[H1res@BestModel]
    H1Ll <- H1res@logLiks[H1res@BestModel]
    ChiSq <- 2*(H1Ll-H0Ll)
    names(H0Ll) <- names(H1Ll) <- names(ChiSq) <- NULL

    if (Mxt=="Hom" || Mxt=="Loc")  {
      df <- p*(k-1)
    }  else if (Mxt=="Het" || Mxt=="Gen") {
      if ( Model=="Normal" || (Model=="NrmandSKN" && bestmod <= 5) )
      {  
        BestConf <- H1res@ModelConfig[H1res@BestModel]
        if (CovCaseArg && BestConf>1)  BestConf <- BestConf+1  # Convert CovCase indices to Config indices
        df <- npar(BestConf,p,q,Ngrps=k,Mxt="Het") - npar(BestConf,p,q,Ngrps=1,Mxt="Hom")
      }  else {
        if (Model=="SKNormal")  {
          BestConf <- H1res@ModelConfig[H1res@BestModel]
        }  else {
          BestConf <- SNH1res@ModelConfig[SNH1res@BestModel]
        }
        if (CovCaseArg && BestConf>1)  BestConf <- BestConf+1  # Convert CovCase indices to Config indices
        df <- SKnpar(BestConf,p,q,Ngrps=k,Mxt="GenMod") - SKnpar(BestConf,p,q,Ngrps=1)
      }
    } 
    pvalue <- pchisq(ChiSq,df,lower.tail=FALSE)
    if (Model=="Normal" && (Mxt=="Hom" || Mxt=="Loc") )  { 
      return( new("IdtClMANOVA",NIVar=Idt@NIVar,grouping=grouping,H0res=H0res,H1res=H1res,
        ChiSq=ChiSq,df=df,pvalue=pvalue,H0logLik=H0Ll,H1logLik=H1Ll)  )
    }  else if (Model=="Normal" && (Mxt=="Het" || Mxt=="Gen") )  { 
      return( new("IdtHetNMANOVA",NIVar=Idt@NIVar,grouping=grouping,H0res=H0res,H1res=H1res,
        ChiSq=ChiSq,df=df,pvalue=pvalue,H0logLik=H0Ll,H1logLik=H1Ll)  )
    }  else if (Model=="SKNormal" && (Mxt=="Hom" || Mxt=="Loc") )  { 
      return( new("IdtLocSNMANOVA",NIVar=Idt@NIVar,grouping=grouping,H0res=H0res,H1res=H1res,
        ChiSq=ChiSq,df=df,pvalue=pvalue,H0logLik=H0Ll,H1logLik=H1Ll)  )
    }  else if (Model=="SKNormal" && (Mxt=="Het" || Mxt=="Gen") )  { 
      return( new("IdtGenSNMANOVA",NIVar=Idt@NIVar,grouping=grouping,H0res=H0res,H1res=H1res,
        ChiSq=ChiSq,df=df,pvalue=pvalue,H0logLik=H0Ll,H1logLik=H1Ll)  )
    }  else if (Model=="NrmandSKN" && (Mxt=="Hom" || Mxt=="Loc") )  { 
      return( new("IdtLocNSNMANOVA",NIVar=Idt@NIVar,grouping=grouping,H0res=H0res,H1res=H1res,
        ChiSq=ChiSq,df=df,pvalue=pvalue,H0logLik=H0Ll,H1logLik=H1Ll)  )
    }  else if (Model=="NrmandSKN" && (Mxt=="Het" || Mxt=="Gen") )  { 
      return( new("IdtGenNSNMANOVA",NIVar=Idt@NIVar,grouping=grouping,H0res=H0res,H1res=H1res,
        ChiSq=ChiSq,df=df,pvalue=pvalue,H0logLik=H0Ll,H1logLik=H1Ll)  )
    }
  }
)

MANOVAPermTest <- function(MANOVAres, Idt, grouping, nrep=200,
    Model=c("Normal","SKNormal","NrmandSKN"), CovCase=1:4,
    SelCrit=c("BIC","AIC"), Mxt=c("Hom","Het","Loc","Gen"), CVtol=1.0e-5, k2max=1e6, 
    OptCntrl=list(), onerror=c("stop","warning","silentNull"), ...)
{
   if (inherits(MANOVAres,"IdtMANOVA")==FALSE) stop("Argument MANOVAres is not of class IdtMANOVA\n")
   ChiSq <- MANOVAres@ChiSq
   if (class(Idt)[1]!="IData") stop("Argument Idt is not of class IData\n")
   if (!is.factor(grouping)) stop("Argument rouping is not a factor\n")
   n <- Idt@NObs
   if (length(grouping) != n) stop("The numbers of data and partition observations are different\n")
   Mxt <- match.arg(Mxt)
   Model <- match.arg(Model)
   SelCrit <- match.arg(SelCrit)
   
   grouping <- factor(grouping,exclude=NULL)
   grplvls <- levels(grouping)
   nk <- as.numeric(table(grouping))
   k <- length(nk)
   cumsumnk <- cumsum(nk)
   curgrouping <- grouping
   empdist <- numeric(nrep)
   for (rep in 1:nrep) {
     permut <- sort.int(runif(n),index.return=TRUE)$ix
     lstind <- 0
     for (g in 1:k) {
       frstind <- lstind+1
       lstind <- cumsumnk[g] 
       curgrouping[permut[frstind:lstind]] <- grplvls[g]
     }
     empdist[rep] <- MANOVA(Idt,curgrouping,Model,CovCase,SelCrit,Mxt,CVtol,k2max,OptCntrl,onerror,...)@ChiSq
   }
   pvalue <- length(which(ChiSq<empdist))/nrep
   cat("Permutation p-value of MANOVA statistic",ChiSq,":",pvalue,"\n")
   pvalue 
} 

setMethod("summary",					
  signature(object = "IdtE"),
  function(object)
  {
    cat("Log likelihoods:\n")
    print(object@logLiks)
    if (object@SelCrit=="AIC") { 
      cat("Akaike Information Criteria:\n")
      print(object@AICs) 
    }  else if (object@SelCrit=="BIC") { 
      cat("Bayesian (Schwartz) Information Criteria:\n")
      print(object@BICs)
    }
    cat("Selected model:\n")
    print(names(object@BestModel))
    invisible()
  }
)

setMethod("show",					
  signature(object = "IdtE"),
  function(object)
  {
    summary(object)
    cat("\nSelected model parameter estimates:\n")
    print(coef(object))
    invisible()
  }
)

setMethod("testMod",					
  signature(ModE = "IdtE"),
  function(ModE, RestMod=ModE@ModelConfig[2]:length(ModE@ModelConfig), FullMod="Next")
  {
    if ( substr(ModE@ModelNames[1],2,8) == "ModCovC" || substr(ModE@ModelNames[1],3,9) == "ModCovC")
    {
      CovCaseArg <- TRUE  
      nCovCases <- 4
    }  else {
      CovCaseArg <- FALSE  
      nCovCases <- 5
    }     
    if (is.character(RestMod))  {  
      RestMod <- sapply(RestMod,function(Rmd) which(Rmd==ModE@ModelNames))
    }
    if (FullMod[1]!="Next" && FullMod[1]!="All" && is.character(FullMod))  { 
      FullMod <- sapply(FullMod,function(Fmd) which(Fmd==ModE@ModelNames))
    }
    if (is.element(ModE@ModelConfig[1],RestMod))  {
      stop("Model",ModE@ModelNames,"can not be sepecified as a restricted model\n",
      "since it is the most general model that has been estimated\n")
    }
    if (length(ModE@ModelType)>nCovCases)  {
      EType <- "NrmandSKN"
    }  else if (ModE@ModelType[1]=="Normal") {
      EType <- "Normal"
    }  else if (ModE@ModelType[1]=="SkewNormal") {
      EType <- "SKNormal"
    }
    if (is.character(RestMod))  {  
      RestMod <- sapply(RestMod,function(Rmd) which(Rmd==ModE@ModelNames))
    }
    if (FullMod[1]!="Next" && FullMod[1]!="All" && is.character(FullMod))  { 
      FullMod <- sapply(FullMod,function(Fmd) which(Fmd==ModE@ModelNames))
    }
    if (EType != "NrmandSKN" && is.element(ModE@ModelConfig[1],RestMod))  { 
      stop("Model",ModE@ModelNames[1],"can not be sepecified as a restricted model\n",
        "since it is the most general model that has been estimated\n")
    }
    if (EType == "NrmandSKN" && is.element(nCovCases+1,RestMod))  { 
      MostGSNInd <- which(RestMod==nCovCases+1)
      if (is.element(ModE@ModelConfig[1],RestMod)) {
        RestMod <- RestMod[-MostGSNInd]
      } else {
        RestMod[MostGSNInd] <- ModE@ModelConfig[1]
      }  
    }
    TestRes <- list()
    RestModels <- character()
    FullModels <- character()
    for (RMind in RestMod)  {
      if (FullMod[1]=="Next") {
        FMindices <- NextModel(RMind,Model=EType,CovCaseArg=CovCaseArg)
      }  else if (FullMod[1]=="All")  {
        FMindices <- NestedBy(RMind,Model=EType,CovCaseArg=CovCaseArg)
      }  else  {
        FMindices <- intersect(FullMod,NestedBy(RMind,Model=EType,CovCaseArg=CovCaseArg))
      }
      for (FMind in FMindices)  {
        H0Ll <- ModE@logLiks[RMind]
        H1Ll <- ModE@logLiks[FMind]
        ChiSq <- 2*(H1Ll-H0Ll)
        q <- ModE@NIVar

        RConfig <- ModE@ModelConfig[RMind]
        FConfig <- ModE@ModelConfig[FMind]
        if (CovCaseArg)  {  # Convert CovCase indices to Config indices
          if (RConfig>1)  RConfig <- RConfig+1
          if (FConfig>1)  FConfig <- FConfig+1
        }
        if (ModE@SngD)  {
          if (ModE@ModelType[RMind]=="Normal")  {
            RMnpar <- npar(RConfig,2*q,q)
          }  else  {
            RMnpar <- SKnpar(RConfig,2*q,q)
          }
          if (ModE@ModelType[FMind]=="Normal")  {
            FMnpar <- npar(FConfig,2*q,q)
          }  else  {
            FMnpar <- SKnpar(FConfig,2*q,q)
          }
        } else {
          if (ModE@ModelType[RMind]=="Normal")  {
            RMnpar <- npar(RConfig,2*q,q,Ngrps=ModE@Ngrps,Mxt="Het")
          }  else  {
            RMnpar <- SKnpar(RConfig,2*q,q,Ngrps=ModE@Ngrps,Mxt="GenMod")
          }
          if (ModE@ModelType[FMind]=="Normal")  {
            FMnpar <- npar(FConfig,2*q,q,Ngrps=ModE@Ngrps,Mxt="Het")
          }  else  {
            FMnpar <- SKnpar(FConfig,2*q,q,Ngrps=ModE@Ngrps,Mxt="GenMod")
          }
        }

        df <-  FMnpar - RMnpar
        pvalue <- pchisq(ChiSq,df,lower.tail=FALSE)
        resi <- new("LRTest",H0logLik=H0Ll,H1logLik=H1Ll,ChiSq=ChiSq,df=df,pvalue=pvalue)
        TestRes <- c(TestRes,resi)
        RestModels <- c(RestModels,ModE@ModelNames[RMind])
        FullModels <- c(FullModels,ModE@ModelNames[FMind])
      }
    }
    new("ConfTests",TestRes=TestRes,RestModels=RestModels,FullModels=FullModels)
  }
)

NestedBy <- function(ModelInd,Model=c("Normal","SKNormal","NrmandSKN"),CovCaseArg)
{
  Model <- match.arg(Model)
  if (!CovCaseArg)  { return(ConfNestedBy(ModelInd,Model)) }
  if (Model=="Normal" || Model=="SKNormal")
  {
    if (ModelInd==1)  { stop("Configuration 1 is the most general model in this analysis\n") }
    if (ModelInd==2 || ModelInd==3)  { return(1) }
    if (ModelInd==4)  { return(1:3) }
  }  else if (Model=="NrmandSKN")  {
    if (ModelInd==1)  { return(5) }
    if (ModelInd==2)  { return(c(1,5,6)) }
    if (ModelInd==3)  { return(c(1,5:7)) }
    if (ModelInd==4)  { return(c(1:3,5:8)) }
    if (ModelInd==5)  { stop("A Skew-Normal distribution with Configuration 1 is the most general model in this analysis\n") }
    if (ModelInd==6 || ModelInd==7)  { return(5) }
    if (ModelInd==8)  { return(5:7) }
  }
}

NextModel <- function(ModelInd,Model=c("Normal","SKNormal","NrmandSKN"),CovCaseArg)
{
  Model <- match.arg(Model)
  if (!CovCaseArg)  { return(ConfNestedBy(ModelInd,Model)) }
  if (Model=="Normal" || Model=="SKNormal")
  {
    if (ModelInd==1)  { stop("Configuration 1 is the most general model in this analysis\n") }
    if (ModelInd==2 || ModelInd==3)  { return(1) }
    if (ModelInd==4)  { return(2:3) }
  }  else if (Model=="NrmandSKN")  {
    if (ModelInd==1)  { return(5) }
    if (ModelInd==2)  { return(c(1,6)) }
    if (ModelInd==3)  { return(c(1,7)) }
    if (ModelInd==4)  { return(c(2:3,8)) }
    if (ModelInd==5)  { stop("A Skew-Normal distribution with Configuration 1 is the most general model in this analysis\n") }
    if (ModelInd==6 || ModelInd==7)  { return(5) }
    if (ModelInd==8)  { return(6:7) }
  }
}

ConfNestedBy <- function(ModelInd,Model=c("Normal","SKNormal","NrmandSKN"))
{
  Model <- match.arg(Model)
  if (Model=="Normal" || Model=="SKNormal")
  {
    if (ModelInd==1)  { stop("Configuration 1 is the most general model in this analysis\n") }
    if (ModelInd==2)  { return(1) }
    if (ModelInd==3 || ModelInd==4)  { return(1:2) }
    if (ModelInd==5)  { return(1:4) }
  }  else if (Model=="NrmandSKN")  {
    if (ModelInd==1)  { return(6) }
    if (ModelInd==2)  { return(c(1,6:7)) }
    if (ModelInd==3)  { return(c(1:2,6:8)) }
    if (ModelInd==4)  { return(c(1:2,6:7,9)) }
    if (ModelInd==5)  { return(c(1:4,6:10)) }
    if (ModelInd==6)  { stop("A Skew-Normal distribution with Configuration 1 is the most general model in this analysis\n") }
    if (ModelInd==7)  { return(6) }
    if (ModelInd==8 || ModelInd==9)  { return(6:7) }
    if (ModelInd==10)  { return(6:9) }
  }
}

ConfNextModel <- function(ModelInd,Model=c("Normal","SKNormal","NrmandSKN"))
{
  Model <- match.arg(Model)
  if (Model=="Normal" || Model=="SKNormal")
  {
    if (ModelInd==1)  { stop("Configuration 1 is the most general model in this analysis\n") }
    if (ModelInd==2)  { return(1) }
    if (ModelInd==3 || ModelInd==4)  { return(2) }
    if (ModelInd==5)  { return(3:4) }
  }  else if (Model=="NrmandSKN")  {
    if (ModelInd==1)  { return(6) }
    if (ModelInd==2)  { return(c(1,7)) }
    if (ModelInd==3)  { return(c(2,8)) }
    if (ModelInd==4)  { return(c(2,9)) }
    if (ModelInd==5)  { return(c(3:4,10)) }
    if (ModelInd==6)  { stop("A Skew-Normal distribution with Configuration 1 is the most general model in this analysis\n") }
    if (ModelInd==7)  { return(6) }
    if (ModelInd==8 || ModelInd==9)  { return(7) }
    if (ModelInd==10)  { return(8:9) }
  }
}

setMethod("show",					
  signature(object = "ConfTests"),
  function(object)
  {
    for (i in 1:length(object@TestRes))  {
      if (!is.na(object@TestRes[[i]]@ChiSq)) 
      {
        cat("Testing Model",object@RestModels[i],"against alternative",object@FullModels[i],":\n")
        print(object@TestRes[[i]])
      }
    }
    invisible()
  }
)

setMethod("BestModel",
  signature(ModE = "IdtE"),
  function(ModE,SelCrit=c("IdtCrt","BIC","AIC"))
  {
    SelCrit <- match.arg(SelCrit)
    if (SelCrit == "IdtCrt")  {
      return(ModE@BestModel)
    }  else if (SelCrit == "AIC")  {
      return(which.min(ModE@AICs))
    }  else if (SelCrit == "BIC")  {
      return(which.min(ModE@BICs))
    }
  }
)

setMethod("coef",
  signature(object = "IdtNDE"),
  function(object,selmodel=BestModel(object),...)
  {
    if (class(object)[1]=="IdtSngNDE" || class(object)[1]=="IdtMxNDE")  {
      return(list(mu=object@mleNmuE,Sigma=object@CovConfCases[[selmodel]]$mleSigE))
    } else if (class(object)[1]=="IdtSngNDRE" || class(object)[1]=="IdtMxNDRE")  {
      return(list(mu=object@RobNmuE,Sigma=object@CovConfCases[[selmodel]]$RobSigE))
    }
  }
)

setMethod("stdEr",
  signature(x = "IdtNDE"),
  function(x,selmodel=BestModel(x),...)
  {
    list(mu=x@mleNmuEse,Sigma=x@CovConfCases[[selmodel]]$mleSigEse)
  }
)

setMethod("vcov",
  signature(object = "IdtNDE"),
  function(object,selmodel=BestModel(object),...)
  {
    object@CovConfCases[[selmodel]]$mlevcov
  }
)

setMethod("mean", signature(x = "IdtNDE"), function(x) coef(x)$mu )
setMethod("var", signature(x ="IdtNDE"), function(x) coef(x)$Sigma )
setMethod("sd", signature(x ="IdtE"), function(x,na.rm=FALSE) sqrt(diag(var(x))) )

setMethod("AIC", signature(object="IdtE"), 
  function(object,...,k=2) {
    if (missing(k)) nbarg <- 2
    else nbarg <- 3
    arguments <- match.call(expand.dots=TRUE)
    if (length(arguments)==nbarg) {
      aic <- object@AICs[object@BestModel]
      if (k==2) return(aic)
      else {
        llik <- object@logLiks[object@BestModel] 
        return(-2*llik + (k/2)*(aic+2*llik))
      }
    }  
    
    nobjs <- length(arguments)-nbarg+1
    aics <- numeric(nobjs)
    objl <- vector("list",nobjs-1)
    for (i in 1:nobjs) {
      if (i==1) {
        aics[1] <- object@AICs[object@BestModel]     
        names(aics)[1] <- names(object@AICs[object@BestModel])
      }  else {
        objl[[i-1]] <- eval(arguments[[1+i]],sys.parent())
        aics[i] <- objl[[i-1]]@AICs[objl[[i-1]]@BestModel]
        names(aics)[i] <- names(objl[[i-1]]@AICs[objl[[i-1]]@BestModel])
      }
    } 

    if (k==2) return(aics)
    else {
      lliks <- numeric(nobjs)
      for (i in 1:nobjs) {
        if (i==1) lliks[1] <- object@logLiks[object@BestModel]     
        else lliks[i] <- objl[[i-1]]@logLiks[objl[[i-1]]@BestModel]
      } 
    }    
    
    return(-2*lliks + (k/2)*(aics+2*lliks))
  }  
)

setMethod("BIC", signature(object="IdtE"), 
  function(object,...) {
    arguments <- match.call(expand.dots=TRUE)
    if (length(arguments)==2) return(object@BICs[object@BestModel])

    nobjs <- length(arguments)-1
    bics <- numeric(nobjs)
    objl <- vector("list",nobjs-1)
    for (i in 1:nobjs) {
      if (i==1) {
        bics[1] <- object@BICs[object@BestModel]     
        names(bics)[1] <- names(object@BICs[object@BestModel])
      }  else {
        objl[[i-1]] <- eval(arguments[[1+i]],sys.parent())
        bics[i] <- objl[[i-1]]@BICs[objl[[i-1]]@BestModel]
        names(bics)[i] <- names(objl[[i-1]]@BICs[objl[[i-1]]@BestModel])
      }
    } 
    return(bics)
  }  
)

setMethod("logLik", signature(object="IdtE"), 
  function(object,...) {
    arguments <- match.call(expand.dots=TRUE)
    if (length(arguments)==2) return(object@logLiks[object@BestModel])
            
    nobjs <- length(arguments)-1
    logLiks <- numeric(nobjs)
    objl <- vector("list",nobjs-1)
    for (i in 1:nobjs) {
      if (i==1) {
        logLiks[1] <- object@logLiks[object@BestModel]     
        names(logLiks)[1] <- names(object@logLiks[object@BestModel])
      }  else {
        objl[[i-1]] <- eval(arguments[[1+i]],sys.parent())
        logLiks[i] <- objl[[i-1]]@logLiks[objl[[i-1]]@BestModel]
        names(logLiks)[i] <- names(objl[[i-1]]@logLiks[objl[[i-1]]@BestModel])
      }
    } 
    return(logLiks)
  }  
)

setMethod("cor",
  signature(x ="IdtNDE"),
  function(x)
  { 
    Sig <- coef(x)$Sigma
    if (length(dim(Sig))==2) {
      return(cov2cor(Sig))
    } else if (length(dim(Sig))==3) {
      return(array(apply(Sig,3,cov2cor),dim=dim(Sig),dimnames=dimnames(Sig)))
    }
  }
)

setMethod("vcov",
  signature(object = "IdtMxNDE"),
  function(object,selmodel=BestModel(object),group=NULL,...)
  {
    if (object@Hmcdt)  {
      object@CovConfCases[[selmodel]]$mlevcov
    } else {
      if (is.null(group))
      {
        warning(paste("vcov returned as three-dimensional array with a different var-cov matrix for each group,\n",
          "which was identified by the level of the third array dimension\n")) 
        object@CovConfCases[[selmodel]]$mlevcov
      } else {
        object@CovConfCases[[selmodel]]$mlevcov[[,,group]]
      }
    }
  }
)

setMethod("ObsLogLiks",					
  signature(object = "IdtSngNDE"),
  function(object,Idt,Conf=object@BestModel)
  {
    p <- 2*Idt@NIVar
    c0 <- -0.5*(p*log(2*pi))
    Xdev <- scale(cbind.data.frame(Idt@MidP,Idt@LogR),center=object@mleNmuE,scale=FALSE)
    if (Conf!=5) {
      SigISr <- t(backsolve(chol(object@CovConfCases[[Conf]]$mleSigE),diag(p)))
      apply(Xdev,1,ILogLikNC1,SigmaSrInv=SigISr,const=c0+sum(log(diag(SigISr))))
    }  else {
      IVar <- 1./diag(object@CovConfCases[[4]]$mleSigE)
      apply(Xdev,1,ILogLikDNC,IVar=IVar,const=c0-0.5*sum(log(IVar)))
    }
  }
)

setMethod("show",					
  signature(object = "LRTest"),
  function(object)
  {
    cat("Null Model log-likelihood:",object@H0logLik,"\n")
    cat("Full Model log-likelihood:",object@H1logLik,"\n")
    cat("Chi-squared statistic:",object@ChiSq,"\n")
    cat("degrees of freedom:",object@df,"\n")
    cat("p-value:",object@pvalue,"\n\n")
    invisible()
  }
)

setMethod("H1res", signature(object = "IdtMANOVA"), function(object) object@H1res)

setMethod("H0res", signature(object = "IdtMANOVA"), function(object) object@H0res)

setMethod("summary",					
  signature(object = "IdtMANOVA"),
  function(object)
  {
    cat("Null Model Log likelihoods:\n")
    print(object@H0res@logLiks)
    cat("Full Model Log likelihoods:\n")
    print(object@H1res@logLiks)
    if (object@H1res@SelCrit=="AIC") { 
      cat("Full Model Akaike Information Criteria:\n")
      print(object@H1res@AICs)
    }  else if (object@H1res@SelCrit=="BIC") { 
      cat("Full Model Bayesian (Schwartz) Information Criteria:\n")
      print(object@H1res@BICs)
    }
    cat("Selected Model:\n")
    print(names(object@H1res@BestModel))
    cat("\n")
    cat("Chi-squared statistic:",object@ChiSq,"\n")
    cat("degrees of freedom:",object@df,"\n")
    cat("p-value:",object@pvalue,"\n\n")
    if ( length(object@grouping)<=30 ) { 
      cat("Note: Given the small sample size, the use of the Chi-square distribution may not be appropriate.\n",
          "Alternatively, consider using the permutation test implemented in function MANOVAPermTest.\n",
          "Note that this may take a long time.\n")
    } 
    invisible()
  }
)

setMethod("show",					
  signature(object = "IdtMANOVA"),
  function(object) {
    summary(object)
    invisible()
  } 
)

