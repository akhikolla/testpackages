setClassUnion("extmatrix",c("matrix","NULL"))
setClassUnion("extlogical",c("logical","NULL"))
setClassUnion("extnumeric",c("numeric","NULL"))
setClassUnion("extcharacter",c("character","NULL"))
setClassUnion("extinteger",c("integer","NULL"))
setClass("IData",slots=c(MidP="data.frame",LogR="data.frame",ObsNames="character",VarNames="character",
  NObs="numeric",NIVar="numeric",NbMicroUnits="integer"))
setClass("IdtE",slots=c(ModelNames="character",ModelType="character",ModelConfig="numeric",NIVar="numeric",SelCrit="character",
  logLiks="numeric",BICs="numeric",AICs="numeric",BestModel="numeric",SngD="logical"),contains="VIRTUAL")
setClass("IdtMclust",slots=c(call="call",data="IData",NObs="numeric",NIVar="numeric",SelCrit="character",Hmcdt="logical",
   BestG="integer",BestC="integer",logLiks="numeric",logLik="numeric",BICs="numeric",bic="numeric",AICs="numeric",aic="numeric",
   parameters="list",z="extmatrix",classification="extcharacter",allres="list"))
setClass("IdtMclustEl",slots=c(NObs="numeric",NIVar="numeric",SelCrit="character",Hmcdt="logical",Conf="integer",nG="integer",
  logLik="numeric",alllnLik="numeric",bic="numeric",aic="numeric",parameters="list",z="extmatrix",classification="extcharacter"))
setClass("IdtSngDE",contains=c("IdtE","VIRTUAL"))
setClass("IdtMxE",slots=c(grouping="factor",Ngrps="numeric"),contains=c("IdtE","VIRTUAL"))
setClass("IdtSngNDE",slots=c(mleNmuE="numeric",mleNmuEse="numeric",CovConfCases="list"),contains="IdtSngDE")
setClass("IdtSngNDRE",slots=c(RobNmuE="numeric",CovConfCases="list",rawSet="numeric",RewghtdSet="numeric",
  RobMD2="numeric",cnp2="numeric",raw.cov="matrix",raw.cnp2="numeric",PerfSt="list"),contains="IdtSngDE")
setClass("IdtMxNDE",slots=c(Hmcdt="logical",mleNmuE="matrix",mleNmuEse="extmatrix",CovConfCases="list"),contains="IdtMxE")
setClass("IdtMxNDRE",slots=c(Hmcdt="logical",RobNmuE="matrix",CovConfCases="list",rawSet="numeric",RewghtdSet="numeric",
  RobMD2="numeric",cnp2="matrix",raw.cov="array",raw.cnp2="matrix",PerfSt="list"),contains="IdtMxE")
setClassUnion("IdtMxtNDE",c("IdtMxNDE","IdtMxNDRE"))
setClassUnion("IdtNDE",c("IdtSngNDE","IdtSngNDRE","IdtMxNDE","IdtMxNDRE"))
setClass("LRTest",slots=c(ChiSq="numeric",df="numeric",pvalue="numeric",H0logLik="numeric",H1logLik="numeric"))
setClass("ConfTests",slots=c(TestRes="list",RestModels="character",FullModels="character"))
setClass("IdtMANOVA",slots=c(NIVar="numeric",grouping="factor",H0res="IdtSngDE",H1res="IdtMxE"),contains="LRTest")
setClass("IdtClMANOVA",contains="IdtMANOVA")
setClass("IdtHetNMANOVA",contains="IdtMANOVA")
setClass("IdtLocSNMANOVA",contains="IdtMANOVA")
setClass("IdtGenSNMANOVA",contains="IdtMANOVA")
setClass("IdtLocNSNMANOVA",contains="IdtMANOVA")
setClass("IdtGenNSNMANOVA",contains="IdtMANOVA")
setClass("Idtlda",slots=c(prior="numeric",means="matrix",scaling="matrix",N="numeric",CovCase="numeric"))
setClass("Idtqda",slots=c(prior="numeric",means="matrix",scaling="array",ldet="numeric",lev="character",CovCase="numeric"))
setClass("IdtSNlocda",slots=c(prior="numeric",ksi="matrix",eta="numeric",scaling="matrix",
  mu="matrix",gamma1="numeric",N="numeric",CovCase="numeric"))
setClass("IdtSNgenda",slots=c(prior="numeric",ksi="matrix",eta="matrix",scaling="array",ldet="numeric",lev="character",
  mu="matrix",gamma1="matrix",CovCase="numeric"))
setClass("IdtSngSNDE",slots=c(CovConfCases="list"),contains="IdtSngDE")
setClass("IdtMxSNDE",slots=c(Hmcdt="logical",CovConfCases="list"),contains="IdtMxE")
setClassUnion("IdtSNDE",c("IdtSngSNDE","IdtMxSNDE"))
setClass("IdtSngNandSNDE",slots=c(NMod="IdtSngNDE",SNMod="IdtSngSNDE"),contains="IdtSngDE")
setClass("IdtMxNandSNDE",slots=c(NMod="IdtMxNDE",SNMod="IdtMxSNDE"),contains="IdtMxE")
setClassUnion("IdtNandSNDE",c("IdtSngNandSNDE","IdtMxNandSNDE"))
setClassUnion("IdtSNda",c("IdtSNlocda","IdtSNgenda"))

setClass("RobEstControl",
  slots=c(
    ncsteps="numeric",
    getalpha="character",
    rawMD2Dist="character",				
    MD2Dist="character",
    eta="numeric",
    multiCmpCor="character",				
    getkdblstar="character",
    outlin="character",
    trialmethod="character",
    m="numeric",
    reweighted="logical",
    k2max="numeric",
    otpType="character"
  ),
  prototype = list(
    ncsteps=200,
    getalpha = "TwoStep75",
    rawMD2Dist="ChiSq",				
    MD2Dist="ChiSq",
    eta=0.025,
    multiCmpCor="never",				
    getkdblstar="Twopplusone",
    outlin="MidPandLogR",
    trialmethod="simple",
    m=1,
    reweighted=TRUE,
    k2max=1e6,
    otpType="SetMD2andEst"
  ),
  contains="CovControlMcd"
)

setClass("EMControl",
  slots=c(nrep="numeric", maxiter="numeric", convtol="numeric", protol="numeric", 
    seed="extnumeric", pertubfct="numeric", k2max="numeric", MaxVarGRt="numeric"),
  prototype = list(nrep = 0, maxiter=1000, convtol=0.01, protol=1e-3, seed=NULL, pertubfct=1, k2max=1e6, MaxVarGRt=1e6)
)
setClass("IdtOutl",slots=c(outliers="extinteger", MD2="numeric",eta="numeric",RefDist="character",
  multiCmpCor="character",NObs="numeric",p="numeric",h="numeric",boolRewind="extlogical"))

setClass("summaryIData",slots=c(MidPsumar="table", Rngsumar="table",LogRsumar="table"))


setGeneric("nrow")
setGeneric("ncol")
setGeneric("rownames")
setGeneric("colnames")
setGeneric("names")
setGeneric("var")
setGeneric("sd",function(x,na.rm=FALSE) standardGeneric("sd"))                      # Check if I really need this ... !!! 
setGeneric("cor")
setGeneric("mean",signature="x")
setGeneric("plot",signature=c("x","y"))
setGeneric("summary",signature="object")
setGeneric("head",package="utils",signature="x")
setGeneric("tail",package="utils",signature="x")
setGeneric("coef",package="stats",signature="object")
setGeneric("stdEr",package="miscTools",signature="x")
setGeneric("vcov",package="stats",signature="object")
setGeneric("predict",package="stats",signature="object")
setGeneric("lda",package="MASS",signature="x")
setGeneric("qda",package="MASS",signature="x")
setGeneric("MidPoints",function(Idt) standardGeneric("MidPoints"))
setGeneric("LogRanges",function(Idt) standardGeneric("LogRanges"))
setGeneric("Ranges",function(Idt) standardGeneric("Ranges"))
setGeneric("mle",
  function(Idt, Model=c("Normal","SKNormal","NrmandSKN"),CovCase=1:4, SelCrit=c("BIC","AIC"), k2max=1e6, OptCntrl=list(),...)
  standardGeneric("mle"))
setGeneric("MANOVA",function(Idt, grouping, Model=c("Normal","SKNormal","NrmandSKN"), CovCase=1:4, SelCrit=c("BIC","AIC"), 
      Mxt=c("Hom","Het","Loc","Gen"), CVtol=1.0e-5, k2max=1e6, OptCntrl=list(), onerror=c("stop","warning","silentNull"), ...)
  standardGeneric("MANOVA"))
setGeneric("BestModel",function(ModE,SelCrit=c("IdtCrt","BIC","AIC"))  standardGeneric("BestModel"))
setGeneric("CovCase",function(object)  standardGeneric("CovCase"))
setGeneric("testMod", function(ModE, RestMod=ModE@ModelConfig[2]:length(ModE@ModelConfig), FullMod="Next")
  standardGeneric("testMod"))
setGeneric("H1res",  function(object) standardGeneric("H1res"))
setGeneric("H0res",  function(object) standardGeneric("H0res"))
setGeneric("snda",function(x, grouping, prior="proportions", ...) standardGeneric("snda"))
setGeneric("ObsLogLiks",function(object,Idt,Conf=object@BestModel) standardGeneric("ObsLogLiks"))
setGeneric("rbind",function(x, y, ...) standardGeneric("rbind"))
setGeneric("cbind",function(x, y, ...) standardGeneric("cbind"))

setClass("summaryIdtMclust",
  slots=c(title="character",modelName="character",Hmcdt="logical",
    NObs="numeric",NIVar="numeric",G="numeric",loglik="numeric",bic="numeric",
    pro="numeric",mean="matrix",covariance="array",classification="extcharacter",
    printParameters ="logical", printClassification="logical") 
)

setGeneric("fulltle",
  function(Idt, CovCase=1:4, SelCrit=c("BIC","AIC"), alpha=0.75, use.correction=TRUE, getalpha="TwoStep", 
    rawMD2Dist=c("ChiSq","HardRockeAsF","HardRockeAdjF"), MD2Dist=c("ChiSq","CerioliBetaF"),
    eta=0.025,multiCmpCor=c("never","always","iterstep"), outlin=c("MidPandLogR","MidP","LogR"), reweighted=TRUE, k2max=1e6,
    force=FALSE, ...)
  standardGeneric("fulltle"))

setGeneric("fasttle",
  function(Idt,
    CovCase=1:4,
    SelCrit=c("BIC","AIC"),
    alpha=control@alpha,
    nsamp = control@nsamp,
    seed=control@seed,
    trace=control@trace,
    use.correction=control@use.correction,
    ncsteps=control@ncsteps,
    getalpha=control@getalpha,
    rawMD2Dist=control@rawMD2Dist,				
    MD2Dist=control@MD2Dist,
    eta=control@eta,
    multiCmpCor=control@multiCmpCor,				
    getkdblstar=control@getkdblstar,
    outlin=control@outlin,
    trialmethod=control@trialmethod,
    m=control@m,
    reweighted = control@reweighted,
    k2max = control@k2max,
    otpType=control@otpType,
    control=RobEstControl(), ...)
  standardGeneric("fasttle"))

setGeneric("getMahaD2",function(IdtOtl) standardGeneric("getMahaD2")) 
setGeneric("geteta",function(IdtOtl) standardGeneric("geteta")) 
setGeneric("getRefDist",function(IdtOtl) standardGeneric("getRefDist")) 
setGeneric("getmultiCmpCor",function(IdtOtl) standardGeneric("getmultiCmpCor")) 

setGeneric("RobMxtDEst",
  function(Idt,grouping,Mxt=c("Hom","Het"),CovEstMet=c("Pooled","Globdev"),
    CovCase=1:4,SelCrit=c("BIC","AIC"),Robcontrol=RobEstControl(), l1medpar=NULL, ...)
  standardGeneric("RobMxtDEst"))
setGeneric("Roblda",
  function(x, grouping, prior="proportions", CVtol=1.0e-5, egvtol=1.0e-10, subset=1:nrow(x),
    CovCase=1:4, SelCrit=c("BIC","AIC"), silent=FALSE,  CovEstMet=c("Pooled","Globdev"), SngDMet=c("fasttle","fulltle"),
    k2max=1e6, Robcontrol=RobEstControl(), ...)
  standardGeneric("Roblda"))
setGeneric("Robqda",
  function(x, grouping, prior="proportions", CVtol=1.0e-5, subset=1:nrow(x),
    CovCase=1:4, SelCrit=c("BIC","AIC"), silent=FALSE, SngDMet=c("fasttle","fulltle"),
    k2max=1e6, Robcontrol=RobEstControl(), ...) 
  standardGeneric("Robqda"))

setGeneric("Idtmclust",
  function(Idt, G=1:9, CovCase=1:4, SelCrit=c("BIC","AIC"), Mxt=c("Hom","Het","HomandHet"), control=EMControl())
  standardGeneric("Idtmclust"))

setGeneric("pcoordplot",
  function(x,title="Parallel Coordinate Plot",Seq=c("AllMidP_AllLogR","MidPLogR_VarbyVar"),G=BestG(x),...)
  standardGeneric("pcoordplot"))

setGeneric("plotInfCrt", function(object,crt=object@SelCrit,legpos="bottomleft",...) standardGeneric("plotInfCrt"))

setGeneric("SelCrit",function (x) standardGeneric("SelCrit"))
setGeneric("Hmcdt",function (x) standardGeneric("Hmcdt"))
setGeneric("BestG",function (x) standardGeneric("BestG"))
setGeneric("BestC",function (x) standardGeneric("BestC"))
setGeneric("PostProb",function(x) standardGeneric("PostProb"))
setGeneric("logLik",function(object,...) standardGeneric("logLik"))                # Check if I really need this ... !!!
setGeneric("BIC",function(object,...) standardGeneric("BIC"))                      # Check if I really need this ... !!!
setGeneric("AIC",function(object,...,k=2) standardGeneric("AIC"))                      # Check if I really need this ... !!! 

setGeneric("parameters",function(x,model="BestModel") standardGeneric("parameters"))
setGeneric("pro",function(x,model="BestModel") standardGeneric("pro"))
setGeneric("classification",function(x,model="BestModel") standardGeneric("classification"))
setGeneric("NbMicroUnits",function(x) standardGeneric("NbMicroUnits"))




