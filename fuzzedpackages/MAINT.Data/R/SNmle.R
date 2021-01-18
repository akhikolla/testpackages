setMethod("coef",
  signature(object = "IdtSNDE"),
  function(object,selmodel=BestModel(object),ParType=c("Centr","Direct","All"),...)
  {
    ParType <- match.arg(ParType)
     if (ParType=="Centr")  {
       return( list(mu=object@CovConfCases[[selmodel]]$muE,Sigma=object@CovConfCases[[selmodel]]$SigmaE,
         gamma1=object@CovConfCases[[selmodel]]$gamma1E) )
     }  else if (ParType=="Direct")  {
       return( list(ksi=object@CovConfCases[[selmodel]]$ksiE,Omega=object@CovConfCases[[selmodel]]$OmegaE,
         alpha=object@CovConfCases[[selmodel]]$alphaE) )
     }  else if (ParType=="All")  {
       return( list(mu=object@CovConfCases[[selmodel]]$muE,Sigma=object@CovConfCases[[selmodel]]$SigmaE,
         gamma1=object@CovConfCases[[selmodel]]$gamma1E,ksi=object@CovConfCases[[selmodel]]$ksiE,
         Omega=object@CovConfCases[[selmodel]]$OmegaE,alpha=object@CovConfCases[[selmodel]]$alphaE) )
     }
  }
)

setMethod("stdEr",
  signature(x = "IdtSNDE"),
  function(x,selmodel=BestModel(x),...)
  {
    modres <- x@CovConfCases[[selmodel]]
    if (modres$status!="Regular")
    {
      if (modres$status=="Invalid") {
        warning("Standard errors were not computed for model ",x@ModelNames[selmodel],"\n",
         "because the computation of the observed Information Matrix was not successful.\n")
#      } else if (modres$status=="SingInf") {
#        warning("Standard errors were not computed for model ",x@ModelNames[selmodel],"\n",
#         "because the mle estimates were too close to the singularity point where all gamma1 equal 0.\n")
#      }  else if (modres$status=="InstInf") {
#        warning("Standard errors were not computed for model ",x@ModelNames[selmodel],"\n",
#        "because the mle estimates resulted in a numerically unstable information matrix.\n")
      }  else if (modres$status=="Onborder") {
        warning("Standard errors were not computed for model ",x@ModelNames[selmodel]," because the mle estimates are too close \n",
        "to the frontier of the parameter space, where classical maximum likeklihood theory does not apply.\n")
#      }  else if (modres$status=="PositiveScore") {
#        warning("Standard errors were not computed for model ",x@ModelNames[selmodel]," because at the current estimate the score\n", 
#         "function (i.e., the derivative of the log-likelihhod) is stricty positive\n")
#      }  else if (modres$status=="SingOmega") {
#        warning("Standard errors were not computed for model ",x@ModelNames[selmodel],"\n",
#        "because the mle estimates resulted in a singular Omega matrix.\n")
      }
      return(NULL)
    }   
    list(mu=modres$muEse,Sigma=modres$SigmaEse,gamma1=modres$gamma1Ese)
  }
)

setMethod("vcov",
  signature(object = "IdtSNDE"),
  function(object,selmodel=BestModel(object),...)
  {
    modres <- object@CovConfCases[[selmodel]]
    if (modres$status!="Regular")
    {
      if (modres$status=="Invalid") {
        warning("The asymptotic covariance matrix was not computed for model ",object@ModelNames[selmodel],"\n",
         "because the computation of the observed Information Matrix was not successful.\n")
#      } else if (modres$status=="SingInf") {
#        warning("The asymptotic covariance matrix was not computed for model ",object@ModelNames[selmodel],"\n",
#         "because the mle estimates were too close to the singularity point where all gamma1 equal 0.\n")
#      }  else if (modres$status=="InstInf") {
#        warning("The asymptotic covariance matrix was not computed for model ",object@ModelNames[selmodel],"\n",
#        "because the mle estimates resulted in a numerically unstable information matrix.\n")
      }  else if (modres$status=="Onborder") {
        warning("The asymptotic covariance matrix was not computed for model ",object@ModelNames[selmodel],
         " because the mle estimates are too close\nto the frontier of the parameter space, where",
         " classical maximum likeklihood theory does not apply.\n")
      }  else if (modres$status=="PositiveScore") {
        warning("The asymptotic covariance matrix was not computed for model ",object@ModelNames[selmodel],
         " because at the current estimate the score\n", 
         "function (i.e., the derivative of the log-likelihhod) is stricty positive\n")
      }  else if (modres$status=="SingOmega") {
        warning("Standard errors were not computed for model ",object@ModelNames[selmodel],"\n",
        "because the mle estimates resulted in a singular Omega matrix.\n")
      }
      return(NULL)
    }
    object@CovConfCases[[selmodel]]$mleCPvcov
  }
)

setMethod("vcov",
  signature(object = "IdtMxSNDE"),
  function(object,selmodel=BestModel(object),group=NULL,...)
  {
    modres <- object@CovConfCases[[selmodel]]
    if (modres$status!="Regular")
    {
      if (modres$status=="Invalid") {
        warning("The asymptotic covariance matrix was not computed for model ",object@ModelNames[selmodel],"\n",
         "because the computation of the observed Information Matrix was not successful.\n")
#      } else if (modres$status=="SingInf") {
#        warning("The asymptotic covariance matrix was not computed for model ",object@ModelNames[selmodel],"\n",
#         "because the mle estimates were too close to the singularity point where all gamma1 equal 0.\n")
#      }  else if (modres$status=="InstInf") {
#        warning("The asymptotic covariance matrix was not computed for model ",object@ModelNames[selmodel],"\n",
#        "because the mle estimates resulted in a numerically unstable information matrix.\n")
      }  else if (modres$status=="Onborder") {
        warning("The asymptotic covariance matrix was not computed for model ",object@ModelNames[selmodel],
         " because the mle estimates are too close\nto the frontier of the parameter space, where",
         " classical maximum likeklihood theory does not apply.\n")
      }  else if (modres$status=="PositiveScore") {
        warning("The asymptotic covariance matrix was not computed for model ",object@ModelNames[selmodel],
         " because at the current estimate the score\n", 
         "function (i.e., the derivative of the log-likelihhod) is stricty positive\n")
      }  else if (modres$status=="SingOmega") {
        warning("Standard errors were not computed for model ",object@ModelNames[selmodel],"\n",
        "because the mle estimates resulted in a singular Omega matrix.\n")
      }
      return(NULL)
    }
    if (object@Hmcdt)  {
      object@CovConfCases[[selmodel]]$mleCPvcov
    } else {
      if (is.null(group))
      {
        warning(paste("vcov returned as three-dimensional array with a different var-cov matrix for each group,\n",
          "which was identified by the level of the third array dimension\n")) 
        object@CovConfCases[[selmodel]]$mleCPvcov
      } else {
        object@CovConfCases[[selmodel]]$mleCPvcov[[,,group]]
      }
    }
  }
)

setMethod("mean", signature(x = "IdtSNDE"), function(x) { coef(x)$mu } )
setMethod("var", signature(x ="IdtSNDE"), function(x) { coef(x)$Sigma } )

setMethod("cor",
  signature(x ="IdtSNDE"),
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

setMethod("coef",
  signature(object = "IdtNandSNDE"),
  function(object,selmodel=BestModel(object),ParType=c("Centr","Direct","All"),...)
  {
    ParType <- match.arg(ParType)
    if (object@ModelType[selmodel]=="Normal")  {
      return(coef(object@NMod))
    } else if (object@ModelType[selmodel]=="SkewNormal")  {
      return(coef(object@SNMod,ParType=ParType))
    }
  }
)

setMethod("stdEr",
  signature(x = "IdtNandSNDE"),
  function(x,selmodel=BestModel(x),...)
  {
    if (x@ModelType[selmodel]=="Normal")  {
      stdEr(x@NMod,selmodel=selmodel)
    }  else if (x@ModelType[selmodel]=="SkewNormal")  {
      stdEr(x@SNMod,selmodel=selmodel-length(x@NMod@ModelConfig))
    }
  }
)

setMethod("vcov",
  signature(object = "IdtNandSNDE"),
  function(object,selmodel=BestModel(object),...)
  {
    if (object@ModelType[selmodel]=="Normal")  {
      return(vcov(object@NMod,selmodel=selmodel))
    }  else if (object@ModelType[selmodel]=="SkewNormal")  {
      return(vcov(object@SNMod,selmodel=selmodel-length(object@NMod@ModelConfig)))
    }
  }
)

setMethod("mean", signature(x = "IdtNandSNDE"), function(x) { coef(x)$mu } )
setMethod("var", signature(x ="IdtNandSNDE"), function(x) { coef(x)$Sigma } )

setMethod("cor",
  signature(x ="IdtNandSNDE"),
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




