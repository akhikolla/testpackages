summary.bess=function(object, ...){
  beta = object$beta
  if(object$method == "sequential"){
   K.opt.aic = which.min(object$AIC)
   K.opt.bic = which.min(object$BIC)
   K.opt.ebic = which.min(object$EBIC)
   predictors.aic = beta[,K.opt.aic]
   predictors.bic = beta[,K.opt.bic]
   predictors.ebic = beta[,K.opt.ebic]
   if(sum(predictors.aic!=0)>1) predictor.a = "predictors" else predictor.a = "predictor"
   if(sum(predictors.bic!=0)>1) predictor.b = "predictors" else predictor.b = "predictor"
   if(sum(predictors.ebic!=0)>1) predictor.e = "predictors" else predictor.e = "predictor"
   cat("-------------------------------------------------------------------------------\n")
   cat("    Primal-dual active algorithm with tuning parameter determined by sequential method", "\n\n")
   cat("    Best model determined by AIC includes" , sum(predictors.aic!=0), predictor.a, "with AIC =", 
       object$AIC[K.opt.aic], "\n\n")
   cat("    Best model determined by BIC includes" , sum(predictors.bic!=0), predictor.b, "with BIC =", 
       object$BIC[K.opt.bic], "\n\n")
   cat("    Best model determined by EBIC includes" , sum(predictors.ebic!=0), predictor.e, "with EBIC =", 
       object$EBIC[K.opt.ebic], "\n")
   cat("-------------------------------------------------------------------------------\n")
  } else {
    cat("------------------------------------------------------------------------------\n")
    cat("    Primal-dual active algorithm with tuning parameter determined by gsection method", "\n\n")
    if(sum(beta[,ncol(beta)]!=0)>0) predictor = "predictors" else predictor = "predictor"
    cat("    Best model includes", sum(beta[,ncol(beta)]!=0), predictor, "with", "\n\n")
    if(logLik(object)[length(logLik(object))]>=0)
      cat("    log-likelihood:   ", logLik(object)[length(logLik(object))],"\n") else cat("    log-likelihood:  ", logLik(object)[length(logLik(object))],"\n")
    
    if(deviance(object)[length(deviance(object))]>=0)
      cat("    deviance:         ", deviance(object)[length(deviance(object))],"\n") else  cat("    deviance:       ", deviance(object)[length(deviance(object))],"\n") 
    
    if(object$AIC[length(object$AIC)]>=0)
      cat("    AIC:              ", object$AIC[length(object$AIC)],"\n") else   cat("    AIC:             ", object$AIC[length(object$AIC)],"\n")
    
    if(object$BIC[length(object$BIC)]>=0)
      cat("    BIC:              ", object$BIC[length(object$BIC)],"\n") else   cat("    BIC:             ", object$BIC[length(object$BIC)],"\n")
    
    if(object$EBIC[length(object$EBIC)]>=0)
      cat("    EBIC:             ", object$EBIC[length(object$EBIC)],"\n") else    cat("    EBIC:            ", object$EBIC[length(object$EBIC)],"\n")
    cat("------------------------------------------------------------------------------\n")
  }
 
}
