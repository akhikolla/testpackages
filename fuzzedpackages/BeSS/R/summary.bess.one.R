summary.bess.one=function(object, ...){
   max.steps = object$max.steps
   df = sum(object$beta!=0)
   predictors = names(which(object$beta!=0))
   a=rbind(predictors, object$beta[predictors])
   cat("----------------------------------------------------------------------\n")
   cat("    Primal-dual active algorithm with maximum iteration being", max.steps, "\n\n")
   cat("    Best model with k =", df, "includes predictors:", "\n\n")
   print(object$beta[predictors])
   cat("\n")
   if(logLik(object)[2]>=0)
   cat("    log-likelihood:   ", logLik(object)[2],"\n") else cat("    log-likelihood:  ", logLik(object)[2],"\n")
   
   if(deviance(object)[2]>=0)
   cat("    deviance:         ", deviance(object)[2],"\n") else cat("    deviance:        ", deviance(object)[2],"\n")
   
   if(object$AIC>=0)
   cat("    AIC:              ", object$AIC,"\n") else cat("    AIC:             ", object$AIC,"\n")
   
   if(object$BIC>=0)
   cat("    BIC:              ", object$BIC,"\n") else cat("    BIC:             ", object$BIC,"\n")
   
   if(object$EBIC>=0)
   cat("    EBIC:             ", object$EBIC,"\n") else cat("    EBIC:            ", object$EBIC,"\n")
   cat("----------------------------------------------------------------------\n")
}

