#' @title SPLMM Tuning Function
#'
#' @description This function fits SPLMM over grids of lambda1 and/or lambda2 and determine the best fit model based on model selection information criterion. 
#'
#' @import penalized
#' 
#' @import methods
#' 
#' @import emulator 
#' 
#' @import miscTools
#' 
#' @import penalized
#' 
#' @param x
#'
#' @param y
#'
#' @param z
#'
#' @param grp
#'
#' @param lam1.seq
#'
#' @param lam2.seq
#'
#' @param nonpen.b
#'
#' @param nonpen.L
#'
#' @param penalty.b
#'
#' @param penalty.L
#'
#' @param CovOpt
#'
#' @param standardize
#'
#' @param control
#'
#' @return splmm.tuning
#'
#' @examples
#'
#' @export


splmmTuning <- function(x,y,z,grp,lam1.seq,lam2.seq,nonpen.b=1,nonpen.L=1,penalty.b=c("lasso","scad"),
                          penalty.L=c("lasso","scad"),CovOpt=c("nlminb","optimize"),standardize=TRUE,control=splmmControl()){
  
  
  if((length(lam1.seq)==1)&(length(lam2.seq)==1)){
    print("single lambda1 and single lambda2 detected, running splmm function.")
    fit = splmm(x=x,y=y,z=z,grp=grp,lam1=lam1.seq,lam2=lam2.seq,nonpen.b=nonpen.b,nonpen.L=nonpen.L,penalty.b=penalty.b,penalty.L=penalty.L,CovOpt=CovOpt,
          standardize=standardize,control=control)
    fit
  }else if((length(lam1.seq)>1)&(length(lam2.seq)==1)){
    print("single lambda2 detected, tuning for lambda1.")
    fit.list = list()
    BIC.vec = vector()
    AIC.vec = vector()
    BBIC.vec = vector()
    EBIC.vec = vector()
    lam1.tuning = TRUE
    lam2.tuning = FALSE
    
    for (i in 1:length(lam1.seq)) {
      fit = splmm(x=x,y=y,z=z,grp=grp,lam1=lam1.seq[i],lam2=lam2.seq,nonpen.b=nonpen.b,nonpen.L=nonpen.L,penalty.b=penalty.b,penalty.L=penalty.L,CovOpt=CovOpt,
                  standardize=standardize,control=control)
      fit.list[[i]] = fit
      BIC.vec[i] = fit$bic
      AIC.vec[i] = fit$aic
      BBIC.vec[i] = fit$bbic
      EBIC.vec[i] = fit$ebic
      
    }
    
    
    min.BIC = min(BIC.vec)
    min.AIC = min(AIC.vec)
    min.BBIC = min(BBIC.vec)
    min.EBIC = min(EBIC.vec)
    
    best.model = which.min(BIC.vec)
    
    best.lam1 = lam1.seq[which.min(BIC.vec)]
    best.fit = fit.list[[which.min(BIC.vec)]]
    
    out = list(lam1.seq=lam1.seq,
               BIC.lam1=BIC.vec,AIC.lam1=AIC.vec,BBIC.lam1=BBIC.vec,EBIC.lam1=EBIC.vec,
               min.BIC=min.BIC,min.AIC=min.AIC,min.BBIC=min.BBIC,min.EBIC=min.EBIC,
               best.model=best.model,best.fit=best.fit,min.lam1=best.lam1,lam1.tuning=lam1.tuning,lam2.tuning=lam2.tuning)
    
    out
    structure(out,class="splmm.tuning")
  }else if((length(lam1.seq)==1)&(length(lam2.seq)>1)){
    print("single lambda1 detected, tuning for lambda2.")
    fit.list = list()
    BIC.vec = vector()
    AIC.vec = vector()
    BBIC.vec = vector()
    EBIC.vec = vector()
    lam1.tuning = FALSE
    lam2.tuning = TRUE
    
    for (i in 1:length(lam2.seq)) {
      fit = splmm(x=x,y=y,z=z,grp=grp,lam1=lam1.seq,lam2=lam2.seq[i],nonpen.b=nonpen.b,nonpen.L=nonpen.L,penalty.b=penalty.b,penalty.L=penalty.L,CovOpt=CovOpt,
                  standardize=standardize,control=control)
      fit.list[[i]] = fit
      BIC.vec[i] = fit$bic
      AIC.vec[i] = fit$aic
      BBIC.vec[i] = fit$bbic
      EBIC.vec[i] = fit$ebic
      
    }
    
    
    min.BIC = min(BIC.vec)
    min.AIC = min(AIC.vec)
    min.BBIC = min(BBIC.vec)
    min.EBIC = min(EBIC.vec)
    
    best.model = which.min(BIC.vec)
    
    best.lam2 = lam2.seq[which.min(BIC.vec)]
    best.fit = fit.list[[which.min(BIC.vec)]]
    
    out = list(lam2.seq=lam2.seq,
               BIC.lam2=BIC.vec,AIC.lam2=AIC.vec,BBIC.lam2=BBIC.vec,EBIC.lam2=EBIC.vec,
               min.BIC=min.BIC,min.AIC=min.AIC,min.BBIC=min.BBIC,min.EBIC=min.EBIC,
               best.model=best.model,best.fit=best.fit,min.lam2=best.lam2,lam1.tuning=lam1.tuning,lam2.tuning=lam2.tuning)
    
    out
    structure(out,class="splmm.tuning")
    
  }else if((length(lam1.seq)>1)&(length(lam2.seq)>1)){
    print("tuning for lambda1 and lambda2.")
    fit.list = list()
    BIC.vec = vector()
    AIC.vec = vector()
    BBIC.vec = vector()
    EBIC.vec = vector()
    
    lam1.tuning = TRUE
    lam2.tuning = TRUE
    
    fit.BIC = matrix(nrow = length(lam1.seq), ncol = length(lam2.seq))
    fit.AIC = matrix(nrow = length(lam1.seq), ncol = length(lam2.seq))
    fit.BBIC = matrix(nrow = length(lam1.seq), ncol = length(lam2.seq))
    fit.EBIC = matrix(nrow = length(lam1.seq), ncol = length(lam2.seq))
    
    for (i in 1:length(lam1.seq)) {
      for (j in 1:length(lam2.seq)) {
        tryCatch({
          fit = splmm(x=x,y=y,z=z,grp=grp,lam1=lam1.seq[i],lam2=lam2.seq[j],nonpen.b=nonpen.b,nonpen.L=nonpen.L,penalty.b=penalty.b,penalty.L=penalty.L,CovOpt=CovOpt,
                      standardize=standardize,control=control)
          fit.BIC[i,j] = fit$bic
          fit.AIC[i,j] = fit$aic
          fit.BBIC[i,j] = fit$bbic
          fit.EBIC[i,j] = fit$ebic},error=function(e){})
        
      }
    }
    
    
    min.bic = which(fit.BIC==min(fit.BIC,na.rm = TRUE), arr.ind = TRUE)
    lam1 = lam1.seq[min.bic[1]]
    lam2 = lam2.seq[min.bic[2]]
    BIC.lam1 = fit.BIC[,min.bic[2]]
    BIC.lam2 = fit.BIC[min.bic[1],]
    
    min.aic = which(fit.AIC==min(fit.AIC,na.rm = TRUE), arr.ind = TRUE)
    AIC.lam2 = fit.AIC[min.aic[1],]
    AIC.lam1 = fit.AIC[,min.aic[2]]
    
    min.bbic = which(fit.BBIC==min(fit.BBIC,na.rm = TRUE), arr.ind = TRUE)
    BBIC.lam2 = fit.BBIC[min.bbic[1],]
    BBIC.lam1 = fit.BBIC[,min.bbic[2]]
    
    min.ebic = which(fit.EBIC==min(fit.EBIC,na.rm = TRUE), arr.ind = TRUE)
    EBIC.lam2 = fit.EBIC[min.ebic[1],]
    EBIC.lam1 = fit.EBIC[,min.ebic[2]]
    
    best.fit = splmm(x=x,y=y,z=z,grp=grp,lam1=lam1,lam2=lam2,nonpen.b=nonpen.b,nonpen.L=nonpen.L,penalty.b=penalty.b,penalty.L=penalty.L,CovOpt=CovOpt,
                     standardize=standardize,control=control)
    
    out = list(lam1.seq=lam1.seq,lam2.seq=lam2.seq,
               BIC.lam1=BIC.lam1,BIC.lam2=BIC.lam2,fit.BIC=fit.BIC,
               AIC.lam1=AIC.lam1,AIC.lam2=AIC.lam2,fit.AIC=fit.AIC,
               BBIC.lam1=BBIC.lam1,BBIC.lam2=BBIC.lam2,fit.BBIC=fit.BBIC,
               EBIC.lam1=EBIC.lam1,EBIC.lam2=EBIC.lam2,fit.EBIC=fit.EBIC,
               best.fit=best.fit,min.lam1=lam1,min.lam2=lam2,lam1.tuning=lam1.tuning,lam2.tuning=lam2.tuning)
    out
    structure(out,class="splmm.tuning")

  }
  
  
  
}