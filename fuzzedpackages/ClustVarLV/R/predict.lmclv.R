#' @title prediction for lmCLV models. 
#' 
#' @description To get the predicted response values based on lmCLV model. 
#' 
#' @param object : result of class \code{lmclv}
#' @param newdata : Data frame of observations for which to make predictions
#' @param shrinkp : shrinkage  parameter.
#' @param ... : further arguments passed to or from other methods
#'  
#' @return a matrix of the predicted values, \cr
#'         each column with an increasing number of CLV component included \cr 
#'         the first column being for the null model \cr
#'         if the response if a binary factor, two additional matrices are provided :\cr
#'         the probabilities of belonging to class 1 and the response values (0 or 1).
#'         
#'   
#' @export
  
predict.lmclv <-    function(object,newdata,shrinkp,...) 
{
  reslmclv<-object
  if (!inherits(reslmclv, "lmclv")) 
      stop("non convenient objects")
      
  
  Nmod=length(reslmclv)
  cas=NULL
  for (mod in 1:Nmod) {
   if (reslmclv[[mod]]$shrinkp==shrinkp){
     cas=mod
     break
   }
  }
  if (is.null(cas))    stop(paste("lm_CLV result not performed with this shrinkage parameter of ",shrinkp," => run lm_CLV with this value"))
  
  res=reslmclv[[cas]]  
  yisfact=FALSE
  if (!is.null(res$yfact)) {
    yisfact=TRUE
    yref=as.numeric(res$yfact)-1
    tabresp=NULL
    tabprob=NULL
    n1=sum(yref==1)
    n0=sum(yref==0)
  }
  
  mX=apply(res$X,2,mean)
  if (res$sX==TRUE) {
    sX=apply(res$X,2,sd)
  } else {
    sX=rep(1,ncol(res$X))
  }
  Xs<-scale(newdata,center=mX,scale=sX)
  Xrefs<-scale(res$X,center=mX,scale=sX)
  # step "cst model"
  ypred=rep(res$Cst,nrow(newdata))
  ypredref=rep(res$Cst,nrow(res$X))
  tabpred=ypred
  if (yisfact) {
    m1=mean(ypredref[yref==1]);  s1=sd(ypredref[yref==1]);
    m0=mean(ypredref[yref==0]);  s0=sd(ypredref[yref==0]);     
    prob=Bayes_classif(ypred,n1,m1,s1,n0,m0,s0)$pclass
    tabprob=prob
    resp=Bayes_classif(ypred,n1,m1,s1,n0,m0,s0)$yclass
    tabresp=resp
  }
 
  # step by step
  nblv=length(res$Group)
  for (lv in 1:nblv) {
    ypred=ypred+Xs%*%res$Beta[[lv]]
    ypredref=ypredref+Xrefs%*%res$Beta[[lv]]
    tabpred=cbind(tabpred,ypred)
    if (yisfact) {
      m1=mean(ypredref[yref==1]);  s1=sd(ypredref[yref==1]);
      m0=mean(ypredref[yref==0]);  s0=sd(ypredref[yref==0]);     
      prob=Bayes_classif(ypred,n1,m1,s1,n0,m0,s0)$pclass
      tabprob=cbind(tabprob,prob)
      resp=Bayes_classif(ypred,n1,m1,s1,n0,m0,s0)$yclass
      tabresp=cbind(tabresp,resp)
    }
  }

  if (!yisfact) listres=list(predictions=tabpred)
  if (yisfact) listres=list(predictions=tabpred,probabilities=tabprob,responses=tabresp)
  return(listres)  
}        
   
 
Bayes_classif=function(ypred,n1,m1,s1,n0,m0,s0) {
  # Bayes classification rule
  # ypred : the predicted values
  # n1,m1,s1,n0,m0,s0 : the nb of obs, mean and sd when y=1 and y=0
  n=n1+n0
  f1=n1/n
  f0=n0/n
  s=sqrt(  (((n1-1)*s1^2)+ ((n0-1)*s0^2)) / (n-2))
  
  if ( (s1>0) & (s0>0) ) {
    a1=f1*dnorm(ypred,mean=m1,sd=s1)
    a0=f0*dnorm(ypred,mean=m0,sd=s0)
  } else {
    a1=rep(f1,length(ypred))
    a0=rep(f0,length(ypred))
  }
  pclass=a1/(a0+a1)
  yclass=apply(cbind(a0,a1),1,which.max)-1
  return(list(yclass=yclass,pclass=pclass))
}