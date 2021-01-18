

####################################################
#####  Regularized CRM                         #####
#####  Adaptive L2 Penalthy                    #####
#####  2P - logistic model                     #####
#####  Algorithm: one-step coordinate descent  #####
####################################################


rCRM = function(x, y, dose0, tp=0.3, family="2P", mlambda=1, nlambda=50, rlambda=NULL, wldose=NULL, nfolds=length(y), foldid=NULL, keep.beta=FALSE, thresh=1e-7, maxit=1e+4, threshP=1e-6, threshB=100) {

  family=match.arg(family)
  if (!all(x %in% dose0)) {stop("dose0 should include all the dose levels in x")}


  fit=switch(family,
             "2P"=P2L2(x,y,dose0,tp,mlambda,nlambda,rlambda,wldose,nfolds,foldid,keep.beta,thresh,maxit,threshP,threshB))
  fit$family=family
  fit$tp=tp
  fit$dose0=dose0

  class(fit)="rCRM"
  return(fit)
}


