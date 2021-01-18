#' Simple MKL 
#'
#' This function conducts Simple MKL for precomputed gramm matrices
#' @param k list of Gramm matrices
#' @param outcome vector of binary outcome -1 and 1
#' @param penalty ppenalty of the smoothness of the resulting desicion rules
#' @param tol change between to iterations is smaller than this, algorithms is considered to have converged
#' @param max.iters maximum number of allowed iteratons
#' @return gamma weight vector for the importnace of each kernel 
#' @return alpha coeffiencents of the dual of MKL
#' @return time total amount of time to train model
#' @return max.iters Numvber of iterations to reach convergence criteria
#' @export
#' @examples
#' library(kernlab)
#' library(caret)
#' library(RMKL)
#' #Load data
#' data(benchmark.data)
#' example.data=benchmark.data[[1]]
#' # Split samples into training and test sets 
#' training.samples=sample(1:dim(example.data)[1],floor(0.7*dim(example.data)[1]),replace=FALSE)
#' # Set up cost parameters and kernels 
#' C=100
#' kernels=rep('radial',3)
#' degree=rep(0,3)
#' scale=rep(0,3)
#' sigma=c(0,2^seq(-3:0))
#' K=kernels.gen(example.data[,1:2], training.samples, kernels, degree, scale, sigma)
#' K.train=K$K.train
#' SimpleMKL.classification(K.train,example.data[training.samples,3], C)

SimpleMKL.classification=function(k,outcome,penalty,tol=10^(-4),max.iters=1000){
  inner=function(Ddag,gammadag,k,J,penalty){
    Jdag=0
    while(Jdag<J&sum(Ddag<0)>0){
      mu=min(which(gammadag==max(gammadag)))
      D=Ddag/sum(abs(Ddag))
      gamma=gammadag
      neg=which(D<0&gamma>0)
      v=neg[which(-gamma[neg]/D[neg]==min(-gamma[neg]/D[neg],na.rm=TRUE))]
      stepsize=-gamma[v]/D[v]
      gammadag=gamma+stepsize*D
      Ddag[mu]=D[mu]-D[v]
      Ddag[v]=0
      kk=Reduce('+',mapply("*", k, gammadag,SIMPLIFY = FALSE))
      h=kk*(outcome%*%t(outcome))
      model=kernlab::ipop(rep(-1,length(outcome)),h,t(outcome),0,rep(0,length(outcome)),rep(penalty,length(outcome)),0)
      Jdag=-1/2*sum(h*kernlab::primal(model)%*%t(kernlab::primal(model)))+sum(kernlab::primal(model))
    }
    return(list('gamma'=gamma,'objval'=Jdag,'direction'=D,'stepsize'=stepsize))
  }
  
  gamma_all=list()
  gamma=rep(1/length(k),length(k))
  epsilon=1
  iters=0
  while(epsilon>tol&iters<max.iters&sum(gamma==1)==0){
   # tic()
    iters=iters+1
    gamma_all[[iters]]=gamma
    kk=Reduce('+',mapply("*", k, gamma,SIMPLIFY = FALSE))
    h=kk*(outcome%*%t(outcome))
    model=kernlab::ipop(rep(-1,length(outcome)),h,t(outcome),0,rep(0,length(outcome)),rep(penalty,length(outcome)),0)
    J=-1/2*sum(h*kernlab::primal(model)%*%t(kernlab::primal(model)))+sum(kernlab::primal(model))
    dJ=sapply(1:length(k), function(a) -1/2*sum(k[[a]]*outcome%*%t(outcome)*kernlab::primal(model)%*%t(kernlab::primal(model))))
    mu=min(which(gamma==max(gamma)))
    gradJ=rep(0,length(k))
    gradJ[-mu]=dJ[-mu]-dJ[mu]
    gradJ[mu]=-sum(gradJ)
    cond1=which(gamma==0&gradJ>0)
    cond2=setdiff(which(gamma>0),mu)
    cond3=mu
    D=rep(0,length(k))
    D[cond1]=0
    D[cond2]=-gradJ[cond2]
    D[cond3]=sum(gradJ[cond2])
    Ddag=D
    innersol=inner(Ddag=D,gammadag = gamma,k=k,J,penalty)
    obj=function(gamma,outcome,penalty,direction,stepsize){
      gammanew=gamma+direction*stepsize
      kk=Reduce('+',mapply("*", k, gammanew,SIMPLIFY = FALSE))
      h=kk*(outcome%*%t(outcome))
      model=kernlab::ipop(rep(-1,length(outcome)),h,t(outcome),0,rep(0,length(outcome)),rep(penalty,length(outcome)),0)
      J=-1/2*sum(h*kernlab::primal(model)%*%t(kernlab::primal(model)))+sum(kernlab::primal(model))
      return(J)
    }
    step_opt=stats::optim(innersol$stepsize/2,obj,outcome=outcome,penalty=penalty,direction=innersol$direction,gamma=innersol$gamma,lower=0,upper=innersol$stepsize,method='L-BFGS-B')
    epsilon=max(gamma-innersol$gamma+step_opt$par*innersol$direction)
    gamma=innersol$gamma+step_opt$par*innersol$direction
    gamma=gamma/sum(gamma)
    #toc(log = TRUE, quiet = TRUE)
  }
  gamma_all[[iters+1]]=gamma
 # log.lst <- tic.log(format = FALSE)
  j=match(kernlab::primal(model)[(kernlab::primal(model)>0)&(kernlab::primal(model)<penalty)][1],kernlab::primal(model))
  b=outcome[j]-sum(kernlab::primal(model)*outcome*kk[,j])
  #tic.clearlog()
  return(list('gamma'=gamma,'iters'=iters,'alpha'=kernlab::primal(model),'b'=b,'gamma_all'=gamma_all))
}
