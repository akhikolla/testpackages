#' Simple and Efficient MKL 
#'
#' This function conducts Simple and Efficnent MKL for precomputed gramm matrices
#' @param k list of Gramm matrices
#' @param outcome vector of binary outcome -1 and 1
#' @param penalty penalty of the smoothness of the resulting desicion rules
#' @param tol change between to iterations is smaller than this, algorithms is considered to have converged
#' @param max.iters maximum number of allowed iteratons
#' @return gamma weight vector for the importnace of each kernel 
#' @return alpha coeffiencents of the dual of MKL
#' @return time total amount of time to train model
#' @return iters Numvber of iterations to reach convergence criteria
#' @return gamma_all Kernel weights for each interation of SEMKL
#' @export
#' @examples
#' data(benchmark.data)
#' example.data=benchmark.data[[1]]
#' #Load data
#' training.samples=sample(1:dim(example.data)[1],floor(0.7*dim(example.data)[1]),replace=FALSE)
#' # Split samples into training and test sets 
#' C=1
#' kernels=c('radial','polynomial')
#' degree=c(0,2)
#' scale=c(0,2)
#' sigma=c(2,0)
#' K=kernels.gen(example.data[,1:2], training.samples, kernels, degree, scale, sigma)
#' K.train=K$K.train
#' SEMKL.classification(K.train,example.data[training.samples,3], C)

SEMKL.classification=function(k,outcome,penalty,tol=0.0001,max.iters=1000){
  delta=rep(1,length(k))
  iters=0
  n=length(outcome)
  m=length(k)
  gamma=rep(1/length(k),length(k))
  gamma_all=list()
  #tic()
  while (max(delta)>tol && iters<max.iters){
    iters=iters+1
    gamma_all[[iters]]=gamma
    kg=lapply(1:length(k), function(a) k[[a]]*gamma[a])
    kk=Reduce('+',kg)
    h=kk*(outcome%*%t(outcome))
    model=kernlab::ipop(rep(-1,n),h,t(outcome),0,rep(0,n),rep(penalty,n),0)
    alpha=kernlab::primal(model)
    fnorm=sapply(1:length(k), function(a) sqrt(gamma[a]^2*(alpha*outcome)%*%k[[a]]%*%(alpha*outcome)))
    temp=gamma
    gamma=fnorm/sum(fnorm)
    delta=abs(temp-gamma)
  }
  #toc(log=TRUE,quiet=TRUE)
  #time=tic.log(format=FALSE)[[1]]$toc-tic.log(format=FALSE)[[1]]$tic
  j=match(alpha[(alpha>0)&(alpha<penalty*0.9999)][1],alpha)
  b=outcome[j]-sum(alpha*outcome*kk[,j])
  gamma_all[[iters+1]]=gamma
  results=list("alpha"=alpha,"b"=b,"gamma"=temp,"iters"=iters,'gamma_all'=gamma_all)
  return(results)
}
