#' Converts Cost from DALMKL to SEMKL or SimpleMKL
#'
#' This function estimates an a comparable cost for SEMKL or SimpleMKL  from DALMKL.
#' @param K.train Gramm matrix of training data 
#' @param DALMKL.model DAL MKL model 
#' @param C.DALMKL Cost used in DAMKL model
#' @return C cost SEMKL or SimpleMKL 
#' @export
#' @examples
#' data(benchmark.data)
#' data.mkl=benchmark.data[[1]]
#' kernels=rep('radial',2)
#' sigma=c(2,1/20)
#' train.samples=sample(1:nrow(data.mkl),floor(0.7*nrow(data.mkl)),replace=FALSE)
#' degree=sapply(1:length(kernels), function(a) ifelse(kernels[a]=='p',2,0))
#' #Kernels.gen splts the data into a training and test set, and generates the desired kernel matrices.
#' #Here we generate two gaussisan kernel matrices with sigma hyperparameter 2 and 0.05
#' K=kernels.gen(data=data.mkl[,1:2],train.samples=train.samples,kernels=kernels,sigma=sigma,
#' degree=degree,scale=rep(0,length(kernels)))
#' C=0.05 #Cost parameter for DALMKL
#' K.train=K$K.train
#' K.test=K$K.test
#'
#' # parameters set up
#' ytr=data.mkl[train.samples,3]
#' #Converts list of kernel matrices in to an array with is appropriate for C++ code
#' k.train=simplify2array(K.train) 
#' k.test=simplify2array(K.test)

#'  #Implement DALMKL with the hinge loss function
#'  spicy_svmb1n=SpicyMKL(K=k.train,y=ytr, loss='hinge',C=C)
#'  #Convert cost from DALMKL to be more compatible withSimpleMKL
#'  C.SimpleMKL=C.convert(K.train,spicy_svmb1n,C)



C.convert=function(K.train,DALMKL.model,C.DALMKL){
norm=function(x,y){y%*%x%*%y}
sum.norms=sum(sapply(1:length(K.train), function(a){
                 norm(K.train[[a]],DALMKL.model$alpha[,a])}))
            C=C.DALMKL*sum.norms
            return(C)
}

