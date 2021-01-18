#' Generate both training and test kernel matrices
#'
#' This function creates gramm matrix for traning set baed upon several types of kernel
#' and specified hyper paremeters. This function is essentially a wrappper functions that combines gramm and grammpred.
#' Additionally this function divides each kernel matrix by it's trace, which is a common transformation used in MKL.
#' @param data List of data matrices
#' @param train.samples Vector of indices that will be used as training samples
#' @param kernels Character vector of  kernel types
#' @param degree  Degree of polynomial kernel matrix
#' @param scale Leading coefficient on the polynomial kernel
#' @param sigma Hyperparameter for the radial basis kernel
#' @return K.train Gramm matricesfor training data
#' @return K.test Gramm matrices for test data
#' @export
#' @examples
#' library(kernlab)
#' data(benchmark.data)
#' example.data=benchmark.data[[1]]
#' #Dividing the samples into a train set and test set.
#' training.samples=sample(1:dim(example.data)[1],floor(0.7*dim(example.data)[1]),replace=FALSE)
#' #Specifying the type and hyperparameters for each kernel.
#' kernels=c('linear',rep('radial',3))
#' degree=rep(0,4)
#' scale=rep(0,4)
#' sigma=c(0,2^seq(-3:0))
#' kernels.gen(example.data[,1:2], training.samples, kernels, degree, scale, sigma)


kernels.gen=function(data,train.samples,kernels,degree,scale,sigma){
  unit.trace=function(K){K/sum(diag(K))}
  trace.kernel=function(K){sum(diag(K))}
  train.data=matrix(data[train.samples,],ncol=ncol(data))
  test.data=matrix(data[-train.samples,],ncol=ncol(data))
  if(length(kernels)==1){
  K.train=RMKL::gramm(x=train.data,degree=degree,scale=scale,sigma=sigma, kernel=kernels)
  trace.train=trace.kernel(K.train)
  K.train=unit.trace(K.train)
  K.test=RMKL::grammpred(xtrain=train.data,xtest=test.data,scale=scale,degree=degree,sigma=sigma,kernel=kernels)
  K.test=K.test/trace.train
 }else{
  K.train=lapply(1:length(kernels),function(i){
                  RMKL:: gramm(x=train.data,degree=degree[i],scale=scale[i],sigma=sigma[i], kernel=kernels[i])})
  trace.train=unlist(lapply(K.train, trace.kernel))
  K.train=lapply(K.train,unit.trace)
  K.test=lapply(1:length(kernels),function(i){
                  RMKL::grammpred(xtrain=train.data,xtest=test.data,scale=scale[i],degree=degree[i],
                  sigma=sigma[i],kernel=kernels[i])})
  K.test=lapply(1:length(K.test), function(a) K.test[[a]]/trace.train[a])
}
  return(list('K.train'=K.train,'K.test'=K.test))
}
