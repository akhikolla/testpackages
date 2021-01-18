#' Prediction from MKL model
#'
#' This function creates gramm matrix for traning set baed upon several types of kernel
#' and specified hyper paremeters. Matrix corresponds to similarity betwween 
#' each sample in the training set.
#' @param ktest Gramm matrix of training data and test data
#' @param model MKL model 
#' @param train.outcome Outcome for the training data
#' @return yhat Predicted value for each test point
#' @return predicted Sign of yhat, which is the final predicted outcome 
#' @export
#' @examples
#' library(kernlab)
#' library(caret)
#' data(benchmark.data)
#' example.data=benchmark.data[[1]]
#' training.samples=sample(1:dim(example.data)[1],floor(0.7*dim(example.data)[1]),replace=FALSE)
#' C=100
#' kernels=rep('radial',3)
#' degree=rep(0,3)
#' scale=rep(0,3)
#' sigma=c(0,2^seq(-3:0))
#' K=kernels.gen(example.data[,1:2], training.samples, kernels, degree, scale, sigma)
#' K.train=K$K.train
#' K.test=K$K.test
#' SEMKL.model=SEMKL.classification(K.train,example.data[training.samples,3], C)
#' predicted=prediction.Classification(SEMKL.model, K.test, example.data[training.samples,3])
#' confusionMatrix(factor(predicted$predict, levels=c(-1,1)),
#'                 factor(example.data[-training.samples,3],levels=c(-1,1)))

prediction.Classification=function(model,ktest,train.outcome) {
  N=length(train.outcome)
  n=dim(ktest[[1]])[1]
  m=length(ktest)
  product=list()
  fushion=matrix(,n,N) #rep(0,n)
  yhat=rep(0,N)
  predict=rep(0,N)
  product=lapply(1:length(ktest), function(a) ktest[[a]]*model$gamma[a])
  fushion=Reduce('+',product)
  yhat=fushion%*%(model$alpha*train.outcome)+model$b
  result=list("yhat"=yhat,"predict"=sign(yhat))
  return(result)
}
