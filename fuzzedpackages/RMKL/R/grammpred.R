#' Gramm Matrix for Test dataset
#'
#' This function creates gramm matrix for test dataset based upon several types of kernel.
#' @param xtrain Matrix of predictors for the training set 
#' @param xtest Matrix of predictors for the test set 
#' @param kernel Type of kernel used to compute a gramm matrix
#' @param sigma Hyperparameters for radial kernels
#' @param degree,scale Hyperparameter for polynomial kernel
#' @return Gramm matrix for test set
#' @export
#' @examples
#' library(kernlab)
#' data(benchmark.data)
#' example.data=benchmark.data[[1]]
#' #Create split between training samples and test samples
#' training.samples=sample(1:dim(example.data)[1],floor(0.7*dim(example.data)[1]),replace=FALSE)
#' xtrain=example.data[training.samples,1:2]
#' xtest=example.data[-training.samples,1:2]
#' #Generate linear kernel
#' grammpred(xtrain,xtest,'linear',0,0,0)
#' #Generate radial kernels with different values for the hyperparameter.
#' grammpred(xtrain,xtest,'radial',2^seq(-3:0),0,0)

grammpred=function(xtrain,xtest,kernel,sigma,degree,scale) {
  {if (kernel=="linear") {
    K=kernlab::kernelMatrix(kernlab::vanilladot(),as.matrix(xtest),as.matrix(xtrain))
  } else if (kernel=="polynomial") {
    K=kernlab::kernelMatrix(kernlab::polydot(degree = degree,scale=1, offset=0),as.matrix(xtest),as.matrix(xtrain))
  } else if (kernel=="radial") {
    K=kernlab::kernelMatrix(kernlab::rbfdot(sigma = sigma),as.matrix(xtest),as.matrix(xtrain))
  } else if (kernel=="sigmoid") {
    K=kernlab::kernelMatrix(kernlab::tanhdot(scale=scale,offset=1),as.matrix(xtest),as.matrix(xtrain))
  }
  else if (kernel=='clinical_cat'){
    int=lapply(1:dim(xtrain)[2], function(a) outer(xtest[,a],xtrain[,a],'=='))
    clinical_mat=Reduce('+',int)/length(int)
    K=kernlab::as.kernelMatrix(clinical_mat)
  }
  else if (kernel=='clinical_cont'){
    int=lapply(1:dim(xtrain)[2], function(a){
      range=max(xtrain[,a])-min(xtrain[,a])
      diff=(range-abs(outer(xtest[,a],xtrain[,a],'-')))/range
      return(diff)
    })
    clinical_mat=Reduce('+',int)/length(int)
    K=kernlab::as.kernelMatrix(clinical_mat)
  }}
   return(K)
        }
    
