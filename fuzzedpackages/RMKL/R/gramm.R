#' Gramm Matrices
#'
#' This function creates a single gramm matrix for traning set based upon several types of kernels.
#' @param x Matrix of predictors 
#' @param kernel Type of kernel used to compute a gramm matrix
#' @param sigma Hyperparameters for radial kernels
#' @param degree,scale Hyperparameter for polynomial kernel
#' @return Gramm matrix
#' @export 
#' @examples
#' library(kernlab)
#' data(benchmark.data)
#' example.data=benchmark.data[[1]]
#' #Generate linear kernel matrix
#' gramm(example.data[,1:2],'linear',0,0,0)
#' #Generate radial kernel matrices with different values for the hyperparameter.
#' gramm(example.data[,1:2],'radial',2^seq(-3:0),0,0)


gramm=function(x,kernel,sigma,degree,scale){
  {if (kernel=="linear") {
    K=kernlab::kernelMatrix(kernlab::vanilladot(),as.matrix(x))
  } else if (kernel=="polynomial") {
    K=kernlab::kernelMatrix(kernlab::polydot(degree = degree, offset=0, scale = 1),as.matrix(x))
  } else if (kernel=="radial") {
    K=kernlab::kernelMatrix(kernlab::rbfdot(sigma = sigma),as.matrix(x))
  } else if (kernel=="sigmoid") {
    K=kernlab::kernelMatrix(kernlab::tanhdot(scale=scale,offset=1),as.matrix(x))
  }
  else if (kernel=='clinical_cat'){
    int=lapply(1:dim(x)[2], function(a) outer(x[,a],x[,a],'=='))
    clinical_mat=Reduce('+',int)/length(int)
    K=kernlab::as.kernelMatrix(clinical_mat)
  }
  else if (kernel=='clinical_cont'){
    int=lapply(1:dim(x)[2], function(a){
      range=max(x[,a])-min(x[,a])
      diff=(range-abs(outer(x[,a],x[,a],'-')))/range
      return(diff)
    })
    clinical_mat=Reduce('+',int)/length(int)
    K=kernlab::as.kernelMatrix(clinical_mat)
  }}
return(K)
}
