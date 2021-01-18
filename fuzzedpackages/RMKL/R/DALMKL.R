#' DALMKL 
#'
#' This function conducts DALMKL for precomputed gramm matrices
#' @param K The multiple kernel cube (3-d array)
#' @param y The outome variable, must be -1/1
#' @param loss The loss function to be used, must be either 'hinge' or 'logistic', default to be 'hinge'
#' @param C tuning parameter for block one norm, default to be .5
#' @param tolOuter change between to iterations is smaller than this, algorithms is considered to have converged for outer loop, default to be .01
#' @param tolInner change between to iterations is smaller than this, algorithms is considered to have converged for inner loop, default to be .000001
#' @param OuterMaxiter maximum number of allowed iteratons for outer loop, default to be 500
#' @param InnerMaxiter maximum number of allowed iteratons for inner loop, default to be 500
#' @param calpha Lagrangian parameter, default to be 10
#' @return b Estimated Intercept 
#' @return alpha coeffiencents of the dual of MKL
#' @return weight Estimated between kernel weight
#' @return rho Estimated within kernel weight
#' @useDynLib RMKL, .registration=TRUE
#' @importFrom Rcpp evalCpp 
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

#' # parameters set up
#'  ytr=data.mkl[train.samples,3]
#' #Converts list of kernel matrices in to an array with is appropriate for C++ code
#'  k.train=simplify2array(K.train)  
#'  k.test=simplify2array(K.test)

#' #Implement DALMKL with the hinge loss function
#'  spicy_svmb1n=SpicyMKL(K=k.train,y=ytr, loss='hinge',C=C)
#'  #Implement DALMKL with the hinge loss function
#'  spicy_logistic=SpicyMKL(K=k.train,y=ytr, loss='logistic',C=C)#' 

SpicyMKL <- function(K, y, loss = 'hinge', C = .5, tolOuter = .01, tolInner = .000001, OuterMaxiter = 500, InnerMaxiter = 500, calpha = 10) {
  if (loss == 'hinge') {
    SpicySVM(K, y, C, tolOuter, tolInner, OuterMaxiter, InnerMaxiter, calpha)
  } else {
    SpicyLogit(K, y, C, tolOuter, tolInner, OuterMaxiter, InnerMaxiter, calpha)
  }
}

#' Predict SpicyMKL
#' 
#' This function is used to predict SpicyMKL models
#' @param alpha coefficient
#' @param b intercept
#' @param k0 the kernel cube needs prediction
#' @return The predicted score
#' @useDynLib RMKL
#' @importFrom Rcpp evalCpp
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
#' K=kernels.gen(data=data.mkl[,1:2],train.samples=train.samples,
#'                 kernels=kernels,sigma=sigma,degree=degree,scale=rep(0,length(kernels)))
#' C=0.05 #Cost parameter for DALMKL
#' K.train=K$K.train
#' K.test=K$K.test

  # parameters set up
#'  ytr=data.mkl[train.samples,3]
#'  #Converts list of kernel matrices in to an array with is appropriate for C++ code
#'  k.train=simplify2array(K.train) 
#'  k.test=simplify2array(K.test)
#'
#'  #Implement DALMKL with the hinge loss function
#'  spicy_svmb1n=SpicyMKL(K=k.train,y=ytr, loss='hinge',C=C)
#'  prediction_logistic=predict_Spicy(spicy_svmb1n$alpha,spicy_svmb1n$b,k.test)
#'  #Implement DALMKL with the hinge loss function
#'  spicy_logistic=SpicyMKL(K=k.train,y=ytr, loss='logistic',C=C)
#'  prediction_logistic=predict_Spicy(spicy_logistic$alpha,spicy_logistic$b,k.test)

predict_Spicy <- function(alpha, b, k0) {
  predictspicy(alpha, b, k0)
}
