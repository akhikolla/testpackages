#' Cross validation for SVM
#'
#' This function performs cross validation to find best combinatiuon of hyper perameters and cost
#' and uses this model to provide prediction performace results. 
#' @param data List of data matrices for each pathways for each pathway    
#' @param outcome Binary outcome variable for MKL 
#' @param train.samples Vector of indices that will be used as training samples
#' @param kernels vector of kernel types
#' @param degree  Degree of polynomial kernel matrix
#' @param scale Leading coefficient on the polynomial kernel
#' @param sigma Hyperparameter for the radial basis kernel
#' @param C Vector of canduidate cost parameters
#' @param nfold Number of folds used in cross validation
#' @return cm Confustion matrix along with a variety of accuracy statistics 
#' @return best.model Model that had the highest accuracy with nfold cross validation
#' @export
#' @examples
#' data(benchmark.data)
#' example.data=benchmark.data[[1]]
#' training.samples=sample(1:nrow(example.data), floor(0.7*nrow(example.data)),replace=FALSE)
#' C=2^c(-1,1)
#' #Find the best cost parameter within the provided range if a linear kernel is used
#' SVM.nfoldcv(example.data[,1:2], as.factor(example.data[,3]),training.samples,C,
#'             'linear',0,0,0,nfold=3)


SVM.nfoldcv=function(data,outcome,train.samples,C,kernels,degree,scale,sigma,nfold=10){
  data=data.frame(data,outcome)
  train.data=data[train.samples,]
  test.data=data[-train.samples,]
  if(kernels=='linear'){
  cv=caret::train(outcome~.,data=train.data, method='svmLinear',trControl=caret::trainControl(method='cv', number=nfold),
                  tuneGrid=expand.grid('C'=C), tuneLength=10)
  cm=caret::confusionMatrix(kernlab::predict(cv$finalModel,test.data[,-dim(test.data)[2]]),test.data$outcome)
  }
  else if(kernels=='polynomial'){
  cv=caret::train(outcome~.,data=train.data, method='svmPoly',trControl=caret::trainControl(method='cv', number=nfold),
                  tuneGrid=expand.grid('C'=C,'degree'=degree,'scale'=scale), tuneLength=10)
  cm=caret::confusionMatrix(kernlab::predict(cv$finalModel,test.data[,-dim(test.data)[2]]),test.data$outcome)
  }
  else if(kernels=='radial'){
  cv=caret::train(outcome~.,data=train.data, method='svmRadial',trControl=caret::trainControl(method='cv', number=nfold),
                  tuneGrid=expand.grid('C'=C,'sigma'=sigma),tuneLength=10)
  cm=caret::confusionMatrix(kernlab::predict(cv$finalModel,test.data[,-dim(test.data)[2]]),test.data$outcome)
  }
  return(list(cm=cm, best.model=cv$finalModel))
}
