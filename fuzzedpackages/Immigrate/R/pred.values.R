
#' pred.values
#'
#' This function performs some statistical value prediction
#' @param y_train label vector for training data
#' @param y_test label vector for test data 
#' @param pred_train predicted probabilities for training data
#' @param pred_test predicted probabilities for test data
#' @keywords predict
#' @return \item{AUC_train}{AUC for training data}
#' @return \item{AUC_test}{ AUC for test data}
#' @return \item{accuracy_test}{ accuracy for test data}
#' @return \item{precision_test}{ precision for test data}
#' @return \item{recall_test}{ recall for test data}
#' @return \item{F1_test}{ F1 score for test data}
#' @return \item{thre}{ threshold to separate two labels, obtained from training data}
#' @importFrom pROC roc
#' 
#' @export
#' @examples
#' y_train<-c(0,1,0,1,0,1)
#' y_test<-c(0,1,0,1)
#' pred_train<-c(0.77,0.89,0.32,0.96,0.10,0.67)
#' pred_test<-c(0.68,0.75,0.50,0.81)
#' re<-pred.values(y_train,y_test,pred_train,pred_test)
#' print(re)
#' 
pred.values<-function(y_train,y_test,pred_train,pred_test){
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Package \"pROC\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (length(unique(y_train)) > 2){
    stop("training data with more than 2 levels")
  }else if (length(unique(y_train)) < 2){
    stop("training data with less than 2 levels")
  }else if (length(unique(y_test)) > 2){
    stop("test data with more than 2 levels")
  }else if(all(sort(unique(y_train))!= sort(unique(y_test)) )){
    stop("training data and test data have different labels")
  }
  
  if ((length(y_train) != length(pred_train))|(length(y_test) != length(pred_test))){
    stop("label and prediction have different lengths")
  }else if ((max(pred_train)>1)|(min(pred_train)<0)|(max(pred_test)>1)|(min(pred_test)<0)){
    stop("Values of prediction out of boundary")
  }
  
  # reshape the label
  y_train[which(y_train == sort(unique(y_train))[1])]<-0
  y_train[which(y_train == sort(unique(y_train))[2])]<-1
  y_test[which(y_test == sort(unique(y_test))[1])]<-0
  y_test[which(y_test == sort(unique(y_test))[2])]<-1
  
  # roc
  score_train<-roc(y_train,pred_train)
  if(length(unique(y_test)) == 2){
    score_test<-roc(y_test,pred_test)
  }
  # threshold 
  score_train$diff<-score_train$sensitivities-(1- score_train$specificities)
  thre<-score_train$thresholds[which.max(score_train$diff)]
  
  ## value for test data
  TC<-sum(as.numeric(pred_test>thre) == y_test)
  TP<-sum((as.numeric(pred_test>thre) == y_test)& (y_test == 1))
  TN<-sum((as.numeric(pred_test>thre) == y_test)& (y_test == 0))
  FP<-sum(y_test==0) - TN
  FN<-sum(y_test==1) - TP
  tc<-TC/max(length(y_test),1)
  tp<-TP/max(sum(y_test==1),1)
  tn<-TN/max(sum(y_test==0),1)
  if((TP+FP)!=0){
    precision<-TP/(TP+FP)
  }else{
    precision<-0
  }
  if((TP+FN)!=0){
    recall<-TP/(TP+FN)
  }else{
    recall<-0
  }
  if((precision!=0)|(recall==0) ){
    F1<-2/(1/precision+1/recall)
  }else{
    F1<-0
  }
  if(length(unique(y_test)) == 1){
    newList<-list("AUC_train"=round(score_train$auc,4),
                  "accuracy_test"=tc,"precision_test"=precision,
                  "recall_test"=recall,"F1_test"=F1,"thre" = thre)
  }else{
    newList<-list("AUC_train"=round(score_train$auc,4),
                  "AUC_test"=round(score_test$auc,4),
                  "accuracy_test"=tc,"precision_test"=precision,
                  "recall_test"=recall,"F1_test"=F1,"thre" = thre) 
  }
  return(newList)
}