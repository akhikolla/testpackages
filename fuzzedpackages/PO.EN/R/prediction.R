#' PO-EN predicting function
#'
#' A prediction function using the linear predictor of PO-EN fitting results.
#'@param X Input design matrix. Should not include the intercept vector.
#'@param beta A coefficients vector from the PO-EN fitting function.
PO.EN.predict<-function(X,beta){
  X<-as.matrix(X)
  X<-cbind(1,X)
  mu<-X%*%beta
  Y.predict=exp(mu)/(1+exp(mu))
  Y.predict[which(is.nan(Y.predict))]=1

  return(Y.predict)
}
