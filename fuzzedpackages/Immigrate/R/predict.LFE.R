#' predict.LFE
#'
#' This function performs predition for LFE(Local Feature Extraction) algorithm.
#' @param object weights obtained from LFE
#' @param xx model matrix of explanatory variables
#' @param yy label vector 
#' @param newx new model matrix to be predicted 
#' @param ... further arguments passed to or from other methods
#' @keywords predict the label of new data based on LFE
#' @return predicted labels for new data (newx)
#' @importFrom stats predict

#' @export
#' @examples
#' data(park)
#' xx<-park$xx
#' yy<-park$yy
#' w<-LFE(xx,yy)
#' pred<-predict(w,xx,yy,xx)
#' print(pred)
#' @references Sun Y, Wu D. A relief based feature extraction algorithm[C]//Proceedings of the 2008 SIAM International Conference on Data Mining. Society for Industrial and Applied Mathematics, 2008: 188-195.

predict.LFE<-function(object,xx,yy,newx,...){
  w<-object
  yy<-as.numeric(yy)
  if (ncol(xx) != ncol(newx)){
    stop("xx and newx have different lengths of explanatory variables")
  }
  pred<-sapply(c(1:nrow(newx)), function(i){
    yy[as.numeric(which.min(rowSums(crossprod(abs(t(xx)-as.numeric(newx[i,])),w)*t(abs(t(xx)-as.numeric(newx[i,])))))[1])]
  })
  return(pred)
}
