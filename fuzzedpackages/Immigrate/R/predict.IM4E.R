#' predict.IM4E
#'
#' This function performs the predition for IM4E(Iterative Margin-Maximization under Max-Min Entropy) algorithm.
#' @param object weight or result of IM4E algorithm
#' @param xx model matrix of explanatory variables
#' @param yy label vector 
#' @param newx new model matrix to be predicted 
#' @param sig sigma used in algorithm, default to be 1
#' @param type the form of final output, default to be "both". One can also choose "response"(predicted probabilities) or "class"(predicted labels). 
#' @param ... further arguments passed to or from other methods
#' @keywords predict the label of new data based on IM4E
#' @return \item{response}{predicted probabilities for new data (newx)}
#' @return \item{class}{predicted class labels for new data (newx)}
#' @importFrom stats predict

#' @export
#' @examples
#' data(park)
#' xx<-park$xx
#' yy<-park$yy
#' index<-c(1:floor(nrow(xx)*0.3))
#' train_xx<-xx[-index,]
#' test_xx<-xx[index,]
#' train_yy<-yy[-index]
#' test_yy<-yy[index]
#' re<-IM4E(train_xx,train_yy)
#' res<-predict(re,train_xx,train_yy,test_xx,type="class")
#' print(res)
#' @references 
#' Bei Y, Hong P. Maximizing margin quality and quantity[C]//Machine Learning for Signal Processing (MLSP), 2015 IEEE 25th International Workshop on. IEEE, 2015: 1-6.

predict.IM4E<-function(object,xx,yy,newx,sig = 1, type = "both",...){
  TYPES <- c("both","response","class") 
  typeIdx <- pmatch(type, TYPES)
  if (is.na(typeIdx)){
    stop("Invalid type")
  }
  yy<-as.numeric(yy)
  if (ncol(xx) != ncol(newx)){
    stop("xx and newx have different lengths of explanatory variables")
  }
  if (is.list(object)){
    w<-matrix(object$w)  
  }else{
    w<-object
  }
  
  label<-unique(yy)
  
  v<-sapply(c(1:length(label)),function(j){
    sapply(c(1:nrow(newx)), function(i){
      yyy<-abs(yy-label[j])
      y<-rep(0,length(yyy))
      y[which(yyy==0)]<-1
      tmp<-matrix((y)*exp(-(t(w)%*%abs(t(xx)-as.numeric(newx[i,]))/sig) -1 ))
      s<-sum(tmp)
      if (s!=0) tmp<-tmp/s
      (t(w)%*%abs(t(xx)-as.numeric(newx[i,])))%*%tmp
    })
  })
  myfun<-sapply(c(1:nrow(newx)), function(i){
    v[i,]<<-v[i,]/sum(v[i,])
  })
  pred<-sapply(c(1:nrow(newx)),function(i){
    label[which.min(v[i,])]
  })
  
  if (type == "both"){
    newList<-list("class"=pred,"prob"=v)
    return(newList) 
  }else if(type == "response"){
    return(v)
  }else if(type == "class"){
    return(pred)
  }else{
    stop("use wrong type")
  }
}

