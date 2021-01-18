#' LFE
#'
#' This function performs LFE(Local Feature Extraction) algorithm.
#' @param xx model matrix of explanatory variables
#' @param yy label vector
#' @param T number of instance used to update weights, default to be 5
#' @keywords LFE
#' @return \item{w}{new weight matrix after LFE algorithm}
#' @export
#' @examples
#' data(park)
#' xx<-park$xx
#' yy<-park$yy
#' re<-LFE(xx,yy)
#' print(re)
#' @references Sun Y, Wu D. A relief based feature extraction algorithm[C]//Proceedings of the 2008 SIAM International Conference on Data Mining. Society for Industrial and Applied Mathematics, 2008: 188-195.
LFE<-function(xx,yy,T=5){
  N<-nrow(xx)
  p<-ncol(xx)
  if ( (T<1)|(T>N)){
    stop("use wrong T")
  }
  M<-matrix(0,nrow = p,ncol = p)
  ins_update<-sample(c(1:N),T)
  weight<-sapply(c(1:T),function(i){
    k<-ins_update[i]
    tmp_yyy<-abs(yy-yy[k])
    tmp_yy<-rep(0,length(tmp_yyy))
    tmp_yy[which(tmp_yyy==0)]<-1
    dis<-colSums((t(xx)-xx[k,])^2)
    dis_h<-tmp_yy*dis
    dis_m<-(1-tmp_yy)*dis
    nh<-which(dis_h==min(dis_h[dis_h>0]))[1]
    nm<-which(dis_m==min(dis_m[dis_m>0]))[1]
    MM<-abs(t(xx)-as.numeric(xx[nm,]))%*%
      (t(abs(t(xx)-as.numeric(xx[nm,]))))-
      abs(t(xx)-as.numeric(xx[nh,]))%*%
      (t(abs(t(xx)-as.numeric(xx[nh,]))))
    M<<-M+MM
  })
  eggvalue<-eigen(-M)
  eggvect<-eggvalue$vectors
  eggvalue<-eggvalue$values
  eggvalue<-pmax(eggvalue,0)
  if (sum(eggvalue)>0){
    eggvalue<-eggvalue/sum(eggvalue)
  } 
  new_w<-eggvect%*%(eggvalue*t(eggvect))
  new_w[which(new_w<0)]=0
  w<-new_w/sqrt(sum(new_w^2))
  class(w)<-"LFE"
  return(w)
}
