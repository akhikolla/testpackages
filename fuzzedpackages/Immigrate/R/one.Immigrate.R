#' one.Immigrate
#'
#' This function performs Immigrate(Iterative Max-Min Entropy Margin-Maximization with Interaction Terms) algorithm for one loop.
#' @param train_xx model matrix of explanatory variables
#' @param train_yy label vector
#' @param W initial weight matrix
#' @param sig sigma used in algorithm, default to be 1
#' @keywords Immigrate one
#' @return \item{W}{new weight matrix after one loop}
#' @return \item{C}{cost after one loop}
#' @export
#' 
#' @examples
#' data(park)
#' xx<-park$xx
#' yy<-park$yy
#' W0 <- diag(rep(1,ncol(xx)),ncol(xx))/sqrt(ncol(xx))
#' re<-one.Immigrate(xx,yy,W0)
#' print(re$w)
#' @seealso Please refer to \url{https://github.com/RuzhangZhao/Immigrate/} for implementation demo.

one.Immigrate<-function(train_xx,train_yy,W,sig=1){
  # compute the sample size and the number of variables
  N<-nrow(train_xx)
  p<-ncol(train_xx)
  # compute coefficients
  entropy <- 0
  # compute new weight
  MM<-matrix(0,nrow = p,ncol = p)
  myfun<-sapply(c(1:N), function(i){
    yyy<-abs(train_yy-train_yy[i])
    yy<-rep(0,length(yyy))
    yy[which(yyy==0)]<-1
    tmp<-exp(-(rowSums(crossprod(abs(t(train_xx)-as.numeric(train_xx[i,])),W)*t(abs(t(train_xx)-as.numeric(train_xx[i,]))))/sig) )
    tmp0 <- yy*tmp
    tmp0[i]<-0
    s0<-sum(tmp0)
    if (s0 != 0 ) tmp0<-tmp0/s0
    tmp1 <- (1-yy)*tmp
    s1<-sum(tmp1)
    if (s1 != 0 ) tmp1<-tmp1/s1
    tmp0[i]<-1
    entropy <<- entropy + sum((tmp0-tmp1)*log(abs(tmp0-tmp1)))
    tmp0[i]<-0
    MM<<-MM+abs(t(train_xx)-as.numeric(train_xx[i,]))%*%
      ((tmp0-tmp1)*t(abs(t(train_xx)-as.numeric(train_xx[i,]))))
    0
  })
  eggvalue<-eigen(-MM)
  eggvect<-eggvalue$vectors
  eggvalue<-eggvalue$values
  eggvalue<-pmax(eggvalue,0)
  if (sum(eggvalue)>0){
    eggvalue<-eggvalue/sum(eggvalue)
  } 
  new_W<-eggvect%*%(eggvalue*t(eggvect))
  new_W[which(new_W<0)]=0
  new_W<-new_W/sqrt(sum(new_W^2))
  C<-sum(crossprod(eggvect,MM)*t(eggvect))+ sig*entropy
  newList<-list("w" = new_W,"C"=C)
  return(newList)
}                          
