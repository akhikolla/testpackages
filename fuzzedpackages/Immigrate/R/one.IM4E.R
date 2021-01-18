#' one.IM4E
#'
#' This function performs (IM4E)Iterative Margin-Maximization under Max-Min Entropy algorithm for one loop.
#' @param train_xx model matrix of explanatory variables
#' @param train_yy label vector
#' @param w initial weight
#' @param sig sigma used in algorithm, default to be 1
#' @param lambda lambda used in algorithm, default to be 1
#' @keywords IM4E
#' @return \item{w}{new weight vector after one loop} 
#' @return \item{C}{cost after one loop} 
#' 
#' 
one.IM4E<-function(train_xx,train_yy,w,sig=1,lambda=1){
  # compute the sample size and the number of variables
  N<-nrow(train_xx)
  p<-ncol(train_xx)
  # compute coefficients
  # coefficients of hit
  entropy <- 0
  vv<-sapply(c(1:N), function(i){
    yyy<-abs(train_yy-train_yy[i])
    yy<-rep(0,length(yyy))
    yy[which(yyy==0)]<-1
    tmp <- exp(-(w%*%abs(t(train_xx)-as.numeric(train_xx[i,]))/sig) )[1,]
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
    abs(t(train_xx)-as.numeric(train_xx[i,]))%*%(tmp0-tmp1)
  })
  v<--rowSums(vv)
  v<-pmax(v,0)
  if (sum(v)>0){
    v<-v/sum(v)
  } 
  C<-v%*%rowSums(vv)+sig*entropy+lambda*sum(v^2)
  newList<-list("w" = v,"C"=C)
  return(newList)
}
