
one.sample.Immigrate<-function(train_xx,train_yy,W,sample_wt,sig=1){
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
    MM<<-MM+sample_wt[i]*abs(t(train_xx)-as.numeric(train_xx[i,]))%*%
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

sample.Immigrate<-function(train_xx,train_yy,sample_wt, W, epsilon=0.01,
                           sig=1, max_iter=10,removesmall=FALSE){
  suppressWarnings(
    return(ImmigrateSampleCpp(onesampleImmigrate = one.sample.Immigrate, train_xx, train_yy, sample_wt, W, 
                              epsilon = epsilon, sig = sig,  max_iter = max_iter, removesmall = removesmall)))
}


calc_error<-function(yy, pred, sample_wt){
  if(length(yy)!=length(pred)){
    stop("The length of labels and pred is not equal")
  }
  pos<-which(yy!=pred)
  return(sum(sample_wt[pos]))
}

update_sample_wt<-function(yy,pred,sample_wt,alp){
  if(length(yy)!=length(pred) || length(yy)!=length(sample_wt))
    stop("three vector sizes should be the same")
  update_wt<-sample_wt
  pos<-which(yy!=pred)
  
  update_wt[pos]<-update_wt[pos]*exp(alp)
  update_wt<-update_wt/sum(update_wt)
  return (update_wt)
}

pred.self.Immigrate<-function(re,xx,yy,sig = 1, type){
  yy<-as.numeric(yy)
  if(is.list(re)){
    w<-re$w
  }else{
    w<-re 
  }
  label<-unique(yy)
  v<-sapply(c(1:length(label)),function(j){
    sapply(c(1:nrow(xx)), function(i){
      yyy<-abs(yy-label[j])
      y<-rep(0,length(yyy))
      y[which(yyy==0)]<-1
      tmp<-matrix((y[-i])*exp(-(rowSums(crossprod(abs(t(xx[-i,])-as.numeric(xx[i,])),w)*
                                          t(abs(t(xx[-i,])-as.numeric(xx[i,]))))/sig) ))[,1]
      s<-sum(tmp)
      if(s!=0) tmp<-tmp/s
      rowSums(crossprod(abs(t(xx[-i,])-as.numeric(xx[i,])),w)*
                t(abs(t(xx[-i,])-as.numeric(xx[i,]))))%*%tmp
    })
  })
  
  
  myfun<-sapply(c(1:nrow(xx)), function(i){
    v[i,]<<-v[i,]/sum(v[i,])
  })
  pred<-sapply(c(1:nrow(xx)),function(i){
    label[which.min(v[i,])]
  })
  
  if (missing(type)){
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

one.boost.Immigrate<-function(train_xx,train_yy,sample_wt,W,sig=1,removesmall = FALSE, max_iter = 10)
{
  p<-ncol(train_xx)
  if (missing(W)){
    W <- matrix(runif(p*p,0,1),p,p)
    W <- W/sqrt(sum(W^2))
  }
  if (removesmall == T){
    re <- sample.Immigrate(train_xx=train_xx,train_yy=train_yy,sample_wt=sample_wt,W = W, sig=sig,max_iter = max_iter,removesmall = removesmall)
  }else{ 
    re <- sample.Immigrate(train_xx=train_xx,train_yy=train_yy,sample_wt=sample_wt,W = W,sig=sig,max_iter = max_iter)
  }
  pred_class<-pred.self.Immigrate(re,train_xx,train_yy,sig = sig,type = "class")
  err = calc_error(train_yy, pred_class, sample_wt)
  alp = 0.0
  if (err>0.5 || err == 0)
  {
    if( err == 0)
      alp = 1.0
    else
    {
      alp = 0.0
    }
  }
  else
  {
    alp = 0.5*log((1.-err)/err)
    sample_wt =  update_sample_wt(train_yy, pred_class, sample_wt, alp)
  }
  
  
  boost_result <- list("w" = re$w, "error" = err,"sample_wt"=sample_wt,"alpha" = alp)
  return (boost_result)
}
