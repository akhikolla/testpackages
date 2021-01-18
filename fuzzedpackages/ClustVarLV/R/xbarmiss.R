xbarmiss<-function (Xmat){
  n<-nrow(Xmat)
  p<-ncol(Xmat)
  sinit=apply(Xmat,2,sd,na.rm=TRUE)
  nbmq<-apply(is.na(Xmat),1,sum)
  
  #initialisation
  xbar = apply(Xmat,1,mean,na.rm=TRUE)
  # if for one observation i, there is 0 or only 1 non-missing value, xbar[i] is set to 0
  xbar[which((p-nbmq)<2)]<-NA
  
 
  # iterative algorithm
  xbar=scale(xbar,center=TRUE,scale=FALSE)  # centering only
  old.xbar<-xbar
  cconv=1
  eps=1e-5
  iter=1
  while (cconv >eps){
    for(i in 1:n) Xmat[i,which(is.na(Xmat[i,]))]<-xbar[i]
    # each column of Xmat is centered and has the same standard deviation than initially
    snow<-apply(Xmat,2,sd,na.rm=TRUE)
    Xmat<-scale(Xmat,center=TRUE,scale=snow/sinit)
    xbar = apply(Xmat,1,mean,na.rm=TRUE) 
    xbar=scale(xbar,center=TRUE,scale=FALSE)     # centering only
    cconv=var(xbar-old.xbar,na.rm=TRUE)/var(old.xbar,na.rm=TRUE)
    if(is.na(cconv)) cconv=0
    old.xbar=xbar
    iter=iter+1
  }
  return(xbar)
}


