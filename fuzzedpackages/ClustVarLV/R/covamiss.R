covamiss<-function (mat,lv,method){
  n<-nrow(mat)
  p<-ncol(mat)
  sinit=apply(mat,2,sd,na.rm=TRUE)
  K<-ncol(lv)
  cova<-matrix(NA,nrow=K,ncol=p)
  cmat<-mat
  lmat<-as.list(data.frame(mat))
  
  for (k in 1:K) {
    if (method==1) compl<-lapply(X=lmat,FUN=compl1j, y=lv[,k])
    if (method==2) compl<-lapply(X=lmat,FUN=compl2j, y=lv[,k])
    cmat<-do.call(cbind,compl)
    snow<-apply(cmat,2,sd,na.rm=TRUE)
    mat<-scale(cmat,center=TRUE,scale=snow/sinit)
    cova[k,]<-cov(lv[,k],cmat,use="pairwise.complete.obs")
  }
 
   return(cova)
}

compl1j<-function(x,y) {x[which(is.na(x))]<-sign(cor(x,y,use="pairwise.complete.obs"))*y[which(is.na(x))]; return(x)}

compl2j<-function(x,y) {x[which(is.na(x))]<-y[which(is.na(x))]; return(x)}


