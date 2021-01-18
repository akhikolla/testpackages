crit_init <-
function(method,X,EXTr,Xr,EXTu,Xu)
### initial clustering criterion (each variable = one group) / hierarchy
{
  n<-nrow(X)
  p<-ncol(X)
  # verification if there are NA values
  valmq=FALSE
  if (sum(is.na(X))>0)  valmq=TRUE
 
  # method 1
  if (method==1) {
    if ((EXTu==0)&(EXTr==0)){ crit<- apply(X,2,var,na.rm=TRUE) } 
    if ((EXTu==0)&(EXTr==1)){ 
      XrX<- t(Xr)%*%X
      crit<- apply(XrX^2/(n-1), 2, sum,na.rm=TRUE)
    }
    if ((EXTu==1)&(EXTr==0)){
      crit=c()
      for (i in 1:p) {
        critk<-var(X[,i],na.rm=TRUE)
        crit=c(crit,critk)
      }
    }
  }
  # method 2                 
  if (method==2) {
     if ((EXTu==0)&(EXTr==0)){ 
      #crit<-apply(X,2,var,na.rm=TRUE)    # version RSA        ck=xbark
      crit<-apply(X,2,sd,na.rm=TRUE)      # version CommStat   ck normalized
    } 
    if ((EXTu==0)&(EXTr==1)){ 
      if (valmq) stop("The matrix X contains missing values. Use a X matrix without missing value for CLV with external data")
        px<-sqrt (diag(tcrossprod(t(X)%*%Xr)))
        crit<- px/(n-1)
    } 
    if ((EXTu==1)&(EXTr==0)){
      if (valmq) stop("The matrix X contains missing values. Use a X matrix without missing value for CLV with external data")
      crit=c()
      for (i in 1:p) {
        critk<- sqrt(crossprod(X[,i])/(n-1))
        crit=c(crit,critk)
      }
    }
  }
  
  return(crit)
}
