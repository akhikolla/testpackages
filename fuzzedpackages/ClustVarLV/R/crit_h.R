crit_h <-
function(method,X12,EXTr,Xr,EXTu,Xu12,tauxNA=0)
### clustering criterion / hierarchy
{
  n<-nrow(X12)
  pk<-ncol(X12)  
  # verification if there are NA values
  valmq=FALSE
  if (sum(is.na(X12))>0)  valmq=TRUE
  
  if (method == 1) {
    if((EXTr==0)&(EXTu==0)) {
      xnew<- X12  
      if(!valmq) {
        if (dim(xnew)[1]>dim(xnew)[2]) vp = powerEigen( crossprod(xnew))
        else vp = powerEigen( tcrossprod(xnew))
        # if (dim(xnew)[1]>dim(xnew)[2]) vp = eigen( t(xnew) %*% xnew)
        # else vp = eigen( xnew %*% t(xnew))
        crit = vp$values[1] /(n-1)  
      }
      if(valmq)  {
        res=cp1miss(xnew)
        if ( ((sum(is.na(res$comp)))/n)> max(tauxNA^2,0.20)   ) {
            crit=NA
         } else {
         crit=var(res$comp,na.rm=TRUE)
       }
      }
    }   
    if((EXTr==1)&(EXTu==0)) {
      xnew<- t(Xr)%*%X12  
      if (dim(xnew)[1]>dim(xnew)[2]) vp = powerEigen( crossprod(xnew))
      else vp = powerEigen( tcrossprod(xnew))
      # if (dim(xnew)[1]>dim(xnew)[2]) vp = eigen( t(xnew) %*% xnew)
      # else vp = eigen( xnew %*% t(xnew))
      crit = vp$values[1] /(n-1)  
    }
    if((EXTr==0)&(EXTu==1))  {
      P<-X12 %*% Xu12
      B<-t(X12)%*%P
      vp <- powerEigen( crossprod(B))
      alpha2<-powerEigen(crossprod(P))$values[1]
      # vp <- eigen(t(B) %*% B)
      # alpha2<-eigen(t(P)%*%P)$values[1]
      crit= vp$values[1]/((n-1)*alpha2)
    }
  }
  
  if (method==2) {
    if ((EXTu==0)&(EXTr==0)){                          # version RSA      ck = xbark
      if(!valmq) {
        xbar = X12 %*% matrix(1,pk,1) /pk        
        crit = pk * sqrt(1/(n-1) * t(xbar) %*% xbar)  
        #crit = pk * 1/(n-1) * t(xbar) %*% xbar        # version CommStat ck normalized
      }
      if(valmq)  {
           xbar= xbarmiss(X12)         # average for each obs is computed only if there are more than 2 non-missing values
           if ( ((sum(is.na(xbar)))/n)> max(tauxNA^2,0.20)   ) {
             crit=NA
           } else {
            crit=pk*sd(xbar,na.rm=TRUE)
           }
      }
    }
    if ((EXTu==0)&(EXTr==1)){                              
      xbar = X12 %*% matrix(1,pk,1) /pk                
      pxbar=sqrt (t(xbar)%*%Xr%*%t(Xr)%*%xbar)
      crit= pk * 1/(n-1) *pxbar 
    }
    if ((EXTu==1)&(EXTr==0)){
      P=X12 %*% Xu12
      alpha2<-sum(diag(crossprod((P))))
      xbar = X12 %*% matrix(1,pk,1) /pk       
      pxbar=sqrt (t(xbar)%*%P%*%t(P)%*%xbar)
                 if (is.nan(pxbar)) { print(Xu12)}
      crit<- pxbar*pk/sqrt(n-1)
      crit<-crit/sqrt(alpha2)
    }
  }
 
return(crit)
}
