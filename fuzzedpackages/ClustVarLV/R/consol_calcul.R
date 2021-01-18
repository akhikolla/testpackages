consol_calcul <-
function(method,X,EXTr,Xr,EXTu,Xu,ind)
#
{
  n<-nrow(X)
  p<-ncol(X)
  pk<-length(ind)
  if (pk==1)   Xk<-matrix(X[,ind],dimnames=list(rownames(X),colnames(X)[ind]))
  if (pk>1)    Xk<-as.matrix(X[,ind])
  # verification if there are NA values
  valmq=FALSE
  if (sum(is.na(Xk))>0)  valmq=TRUE
  
  if (method==1){
    if ((EXTr==0)&(EXTu==0)) { 
      if (!valmq) {
        ressvd = svd(Xk)
        compnorm = ressvd$u[,1]
        comp = ressvd$u[,1] * (ressvd$d[1])
        critere = (ressvd$d[1])^2/(n-1)
      }
      if (valmq) {
        res<-cp1miss(Xk)
        comp<-res$comp
  #     compnorm<-comp/sqrt(t(comp)%*%comp)
        critere=var(res$comp,na.rm=TRUE)
      }
      veccor<-cor(comp,Xk,use="pairwise.complete.obs")
      jj<-which.max(abs(veccor))  # modification of the sign of comp so that it is positively correlated with the closest variable  
      comp=sign(veccor[jj])*comp
      res<-list(comp=comp,critere=critere)
    }
    if((EXTr==1)&(EXTu==0)) {
      if ( dim(Xr)[2] > dim(Xk)[2]   ) {  
        # vp <- eigen( t(t(Xr)%*%Xk)%*%t(Xr)%*%Xk)
        # a<-t(Xr)%*%Xk%*%vp$vectors[,1]%*%(vp$values[1])^(-1/2)
        vp <- powerEigen( crossprod(t(Xr)%*%Xk))
        a<-t(Xr)%*%Xk%*%vp$vectors%*%(vp$values)^(-1/2)
       } else {
        vp <- powerEigen(tcrossprod(t(Xr)%*%Xk))
        a<- vp$vectors
         # vp <- eigen(t(Xr)%*%Xk %*% t(t(Xr)%*%Xk))
         # a<- vp$vectors[,1] 
    }  
                     
      comp<- Xr%*% a        
      jj<-which.max(cor(comp,Xk))  # modification of the sign of comp so that it is positively correlated with the closest variable  
      if (cor(comp,Xk[,jj])<0) comp=(-1)*comp
      critere<-  vp$values /(n-1)
      res<-list(comp=comp,a=a,critere=critere)
    }
    if ((EXTr==0)&(EXTu==1)) {
      P<-Xk %*% Xu[ind,]
      if(sum(P^2)==0) stop("error in P")
      B<-t(Xk)%*% P
      vp = powerEigen( B)
      alpha2<-powerEigen(P)$values
      # vp = eigen(t(B) %*% B)
      # alpha2<-eigen(t(P)%*%P)$values[1] 
      crit<- vp$values/((n-1)*alpha2)      
      u<-vp$vectors
      comp<-P%*%u /sqrt(alpha2)
      jj<-which.max(cor(comp,P))  # modification of the sign of comp so that it is positively correlated with the closest variable  
      if (cor(comp,P[,jj])<0) comp=(-1)*comp
      res<-list(comp=comp,u=u,critere=crit)
    }
  }

  if (method==2){
    if ((EXTu==0)& (EXTr==0)) {
      if (!valmq) {
        comp <- Xk %*% matrix(1,pk,1) /pk                 
        #critere<-pk*var(comp)                                     # version RSA
        #res<-list(comp=comp,critere=critere)                      # version RSA       
        #compnorm<-comp/as.numeric(sqrt(t(comp)%*%comp))*sqrt(n-1)  # version CommStat  ck normalized
        critere<-pk*sd(comp,na.rm=TRUE)  
      }
      if (valmq) {
          comp= xbarmiss(Xk)         
          critere=pk*sd(comp,na.rm=TRUE)
      }
      res<-list(comp=comp,critere=critere) 
    }    
    if ((EXTu==0)& (EXTr==1)){                 
      aa = t(Xr)%*% Xk %*% matrix(1,pk,1) /pk
      a<-aa/as.numeric(sqrt(t(aa)%*%aa))
      comp<-(Xr%*%a)
      critere<-pk*sqrt(crossprod(aa))/(n-1)
      res<-list(comp=comp,a=a,critere=critere)
    }
    if ((EXTu==1)& (EXTr==0)) {
      Xugroupe<-Xu[ind,]
      P=Xk%*%Xugroupe
      if(sum(P^2)==0) stop("error in P")
      alpha2<-sum(diag(t(P)%*%P))
      uu = t(P)%*% Xk %*% matrix(1,pk,1) /pk
      u<-uu/as.numeric(sqrt(t(uu)%*%uu))
      comp<-(P%*%u)/sqrt(alpha2)
      critere<-pk*sqrt(t(uu)%*%uu)/sqrt(n-1)
      critere<-critere/sqrt(alpha2)
      res<-list(comp=comp,u=u,critere=critere)
 
    }   
  }
  
  return(res)
}
