consol_affect <-
function(method,X,Xr,Xu,EXTr,EXTu,comp,a,u)
{   

  p<-ncol(X)
  gtmp<-rep(0,p)
  n<-nrow(X)
  # verification if there are NA values
  valmq=FALSE
  if (sum(is.na(X))>0)  valmq=TRUE
 
  if (method == 1){
    if ((EXTu==0)&(EXTr==0)) {
      if (!valmq) {
      #comp<-comp/matrix(sqrt(diag(t(comp)%*%comp)),nrow=n,ncol=ncol(comp),byrow=T)  #normalization de comp
      comp<-scale(comp)
      cova = t(comp) %*% X/(n-1)
      covaC<-  cova^2
      }
      if (valmq) {
        comp<-scale(comp)
        #cova = cov(comp,X,use="pairwise.complete.obs")
        cova=covamiss(X,comp,method)
        covaC<-  cova^2
      }
      for (j in (1:p)) {   
        vec=covaC[,j] 
        maxj<-which.max(vec)
        gtmp[j] = maxj[1]   
      }
    } 
    if ((EXTu==0)&(EXTr==1)) {
      cova = t(comp) %*% X/(n-1)
      covaC <- cova^2
      for (j in (1:p)) {
        vec=covaC[,j]
        maxj<-which.max(vec)
        gtmp[j] = maxj[1]   
      }
    }  
    if ((EXTu==1)&(EXTr==0)){          
      cova = t(comp) %*% X/(n-1)
      covaC <- cova^2
      for (j in (1:p)) {
        vec=covaC[,j]
        maxj<-which.max(vec)
        gtmp[j] = maxj[1]   
      }
    }   
  }
  
  if (method==2){
    if ((EXTu==0)&(EXTr==0)) {
      if (!valmq) {
      #comp<-comp/matrix(sqrt(diag(t(comp)%*%comp/(n-1))),nrow=n,ncol=ncol(comp),byrow=T)  #normalization de comp
      comp<-scale(comp)
      cova = t(comp) %*% X/(n-1)
      }
      if (valmq) {
        comp<-scale(comp)
        #cova = cov(comp,X,use="pairwise.complete.obs")
        cova=covamiss(X,comp,method)
      }
      for (j in (1:p)) {
        vec=cova[,j]
        maxj<-which.max(vec)
        gtmp[j] = maxj[1]   
      }
    } 
    if ((EXTu==0)&(EXTr==1)) {
      cova = t(comp) %*% X/(n-1)
      covaC <- cova
      for (j in (1:p)) {
        vec=cova[,j]
        maxj<-which.max(vec)
        gtmp[j] = maxj[1]   
      }
    }  
    if ((EXTu==1)&(EXTr==0)){          
      for (j in 1:p) {
        Pj<- as.matrix(X[,j])%*%Xu[j,]
        critind=diag(t(u)%*%t(Pj)%*%comp)  
        vec=critind
        maxj<-which.max(vec)
        gtmp[j] = maxj[1]                    
      }
    }   
  } 
  
  return(gtmp)
}
