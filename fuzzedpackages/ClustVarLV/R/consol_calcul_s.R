consol_calcul_s <-
function(method,X,EXTr,Xr,EXTu,Xu,ind,max.iter=20, eps = 0.001,rlevel)
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
  
  dur <- function(a,para)
  {
    b = a
    b[which((abs(a) - para)<0)] = 0
    b
  }
  soft <- function(a,para)
  {    
    #b <- sort(abs(a))?
    b <- abs(a) - para
    b <- (b + abs(b))/2 
    b <- sign(a) * b
    b 
  }
  normvec <- function(a){     
    b = sqrt(sum(a^2))         
    if (b==0) b = 1
    b
  }
  trimm <- function(a,para)
  {
    b = a
    b[which((a - para)<0)] = 0
    b
  }
  
  
  
  if ((EXTr==0)&(EXTu==0)) { 
  #####################################################################################
  if (method==1){
  
      if (!valmq) {
        ressvd = svd(Xk)
        comp = ressvd$u[,1] * (ressvd$d[1])
        compnorm=scale(comp)
        critere = (ressvd$d[1])^2/(n-1)
      } 
      if (valmq) {
        rescp1<-cp1miss(Xk)
        comp<-rescp1$comp
        compnorm<-scale(comp)
        critere=var(rescp1$comp,na.rm=TRUE)
        matcov=as.matrix(cov(Xk,use="pairwise.complete.obs"))
     }
      Ck=compnorm  
        if (!valmq) cova = cov(Xk,Ck) 
        if (valmq)  cova=t(covamiss(Xk,Ck,method))
      para = rlevel*apply(Xk,2,sd,na.rm=TRUE)
      #beta = dur(cova,para)
                                   
      beta = soft(cova,para)   
      temp <- beta   # in order to check convergence 
      temp <- temp/normvec(temp)
      
      k <- 0
      diff <- 1
                       if ((mean(temp)==0)&(sd(temp)==0))  k=max.iter
      while ((k < max.iter) & (diff > eps)) {
      k <- k + 1 
       if (!valmq) {
           alpha <- t(Xk)%*%Xk%*%beta/normvec(t(Xk)%*%Xk%*%beta)
#           Ck = Xk%*%alpha/normvec(Xk%*%alpha)
           Ck = scale(Xk%*%alpha)
           cova = cov(Xk,Ck) 
       }
       if (valmq) {
            alpha <- matcov%*%beta/normvec(matcov%*%beta)  
            XXk<-Xk
            XXk[which(is.na(Xk))]<-0
            Ck = scale(XXk%*%alpha)
            cova=t(covamiss(Xk,Ck,method))    
       }
                                        
#       beta = dur(cova,para)
        beta = soft(cova,para)
        beta2 = beta/normvec(beta)
        diff <- mean(abs(beta2 - temp))
        temp <- beta2
         
      }  # end of loop of iteration "k"
                            
      beta = beta/normvec(beta)
      colnames(beta) = "loading" 
      loading=beta
      rownames(loading)= colnames(Xk)
                                 
      
      XXk<-Xk
      XXk[which(is.na(Xk))]<-0
      comp<-XXk%*%beta 
       
      if (sd(comp)>0)  {                                 
        veccor<-cor(comp,Xk,use="pairwise.complete.obs")
        jj<-which.max(abs(veccor))  # modification of the sign of comp so that it is positively correlated with the closest variable  
        comp =sign(veccor[jj])*comp 
      
        compnorm=scale(comp)
        critere<-sum(cov(Xk,compnorm,use="pairwise.complete.obs")^2)
      } else {
        critere=0
      }
      
      res<-list(comp=comp,loading=loading,critere=critere)
 
   }
  ###end of method=1 #########################################################################################
  
  
  
  
  
  #############################################################################################################
  if (method==2){
                                                
      if (!valmq) {
         comp <- Xk %*% matrix(1,pk,1) /pk                 
         #critere<-pk*var(comp)                             # version RSA
         #res<-list(comp=comp,critere=critere)              # version RSA       
         compnorm<-scale(comp)                              # version CommStat  ck normalized
         critere<-pk*sd(comp,na.rm=TRUE)                    # version CommStat
      }
      if (valmq) {
        comp= xbarmiss(Xk)         
        compnorm<-scale(comp)     
        critere=pk*sd(comp,na.rm=TRUE)
      }
      
      alpha = matrix(rep(1/pk,pk),pk,1)
    
      Ck = compnorm
      if (!valmq) cova = cov(Xk,Ck) 
      if (valmq)  cova=t(covamiss(Xk,Ck,method))
      para = rlevel*apply(Xk,2,sd,na.rm=TRUE)
    # beta = dur(cova,para)
      beta = trimm(cova,para)
      temp = beta
      
      k <- 0
      diff <- 1
      while ((k < max.iter) & (diff > eps)) {
        k <- k + 1
        
         # beta includes 0 or 1
          beta[which(beta!=0)]=1
          b = sum(beta)
          if (b==0) b = 1
          alpha <- beta/b
          if (!valmq) {
             Cmean = Xk%*%alpha
             Ck = scale(Cmean)  
             cova = cov(Xk,Ck)
          }
          if (valmq) {
            XXk<-Xk
            XXk[which(is.na(Xk))]<-0
            Cmean = XXk%*%alpha
            Ck = scale(Cmean)  
            cova = t(covamiss(Xk,Ck,method))
          }
        # beta = dur(cova,para)
        beta = trimm(cova,para)
        beta2<-beta
        
        beta[which(beta!=0)]=1

        diff <- mean(abs(beta2 - temp))
        temp <- beta2
      }
       
      
      comp = Cmean
   #   loading = beta
      loading =alpha
      rownames(loading)= colnames(Xk)
      colnames(loading) = "loading"  
 
      compnorm=scale(comp)
      critere<-sum(cov(Xk,compnorm,use="pairwise.complete.obs"))
      res<-list(comp=comp,loading=loading,critere=critere)

  }
  ###end of method=2 #########################################################################################
  }
  ############################################################################################################
  
  return(res)
}
