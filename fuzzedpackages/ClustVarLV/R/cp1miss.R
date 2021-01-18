cp1miss<-function (Xmat)   {
  
  # Xmat : matrix (already centered and possibly standardized)
  # H : nb dim to extract
  
  p=ncol(Xmat)
  n=nrow(Xmat)
  sinit=apply(Xmat,2,sd,na.rm=TRUE)
  nbmq<-apply(is.na(Xmat),1,sum)
  
  # initialisation
  res=nipals_1(Xmat,H=1)
  pc1=res$comp               # pc1 may contains some NA
  load=res$load   
 
  
  # iterative algorithm
  pc1=scale(pc1,center=TRUE,scale=FALSE)  # centering only
  old.pc1<-pc1
  cconv=1
  eps=1e-5
  iter=1
  while (cconv >eps){
    # each missing value is imputed by the corresponding value in pc1, taking account of the sign of the correlation between variables and pc1
    for (j in 1:p)  {
          Xmat[which(is.na(Xmat[,j])),j]<-pc1[which(is.na(Xmat[,j]))]*sign(cor(Xmat[,j],pc1,use="pairwise.complete.obs"))
    }
    # each column of Xmat is centered and has the same standard deviation than initially
    snow<-apply(Xmat,2,sd,na.rm=TRUE)
    Xmat<-scale(Xmat,center=TRUE,scale=snow/sinit)
    res=nipals_1(Xmat,H=1)  #ressvd = svd(Xmat,nu=1,nv=1)
    pc1=res$comp            #pc1<-ressvd$u%*%ressvd$d[1]
    load=res$load          #xload<-ressvd$v
    cconv=var(pc1-old.pc1,na.rm=TRUE)/var(old.pc1,na.rm=TRUE)
    if(is.na(cconv)) cconv=0
    old.pc1=pc1
    iter=iter+1
  }
                            
  
  res=list(load=load,comp=pc1)
  return(res)
  
}
  



#-----------------------------------------------------------------------------
#                            PCA NIPALS algorithm for the first component
#                            missing values allowed
#-----------------------------------------------------------------------------
nipals_1=function(X,H=1){
  
  # X : matrix (already centered and possibly standardized)
  # H : nb dim to extract
  
  # ---------------- data centered and normalized  ---------------------------
  #X<-scale(X,center=T,scale=reduc)
  
  p=ncol(X)
  n=nrow(X)
  
  
  # ------------------------- location of missing values ---------------------
  loc<-which(is.na(X),arr.ind=T)
  nbmq<-apply(is.na(X),1,sum)

 #-- algorithm NIPALS ------------------------------------------
 # H=1

    eps=1e-5
    maxiter=50
    iter=0
    xy=1
    # intial component : 
    # one among the last column without missing data 
    # or the last one with imputation of NA with 0   
    # why the last : because with CLV hierarchical algo, the last var already belongs to a group
    choix<-setdiff(1:p,loc[,2])
    if (length(choix)>0) { 
        xpc<-X[,choix[length(choix)]]
    } else {
        xpc=X[,p] 
    }
    old.xpc=xpc
    lX<-as.list(data.frame(X))
    ltX<-as.list(data.frame(t(X)))
    while ((xy >eps)&(iter<maxiter)){ 
      iter=iter+1
      xload=rep(0,p)
      lload<-lapply(X=lX,FUN=fload, y=xpc)
      xload<-as.vector(do.call(cbind,lload))
#       for (j in 1:p){
#         setobs=intersect(which(!is.na(X[, j])),which(!is.na(xpc))) 
#         xload[j] <- sum(X[setobs, j] * xpc[setobs])/sqrt(sum(xpc[setobs]^2))
#       }
      xload=xload/sqrt(sum(xload*xload, na.rm = TRUE))  #normalization
      
      lpc<-lapply(X=ltX,FUN=fpc, y=xload)
      xpc<-as.vector(do.call(cbind,lpc))
#      for (i in 1:n){
#         setvar=intersect(which(!is.na(X[i,])),which(!is.na(xload)))
#         if (length(setvar)==0) {
#           xpc[i]<-NA
#         } else {
#           xpc[i] <- sum(X[i,setvar ] * xload[setvar], na.rm = TRUE)/sqrt(sum(xload[setvar]^2))
#         }
#      }
      if (iter>1){
        xy=var(xpc-old.xpc,na.rm=TRUE)/var(old.xpc,na.rm=TRUE)
        if(is.na(xy)) xy=0
      }
      old.xpc=xpc
    }
                                     
  res.nipals=list(load=xload,comp=xpc)
  return(res.nipals)
}

# functions required for nipals_1 

fload<-function(x,y) {
       setobs=intersect(which(!is.na(x)),which(!is.na(y)))
       z <- sum(x[setobs] * y[setobs])/sqrt(sum(y[setobs]^2))
       return(z)
}                                               
                                            
fpc<-function(x,y) {
       setvar=intersect(which(!is.na(x)),which(!is.na(y)))
       if (length(setvar)==0) {
           z<-NA
       } else {
          z <- sum(x[setvar ] * y[setvar])/sqrt(sum(y[setvar]^2))
       }
       return(z)
}

##################################################################################
