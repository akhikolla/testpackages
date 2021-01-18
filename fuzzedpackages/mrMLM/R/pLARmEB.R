pLARmEB<-function(gen,phe,outATCG,genRaw,kk,psmatrix,CriLOD,lars1,Genformat,Bootstrap,CLO){
  
    lodvalue<-CriLOD
    gene.data<-gen
    
    rm(gen)
    gc()
    
    inputform<-Genformat
    
    if(is.null(psmatrix)){
      flagps<-1
    }else{
      flagps<-0
    }
    
    if(is.null(lodvalue)==TRUE||is.null(lars1)==TRUE){
      warning("Please set parameter!")
    }
    if(lodvalue<0)
    {
      warning("Please input critical LOD score: > 0 !")
    }
    if(lars1<0||lars1>=nrow(phe))
    {
      warning("Please input the number of most relevant variables select by LARS: >0 and less than numbers of sample!")
    }
    if(is.null(gene.data)==TRUE)
    {
      warning("Please input correct genotypic data !")
      
    }
    if(is.null(phe)==TRUE)
    {
      warning("Please input correct phenotypic data !")
    }
    if((is.null(gene.data)==FALSE)&&(is.null(phe)==FALSE)&&(ncol(gene.data)!=(nrow(phe)+2)))
    {
      warning("Sample size in genotypic dataset doesn't equal to the sample size in phenotypic dataset !")
    }
    
    if((is.null(gene.data)==FALSE)&&(is.null(phe)==FALSE)&&((ncol(gene.data)==(nrow(phe)+2)))&&(lodvalue>=0)&&(lars1>0))
    {
      
      wan<-NULL
      result<-NULL
      
      multinormal<-function(y,mean,sigma)
      {
        pdf_value<-(1/sqrt(2*3.14159265358979323846*sigma))*exp(-(y-mean)*(y-mean)/(2*sigma));
        return (pdf_value)
      }
      
      ebayes_EM<-function(x,z,y)
      {
        n<-nrow(z);k<-ncol(z)
        
        if(abs(min(eigen(crossprod(x,x))$values))<1e-6){
          b<-solve(crossprod(x,x)+diag(ncol(x))*1e-8)%*%crossprod(x,y)
        }else{
          b<-solve(crossprod(x,x))%*%(crossprod(x,y))
        }
        
        v0<-as.numeric(crossprod((y-x%*%b),(y-x%*%b))/n)
        u<-matrix(rep(0,k),k,1)
        v<-matrix(rep(0,k),k,1)
        s<-matrix(rep(0,k),k,1)
        for(i in 1:k)
        {
          zz<-z[,i]
          s[i]<-((crossprod(zz,zz)+1e-100)^(-1))*v0
          u[i]<-s[i]*crossprod(zz,(y-x%*%b))/v0
          v[i]<-u[i]^2+s[i]
        }
        
        vv<-matrix(rep(0,n*n),n,n);
        for(i in 1:k)
        {
          zz<-z[,i]
          vv=vv+tcrossprod(zz,zz)*v[i]
        }
        vv<-vv+diag(n)*v0
        
        iter<-0;err<-1000;iter_max<-500;err_max<-1e-8
        tau<-0;omega<-0
        while((iter<iter_max)&&(err>err_max))
        {
          iter<-iter+1
          v01<-v0
          v1<-v
          b1<-b
          vi<-solve(vv)
          xtv<-crossprod(x,vi)
          
          if(ncol(x)==1)
          {
            b<-((xtv%*%x)^(-1))*(xtv%*%y)
          }else{
            if(abs(min(eigen(xtv%*%x)$values))<1e-6){
              b<-solve((xtv%*%x)+diag(ncol(x))*1e-8)%*%(xtv%*%y)
            }else{
              b<-solve(xtv%*%x)%*%(xtv%*%y)
            }
          }
          r<-y-x%*%b
          ss<-matrix(rep(0,n),n,1)
          for(i in 1:k)
          {
            zz<-z[,i]
            zztvi<-crossprod(zz,vi)
            u[i]<-v[i]*zztvi%*%r
            s[i]<-v[i]*(1-zztvi%*%zz*v[i])
            v[i]<-(u[i]^2+s[i]+omega)/(tau+3)
            ss<-ss+zz*u[i]
          }
          v0<-as.numeric(crossprod(r,(r-ss))/n)
          
          vv<-matrix(rep(0,n*n),n,n)
          for(i in 1:k)
          {
            zz<-z[,i]
            vv<-vv+tcrossprod(zz,zz)*v[i]
          }
          vv<-vv+diag(n)*v0
          
          err<-(crossprod((b1-b),(b1-b))+(v01-v0)^2+crossprod((v1-v),(v1-v)))/(2+k)
          beta<-t(b)
          sigma2<-v0
        }
        
        wang<-matrix(rep(0,k),k,1)
        for (i in 1:k){
          stderr<-sqrt(s[i]+1e-20)
          t<-abs(u[i])/stderr
          f<-t*t
          p<-pchisq(f,1,lower.tail = F)
          wang[i]<-p
        }
        
        return(list(u=u,sigma2=sigma2,wang=wang))
      }
      
      likelihood<-function(xxn,xxx,yn,bbo)
      {
        nq<-ncol(xxx)
        ns<-nrow(yn)
        at1<-0
        
        if(is.null(bbo)==TRUE){
          ww1<-1:ncol(xxx)
          ww1<-as.matrix(ww1)
        }else{
          ww1<-as.matrix(which(abs(bbo)>1e-5))
        }
        at1<-dim(ww1)[1]
        lod<-matrix(rep(0,nq),nq,1)
        if(at1>0.5)
          ad<-cbind(xxn,xxx[,ww1])
        else
          ad<-xxn
        if(abs(min(eigen(crossprod(ad,ad))$values))<1e-6)
          bb<-solve(crossprod(ad,ad)+diag(ncol(ad))*0.01)%*%crossprod(ad,yn)
        else
          bb<-solve(crossprod(ad,ad))%*%crossprod(ad,yn)
        vv1<-as.numeric(crossprod((yn-ad%*%bb),(yn-ad%*%bb))/ns);
        ll1<-sum(log(abs(multinormal(yn,ad%*%bb,vv1))))
        sub<-1:ncol(ad);
        if(at1>0.5)
        {
          for(i in 1:at1)
          {
            ij<-which(sub!=sub[i+ncol(xxn)])
            ad1<-ad[,ij]
            if(abs(min(eigen(crossprod(ad1,ad1))$values))<1e-6)
              bb1<-solve(crossprod(ad1,ad1)+diag(ncol(ad1))*0.01)%*%crossprod(ad1,yn)
            else
              bb1<-solve(crossprod(ad1,ad1))%*%crossprod(ad1,yn) 
            vv0<-as.numeric(crossprod((yn-ad1%*%bb1),(yn-ad1%*%bb1))/ns);
            ll0<-sum(log(abs(multinormal(yn,ad1%*%bb1,vv0))))
            lod[ww1[i]]<--2.0*(ll0-ll1)/(2.0*log(10))
          }
        }
        return (lod)
      }
      
      
      emma.eigen.L <- function(Z,K,complete=TRUE) {
        if ( is.null(Z) ) {
          return(emma.eigen.L.wo.Z(K))
        }
        else {
          return(emma.eigen.L.w.Z(Z,K,complete))
        }
      }
      #likelihood
      emma.eigen.L.wo.Z <- function(K) {
        eig <- eigen(K,symmetric=TRUE)
        return(list(values=eig$values,vectors=eig$vectors))
      }
      #likelihood
      emma.eigen.L.w.Z <- function(Z,K,complete=TRUE) {
        if ( complete == FALSE ) {
          vids <- colSums(Z)>0
          Z <- Z[,vids]
          K <- K[vids,vids]
        }
        eig <- eigen(K%*%crossprod(Z,Z),symmetric=FALSE,EISPACK=TRUE)
        return(list(values=eig$values,vectors=qr.Q(qr(Z%*%eig$vectors),complete=TRUE)))
      }
      
      #restricted likelihood
      emma.eigen.R <- function(Z,K,X,complete=TRUE) {
        if ( ncol(X) == 0 ) {
          return(emma.eigen.L(Z,K))
        }
        else if ( is.null(Z) ) {
          return(emma.eigen.R.wo.Z(K,X))
        }
        else {
          return(emma.eigen.R.w.Z(Z,K,X,complete))
        }
      }
      #restricted likelihood
      emma.eigen.R.wo.Z <- function(K, X) {
        n <- nrow(X)
        q <- ncol(X)
        S <- diag(n)-X%*%solve(crossprod(X,X))%*%t(X)
        eig <- eigen(S%*%(K+diag(1,n))%*%S,symmetric=TRUE)
        stopifnot(!is.complex(eig$values))
        return(list(values=eig$values[1:(n-q)]-1,vectors=eig$vectors[,1:(n-q)]))
      }
      
      emma.eigen.R.w.Z <- function(Z, K, X, complete = TRUE) {
        if ( complete == FALSE ) {
          vids <-  colSums(Z) > 0
          Z <- Z[,vids]
          K <- K[vids,vids]
        }
        n <- nrow(Z)
        t <- ncol(Z)
        q <- ncol(X)
        
        SZ <- Z - X%*%solve(crossprod(X,X))%*%crossprod(X,Z)
        eig <- eigen(K%*%crossprod(Z,SZ),symmetric=FALSE)
        if ( is.complex(eig$values) ) {
          eig$values <- Re(eig$values)
          eig$vectors <- Re(eig$vectors)    
        }
        qr.X <- qr.Q(qr(X))
        return(list(values=eig$values[1:(t-q)],
                    vectors=qr.Q(qr(cbind(SZ%*%eig$vectors[,1:(t-q)],qr.X)),
                                 complete=TRUE)[,c(1:(t-q),(t+1):n)]))   
      }
      
      emma.delta.ML.LL.wo.Z <- function(logdelta, lambda, etas, xi) {
        n <- length(xi)
        delta <- exp(logdelta)
        return( 0.5*(n*(log(n/(2*pi))-1-log(sum((etas*etas)/(delta*lambda+1))))-sum(log(delta*xi+1))) )  
      }
      
      emma.delta.ML.LL.w.Z <- function(logdelta, lambda, etas.1, xi.1, n, etas.2.sq ) {
        delta <- exp(logdelta)
        return( 0.5*(n*(log(n/(2*pi))-1-log(sum(etas.1*etas.1/(delta*lambda+1))+etas.2.sq))-sum(log(delta*xi.1+1)) ))
        
      }
      
      emma.delta.ML.dLL.wo.Z <- function(logdelta, lambda, etas, xi) {
        n <- length(xi)
        delta <- exp(logdelta)
        etasq <- etas*etas
        ldelta <- delta*lambda+1
        return( 0.5*(n*sum(etasq*lambda/(ldelta*ldelta))/sum(etasq/ldelta)-sum(xi/(delta*xi+1))) )
      }
      
      emma.delta.ML.dLL.w.Z <- function(logdelta, lambda, etas.1, xi.1, n, etas.2.sq ) {
        delta <- exp(logdelta)
        etasq <- etas.1*etas.1
        ldelta <- delta*lambda+1
        return( 0.5*(n*sum(etasq*lambda/(ldelta*ldelta))/(sum(etasq/ldelta)+etas.2.sq)-sum(xi.1/(delta*xi.1+1))) )
      }
      
      emma.delta.REML.LL.wo.Z <- function(logdelta, lambda, etas) {
        nq <- length(etas)
        delta <-  exp(logdelta)
        return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas*etas/(delta*lambda+1))))-sum(log(delta*lambda+1))) )
      }
      
      emma.delta.REML.LL.w.Z <- function(logdelta, lambda, etas.1, n, t, etas.2.sq ) {
        tq <- length(etas.1)
        nq <- n - t + tq
        delta <-  exp(logdelta)
        return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas.1*etas.1/(delta*lambda+1))+etas.2.sq))-sum(log(delta*lambda+1))) ) 
      }
      
      emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
        nq <- length(etas)
        delta <- exp(logdelta)
        etasq <- etas*etas
        ldelta <- delta*lambda+1
        return( 0.5*(nq*sum(etasq*lambda/(ldelta*ldelta))/sum(etasq/ldelta)-sum(lambda/ldelta)) )
      }
      
      emma.delta.REML.dLL.w.Z <- function(logdelta, lambda, etas.1, n, t1, etas.2.sq ) {
        t <- t1
        tq <- length(etas.1)
        nq <- n - t + tq
        delta <- exp(logdelta)
        etasq <- etas.1*etas.1
        ldelta <- delta*lambda+1
        return( 0.5*(nq*sum(etasq*lambda/(ldelta*ldelta))/(sum(etasq/ldelta)+etas.2.sq)-sum(lambda/ldelta) ))
      }
      
      
      emma.MLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10,
                           esp=1e-10, eig.L = NULL, eig.R = NULL)
      {
        n <- length(y)
        t <- nrow(K)
        q <- ncol(X)
        stopifnot(ncol(K) == t)
        stopifnot(nrow(X) == n)
        if ( det(crossprod(X,X)) == 0 ) {
          warning("X is singular")
          return (list(ML=0,delta=0,ve=0,vg=0))
        }
        
        if ( is.null(Z) ) {
          if ( is.null(eig.L) ) {
            eig.L <- emma.eigen.L.wo.Z(K)
          }
          if ( is.null(eig.R) ) {
            eig.R <- emma.eigen.R.wo.Z(K,X)
          }
          etas <- crossprod(eig.R$vectors,y)
          logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
          m <- length(logdelta)
          delta <- exp(logdelta)
          Lambdas.1<-matrix(eig.R$values,n-q,m)    
          Lambdas <- Lambdas.1 * matrix(delta,n-q,m,byrow=TRUE)+1
          Xis.1<-matrix(eig.L$values,n,m)
          Xis <- Xis.1* matrix(delta,n,m,byrow=TRUE)+1
          Etasq <- matrix(etas*etas,n-q,m)
          dLL <- 0.5*delta*(n*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(Xis.1/Xis))
          optlogdelta <- vector(length=0)
          optLL <- vector(length=0)
          if ( dLL[1] < esp ) {
            optlogdelta <- append(optlogdelta, llim)
            optLL <- append(optLL, emma.delta.ML.LL.wo.Z(llim,eig.R$values,etas,eig.L$values))
          }
          if ( dLL[m-1] > 0-esp ) {
            optlogdelta <- append(optlogdelta, ulim)
            optLL <- append(optLL, emma.delta.ML.LL.wo.Z(ulim,eig.R$values,etas,eig.L$values))
          }
          
          for( i in 1:(m-1) )
          {
            if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
            {
              r <- uniroot(emma.delta.ML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas, xi=eig.L$values)
              optlogdelta <- append(optlogdelta, r$root)
              optLL <- append(optLL, emma.delta.ML.LL.wo.Z(r$root,eig.R$values, etas, eig.L$values))
            }
          }
        }
        else {
          if ( is.null(eig.L) ) {
            eig.L <- emma.eigen.L.w.Z(Z,K)
          }
          if ( is.null(eig.R) ) {
            eig.R <- emma.eigen.R.w.Z(Z,K,X)
          }
          etas <- crossprod(eig.R$vectors,y)
          etas.1 <- etas[1:(t-q)]
          etas.2 <- etas[(t-q+1):(n-q)]
          etas.2.sq <- sum(etas.2*etas.2)
          logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
          m <- length(logdelta)
          delta <- exp(logdelta)
          Lambdas.1<-matrix(eig.R$values,t-q,m)
          Lambdas <- Lambdas.1 * matrix(delta,t-q,m,byrow=TRUE) + 1
          Xis.1<-matrix(eig.L$values,t,m)
          Xis <- Xis.1 * matrix(delta,t,m,byrow=TRUE) + 1
          Etasq <- matrix(etas.1*etas.1,t-q,m)
          dLL <- 0.5*delta*(n*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/(colSums(Etasq/Lambdas)+etas.2.sq)-colSums(Xis.1/Xis))
          optlogdelta <- vector(length=0)
          optLL <- vector(length=0)
          if ( dLL[1] < esp ) {
            optlogdelta <- append(optlogdelta, llim)
            optLL <- append(optLL, emma.delta.ML.LL.w.Z(llim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
          }
          if ( dLL[m-1] > 0-esp ) {
            optlogdelta <- append(optlogdelta, ulim)
            optLL <- append(optLL, emma.delta.ML.LL.w.Z(ulim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
          }
          
          for( i in 1:(m-1) )
          {
            if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
            {
              r <- uniroot(emma.delta.ML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, xi.1=eig.L$values, n=n, etas.2.sq = etas.2.sq )
              optlogdelta <- append(optlogdelta, r$root)
              optLL <- append(optLL, emma.delta.ML.LL.w.Z(r$root,eig.R$values, etas.1, eig.L$values, n, etas.2.sq ))
            }
          }
        }
        maxdelta <- exp(optlogdelta[which.max(optLL)])
        optLL=replaceNaN(optLL)
        maxLL <- max(optLL)
        if ( is.null(Z) ) {
          maxve <- sum(etas*etas/(maxdelta*eig.R$values+1))/n    
        }
        else {
          maxve <- (sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/n
        }
        maxvg <- maxve*maxdelta
        
        return (list(ML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg))
      }
      
      emma.REMLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10,
                             esp=1e-10, eig.L = NULL, eig.R = NULL) {
        n <- length(y)
        t <- nrow(K)
        q <- ncol(X)
        stopifnot(ncol(K) == t)
        stopifnot(nrow(X) == n)
        if ( det(crossprod(X,X)) == 0 ) {
          warning("X is singular")
          return (list(REML=0,delta=0,ve=0,vg=0))
        }
        
        if ( is.null(Z) ) {
          if ( is.null(eig.R) ) {
            eig.R <- emma.eigen.R.wo.Z(K,X)
          }
          etas <- crossprod(eig.R$vectors,y)
          logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
          m <- length(logdelta)
          delta <- exp(logdelta)
          Lambdas.1<-matrix(eig.R$values,n-q,m)
          Lambdas <- Lambdas.1 * matrix(delta,n-q,m,byrow=TRUE) + 1
          Etasq <- matrix(etas*etas,n-q,m)
          dLL <- 0.5*delta*((n-q)*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(Lambdas.1/Lambdas))
          optlogdelta <- vector(length=0)
          optLL <- vector(length=0)
          if ( dLL[1] < esp ) {
            optlogdelta <- append(optlogdelta, llim)
            optLL <- append(optLL, emma.delta.REML.LL.wo.Z(llim,eig.R$values,etas))
          }
          if ( dLL[m-1] > 0-esp ) {
            optlogdelta <- append(optlogdelta, ulim)
            optLL <- append(optLL, emma.delta.REML.LL.wo.Z(ulim,eig.R$values,etas))
          }
          
          for( i in 1:(m-1) )
          {
            if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
            {
              r <- uniroot(emma.delta.REML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas)
              optlogdelta <- append(optlogdelta, r$root)
              optLL <- append(optLL, emma.delta.REML.LL.wo.Z(r$root,eig.R$values, etas))
            }
          }
        }
        else {
          if ( is.null(eig.R) ) {
            eig.R <- emma.eigen.R.w.Z(Z,K,X)
          }
          etas <- crossprod(eig.R$vectors,y)
          etas.1 <- etas[1:(t-q)]
          etas.2 <- etas[(t-q+1):(n-q)]
          etas.2.sq <- sum(etas.2*etas.2)
          logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
          m <- length(logdelta)
          delta <- exp(logdelta)
          Lambdas.1 <- matrix(eig.R$values,t-q,m) 
          Lambdas <- Lambdas.1 * matrix(delta,t-q,m,byrow=TRUE) + 1
          Etasq <- matrix(etas.1*etas.1,t-q,m)
          dLL <- 0.5*delta*((n-q)*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/(colSums(Etasq/Lambdas)+etas.2.sq)-colSums(Lambdas.1/Lambdas))
          optlogdelta <- vector(length=0)
          optLL <- vector(length=0)
          if ( dLL[1] < esp ) {
            optlogdelta <- append(optlogdelta, llim)
            optLL <- append(optLL, emma.delta.REML.LL.w.Z(llim,eig.R$values,etas.1,n,t,etas.2.sq))
          }
          if ( dLL[m-1] > 0-esp ) {
            optlogdelta <- append(optlogdelta, ulim)
            optLL <- append(optLL, emma.delta.REML.LL.w.Z(ulim,eig.R$values,etas.1,n,t,etas.2.sq))
          }
          
          for( i in 1:(m-1) )
          {
            if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
            {
              r <- uniroot(emma.delta.REML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, n=n, t1=t, etas.2.sq = etas.2.sq )
              optlogdelta <- append(optlogdelta, r$root)
              optLL <- append(optLL, emma.delta.REML.LL.w.Z(r$root,eig.R$values, etas.1, n, t, etas.2.sq ))
            }
          }
        }  
        maxdelta <- exp(optlogdelta[which.max(optLL)])
        optLL=replaceNaN(optLL)
        maxLL <- max(optLL)
        if ( is.null(Z) ) {
          maxve <- sum(etas*etas/(maxdelta*eig.R$values+1))/(n-q)    
        }
        else {
          maxve <- (sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/(n-q)
        }
        maxvg <- maxve*maxdelta
        return (list(REML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg))
      }
      
      
      emma.maineffects.B<-function(Z=NULL,K,deltahat.g,complete=TRUE){
        if( is.null(Z) ){
          return(emma.maineffects.B.Zo(K,deltahat.g))
        }
        else{
          return(emma.maineffects.B.Z(Z,K,deltahat.g,complete))
        }
      }
      
      
      emma.maineffects.B.Zo <-function(K,deltahat.g){
        t <- nrow(K)
        stopifnot(ncol(K) == t)
        B<-deltahat.g*K+diag(1,t)
        eig<-eigen(B,symmetric=TRUE)
        qr.B<-qr(B)
        q<-qr.B$rank
        stopifnot(!is.complex(eig$values))
        A<-diag(1/sqrt(eig$values[1:q]))
        Q<-eig$vectors[,1:q]
        C<-Q%*%A%*%t(Q)
        return(list(mC=C,Q=Q,A=A))
      }
      
      emma.maineffects.B.Z <- function(Z,K,deltahat.g,complete=TRUE){
        if ( complete == FALSE ) {
          vids <- colSums(Z)>0
          Z <- Z[,vids]
          K <- K[vids,vids]
        }
        n <- nrow(Z)  
        B <- deltahat.g*Z%*%K%*%t(Z)+diag(1,n)
        eig <- eigen(B,symmetric=TRUE,EISPACK=TRUE)
        qr.B<-qr(B)
        q<-qr.B$rank
        stopifnot(!is.complex(eig$values))
        A<-diag(1/sqrt(eig$values[1:q]))
        Q<-eig$vectors[,1:q]
        C<-Q%*%A%*%t(Q)
        return(list(mC=C,Q=Q,A=A,complete=TRUE))
      }
      emma.MLE0.c <- function(Y_c,W_c){
        n <- length(Y_c)
        stopifnot(nrow(W_c)==n)
        M_c<-diag(1,n)-W_c%*%solve(crossprod(W_c,W_c))%*%t(W_c)
        etas<-crossprod(M_c,Y_c)
        LL <- 0.5*n*(log(n/(2*pi))-1-log(sum(etas*etas)))
        return(list(ML=LL))
      }
      
      emma.REMLE0.c <- function(Y_c,W_c){
        n <- length(Y_c)
        stopifnot(nrow(W_c)==n)
        M_c <-diag(1,n)-W_c%*%solve(crossprod(W_c,W_c))%*%t(W_c)
        eig <-eigen(M_c)
        t <-qr(W_c)$rank
        v <-n-t
        U_R <-eig$vector[,1:v]
        etas<-crossprod(U_R,Y_c)
        LL <- 0.5*v*(log(v/(2*pi))-1-log(sum(etas*etas)))
        return(list(REML=LL))
      }
      
      replaceNaN<-  function(LL) {
        index=(LL=="NaN")
        if(length(index)>0) theMin=min(LL[!index])
        if(length(index)<1) theMin="NaN"
        LL[index]=theMin
        return(LL)    
      }
      
      lars <-  function(x, y, type = c("lasso", "lar", "forward.stagewise","stepwise"), trace = FALSE,
                        normalize=TRUE, intercept=TRUE, Gram, 
                        eps = .Machine$double.eps,  max.steps, use.Gram = TRUE)
      {
        
        call <- match.call()
        type <- match.arg(type)
        TYPE <- switch(type,
                       lasso = "LASSO",
                       lar = "LAR",
                       forward.stagewise = "Forward Stagewise",
                       stepwise = "Forward Stepwise")
        if(trace)
          cat(paste(TYPE, "sequence\n"))
        
        nm <- dim(x)
        n <- nm[1]
        m <- nm[2]
        im <- inactive <- seq(m)
        one <- rep(1, n)
        vn <- dimnames(x)[[2]]  
        ### Center x and y, and scale x, and save the means and sds
        if(intercept){
          meanx <- drop(one %*% x)/n
          x <- scale(x, meanx, FALSE)  # centers x
          mu <- mean(y)
          y <- drop(y - mu)
        }
        else {
          meanx <- rep(0,m)
          mu <- 0
          y <- drop(y)
        }
        if(normalize){
          normx <- sqrt(drop(one %*% (x^2)))
          nosignal<-normx/sqrt(n) < eps
          if(any(nosignal))# ignore variables with too small a variance
          {
            ignores<-im[nosignal]
            inactive<-im[-ignores]
            normx[nosignal]<-eps*sqrt(n)
            if(trace)
              cat("LARS Step 0 :\t", sum(nosignal), "Variables with Variance < eps; dropped for good\n")  #
          }
          else ignores <- NULL #singularities; augmented later as well
          names(normx) <- NULL
          x <- scale(x, FALSE, normx)	# scales x
        }
        else {
          normx <- rep(1,m)
          ignores <- NULL
        }
        if(use.Gram & missing(Gram)) {
          if(m > 500 && n < m)
            cat("There are more than 500 variables and n<m;\nYou may wish to restart and set use.Gram=FALSE\n"
            )
          if(trace)
            cat("Computing X'X .....\n")
          Gram <- t(x) %*% x	#Time saving
        }
        Cvec <- drop(t(y) %*% x)
        ssy <- sum(y^2)	### Some initializations
        residuals <- y
        if(missing(max.steps))
          max.steps <- 8*min(m, n-intercept)
        beta <- matrix(0, max.steps + 1, m)	# beta starts at 0
        lambda=double(max.steps)
        Gamrat <- NULL
        arc.length <- NULL
        R2 <- 1
        RSS <- ssy
        first.in <- integer(m)
        active <- NULL	# maintains active set
        actions <- as.list(seq(max.steps))	
        
        drops <- FALSE
        Sign <- NULL
        R <- NULL	###
        
        k <- 0
        while((k < max.steps) & (length(active) < min(m - length(ignores),n-intercept)) )
        {
          action <- NULL
          C <- Cvec[inactive]	#
          
          Cmax <- max(abs(C))
          if(Cmax<eps*100){ 
            if(trace)cat("Max |corr| = 0; exiting...\n")
            break
          }
          k <- k + 1
          lambda[k]=Cmax
          
          if(!any(drops)) {
            new <- abs(C) >= Cmax - eps
            C <- C[!new]	# for later
            new <- inactive[new]	# Get index numbers
            
            for(inew in new) {
              if(use.Gram) {
                R <- updateR(Gram[inew, inew], R, drop(Gram[
                  inew, active]), Gram = TRUE,eps=eps)
              }
              else {
                R <- updateR(x[, inew], R, x[, active], Gram
                             = FALSE,eps=eps)
              }
              if(attr(R, "rank") == length(active)) {
                
                nR <- seq(length(active))
                R <- R[nR, nR, drop = FALSE]
                attr(R, "rank") <- length(active)
                ignores <- c(ignores, inew)
                action <- c(action,  - inew)
                if(trace)
                  cat("LARS Step", k, ":\t Variable", inew, 
                      "\tcollinear; dropped for good\n")	#
              }
              else {
                if(first.in[inew] == 0)
                  first.in[inew] <- k
                active <- c(active, inew)
                Sign <- c(Sign, sign(Cvec[inew]))
                action <- c(action, inew)
                if(trace)
                  cat("LARS Step", k, ":\t Variable", inew, 
                      "\tadded\n")
              }
            }
          }
          else action <-  - dropid
          Gi1 <- backsolve(R, backsolvet(R, Sign))	
          
          dropouts<-NULL
          if(type == "forward.stagewise") {
            directions <- Gi1 * Sign
            if(!all(directions > 0)) {
              if(use.Gram) {
                nnls.object <- nnls.lars(active, Sign, R, 
                                         directions, Gram[active, active], trace = 
                                           trace, use.Gram = TRUE,eps=eps)
              }
              else {
                nnls.object <- nnls.lars(active, Sign, R, 
                                         directions, x[, active], trace = trace, 
                                         use.Gram = FALSE,eps=eps)
              }
              positive <- nnls.object$positive
              dropouts <-active[-positive]
              action <- c(action, -dropouts)
              active <- nnls.object$active
              Sign <- Sign[positive]
              Gi1 <- nnls.object$beta[positive] * Sign
              R <- nnls.object$R
              C <- Cvec[ - c(active, ignores)]
            }
          }
          A <- 1/sqrt(sum(Gi1 * Sign))
          w <- A * Gi1	# note that w has the right signs
          if(!use.Gram) u <- drop(x[, active, drop = FALSE] %*% w)	###
          
          if( (length(active) >=  min(n-intercept, m - length(ignores) ) )|type=="stepwise") {
            gamhat <- Cmax/A
          }
          else {
            if(use.Gram) {
              a <- drop(w %*% Gram[active,  - c(active,ignores), drop = FALSE])
            }
            else {
              a <- drop(u %*% x[,  - c(active, ignores), drop=FALSE])
            }
            gam <- c((Cmax - C)/(A - a), (Cmax + C)/(A + a))	
            
            gamhat <- min(gam[gam > eps], Cmax/A)	
          }
          if(type == "lasso") {
            dropid <- NULL
            b1 <- beta[k, active]	# beta starts at 0
            z1 <-  - b1/w
            zmin <- min(z1[z1 > eps], gamhat)
            if(zmin < gamhat) {
              gamhat <- zmin
              drops <- z1 == zmin
            }
            else drops <- FALSE
          }
          beta[k + 1,  ] <- beta[k,  ]
          beta[k + 1, active] <- beta[k + 1, active] + gamhat * w
          if(use.Gram) {
            Cvec <- Cvec - gamhat * Gram[, active, drop = FALSE] %*% w
          }
          else {
            residuals <- residuals - gamhat * u
            Cvec <- drop(t(residuals) %*% x)
          }
          Gamrat <- c(Gamrat, gamhat/(Cmax/A))
          arc.length <- c(arc.length, gamhat)	
          if(type == "lasso" && any(drops)) {
            dropid <- seq(drops)[drops]	
            for(id in rev(dropid)) {
              if(trace)
                cat("Lasso Step", k+1, ":\t Variable", active[
                  id], "\tdropped\n")
              R <- downdateR(R, id)
            }
            dropid <- active[drops]
            beta[k+1,dropid]<-0  
            active <- active[!drops]
            Sign <- Sign[!drops]
          }
          if(!is.null(vn))
            names(action) <- vn[abs(action)]
          actions[[k]] <- action
          inactive <- im[ - c(active, ignores)]
          if(type=="stepwise")Sign=Sign*0
        }
        beta <- beta[seq(k + 1), ,drop=FALSE ]	
        lambda=lambda[seq(k)]
        dimnames(beta) <- list(paste(0:k), vn)
        if(trace)
          cat("Computing residuals, RSS etc .....\n")
        residuals <- y - x %*% t(beta)
        beta <- scale(beta, FALSE, normx)
        RSS <- apply(residuals^2, 2, sum)
        R2 <- 1 - RSS/RSS[1]
        actions=actions[seq(k)]
        netdf=sapply(actions,function(x)sum(sign(x)))
        df=cumsum(netdf)### This takes into account drops
        if(intercept)df=c(Intercept=1,df+1)
        else df=c(Null=0,df)
        rss.big=rev(RSS)[1]
        df.big=n-rev(df)[1]
        if(rss.big<eps|df.big<eps)sigma2=NaN
        else
          sigma2=rss.big/df.big
        Cp <- RSS/sigma2 - n + 2 * df
        attr(Cp,"sigma2")=sigma2
        attr(Cp,"n")=n
        object <- list(call = call, type = TYPE, df=df, lambda=lambda,R2 = R2, RSS = RSS, Cp = Cp, 
                       actions = actions[seq(k)], entry = first.in, Gamrat = Gamrat, 
                       arc.length = arc.length, Gram = if(use.Gram) Gram else NULL, 
                       beta = beta, mu = mu, normx = normx, meanx = meanx)
        class(object) <- "lars"
        object
      }
      
      Y.data<-as.matrix(phe)
      if(is.null(psmatrix)==FALSE){
        psmatrix<-as.matrix(psmatrix)
      }
      nsam <-ncol(gene.data)-2
      chrnum<-length(unique(gene.data[,1]))
      
      W.orig<-matrix(1,nsam,1)
      if(is.null(psmatrix)==FALSE){
        W1 <-cbind(W.orig,psmatrix)
      }else{
        W1<-W.orig
      }
      
      kk<-list(NULL)
      cc<-list(NULL)
      kktotal<-matrix(0,nsam,nsam)
      
      for(i in 1:chrnum){
        xxot <- as.matrix(gene.data[gene.data[,1]==i,3:ncol(gene.data)])
        xot <-t(xxot)
        nmarkot <-ncol(xot)
        # kk[[i]]<-matrix(0,nsam,nsam)
        # for(k in 1:nmarkot)
        # {
        #   z<-as.matrix(xot[,k])
        #   kk[[i]]<-kk[[i]]+z%*%t(z)
        # }
        kk[[i]]<-mrMLM::multiplication_speed(xot,t(xot))
        cc[[i]]<-mean(diag(kk[[i]]))
        kktotal<-kktotal+kk[[i]]
      }
      
      rm(xot)
      gc()
      
      
      cl.cores <- detectCores()
      if((cl.cores<=2)||(is.null(CLO)==FALSE)){
        cl.cores<-1
      }else if(cl.cores>2){
        if(cl.cores>10){
          cl.cores<-10
        }else {  
          cl.cores <- detectCores()-1
        }
      }   
      
      cl <- makeCluster(cl.cores)
      registerDoParallel(cl)
      
      
      larsres<-foreach(i=1:chrnum,.multicombine=TRUE,.combine='rbind')%dopar%{
        
        requireNamespace("foreach")
        requireNamespace("lars")
        requireNamespace("sampling") 
        
        xxot <- as.matrix(gene.data[(gene.data[,1])!=i,3:ncol(gene.data)])
        
        xx1 <- as.matrix(gene.data[gene.data[,1]==i,3:ncol(gene.data)])
        YY1 <- matrix(Y.data,,1)
        
        K1 <- (kktotal-kk[[i]])/(sum(unlist(cc))-as.numeric(cc[i]))
        
        repl<-numeric()
        if(Bootstrap==TRUE){
          
          res1<-foreach(repl=1:5,.multicombine=TRUE,.combine='cbind')%do%{
            
            if(repl==1){
              YY<-YY1
              xx<-xx1
              K<-K1
              W<-W1
            }else{
              s<-srswr(nrow(YY1),nrow(YY1))
              ind<-(1:nrow(YY1))[s!=0]
              n<-s[s!=0]
              ind<-rep(ind,times=n)
              YY<-as.matrix(YY1[ind,])
              xx<-xx1[,ind]
              
              xxot2<-xxot[,ind]
              xot2 <-t(xxot2)
              nmarkot2 <-ncol(xot2)
              kk2<-matrix(0,nsam,nsam)
              for(k in 1:nmarkot2)
              {
                z2<-as.matrix(xot2[,k])
                kk2<-kk2+z2%*%t(z2)
              }
              cc2<-mean(diag(kk2))
              K <- numeric()
              K <- kk2/cc2
              W<-as.matrix(W1[ind,])
            }
            
            remle2<-emma.REMLE(YY, W, K, Z=NULL, ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL)
            remle1.B1<-emma.maineffects.B(Z=NULL,K,remle2$delta)
            rm(K)
            gc()
            C2<-remle1.B1$mC
            Y_c <- C2%*%YY
            W_c <- C2%*%W
            G_c <- C2%*%t(xx)
            GGG <- t(G_c)
            rm(G_c)
            gc()
            ylars <- as.matrix(Y_c)
            xlars <- cbind(W_c,t(GGG))
            rm(GGG)
            gc()
            LAR <- lars(xlars,ylars,type="lar",use.Gram=F,max.steps=lars1)
            rm(xlars)
            gc()
            LAR$beta[nrow(LAR$beta),]
          }
          
        }else if(Bootstrap==FALSE){
          res1 <- numeric()
          remle2<-emma.REMLE(YY1, W1, K1, Z=NULL, ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL)
          remle1.B1<-emma.maineffects.B(Z=NULL,K1,remle2$delta)
          rm(K1)
          gc()
          C2<-remle1.B1$mC
          Y_c <- C2%*%YY1
          W_c <- C2%*%W1
          G_c <- C2%*%t(xx1)
          rm(xx1)
          gc()
          GGG <- t(G_c)
          rm(G_c)
          gc()
          ylars <- as.matrix(Y_c)
          xlars <- cbind(W_c,t(GGG))
          rm(GGG)
          gc()
          LAR <- lars(xlars,ylars,type="lar",use.Gram=F,max.steps=lars1)
          rm(xlars)
          gc()
          res1<-cbind(res1,LAR$beta[nrow(LAR$beta),])
        }  
        
        if(is.null(psmatrix)==FALSE){
          rr <- as.matrix(res1[-c(1:(ncol(psmatrix)+1)),])
        }else{
          rr <- as.matrix(res1[-1,])
        }
      }
      
      stopCluster(cl)
      
      rm(kk,kktotal)
      gc()
      
     if(Bootstrap==TRUE){
        count <- matrix(rep(0,nrow(larsres)),nrow(larsres),1)
        
        ttt <- numeric()
        for(ii in 1:nrow(larsres))
        {
          tt <- 0
          for(jj in 1:ncol(larsres))
          {
            if ((abs(larsres[ii,jj]))>0){tt <- tt+1}
          }
          count[ii] <-tt
        }
        larsres <-cbind(larsres,count)
        
        for(ii in 1:nrow(larsres))
        {
          if(larsres[ii,ncol(larsres)]>=3){ttt <- cbind(ttt,ii)}
        }
        
        countnum <- ttt
        
      }else{
        
        countnum <- numeric()
        for(ii in 1:nrow(larsres))
        {
          if ((abs(larsres[ii]))>0){countnum <- cbind(countnum,ii)}
        }
      }
      if(ncol(countnum)>nrow(phe)){
        
        if(length(countnum)==1){
          xx2 <- matrix(gene.data[c(countnum),3:ncol(gene.data)],1,)
          
        }else{
          xx2 <- as.matrix(gene.data[c(countnum),3:ncol(gene.data)])
        }
        YY2 <- matrix(Y.data,,1)
        
        ylars <- as.matrix(YY2)
        xlars <- cbind(W1,t(xx2))
        LAR <- lars(xlars,ylars,type="lar",use.Gram=F)
        
        res1<-as.matrix(LAR$beta[nrow(LAR$beta),])
        
        rm(xlars,xx2)
        gc()
        
        if(is.null(psmatrix)==FALSE){
          rr <- as.matrix(res1[-c(1:(ncol(psmatrix)+1)),])
        }else{
          rr <- as.matrix(res1[-1,])
        }
        
        ct <- numeric()
        for(ii in 1:nrow(rr))
        {
          if ((abs(rr[ii]))>0){ct <- cbind(ct,ii)}
        }
        
        inct<-c(ct)
        countnum<-countnum[,inct]
        
      }
      
      if(length(countnum)==1){
        xeb <- matrix(gene.data[c(countnum),],1,)
        ebrow <-matrix(xeb[,1:2],,2)
        xeb1<-matrix(xeb[,3:ncol(xeb)],1,)
        xxeb <- as.matrix(t(xeb1))
        nmak <- ncol(xxeb)
        
      }else{
        xeb <- as.matrix(gene.data[c(countnum),])
        ebrow <-as.matrix(xeb[,1:2])
        xeb1<-as.matrix(xeb[,3:ncol(xeb)])
        xxeb <- as.matrix(t(xeb1))
        nmak <- ncol(xxeb)
      }
      bayeslodres <- numeric()
      
      genname<-gene.data[,1:2] 
      
      rm(xeb,gene.data)
      gc()
      
      yeb <- as.matrix(phe)
      
      if(is.null(psmatrix)==FALSE){
        u1<-ebayes_EM(cbind(matrix(1,nrow(xxeb),1),psmatrix),xxeb,yeb)
        xb<-u1$u
      }else{
        u1<-ebayes_EM(matrix(1,nrow(xxeb),1),xxeb,yeb)
        xb<-u1$u
      }
      xb<-as.matrix(xb)
      if(is.null(psmatrix)==FALSE){
        temp<-cbind(matrix(1,nrow(xxeb),1),psmatrix)
      }else{
        temp<-matrix(1,nrow(xxeb),1)
      }
      
      lodres<-likelihood(temp,xxeb,yeb,xb)
      lodres<-as.matrix(lodres)
      #### compute heredity#######
      ch_er <- as.numeric()
      ch_x <- cbind(matrix(1,nrow(xxeb),1),xxeb)
      
      ch_bb <- rbind(mean(yeb),as.matrix(xb))
      
      rm(xxeb)
      gc()
      
      
      for(i in 1:(ncol(ch_x)-1))
      {
        ch_xi <- ch_x[,(1+i)]
        as1 <- length(which(ch_xi==1))/nrow(ch_x)
        as2 <- 1-as1
        ch_er <- rbind(ch_er,(1-(as1-as2)*(as1-as2))*ch_bb[i+1]*ch_bb[i+1])
      }
      ch_v0 <- (1/(nrow(ch_x)-1))*(t(yeb-ch_x%*%ch_bb)%*%(yeb-ch_x%*%ch_bb))
      
      rm(ch_x)
      gc()
      
      
      if(var(yeb)>=sum(ch_er)+ch_v0){
        hered <- (ch_er/as.vector(var(yeb)))*100 
      }else{
        hered <- (ch_er/as.numeric(sum(ch_er)+ch_v0))*100
      }
      
      bayeslodres<-cbind(ebrow,xb,lodres,hered)
      
      
      lodid<-which(bayeslodres[,4]>lodvalue)
      if(length(lodid)!=0){
        
        if(length(lodid)==1){
          lastres<-matrix(bayeslodres[lodid,],1,)
          xeb2<-matrix(xeb1[lodid,],1,)
        }else{
          lastres<-bayeslodres[lodid,]
          xeb2<-as.matrix(xeb1[lodid,])
        }
        
        rm(xeb1)
        gc()
        
        xxmaf<- xeb2
        leng.maf<-dim(xxmaf)[2]
        
        maf.fun<-function(snp){
          leng<-length(snp)
          snp1<-length(which(snp==1))
          snp11<-length(which(snp==-1))
          snp0<-length(which(snp==0))
          ma1<-(2*snp1+snp0)/(2*leng)
          ma2<-(2*snp11+snp0)/(2*leng)
          maf<-min(ma1,ma2)
          return(maf)
        }
        
        maf<-apply(xxmaf,1,maf.fun)
        maf<-as.matrix(round(maf,4))
        vee<-round(u1$sigma2,4)
        pee<-round(var(yeb),4)
        
        vees<-matrix("",nrow = nrow(lastres),1)
        pees<-matrix("",nrow = nrow(lastres),1)
        pees[1,1]<-pee
        vees[1,1]<-vee
        result<-lastres
        result<-result
        
        if(nrow(result)>1){
          temp<-as.matrix(result[,3:5])
          temp[which(abs(temp)>=1e-4)]<-round(temp[abs(temp)>=1e-4],4)
          temp[which(abs(temp)<1e-4)]<-as.numeric(sprintf("%.4e",temp[abs(temp)<1e-4]))
          wan<-cbind(result[,1:2],temp)
        }else{
          temp<-t(as.matrix(result[,3:5]))
          temp[which(abs(temp)>=1e-4)]<-round(temp[abs(temp)>=1e-4],4)
          temp[which(abs(temp)<1e-4)]<-as.numeric(sprintf("%.4e",temp[abs(temp)<1e-4]))
          wan<-cbind(t(as.matrix(result[,1:2])),temp)  
        }
        
        if(inputform==1){
          genRaw<-as.data.frame(genRaw)
          genraw<-genRaw[-1,1:4]
          
          wan_len<-dim(wan)[1]
          marker<-character()
          snp<-character()
          
          for(i in 1:wan_len){
            chr_pos<-which(genraw[,2]==wan[i,1])
            new_matrix<-genraw[chr_pos,]
            posi_pos<-which(new_matrix[,3]==wan[i,2])
            mark<-matrix(new_matrix[posi_pos,1],1,)
            marker<-rbind(marker,mark)
            sn<-matrix(new_matrix[posi_pos,4],1,)
            snp<-rbind(snp,sn)
          }
        }
        if(inputform==2){
          
          genRaw<-as.data.frame(genRaw)
          genraw<-genRaw[-1,1:4]
          
          wan_len<-dim(wan)[1]
          marker<-character()
          snp<-character()
          for(i in 1:wan_len){
            chr_pos<-which(genraw[,2]==wan[i,1])
            new_matrix<-genraw[chr_pos,]
            posi_pos<-which(new_matrix[,3]==wan[i,2])
            mark<-matrix(new_matrix[posi_pos,1],1,)
            marker<-rbind(marker,mark)
            sn<-matrix(new_matrix[posi_pos,4],1,)
            snp<-rbind(snp,sn)
          }
        }
        if(inputform==3){
          genRaw<-as.data.frame(genRaw)
          genraw<-genRaw[-1,c(1,3,4,12)]
          
          wan_len<-dim(wan)[1]
          marker<-character()
          snp<-character()
          for(i in 1:wan_len){
            chr_pos<-which(genraw[,2]==wan[i,1])
            new_matrix<-genraw[chr_pos,]
            posi_pos<-which(new_matrix[,3]==wan[i,2])
            mark<-matrix(new_matrix[posi_pos,1],1,)
            marker<-rbind(marker,mark)
            sn<-matrix(new_matrix[posi_pos,4],1,)
            snp<-rbind(snp,sn)
          }
        }
        
        wan<-cbind(marker,wan,maf,snp,vees,pees)
        tempwan <- wan
        lodscore1 <- as.numeric(tempwan[,5])
        log10P <- as.matrix(round(-log10(pchisq(lodscore1*4.605,1,lower.tail = F)),4))
        if(nrow(tempwan)>1){
          tempwan1 <- cbind(tempwan[,1:5],log10P,tempwan[,6:10])
        }else{
          tempwan1 <- cbind(t(as.matrix(tempwan[,1:5])),log10P,t(as.matrix(tempwan[,6:10])))  
        }
        wan <- tempwan1
        colnames(wan)<-c("RS#","Chromosome","Marker position (bp)","QTN effect","LOD score","'-log10(P)'","r2 (%)","MAF","Genotype for code 1","Var_error","Var_phen (total)")
        wan<-as.data.frame(wan)
      }
      output<-list(result=wan)
    } 
    return(output) 
  }