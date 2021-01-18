pKWmEB<-function(gen,phe,outATCG,genRaw,kk,psmatrix,svpal,svmlod,Genformat,CLO){
  
    inputform<-Genformat
    pheRAW<-phe
    
    if(is.null(kk)){
      if(is.null(gen)==TRUE)
      {
        warning("Please input correct genotypic dataset !")
      }else{
        envgen<-gen[,3:ncol(gen)]
        envgen<-t(envgen)
        m<-ncol(envgen)
        n<-nrow(envgen)
        #kk1<-matrix(0,n,n)
        # for(k in 1:m){
        #   z<-as.matrix(envgen[,k])
        #   kk1<-kk1+z%*%t(z)
        # }
        kk1<-mrMLM::multiplication_speed(envgen,t(envgen))
        cc<-mean(diag(kk1))
        kk1<-kk1/cc
        kk<-as.matrix(kk1)
      }
      rm(envgen,kk1)
      gc()
      
    } 
    
    if(is.null(psmatrix)){
      flagps<-1
    }else{
      flagps<-0
    }
    
    if((flagps==1)||(exists("psmatrix")==FALSE))
    {
      phe<-phe
    }else if(flagps==0)
    {
      phe<-phe
      fixps <- cbind(matrix(1,nrow(phe),1),psmatrix)
      cui<-det(t(fixps)%*%fixps)
      p1<-rep(1,ncol(fixps))
      p2<-diag(p1)
      if (cui<1e-6){bbps<-solve(t(fixps)%*%fixps+p2*0.01)%*%t(fixps)%*%phe}
      if (cui>=1e-6){ bbps<-solve(t(fixps)%*%fixps)%*%t(fixps)%*%phe }
      bbps <- bbps[2:(nrow(bbps)),1]
      phe <- as.matrix(phe) - as.matrix(psmatrix)%*%as.matrix(bbps)
    }
    
    
    if(is.null(svpal)==TRUE||is.null(svmlod)==TRUE){
      warning("Please set parameter!")
    }
    
    if((svpal<0)||(svpal>1))
    {
      warning("Please input critical P-value between 0 and 1!")
    }
    if(svmlod<0)
    {
      warning("Please input critical LOD score: > 0 !")
    }
    
    if(exists("gen")==FALSE)
    {
      warning("Please input correct genotypic dataset !")
    }
    if(exists("phe")==FALSE)
    {
      warning("Please input correct phenotypic dataset !")
    }
    if(exists("kk")==FALSE)
    {
      warning("Please input correct kinship (K) dataset !")
    }
    if((exists("gen")==TRUE)&&(exists("phe")==TRUE)&&(ncol(gen)!=(nrow(phe)+2)))
    {
      warning("Sample size in genotypic dataset doesn't equal to the sample size in phenotypic dataset !")
    }
    
    if((exists("gen")==TRUE)&&(exists("phe")==TRUE)&&(exists("kk")==TRUE)&&((ncol(gen)==(nrow(phe)+2)))&&(svpal>=0)&&(svpal<=1)&&(svmlod>=0))
    {
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
      
      parmsShow<-NULL
      wan<-NULL
      parms.pchange<-NULL
      parmsm<-NULL
      
      K.data <- kk
      Y.data <- phe
      rawgen <- gen
      rawphe <- Y.data

      gene.data<-rawgen[,3:ncol(rawgen)]
      nsample <- ncol(gene.data)
      
      fix <- matrix(1,nsample,1)
      sam <- nsample
      Y.data <- matrix(Y.data,nsample,1)
      n<-dim(Y.data)[1]
      W.orig<-matrix(1,n,1)
      W <- W.orig
      K <- K.data
      YY <- Y.data
      rm(K.data)
      gc()
      p_value <- svpal
      ffpptotal <- numeric()
      gglartotal <- numeric()
      pvaluetotal <- numeric()
      
      #for(ii in 1:1){
      ii<-1
      remle2<-emma.REMLE(YY[,ii], W, K, Z=NULL, ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL)
      remle1.B1<-emma.maineffects.B(Z=NULL,K,remle2$delta)
      C2<-remle1.B1$mC
      
      rm(K,remle1.B1)
      gc()
      
      Y_c <- C2%*%YY[,ii]
      W_c <- C2%*%W
      G_c <- C2%*%t(gene.data)
      
      GGG <- t(G_c)
      
      rm(C2,G_c)
      gc()
      
      allrowmean <- rowMeans(GGG)
      nnG <- nrow(GGG)
      
      for(jj in 1:nnG)
      {
        GGG[jj,which(GGG[jj,]>=allrowmean[jj])] <- 1
        GGG[jj,which(GGG[jj,]<allrowmean[jj])] <- -1
      }
      
      gentran <- GGG
      
      rm(GGG)
      gc()
      
      phetran <- Y_c
      nn <- dim(gentran)[1]
      bb<-numeric()
      cc <- numeric()
      ff <- numeric()
      
      newphe <- cbind(matrix(c(1:sam),,1),phetran)
      ph <- unique(newphe[,2])
      newph <- newphe[match(ph,newphe[,2],0L),]
      newy <- newph[,2]
      sob <- newph[,1]
      
      ff<- foreach(i=1:nn)%do%
      {
        temp <- as.matrix(gentran[i,sob])
        temp<-factor(temp)
        loc<-which(as.numeric(levels(temp))==1)
        
      }
      
      fff<-unlist(ff)
      sameloc<-which(fff==1)
      
      if(length(sameloc)!=0){
        gentran1<-gentran[-c(sameloc),]
      }else if(length(sameloc)==0){
        gentran1<-gentran 
      }
      
      nnn<-dim(gentran1)[1]
      
      rm(gentran)
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
      
      unsameloc=foreach(i=1:nnn, .combine = 'rbind')%dopar%
      {
        requireNamespace("coin")
        requireNamespace("lars")
        temp <- as.matrix(gentran1[i,sob])
        xy <- cbind(temp,newy)
        b <- unique(xy[,1])
        
        temp <- factor(temp)
        snp <- data.frame(newy,temp)
        kw <- kruskal_test(newy~temp, data = snp,distribution = "asymptotic")
        kw <- pvalue(kw)
        aa <- kw[1]
      }
      
      stopCluster(cl)
      a<-matrix(0,nrow = nn,ncol=1)
      a[c(sameloc)]<-1
      a[which(a[]==0)]<-unsameloc
      bb<-a
      
      rm(gentran1)
      gc()
      
      kk <- matrix(seq(1:nn),nn,1)
      bb <- matrix(bb,nn,1)
      cc <- cbind(ii,kk,bb)
      pvaluetotal <- cc[,2:3]
      ff <- cc[which(cc[,3] < p_value),]
      ffpptotal <- ff
      pvaluetotal <- pvaluetotal
      
      ############lars###########################
      gg <- numeric()
      nchoice <- ff[,2]   
      genchoice <- gene.data[nchoice,]
      newpheno <- as.matrix(rawphe[((ii-1)*sam+1):(ii*sam),1])
      
      aall <- lars(t(genchoice),newpheno,type="lar",use.Gram=FALSE)
      bb2 <- aall$beta[nrow(aall$beta),]
      var <- unlist(aall[[8]])
      
      tempnn <- dim(ff)[1]
      if(tempnn<=150)
      {
        if(tempnn>=nsample)
        {
          tempnn <- nsample - 1
        }else if(tempnn <nsample)
        {
          tempnn <- dim(ff)[1]
        }
        var1 <- var[1:tempnn]
        bb2 <- bb2[abs(var1)]
        gg <- as.matrix(nchoice[abs(var1)])
        ############Empirical Bayes##################
        ggbayes <- numeric()
        optloci <- gg
        optgen <- gene.data[c(optloci),]
        newphebayes <- as.matrix(rawphe[((ii-1)*sam+1):(ii*sam),1])
        bbeff <- ebayes_EM(fix,t(optgen),newphebayes)
        lod <- likelihood(fix,t(optgen),newphebayes,bbeff$u)
        optlod <- which(lod>svmlod)
        if(length(optlod)>0){
          locich <- optloci[optlod]
          ggbayes <- cbind(ii,locich,matrix(rawgen[locich,1:2],,2),bbeff$u[optlod],lod[optlod],bbeff$sigma2)
        }
        gglartotal <- ggbayes
        rm(rawgen)
        gc()
        
        
      }else if((tempnn > 150)&&(nsample > 150))
      {
        if(tempnn>=nsample)
        {
          tempnn <- nsample - 1
        }else if(tempnn <nsample)
        {
          tempnn <- dim(ff)[1]
        }
        var1 <- var[1:tempnn]
        bb2 <- bb2[abs(var1)]
        gg <- as.matrix(nchoice[abs(var1)])
        
        aic <- numeric()
        hhbayes50 <- numeric()
        hhbayes100 <- numeric()
        hhbayes150 <- numeric()
        
        ggbayes <- numeric()
        ggbayes50 <- numeric()
        ggbayes100 <- numeric()
        ggbayes150 <- numeric()
        
        optloci <- gg
        optloci50 <- as.matrix(optloci[1:50])
        optloci100 <- as.matrix(optloci[1:100])
        optloci150 <- as.matrix(optloci[1:150])
        
        ##################choose 50 number variable from lars######################
        optgen50 <- gene.data[c(optloci50),]
        phebayes <- as.matrix(rawphe[((ii-1)*sam+1):(ii*sam),1])
        bbeff50 <- ebayes_EM(fix,t(optgen50),phebayes)
        lod50 <- likelihood(fix,t(optgen50),phebayes,bbeff50$u)
        
        optlod50 <- which(lod50>svmlod)
        if(length(optlod50)>0){
          locich50 <- optloci50[c(optlod50)]
          ggbayes50 <- cbind(ii,locich50,matrix(rawgen[locich50,1:2],,2),bbeff50$u[optlod50],lod50[optlod50],bbeff50$sigma2)
          hhbayes50 <- rbind(hhbayes50,ggbayes50)  
        }
        
        ##################choose 100 number variable from lars#####################
        optgen100 <- gene.data[c(optloci100),]
        phebayes <- as.matrix(rawphe[((ii-1)*sam+1):(ii*sam),1])
        bbeff100 <- ebayes_EM(fix,t(optgen100),phebayes)
        lod100 <- likelihood(fix,t(optgen100),phebayes,bbeff100$u)
        
        optlod100 <- which(lod100>svmlod)
        if(length(optlod100)>0){
          locich100 <- optloci100[optlod100]
          ggbayes100 <- cbind(ii,locich100,matrix(rawgen[locich100,1:2],,2),bbeff100$u[optlod100],lod100[optlod100],bbeff100$sigma2)
          hhbayes100 <- rbind(hhbayes100,ggbayes100)  
        }
        
        ##################choose 150 number variable from lars#####################
        optgen150 <- gene.data[c(optloci150),]
        phebayes <- as.matrix(rawphe[((ii-1)*sam+1):(ii*sam),1])
        bbeff150 <- ebayes_EM(fix,t(optgen150),phebayes)
        lod150 <- likelihood(fix,t(optgen150),phebayes,bbeff150$u)
        
        optlod150 <- which(lod150>svmlod)
        if(length(optlod150)>0){
          locich150 <- optloci150[optlod150]
          ggbayes150 <- cbind(ii,locich150,matrix(rawgen[locich150,1:2],,2),bbeff150$u[optlod150],lod150[optlod150],bbeff150$sigma2)
          hhbayes150 <- rbind(hhbayes150,ggbayes150)  
        }
        
        rm(rawgen)
        gc()
        
        ####################################AIC#####################################
        if(length(optlod50)==0)
        {
          lmres1 <- lm(phebayes~fix)
          aic1 <- AIC(lmres1)
        }
        if(length(optlod100)==0)
        {
          lmres2 <- lm(phebayes~fix)
          aic2 <- AIC(lmres2)
        }
        if(length(optlod150)==0)
        {
          lmres3 <- lm(phebayes~fix)
          aic3 <- AIC(lmres3)
        }
        
        if(length(optlod50)==1)
        {
          xx1 <- as.matrix(gene.data[ggbayes50[,2],])
          lmres1 <- lm(phebayes~xx1)
          aic1 <- AIC(lmres1)
        }
        if(length(optlod100)==1)
        {
          xx2 <- as.matrix(gene.data[ggbayes100[,2],])
          lmres2 <- lm(phebayes~xx2)
          aic2 <- AIC(lmres2)
        }
        if(length(optlod150)==1)
        {
          xx3 <- as.matrix(gene.data[ggbayes150[,2],])
          lmres3 <- lm(phebayes~xx3)
          aic3 <- AIC(lmres3)
        }
        
        if(length(optlod50)>1)
        {
          xx1 <- t(gene.data[ggbayes50[,2],])
          lmres1 <- lm(phebayes~xx1)
          aic1 <- AIC(lmres1)
        }
        if(length(optlod100)>1)
        {
          xx2 <- t(gene.data[ggbayes100[,2],])
          lmres2 <- lm(phebayes~xx2)
          aic2 <- AIC(lmres2)
        }
        if(length(optlod150)>1)
        {
          xx3 <- t(gene.data[ggbayes150[,2],])
          lmres3 <- lm(phebayes~xx3)
          aic3 <- AIC(lmres3)
        }
        
        aic <- rbind(aic,matrix(c(ii,aic1,aic2,aic3),1,4))
        
        ############################################################################
        if(aic1==min(aic1,aic2,aic3))
        {
          ggbayes <- ggbayes50
        }else if(aic2==min(aic1,aic2,aic3)){
          ggbayes <- ggbayes100
        }else if(aic3==min(aic1,aic2,aic3)){
          ggbayes <- ggbayes150
        }
        gglartotal <- ggbayes
      }
      #}
      
      gglartotal <- gglartotal
      
      if(inputform==1){
        #output result1 using mrMLM numeric format
        parmsShow<-as.matrix(-log10(pvaluetotal[,2]))
        tempparms<-parmsShow
        tempparms[which(abs(tempparms)>=1e-4)]<-round(tempparms[which(abs(tempparms)>=1e-4)],4)
        tempparms[which(abs(tempparms)<1e-4)]<-as.numeric(sprintf("%.4e",tempparms[which(abs(tempparms)<1e-4)]))
        kong<-matrix("",nrow(tempparms),1)
        parmsShow<-data.frame(genRaw[-1,1],gen[,1:2],kong,tempparms,genRaw[-1,4])
        colnames(parmsShow)<-c("RS#","Chromosome","Marker position (bp)","SNP effect (pKWmEB)","'-log10(P) (pKWmEB)'","Genotype for code 1")
        
      }
      if(inputform==2){
        #output result1 using mrMLM character format
        parmsShow<-as.matrix(-log10(pvaluetotal[,2]))
        outATCG<-matrix(outATCG,,1)
        tempparms<-parmsShow
        tempparms[which(abs(tempparms)>=1e-4)]<-round(tempparms[which(abs(tempparms)>=1e-4)],4)
        tempparms[which(abs(tempparms)<1e-4)]<-as.numeric(sprintf("%.4e",tempparms[which(abs(tempparms)<1e-4)]))
        kong<-matrix("",nrow(tempparms),1)
        parmsShow<-data.frame(genRaw[-1,1],gen[,1:2],kong,tempparms,outATCG)
        colnames(parmsShow)<-c("RS#","Chromosome","Marker position (bp)","SNP effect (pKWmEB)","'-log10(P) (pKWmEB)'","Genotype for code 1")
        
      }
      if(inputform==3){
        #output result1 using TASSEL format
        parmsShow<-as.matrix(-log10(pvaluetotal[,2]))
        outATCG<-matrix(outATCG,,1)
        #outATCG<-unlist(strsplit(outATCG,""))
        #outATCG<-matrix(outATCG[c(TRUE,FALSE)],,1)
        tempparms<-parmsShow
        tempparms[which(abs(tempparms)>=1e-4)]<-round(tempparms[which(abs(tempparms)>=1e-4)],4)
        tempparms[which(abs(tempparms)<1e-4)]<-as.numeric(sprintf("%.4e",tempparms[which(abs(tempparms)<1e-4)]))
        kong<-matrix("",nrow(tempparms),1)
        parmsShow<-data.frame(genRaw[-1,1],gen[,1:2],kong,tempparms,outATCG)
        colnames(parmsShow)<-c("RS#","Chromosome","Marker position (bp)","SNP effect (pKWmEB)","'-log10(P) (pKWmEB)'","Genotype for code 1")
      }
      
      finalres <- gglartotal
      
      if(length(finalres)!=0){
      
      if(length(finalres[,2])>1){
        
        if((flagps==1)||(exists("psmatrix")==FALSE))
        {
          ex<-cbind(fix,t(gene.data[finalres[,2],]))
        }else if(flagps==0)
        {
          ex<-cbind(cbind(fix,psmatrix),t(gene.data[finalres[,2],]))
        }
        
      }else{
        
        if((flagps==1)||(exists("psmatrix")==FALSE))
        {
          ex<-cbind(fix,as.matrix(gene.data[finalres[,2],]))
        }else if(flagps==0)
        {
          ex<-cbind(cbind(fix,psmatrix),as.matrix(gene.data[finalres[,2],]))
        }  
      }
      
      ex<-as.matrix(ex)
      cui<-det(t(ex)%*%ex)
      p1<-rep(1,ncol(ex))
      p2<-diag(p1)
      if (cui<1e-6){bbbb<-solve(t(ex)%*%ex+p2*0.01)%*%t(ex)%*%phe}
      if (cui>=1e-6){ bbbb<-solve(t(ex)%*%ex)%*%t(ex)%*%phe }
      if((flagps==1)||(exists("psmatrix")==FALSE))
      {
        eeff<-bbbb[2:(nrow(bbbb)),1]
      }else if(flagps==0)
      {
        eeff<-bbbb[(2+ncol(psmatrix)):(nrow(bbbb)),1]
      }
      
      eeff<-as.matrix(eeff)
      er<-as.numeric()
      her<-as.numeric()
      if((flagps==1)||(exists("psmatrix")==FALSE))
      {
        excol<-ncol(ex)
        for(i in 1:(excol-1))
        {
          em<-ex[,(1+i)]
          as1<-length(which(em==1))/nrow(ex)
          as2<-1-as1
          er<-rbind(er,(1-(as1-as2)*(as1-as2))*eeff[i]*eeff[i])
        }
        v0<-(1/(nrow(ex)-1))*(t(phe-ex%*%bbbb)%*%(phe-ex%*%bbbb))
        
        if(var(phe)>=sum(er)+v0){
          her<-(er/as.vector(var(phe)))*100 
        }else{
          
          her<-(er/as.numeric(sum(er)+v0))*100 
        }
        
      }else if(flagps==0)
      {
        excol<-ncol(ex)
        for(i in 1:(excol-1-ncol(psmatrix)))
        {
          em<-ex[,(1+ncol(psmatrix)+i)]
          as1<-length(which(em==1))/nrow(ex)
          as2<-1-as1
          er<-rbind(er,(1-(as1-as2)*(as1-as2))*eeff[i]*eeff[i])
        }
        v0<-(1/(nrow(ex)-1))*(t(phe-ex%*%bbbb)%*%(phe-ex%*%bbbb))
        
        if(var(phe)>=sum(er)+v0){
          her<-(er/as.vector(var(phe)))*100 
        }else{
          
          her<-(er/as.numeric(sum(er)+v0))*100 
        } 
      }
      
      X1<-t(gene.data[,])[3:ncol(gene.data),]
      rm(gene.data)
      gc()
  
      xxxx<-as.matrix(X1[,finalres[,2]])
      
      rm(X1)
      gc()
      
      xxmaf<-t(xxxx)
      
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
      
      eeff <- finalres[,5]
      lo <- finalres[,6]
      eeff[which(abs(eeff)>=1e-4)] <- round(eeff[which(abs(eeff)>=1e-4)],4)
      eeff[which(abs(eeff)<1e-4)] <- as.numeric(sprintf("%.4e",eeff[which(abs(eeff)<1e-4)]))
      lo[which(abs(lo)>=1e-4)] <- round(lo[which(abs(lo)>=1e-4)],4)
      lo[which(abs(lo)<1e-4)] <- as.numeric(sprintf("%.4e",lo[which(abs(lo)<1e-4)]))
      her[which(abs(her)>=1e-4)] <- round(her[which(abs(her)>=1e-4)],4)
      her[which(abs(her)<1e-4)] <- as.numeric(sprintf("%.4e",her[which(abs(her)<1e-4)]))
      needrs <- genRaw[-1,1]
      needrs <- as.matrix(needrs[finalres[,2]])
      needgenofor <- as.character()
      if(inputform==1)
      {
        needgenofor <- genRaw[-1,4]
        needgenofor <- as.matrix(needgenofor[finalres[,2]])
      }
      if(inputform==2)
      {
        needgenofor <- outATCG
        needgenofor <- as.matrix(needgenofor[finalres[,2]])
      }
      if(inputform==3)
      {
        needgenofor <- outATCG
        needgenofor <- as.matrix(needgenofor[finalres[,2]])
      }
      
      phevartotal<-var(pheRAW)
      if(finalres[1,7]>=1e-4){finalres[1,7]<-round(finalres[1,7],4)}
      if(finalres[1,7]<1e-4){finalres[1,7]<-as.numeric(sprintf("%.4e",finalres[1,7]))}
      if(phevartotal>=1e-4){phevartotal<-round(phevartotal,4)}
      if(phevartotal<1e-4){phevartotal<-as.numeric(sprintf("%.4e",phevartotal))}
      tempvar <- dim(as.matrix(lo))[1]
      if(tempvar==1)
      {
        wan<-data.frame(needrs,t(as.matrix(gen[finalres[,2],1:2])),as.matrix(eeff),as.matrix(lo),her,maf,needgenofor,as.matrix(finalres[,7]),phevartotal) 
      }else if(tempvar>1)
      {
        wan<-data.frame(needrs,gen[finalres[,2],1:2],eeff,lo,her,maf,needgenofor) 
        wan<-wan[order(wan[,2]),]
        wan<-data.frame(wan,rbind(finalres[1,7],as.matrix(rep("",(tempvar-1)))),rbind(phevartotal,as.matrix(rep("",(tempvar-1)))))
      }
      
      tempwan <- wan
      lodscore1 <- as.numeric(tempwan[,5])
      log10P <- as.matrix(round(-log10(pchisq(lodscore1*4.605,1,lower.tail = F)),4))
      tempwan1 <- cbind(tempwan[,1:5],log10P,tempwan[,6:10])
      wan <- tempwan1
      
      colnames(wan)<-c("RS#","Chromosome","Marker position (bp)","QTN effect","LOD score","'-log10(P)'","r2 (%)","MAF","Genotype for code 1","Var_error","Var_phen(total)")
      wan<-as.data.frame(wan)
      }#change20190125
      output<-list(result1=parmsShow,result2=wan)
      return(output)
    } 
   }