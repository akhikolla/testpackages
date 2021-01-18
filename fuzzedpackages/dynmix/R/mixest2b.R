

### Nagy, I., Suzdaleva, E., Karny, M., Mlynarova, T., 2011,
### Bayesian estimation of dynamic finite mixtures,
### International Journal of Adaptive Control and Signal Processing 25, 765-787
### DOI:10.1002/acs.1239


.mixest2b <- function(y,x,mods=NULL,ftype=NULL,V=NULL,W=NULL,Tvar=NULL)
  {
    if (is.null(ftype)) { ftype <- 1 }
    
    if (is.null(colnames(x)))
      {
        colnames(x) <- colnames(x,do.NULL=FALSE,prefix="X")
      }
    
    if (is.null(mods))
      {
        mods <- expand.grid(rep.int(list(0:1),ncol(x)))
        mods <- as.matrix(cbind(rep.int(1,nrow(mods)),mods))
      }
    if (is.null(V)) { V <- 1 }
    if (is.null(W)) { W <- 1 }

    if (is.null(Tvar)) { Tvar <- 30 }

    x <- cbind(1,x)
    colnames(x)[1] <- "const"
    
    T <- nrow(y)
    nc <- nrow(mods)
    thetas <- list()
    predpdfs <- matrix(0,nrow=1,ncol=nc)
    y.pred <- matrix(0,ncol=nc,nrow=T+1)
    
    w <- matrix(1/nc,nrow=T+1,ncol=nc)
    a <- matrix(1/nc,nrow=nc,ncol=nc)
    v <- matrix(1,nrow=nc,ncol=nc)
    theta.av <- matrix(0,nrow=T+1,ncol=ncol(x))
    theta.out <- matrix(0,nrow=T+1,ncol=ncol(x))
    kappa <- matrix(ncol(x),nrow=T+1,ncol=nc)
    kappa[1,] <- rowSums(mods)
    Vs <- list()
    Rs <- vector()
    
    fd <- matrix(0,nrow=1,ncol=nc)
    D.old <- list()
    L.old <- list()
    D.new <- list()
    L.new <- list()
    C <- list()
    
    for (i in 1:nc)
      {
        thetas[[i]] <- matrix(0,ncol=ncol(x),nrow=nrow(y)+1)
        Vs[[i]] <- diag(W^(-1),ncol(x)+1)
        Vs[[i]][1,1] <- V
        Rs[i] <- V
        C[[i]] <- diag(W,ncol(x))
      }
    V.out <- matrix(V,nrow=T+1,ncol=1)
    R.out <- matrix(W,ncol=ncol(x),nrow=T+1)
    
    for (t in 1:T)
      {
        for (i in 1:nc)
          {
            x.mod <- x
            x.mod[,which(mods[i,]==0)] <- 0
            xt <- x.mod[t,,drop=FALSE]
            
            D <- cbind(y[t,,drop=FALSE],xt)
            D <- t(D) %*% D
            
            L.old[[i]] <- ltdl(Vs[[i]])$L
            D.old[[i]] <- ltdl(Vs[[i]])$D

            Vs[[i]] <- Vs[[i]] + D
            kappa[t+1,i] <- kappa[t,i] + 1
            
            L.new[[i]] <- ltdl(Vs[[i]])$L
            D.new[[i]] <- ltdl(Vs[[i]])$D
            
            thetas[[i]][t+1,] <- ginv(L.new[[i]][-1,-1,drop=FALSE]) %*% L.new[[i]][2:nrow(L.new[[i]]),1,drop=FALSE]
            
            C[[i]] <- ginv(L.new[[i]][-1,-1,drop=FALSE]) %*% ginv(D.new[[i]][-1,-1,drop=FALSE]) %*% t(ginv(L.new[[i]][-1,-1,drop=FALSE]))
            
            I.new <- gamma(0.5 * kappa[t+1,i]) * (D.new[[i]][1,1,drop=FALSE])^(-0.5 * kappa[t+1,i]) * (det(D.new[[i]][-1,-1,drop=FALSE]))^(-0.5) * 2^(0.5 * kappa[t+1,i]) * (2 * pi)^(0.5 * ncol(Vs[[i]]))
            I.old <- gamma(0.5 * kappa[t,i]) * (D.old[[i]][1,1,drop=FALSE])^(-0.5 * kappa[t,i]) * (det(D.old[[i]][-1,-1,drop=FALSE]))^(-0.5) * 2^(0.5 * kappa[t,i]) * (2 * pi)^(0.5 * ncol(Vs[[i]]))
            
            fd[1,i] <- I.new / I.old
          }

        v.bar <- colSums(v)
        for (i in 1:ncol(v))
          {
            a[,i] <- v[,i] / v.bar[i]
          }
          
        w.bar <- t(fd) %*% (w[t,,drop=FALSE] %*% t(a)) 
        w.bar <- w.bar / sum(w.bar)
        w[t+1,] <- rowSums(w.bar)

        v <- mKIapprox(w.bar,v)
       
        A <- 0
        B <- 0
        for (i in 1:nc)
          {
            A <- A + w[t+1,i] * kappa[t+1,i] / D.new[[i]][1,1]
            B <- B + w[t+1,i] * log(D.new[[i]][1,1]) - w[t+1,i] * digamma(0.5 * kappa[t+1,i])
          }
        B <- B + log(A)  
        
        kappa.old <- kappa[t+1,,drop=FALSE]
        for (i in 1:nc)
          {
            kappa[t+1,i] <- (1 + sqrt(1 + (4/3) * (B - log(2)))) / (2 * (B - log(2)))
          }
          
        D.d <- kappa[t+1,1] / A
       
        theta.temp <- matrix(0,nrow=1,ncol=ncol(x))
        for (i in 1:nc)
          {
            theta.temp <- theta.temp + w[t+1,i] * kappa.old[1,i] / as.numeric(D.new[[i]][1,1]) * thetas[[i]][t+1,,drop=FALSE]
          }
        theta.temp <- theta.temp / A 
        theta.av[t+1,] <- theta.temp
        
        if (t >= Tvar)
          {
            C.temp <- matrix(0,ncol=ncol(x),nrow=ncol(x))
            for (i in 1:nc)
              {
                thetas.d <- thetas[[i]][t+1,,drop=FALSE] - theta.av[t+1,,drop=FALSE]
                C.temp <- C.temp + w[t+1,i] * (C[[i]] + kappa.old[1,i] / as.numeric(D.new[[i]][1,1]) * t(thetas.d) %*% thetas.d)              
              }
          }
        else
          {
            C.temp <- diag(W,ncol(x))
          }
        
        D.psi.t <- ldlt(C.temp)
        D.psi <- ginv(D.psi.t$D)
        L.psi <- ginv(D.psi.t$L)
        L.dpsi <- L.psi %*% t(theta.av[t+1,,drop=FALSE])
            
        D <- matrix(0,ncol(Vs[[1]]),ncol(Vs[[1]]))
        D[1,1] <- D.d
        D[-1,-1] <- D.psi
            
        L <- diag(1,ncol(Vs[[1]]),ncol(Vs[[1]]))
        L[2:nrow(L),1] <- L.dpsi
        L[-1,-1] <- L.psi
            
        V.av <- t(L) %*% D %*% L

        for (i in 1:nc)
          {
            x.mod <- x
            x.mod[,which(mods[i,]==0)] <- 0
            xt <- x.mod[t,,drop=FALSE]
           
            yhat <- as.numeric(xt %*% t(thetas[[i]][t,,drop=FALSE]))
            y.pred[t,i] <- yhat
            
            if (t >= Tvar)
              {
                Vs[[i]] <- V.av
                Rs[i] <- D.d^2 * 2 / ((kappa[t+1,1] - 2)^2 * (kappa[t+1,1] - 4))
              }
            else
              {
                Vs[[i]] <- diag(W^(-1),ncol(x)+1)
                Vs[[i]][1,1] <- V
                Rs[i] <- V
                C[[i]] <- diag(W,ncol(x))
              }
          }
              
        thetas.mods <- matrix(0,nrow=nc,ncol=ncol(x))
        for (i in 1:nc)
          {
            thetas.mods[i,] <- thetas[[i]][t+1,]
          }
        
        if (ftype==1)
          {
            theta.out[t+1,] <- w[t+1,,drop=FALSE] %*% thetas.mods
            V.out[t+1,1] <- Rs[1]
            R.out[t+1,] <- diag(C.temp)
          }
        if (ftype==2)
          {
            theta.out[t+1,] <- thetas.mods[which.max(w[t+1,]),]
            V.out[t+1,1] <- Rs[which.max(w[t+1,])]
            R.out[t+1,] <- diag(C[[which.max(w[t+1,])]])
          }
        if (ftype==3)
          {
            j.var <- as.vector(w[t+1,] %*% mods)
            j.var <- which(j.var >= 0.5)
            j <- matrix(0,nrow=1,ncol=ncol(mods))
            j[1,j.var] <- 1
            j.mod <- which(apply(mods,1,function(x) all(x == j[1,])))            
            theta.out[t+1,] <- thetas.mods[j.mod,]
            V.out[t+1,1] <- Rs[j.mod]
            R.out[t+1,] <- diag(C[[j.mod]])
          }
      }

    if (ftype==1)
      {
        y.hat <- as.numeric(diag(y.pred %*% t(w)))
        pip <- w %*% mods
      }
    if (ftype==2)
      {
        y.hat <- vector()
        pip <- matrix(0,nrow=nrow(w),ncol=ncol(x))
        for (i in 1:nrow(w))
          {
            y.hat[i] <- y.pred[i,which.max(w[i,])]
            pip[i,] <- mods[which.max(w[i,]),]
          }
       }
    if (ftype==3)
      {
        y.hat <- vector()
        pip <- matrix(0,nrow=nrow(w),ncol=ncol(x))
        j <- matrix(0,ncol=ncol(x),nrow=nrow(w))
        for (i in 1:nrow(w))
          {
            j.var <- as.vector(w[i,] %*% mods)
            j.var <- which(j.var >= 0.5)
            j[i,j.var] <- 1
            j.mod <- which(apply(mods,1,function(x) all(x == j[i,])))
            y.hat[i] <- y.pred[i,j.mod]
            pip[i,] <- mods[j.mod,]
          }
       }
       
    y.hat <- y.hat[-length(y.hat)] 
    w <- w[-nrow(w),,drop=FALSE]
    pip <- pip[-nrow(pip),,drop=FALSE]
    theta.out <- theta.out[-nrow(theta.out),,drop=FALSE]
    
    V.out <- as.vector(V.out)
    V.out <- V.out[-length(V.out)]
    R.out <- R.out[-nrow(R.out),]

    colnames(pip) <- colnames(x)
    colnames(theta.out) <- colnames(x)
    colnames(R.out) <- colnames(x)
    rownames(theta.out) <- rownames(x)

    params <- c("normal regression components",ftype,V,W,Tvar,1)
    names(params) <- c("mixture type","forecasting method","V0","W0","Tvar","approximation method")
    colnames(mods) <- colnames(x)
    out <- list(y.hat,pip,theta.out,w,V.out,R.out,mods,params)
    names(out) <- c("y.hat","rvi","coef","weights","V","R","components","parameters")
    class(out) <- "mixest"
    
    return(out)
  }

  