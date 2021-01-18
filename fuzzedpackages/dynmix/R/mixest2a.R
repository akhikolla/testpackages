

### Nagy, I., Suzdaleva, E., 2017, 
### Algorithms and Programs of Dynamic Mixture Estimation, Springer.


.mixest2a <- function(y,x,mods=NULL,ftype=NULL,V=NULL,W=NULL,Tvar=NULL)
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
    kappa <- matrix(ncol(x),nrow=T+1,ncol=nc)
    kappa[1,] <- rowSums(mods)
    Vs <- list()
    Rs <- vector()
    Cs <- list()

    for (i in 1:nc)
      {
        thetas[[i]] <- matrix(0,ncol=ncol(x),nrow=nrow(y)+1)
        Vs[[i]] <- diag(W^(-1),ncol(x)+1)
        Vs[[i]][1,1] <- V
        Rs[i] <- V
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
           
            yhat <- as.numeric(xt %*% t(thetas[[i]][t,,drop=FALSE]))
            e <- as.numeric(y[t,,drop=FALSE]) - yhat

            predpdfs[1,i] <- exp(-0.5 * e^2 / Rs[i]) / sqrt(2 * pi * Rs[i])
            y.pred[t,i] <- yhat
          }
          
        w.bar <- (t(predpdfs) %*% w[t,,drop=FALSE]) * a
        w.bar <- w.bar / sum(w.bar)
        w[t+1,] <- rowSums(w.bar)
           
        for (i in 1:nc)
          {
            x.mod <- x
            x.mod[,which(mods[i,]==0)] <- 0
            xt <- x.mod[t,,drop=FALSE]
            
            D <- cbind(y[t,,drop=FALSE],xt)
            D <- t(D) %*% D

            Vs[[i]] <- Vs[[i]] + w[t,i] * D
            kappa[t+1,i] <- kappa[t,i] + w[t,i]
          }
          
        v <- v + w.bar
            
        for (i in 1:nc)
          {
            V.temp <- Vs[[i]]
            Vy <- V.temp[1,1,drop=FALSE]
            Vypsi <- V.temp[2:nrow(V.temp),1,drop=FALSE]
            Vpsi <- V.temp[2:nrow(V.temp),2:ncol(V.temp),drop=FALSE]
            thetas[[i]][t+1,] <- ginv(Vpsi) %*% Vypsi
            if (t >= Tvar) 
              {
                R <- (Vy - (t(Vypsi) %*% ginv(Vpsi) %*% Vypsi)) / kappa[t+1,i]
                Cs.temp <- ltdl(V.temp)
                L.psi <- Cs.temp$L[-1,-1,drop=FALSE]
                D.psi <- Cs.temp$D[-1,-1,drop=FALSE]
                Cs[[i]] <- (ginv(L.psi) %*% ginv(D.psi) %*% t(ginv(L.psi))) / kappa[t+1,i]
              }
            else
              {
                R <- V
                Cs[[i]] <- diag(W,ncol(x))
              }
            Rs[i] <- as.numeric(R) 
          }
           
        a.bar <- colSums(v)
        for (i in 1:ncol(v))
          {
            a[,i] <- v[,i] / a.bar[i]
          }
           
        thetas.mods <- matrix(0,nrow=nc,ncol=ncol(x))
        for (i in 1:nc)
          {
            thetas.mods[i,] <- thetas[[i]][t+1,]
          }
        
        if (ftype==1)
          {
            theta.av[t+1,] <- w[t+1,,drop=FALSE] %*% thetas.mods
            V.out[t+1,1] <- sum(as.vector(w[t+1,,drop=FALSE]) * Rs)
            C.temp <- matrix(0,ncol(Cs[[1]]),ncol(Cs[[1]]))
            for (i in 1:nc)
              {
                C.temp <- C.temp + w[t+1,i] * Cs[[i]]
              }
            R.out[t+1,] <- diag(C.temp)
          }
        if (ftype==2)
          {
            theta.av[t+1,] <- thetas.mods[which.max(w[t+1,]),]
            V.out[t+1,1] <- Rs[which.max(w[t+1,])]
            R.out[t+1,] <- diag(Cs[[which.max(w[t+1,])]])
          }
        if (ftype==3)
          {
            j.var <- as.vector(w[t+1,] %*% mods)
            j.var <- which(j.var >= 0.5)
            j <- matrix(0,nrow=1,ncol=ncol(mods))
            j[1,j.var] <- 1
            j.mod <- which(apply(mods,1,function(x) all(x == j[1,])))            
            theta.av[t+1,] <- thetas.mods[j.mod,]
            V.out[t+1,1] <- Rs[j.mod]
            R.out[t+1,] <- diag(Cs[[j.mod]])
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
    theta.av <- theta.av[-nrow(theta.av),,drop=FALSE]
    
    V.out <- as.vector(V.out)
    V.out <- V.out[-length(V.out)]
    R.out <- R.out[-nrow(R.out),]

    colnames(pip) <- colnames(x)
    colnames(theta.av) <- colnames(x)
    colnames(R.out) <- colnames(x)
    rownames(theta.av) <- rownames(x)

    params <- c("normal regression components",ftype,V,W,Tvar,0)
    names(params) <- c("mixture type","forecasting method","V0","W0","Tvar","approximation method")
    colnames(mods) <- colnames(x)
    out <- list(y.hat,pip,theta.av,w,V.out,R.out,mods,params)
    names(out) <- c("y.hat","rvi","coef","weights","V","R","components","parameters")
    class(out) <- "mixest"
    
    return(out)
  }

  