

### Nagy, I., Suzdaleva, E., 2013,
### Mixture estimation with state-space components and Markov model of switching,
### Applied Mathematical Modelling 37, 9970-9984
### DOI:10.1016/j.apm.2013.05.038


mixest1 <- function(y,x,mods=NULL,ftype=NULL,lambda=NULL,kappa=NULL,V=NULL,W=NULL,atype=NULL)
  {
    if (is.null(lambda) && is.null(kappa))
      {
        kft <- 2
      }
    else
      {
        kft <- 1
        if (is.null(lambda)) { lambda <- 1 }
      }
    
    if (is.null(ftype)) { ftype <- 0 }
    
    if (is.null(atype)) { atype <- 0 }
    
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

    x <- cbind(1,x)
    colnames(x)[1] <- "const"
    
    T <- nrow(y)
    nc <- nrow(mods)
    thetas <- list()
    Es <- list()
    Vs <- list()
    predpdfs <- matrix(0,nrow=1,ncol=nc)
    y.pred <- matrix(0,ncol=nc,nrow=T+1)
    
    w <- matrix(1/nc,nrow=T+1,ncol=nc)
    a <- matrix(1/nc,nrow=nc,ncol=nc)
    v <- matrix(1,nrow=nc,ncol=nc)
    theta.av <- matrix(0,nrow=T+1,ncol=ncol(x))
    theta.out <- theta.av
    R <- list()
    R[[1]] <- diag(W,ncol(x))

    for (i in 1:nc)
      {
        thetas[[i]] <- matrix(0,ncol=ncol(x),nrow=nrow(y)+1)
        Es[[i]] <- diag(W,ncol(x))
        Vs[[i]] <- V
      }
    R.out <- matrix(W,ncol=ncol(x),nrow=T+1)
    V.out <- matrix(V,ncol=1,nrow=T+1)

    for (t in 1:T)
      {
        for (i in 1:nc)
          {
            x.mod <- x
            x.mod[,which(mods[i,]==0)] <- 0
            if (kft==1)
              {
                kf <- .kalman(y=as.numeric(y[t,,drop=FALSE]),x=x.mod[t,,drop=FALSE],
                              theta=theta.av[t,,drop=FALSE],E=R[[t]],
                              V=Vs[[i]],lambda=lambda,kappa=kappa,t=t) 
              }
            else
              {
                kf <- .kalman2(y=as.numeric(y[t,,drop=FALSE]),x=x.mod[t,,drop=FALSE],
                              theta=theta.av[t,,drop=FALSE],R=R[[t]],
                              t=t,Rw=R[[1]],Vv=V) 
              }
            y.pred[t,i] <- kf$y.hat
            thetas[[i]][t+1,] <- kf$theta
            Es[[i]] <- kf$E
            Vs[[i]] <- kf$V
            predpdfs[1,i] <- kf$pdens
          }
       
        w.bar <- (t(predpdfs) %*% w[t,,drop=FALSE]) * a
        w.bar <- w.bar / sum(w.bar)
        w[t+1,] <- rowSums(w.bar)
        
        if (atype==0)
          {
            v <- v + w.bar
          }
        else
          {
            v <- mKIapprox(w.bar,v)
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
        theta.av[t+1,] <- w[t+1,,drop=FALSE] %*% thetas.mods
        
        R.temp <- matrix(0,nrow=ncol(x),ncol=ncol(x))
        for (i in 1:nc)
          {
            e.theta <- theta.av[t+1,,drop=FALSE] - thetas[[i]][t+1,,drop=FALSE]
            R.temp <- R.temp + as.numeric(w[t+1,i]) * Es[[i]] + as.numeric(w[t+1,i]) * (t(e.theta) %*% e.theta)
          }
        R[[t+1]] <- R.temp
        
        if (ftype==0 || ftype==1)
          {
            theta.out[t+1,] <- theta.av[t+1,]
            R.out[t+1,] <- diag(R[[t+1]])
            
            V.temp <- 0
            for (i in 1:nc)
              {
                V.temp <- V.temp + w[t+1,i] * Vs[[i]]
              }
            V.out[t+1,] <- V.temp
          }
        if (ftype==2)
          {
            theta.out[t+1,] <- thetas.mods[which.max(w[t+1,]),]
            R.out[t+1,] <- diag(Es[[which.max(w[t+1,])]])
            V.out[t+1,] <- Vs[[which.max(w[t+1,])]]
          }
        if (ftype==3)
          {
            j <- matrix(0,ncol=ncol(x),nrow=1)
            j.var <- as.vector(w[t+1,] %*% mods)
            j.var <- which(j.var >= 0.5)
            j[1,j.var] <- 1
            j.mod <- which(apply(mods,1,function(x) all(x == j[1,])))
            theta.out[t+1,] <- thetas.mods[j.mod,]
            R.out[t+1,] <- diag(Es[[j.mod]])
            V.out[t+1,] <- Vs[[j.mod]]
          }
      }

    for (i in 1:nc)
      {
        thetas[[i]] <- thetas[[i]][-nrow(thetas[[i]]),,drop=FALSE]
      }
    
    if (ftype==0)
      {
        y.hat <- as.numeric(diag(x %*% t(theta.av)))
        pip <- w %*% mods
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
    
    if (!ftype==0) { y.hat <- y.hat[-length(y.hat)] }
    w <- w[-nrow(w),,drop=FALSE]
    pip <- pip[-nrow(pip),,drop=FALSE]
    theta.out <- theta.out[-nrow(theta.out),,drop=FALSE]
    
    R.out <- R.out[-nrow(R.out),]
    V.out <- as.vector(V.out[-nrow(V.out),])

    colnames(pip) <- colnames(x)
    colnames(theta.out) <- colnames(x)
    colnames(R.out) <- colnames(x)
    rownames(theta.out) <- rownames(x)

    if (is.null(lambda)) { lambda <- NA }
    if (is.null(kappa)) { kappa <- NA }
    params <- c("state-space components",ftype,lambda,kappa,V,W,atype)
    names(params) <- c("mixture type","forecasting method","lambda","kappa","V0","W0","approximation method")
    colnames(mods) <- colnames(x)
    out <- list(y.hat,pip,theta.out,w,V.out,R.out,mods,params)
    names(out) <- c("y.hat","rvi","coef","weights","V","R","components","parameters")
    class(out) <- "mixest"
    
    return(out)
  }

