

tvp.reg <- function(y,x,lambda=NULL,kappa=NULL,V=NULL,W=NULL)
  {
    if (is.null(colnames(x)))
      {
        colnames(x) <- colnames(x,do.NULL=FALSE,prefix="X")
      }
    
    x <- cbind(1,x)
    colnames(x)[1] <- "const"
    
    if (is.null(V)) { V <- 1 }
    if (is.null(W)) { W <- 1 }
    E <- diag(W,ncol(x))
    E.init <- E
    y.pred <- matrix(NA,ncol=1,nrow=nrow(y))
    thetas <- matrix(NA,ncol=ncol(x),nrow=nrow(y)+1)
    y.pred[1,] <- rep(0,ncol(y.pred))
    thetas[1,] <- rep(0,ncol(thetas))
    R.out <- matrix(NA,ncol=ncol(x),nrow=nrow(y)+1)
    R.out[1,] <- diag(E)
    V.out <- vector()
    V.out[1] <- V

    if (is.null(lambda) && is.null(kappa))
      {
        kft <- 2
      }
    else
      {
        kft <- 1
        if (is.null(lambda)) { lambda <- 1 }
      }

    for (t in 1:nrow(y))
      {
        if (kft==2)
          {
            kf <- .kalman2(y=as.numeric(y[t,,drop=FALSE]),x=x[t,,drop=FALSE],
                          theta=thetas[t,,drop=FALSE],R=E,Rw=E.init,Vv=V.out[1],t=t)
          }  
        else
          {
            kf <- .kalman(y=as.numeric(y[t,,drop=FALSE]),x=x[t,,drop=FALSE],
                          theta=thetas[t,,drop=FALSE],E=E,V=V,lambda=lambda,kappa=kappa,t=t)
          }
        y.pred[t,] <- kf$y.hat
        thetas[t+1,] <- kf$theta
        E <- kf$E
        V <- kf$V
        
        V.out[t+1] <- V
        R.out[t+1,] <- diag(E)
      }
    
    thetas <- thetas[-nrow(thetas),,drop=FALSE]
    colnames(thetas) <- colnames(x)
    colnames(R.out) <- colnames(x)
    R.out <- R.out[-nrow(R.out),,drop=FALSE]
    V.out <- V.out[-length(V.out)]
    out <- list(as.vector(y.pred),thetas,R.out,V.out)
    class(out) <- "tvpreg"
    names(out) <- c("y.hat","coef","R","V")
    return(out)
  }

