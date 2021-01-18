

### Karny, M., Kadlec, J., Sutanto, E.L., 1998,
### Quasi-Bayes estimation applied to normal mixture,
### Preprints of The 3rd European IEEE Workshop on Computer-Intensive Methods in Control and Data Processing, 
### (eds.) Rojicek, J., Valeckova, M., Karny, M., Warwick K.,
### UTIA AV CR, Prague, 77-82  

### y   - data matrix of nrow(y) observations of ncol(y)-tuples
### m   - number of components (clusters), if not specified m = 2 is taken
### mu0 - initial means, should be a list of m matrices of dimension 1 x ncol(y),
###       if not specified random values are taken
### R0  - initial covariance matrices, should be a list of m matrices of dimension ncol(y) x ncol(y),
###       if not specified identity matrices are taken


qbnmix <- function(y,m=2,mu0=NULL,R0=NULL)
  {
    y <- as.matrix(y)
    
    T <- nrow(y)
    K <- ncol(y)
    
    if (is.null(mu0))
      {
        mu <- list()
        for (i in 1:m)
          {
            mu[[i]] <- matrix(rnorm(K),nrow=1,ncol=K)
          }
      }
    else
      {
        mu <- mu0
      }
    
    if (is.null(R0))
      {
        R <- list()
        for (i in 1:m)
          {
            R[[i]] <- diag(1,K)
          }
      }
    else
      {
        R <- R0
      }
    
    kappa <- rep.int(1/m,m)
    alpha <- kappa
    w.out <- kappa
    w <- vector()
    mu.out <- list()
    for (i in 1:m)
      {
        mu.out[[i]] <- mu[[i]]
      }
    g <- vector()
    ee <- list()
    
    for (t in 1:T)
      {
        yt <- y[t,,drop=FALSE]
        for (i in 1:m)
          {
            ee[[i]] <- yt - mu[[i]]
            ee[[i]] <- matrix(ee[[i]],ncol=1,nrow=K)
            w[i] <- ((det(2*pi*R[[i]]))^(-0.5))*exp(as.numeric(-0.5*t(ee[[i]])%*%solve(R[[i]])%*%ee[[i]]))*(kappa[i]/(t+1))
          }
        w <- w / sum(w) 
        for (i in 1:m)
          {
            kappa[i] <- kappa[i] + w[i]
          }
        for (i in 1:m)
          {
            g[i] <- as.numeric(w[i] / kappa[i])
          }
        for (i in 1:m)
          {
            mu[[i]] <- mu[[i]] + as.numeric(g[i])*t(ee[[i]]) 
            R[[i]] <- R[[i]] + as.numeric(g[i])*(ee[[i]]%*%t(ee[[i]])*(1-as.numeric(g[i])) - R[[i]])
            mu.out[[i]] <- rbind(mu.out[[i]],mu[[i]])
          }
        alpha.new <- kappa / sum(kappa)
        alpha <- rbind(alpha,alpha.new)
        w.out <- rbind(w.out,w)
      }
    
    alpha <- alpha[-nrow(alpha),,drop=FALSE] 
    rownames(alpha) <- seq(from=1,to=nrow(alpha),by=1)
    colnames(alpha) <- seq(from=1,to=m,by=1)
    w.out <- w.out[-nrow(w.out),,drop=FALSE] 
    rownames(w.out) <- seq(from=1,to=nrow(w.out),by=1)
    colnames(w.out) <- seq(from=1,to=m,by=1)
    for (i in 1:m)
      {
        mu.out[[i]] <- (mu.out[[i]])[-nrow(mu.out[[i]]),,drop=FALSE] 
        rownames(mu.out[[i]]) <- seq(from=1,to=nrow(mu.out[[i]]),by=1)
        colnames(mu.out[[i]]) <- seq(from=1,to=K,by=1)
      }
    
    out <- list(mu.out,R,alpha,w.out,mu0,R0)
    class(out) <- "qbnmix"
    names(out) <- c("mu","R","alpha","w","mu0","R0")
    return(out)
  }


### mu    - list of estimated means
### R     - list of estimated covariance matrices (from last step only)
### alpha - matrix of estimates of mixing weights (components columnwise)
### w     - matrix of posterior probabilities (components columnwise)
### mu0   - list of initial means matrices
### R0    - list of initial covaraince matrices

