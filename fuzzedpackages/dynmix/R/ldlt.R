

ldlt <- function(A)
  {
    if (! is.matrix(A)) { A <- as.matrix(A) }
    
    n <- ncol(A)
    
    L <- diag(1,nrow=n,ncol=n)
    D <- matrix(0,nrow=n,ncol=1)
    
    for (i in 1:n)
      {
        D[i,1] <- A[i,i] - (L[i,1:(i-1),drop=FALSE])^2 %*% D[1:(i-1),1,drop=FALSE]
        if (i<n)
          {
            for (j in (i+1):n)
              {
                L[j,i] <- (A[j,i] - (L[j,1:(i-1),drop=FALSE] * L[i,1:(i-1),drop=FALSE]) %*% D[1:(i-1),1,drop=FALSE]) / D[i,1]
              }
          }
      }
    
    D <- diag(as.vector(D),ncol=n,nrow=n)
    
    out <- list(L,D)
    names(out) <- c("L","D")
    return(out)
  }
    
