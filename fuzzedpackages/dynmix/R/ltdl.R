

ltdl <- function(A)
  {
    if (! is.matrix(A)) { A <- as.matrix(A) }
    
    n <- ncol(A)
    
    L <- matrix(0,n,n)
    D <- L
    
    for (i in (seq(from=n,to=1,by=-1)))
      {
         D[i,i] <- A[i,i]
         L[i,1:i] <- A[i,1:i] / sqrt(as.numeric(A[i,i]))
         for (j in (1:(i-1)))
          {
            A[j,1:j] <- A[j,1:j] - L[i,1:j] * as.numeric(L[i,j])
          }
         L[i,1:i] <- L[i,1:i] / as.numeric(L[i,i])
      }
    rm(n,i,j)
    
    out <- list(L,D)
    names(out) <- c("L","D")
    return(out)
  }
  
