

sqrtmat <- function(A)
  {
    if (ncol(A)>1)
      {
        r <- eigen(A)
        V <- r$vectors
        D <- diag(r$values)
        D <- sqrt(D)
    
        out <- V %*% D %*% t(V)
      }
    else
      {
        out <- sqrt(A)
      }
    
    return(out)
  }

