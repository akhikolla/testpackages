
### produces matrix indicating models with a constant and just one extra variable, 
### and a model with a constant only
### x should be a matrix
### the function is designed to simplify working with one-variable models in fDMA()

onevar <- function(x)
  {
    i <- ncol(x)
    m <- diag(1,i,i)
    m <- cbind(1,m)
    m <- rbind(0,m)
    m[1,1] <- 1
    return(m)
  }
