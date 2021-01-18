orthogonalize <- function(X, rho, corr){

  if (is.matrix(X) == FALSE) {
    X <- as.matrix(X)
  }
  d <- ncol(X)
  n <- nrow(X)
  LM <- matrix(0, n, n)

  if(corr$type == "exponential" & corr$power == 2) { # closed form provided
    erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
    phi <- 10^(-rho/2)
    IM <- 2 * sqrt(pi) * phi * erf(2/phi) - phi^2 * (1 - exp(-4/phi^2))
    ILL <- phi^4 / 6 * (1 - exp(-4/phi^2)) - phi^2 / 3 * (3 - exp(-4/phi^2)) + 2 * sqrt(pi) / 3 * phi * erf(2/phi)
    Phi <- matrix(phi, nrow = n, ncol = d, byrow = TRUE)
    M <- sqrt(pi) / 2 * Phi * (erf((X+1)/Phi) - erf((X-1)/Phi))
    L <- Phi^2/2 * (exp(-(X+1)^2/Phi^2) - exp(-(X-1)^2/Phi^2)) + X * M
    for(i in 1:d){
      L.tmp <- L[,i, drop = FALSE] / sqrt(ILL[i])
      M.tmp <- matrix(apply(M[,-i, drop = FALSE], 1, prod) / sqrt(prod(IM[-i])), ncol = 1)
      LM <- LM + (L.tmp %*% t(L.tmp)) * (M.tmp %*% t(M.tmp))
    }
  }

  return(LM)
}
