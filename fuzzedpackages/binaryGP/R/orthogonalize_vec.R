orthogonalize_vec <- function(x, X, rho, corr){

  if (is.matrix(X) == FALSE) {
    X <- as.matrix(X)
  }
  d <- ncol(X)
  n <- nrow(X)
  LM <- matrix(0, nrow = 1, ncol = n)

  if(corr$type == "exponential" & corr$power == 2) { # closed form provided
    erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
    phi <- 10^(-rho/2)
    IM <- 2 * sqrt(pi) * phi * erf(2/phi) - phi^2 * (1 - exp(-4/phi^2))
    ILL <- phi^4 / 6 * (1 - exp(-4/phi^2)) - phi^2 / 3 * (3 - exp(-4/phi^2)) + 2 * sqrt(pi) / 3 * phi * erf(2/phi)
    Phi <- matrix(phi, nrow = n, ncol = d, byrow = TRUE)
    M <- sqrt(pi) / 2 * Phi * (erf((X+1)/Phi) - erf((X-1)/Phi))
    L <- Phi^2/2 * (exp(-(X+1)^2/Phi^2) - exp(-(X-1)^2/Phi^2)) + X * M
    M.x <- sqrt(pi) / 2 * Phi[1,,drop=FALSE] * (erf((x+1)/Phi[1,,drop=FALSE]) - erf((x-1)/Phi[1,,drop=FALSE]))
    L.x <- Phi[1,,drop=FALSE]^2/2 * (exp(-(x+1)^2/Phi[1,,drop=FALSE]^2) - exp(-(x-1)^2/Phi[1,,drop=FALSE]^2)) + x * M.x
    for(i in 1:d){
      L.tmp <- L[,i, drop = FALSE] / ILL[i]
      M.tmp <- matrix(apply(M[,-i, drop = FALSE], 1, prod) / prod(IM[-i]), ncol = 1)
      LM <- LM + (L.x[i] * t(L.tmp)) * (prod(M.x[-i]) * t(M.tmp))
    }
  }

  return(LM)
}
