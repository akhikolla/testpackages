lvglasso_calc <- function(Sigmahat,lambda1,lambda2,convergence=1e-10,maxiter=1000,rho=2.5){

  p <- nrow(Sigmahat)

  # Variables Initialization
  oldR <- R <- S <- L <- diag(rep(1,p))
  Gamma <- matrix(0,p,p)
  criteria <- 1e10
  i <- 1
  # While loop for the iterations
  while(criteria > convergence && i <= maxiter){

    R <- updateR(S,L,Sigmahat,Gamma,rho)

    S <- updateS(R,L,Gamma,lambda1,rho)

    L <- updateL(R,S,Gamma,lambda2,rho)

    Gamma <- Gamma+rho*(R-S+L)

    R <- S-L
    aa <- sum((oldR)^2) + 1e-7
    criteria <- sum((R-oldR)^2)/aa
    #criteria <- sum((R-oldR)^2)/sum((oldR)^2)
    oldR <- R
    i <- i+1

  }
  R <- S-L

  if(i>maxiter){
    warning("The algorithm has not converged by the specified maximum number of iteration")
  }

  return(list(R=R,S=S,L=L,Gamma=Gamma,iteration=i))

}
