\name{qrfit.lasso}
\alias{qrfit.lasso}
\title{ Quantile Regression (qr) with Adaptive Lasso Penalty (lasso)}
\description{
   High level function for estimating and selecting parameter by quantile regression with adaptive lasso penalty.
}
\usage{
qrfit.lasso(X,y,tau,lambda,beta,method,maxit,toler,rho)
}
\arguments{
  \item{X}{the design matrix}
  \item{y}{response variable}
  \item{tau}{ quantile level}
  \item{method}{"mm" for majorize and minimize method,"cd" for coordinate descent method, "admm" for Alternating method of mulipliers method,"ip" for interior point mehod}
  \item{lambda}{The constant coefficient of penalty function. (default lambda=1)}
  \item{rho}{augmented Lagrangian parameter}
  \item{beta}{initial value of estimate coefficient (default naive guess by least square estimation) }
  \item{maxit}{maxim iteration (default 200)}
  \item{toler}{the tolerance critical for stop the algorithm (default 1e-3)}
}  
\value{
  a \code{\link{list}} structure is with components
  \item{beta}{the vector of estimated coefficient}
  \item{b}{intercept}
}
\note{
  qrfit.lasso(x,y,tau) work properly only if the least square estimation is good. 
  Interior point method is done by quantreg.
}
\keyword{Quantile regression}
