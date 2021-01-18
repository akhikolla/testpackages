\name{qrfit}
\alias{qrfit}
\title{ Quantile Regression (qr) model fitting}
\description{
   High level function for estimating parameters by quantile regression 
}
\usage{
qrfit(X,y,tau,beta,method,maxit,toler,rho)
}
\arguments{
  \item{X}{the design matrix}
  \item{y}{response variable}
  \item{tau}{ quantile level}
  \item{method}{"mm" for majorize and minimize method,"cd" for coordinate descent method, "admm" for Alternating method of mulipliers method,"ip" for interior point mehod}
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
  qrfit(x,y,tau) work properly only if the least square estimation is good. 
  Interior point method is done by quantreg.
}
\keyword{Quantile regression}
