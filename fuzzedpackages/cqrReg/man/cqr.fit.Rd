\name{cqr.fit}
\alias{cqr.fit}
\title{Composite Quantile Regression (cqr) model fitting}
\description{
   Composite quantile regression (cqr) find the estimated coefficient which minimize the absolute error for various quantile level.
  High level function for estimating parameter by composite quantile regression.
}
\usage{
cqr.fit(X,y,tau,beta,method,maxit,toler,rho)
}
\arguments{
  \item{X}{the design matrix}
  \item{y}{response variable}
  \item{tau}{vector of quantile level}
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
  cqr.fit(x,y,tau) work properly only if the least square estimation is good. 
  Interior point method is done by quantreg.
}
\keyword{Composite quantile regression}
