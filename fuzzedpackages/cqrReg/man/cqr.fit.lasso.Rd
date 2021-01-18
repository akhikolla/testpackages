\name{cqr.fit.lasso}
\alias{cqr.fit.lasso}
\title{Composite Quantile Regression (cqr) with Adaptive Lasso Penalty (lasso)}
\description{
   Composite quantile regression (cqr) find the estimated coefficient which minimize the absolute error for various quantile level.
  High level function for estimating and selecting parameter by composite quantile regression with adaptive lasso penalty.
}
\usage{
cqr.fit.lasso(X,y,tau,lambda,beta,method,maxit,toler,rho)
}
\arguments{
  \item{X}{the design matrix}
  \item{y}{response variable}
  \item{tau}{vector of quantile level}
  \item{method}{"mm" for majorize and minimize method,"cd" for coordinate descent method, "admm" for Alternating method of mulipliers method}
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
  cqr.fit.lasso(x,y,tau) work properly only if the least square estimation is good. 
}
\keyword{Composite quantile regression}
