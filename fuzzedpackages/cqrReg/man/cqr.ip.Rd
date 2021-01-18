\name{cqr.ip}
\alias{cqr.ip}
\title{Composite Quantile Regression (cqr) use Interior Point (ip) Method}
\description{
   The function use the interior point method from quantreg to solve the quantile regression problem.
}
\usage{
cqr.ip(X,y,tau)
}
\arguments{
  \item{X}{the design matrix}
  \item{y}{response variable}
  \item{tau}{vector of quantile level}
 } 
\value{
  a \code{\link{list}} structure is with components
  \item{beta}{the vector of estimated coefficient}
  \item{b}{intercept}
}
\note{
  Need to install quantreg package from CRAN. 
}
\references{
Koenker, R. and S. Portnoy (1997).
The Gaussian Hare and the Laplacian Tortoise:
Computability of squared-error vs. absolute-error estimators, with discussion,
\emph{Statistical Science}, \bold{12}, 279-300.

Hui Zou and Ming Yuan(2008). Composite Quantile Regression and the Oracle Model Selection Theory, \emph{The Annals of Statistics}, \bold{36}, Number 3, Page 1108--1126.
}
\examples{
set.seed(1)
n=100
p=2
a=rnorm(n*p, mean = 1, sd =1)
x=matrix(a,n,p)
beta=rnorm(p,1,1)
beta=matrix(beta,p,1)
y=x\%*\%beta-matrix(rnorm(n,0.1,1),n,1)
tau=1:5/6
# x is 1000*10 matrix, y is 1000*1 vector, beta is 10*1 vector
#you should install quantreg first to run following command
#cqr.ip(x,y,tau)
}

\keyword{IP}
