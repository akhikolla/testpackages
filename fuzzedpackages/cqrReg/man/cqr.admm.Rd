\name{cqr.admm}
\alias{cqr.admm}
\title{Composite Quantile regression (cqr) use Alternating Direction Method of Multipliers (ADMM) algorithm.}
\description{
   Composite quantile regression (cqr) find the estimated coefficient which minimize the absolute error for various quantile level.
   The problem is well suited to distributed convex optimization and is based on Alternating Direction Method of Multipliers (ADMM) algorithm .
}
\usage{
cqr.admm(X,y,tau,rho,beta, maxit, toler)
}
\arguments{
  \item{X}{the design matrix}
  \item{y}{response variable}
  \item{tau}{vector of quantile level}
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
  cqr.admm(x,y,tau) work properly only if the least square estimation is good. 
}
\references{
 S. Boyd, N. Parikh, E. Chu, B. Peleato and J. Eckstein.(2010) Distributed Optimization and Statistical Learning via the Alternating Direction. Method of Multipliers \emph{Foundations and Trends in Machine Learning}, \bold{3}, No. 1, 1--122 

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
cqr.admm(x,y,tau)
}
\keyword{ADMM}
