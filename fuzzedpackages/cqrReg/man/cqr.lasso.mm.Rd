\name{cqr.lasso.mm}
\alias{cqr.lasso.mm}
\title{Composite Quantile Regression (cqr) with Adaptive Lasso Penalty (lasso) use Majorize and Minimize (mm) Algorithm}
\description{
   The adaptive lasso penalty parameter base on the estimated coefficient without penalty function.
   Composite quantile regression find the estimated coefficient which minimize the absolute error for various quantile level.
   The algorithm majorizing the objective function by a quadratic function followed by minimizing that quadratic. 
}
\usage{
cqr.lasso.mm(X,y,tau,lambda,beta,maxit,toler)
}
\arguments{
  \item{X}{the design matrix}
  \item{y}{response variable}
  \item{tau}{vector of quantile level}
  \item{lambda}{The constant coefficient of penalty function. (default lambda=1)}
  \item{beta}{initial value of estimate coefficient (default naive guess by least square estimation) }
  \item{maxit}{maxim iteration (default 200)}
  \item{toler}{the tolerance critical for stop the algorithm (default 1e-3)}
}  
\value{
  a \code{\link{list}} structure is with components
  \item{beta}{the vector of estimated coefficient}
  \item{b}{intercept for various quantile level}
}
\note{
  cqr.lasso.mm(x,y,tau) work properly only if the least square estimation is good. 
}
\references{
 David R.Hunter and Runze Li.(2005) Variable Selection Using MM Algorithms,\emph{The Annals of Statistics} \bold{33}, Number 4, Page 1617--1642.

 Hui Zou and Ming Yuan(2008). Composite Quantile Regression and the Oracle Model Selection Theory, \emph{The Annals of Statistics}, \bold{36}, Number 3, Page 1108--1126.
}
\examples{
set.seed(1)
n=100
p=2
a=2*rnorm(n*2*p, mean = 1, sd =1)
x=matrix(a,n,2*p)
beta=2*rnorm(p,1,1)
beta=rbind(matrix(beta,p,1),matrix(0,p,1))
y=x\%*\%beta-matrix(rnorm(n,0.1,1),n,1)
tau=1:5/6
# x is 1000*20 matrix, y is 1000*1 vector, beta is 20*1 vector with last ten zero value elements. 
cqr.lasso.mm(x,y,tau)
}
\keyword{MM}
