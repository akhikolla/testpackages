\name{QR.lasso.mm}
\alias{QR.lasso.mm}
\title{Quantile Regression (QR) with Adaptive Lasso Penalty (lasso) use Majorize and Minimize (mm) algorithm}
\description{
   The adaptive lasso parameter base on the estimated coefficient without penalty function.
   The algorithm majorizing the objective function by a quadratic function followed by minimizing that quadratic.
}
\usage{
QR.lasso.mm(X,y,tau,lambda,beta,maxit,toler)
}
\arguments{
  \item{X}{the design matrix.}
  \item{y}{response variable.}
  \item{tau}{quantile level.}
  \item{lambda}{The constant coefficient of penalty function. (default lambda=1)}
  \item{beta}{initial value of estimate coefficient.(default naive guess by least square estimation) }
  \item{maxit}{maxim iteration. (default 200)}
  \item{toler}{the tolerance critical for stop the algorithm. (default 1e-3)}
 } 
\value{
  a \code{\link{list}} structure is with components
  \item{beta}{the vector of estimated coefficient}
  \item{b}{intercept}
}
\note{
  QR.lasso.mm(x,y,tau) work properly only if the least square estimation is good. 
}
\references{
 David R.Hunter and Runze Li.(2005) Variable Selection Using MM Algorithms,\emph{The Annals of Statistics} \bold{33}, Number 4, Page 1617--1642.

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
# x is 1000*20 matrix, y is 1000*1 vector, beta is 20*1 vector with last ten zero value elements. 
QR.lasso.mm(x,y,0.1)
}
\keyword{MM}
