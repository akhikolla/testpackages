\name{QR.cd}
\alias{QR.cd}
\title{Quantile Regression (QR) use Coordinate Descent (cd) Algorithms}
\description{
   The algorithm base on greedy coordinate descent and Edgeworth's for ordinary \eqn{l_1} regression. 
}
\usage{
QR.cd(X,y,tau,beta,maxit,toler)
}
\arguments{
  \item{X}{the design matrix}
  \item{y}{response variable}
  \item{tau}{quantile level}
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
  QR.cd(x,y,tau) work properly only if the least square estimation is good. 
}
\references{
 Wu, T.T. and Lange, K. (2008). Coordinate Descent Algorithms for Lasso Penalized Regression. \emph{Annals of Applied Statistics}, \bold{2}, No 1, 224--244.

Koenker, Roger. \emph{Quantile Regression}, New York, 2005. Print.
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
# x is 1000*10 matrix, y is 1000*1 vector, beta is 10*1 vector
QR.cd(x,y,0.1)
}
\keyword{CD}
