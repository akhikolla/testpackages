\name{QR.admm}
\alias{QR.admm}
\title{Quantile Regression (QR) use Alternating Direction Method of Multipliers (ADMM) algorithm}
\description{
   The problem is well suited to distributed convex optimization and is based on Alternating Direction Method of Multipliers (ADMM) algorithm .
}
\usage{
QR.admm(X,y,tau,rho,beta, maxit, toler)
}
\arguments{
  \item{X}{the design matrix}
  \item{y}{response variable}
  \item{tau}{quantile level}
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
  QR.admm(x,y,tau) work properly only if the least square estimation is good. 
}
\references{
 S. Boyd, N. Parikh, E. Chu, B. Peleato and J. Eckstein.(2010) Distributed Optimization and Statistical Learning via the Alternating Direction.Method of Multipliers \emph{Foundations and Trends in Machine Learning}, \bold{3}, No.1, 1--122 

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
QR.admm(x,y,0.1)
}
\keyword{ADMM}
