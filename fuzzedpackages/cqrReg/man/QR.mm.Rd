\name{QR.mm}
\alias{QR.mm}
\title{Quantile Regression (QR) use Majorize and Minimize (mm) algorithm}
\description{
   The algorithm majorizing the objective function by a quadratic function followed by minimizing that quadratic.
}
\usage{
QR.mm(X,y,tau,beta,maxit,toler)
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
  QR.mm(x,y,tau) work properly only if the least square estimation is good. 
}
\references{
 David R.Hunter and Kenneth Lange. Quantile Regression via an MM Algorithm, \emph{Journal of Computational and Graphical Statistics}, \bold{9}, Number 1, Page 60--77
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
QR.mm(x,y,0.1)
}
\keyword{MM}
