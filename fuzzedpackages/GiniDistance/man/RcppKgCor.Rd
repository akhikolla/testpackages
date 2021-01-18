\name{RcppKgCor}
\title{ Kernel Gini Distance Correlation Statistics}
\alias{RcppKgCor}
\description{
  Computes Kernel Gini distance correlation statistics,
  in which Xs are quantitative, Y are categorical, sigma is kernel standard deviation and return the kernel Gini mean difference.
}
\usage{
  RcppKgCor(x, y, sigma)
}
\arguments{
  \item{x}{ data }
  \item{y}{ label of data or univariate response variable}
  \item{sigma}{ kernel standard deviation}
}
\details{
  \code{RcppKgCor} compute kernel Gini distance correlation statistics for data.
  It is Rcpp version of \link{KgCor}.  
}
\value{
  \code{RcppKgCor} returns the sample Kernel Gini distance covariance
}
\seealso{
  \code{\link{gCov}}  \code{\link{gCor}}  \code{\link{dCor}}
}
\examples{
  n=100
  x<-runif(n)
  y<-c(rep(1,n/2),rep(2,n/2))
  RcppKgCor(x, y, sigma=1)
}
\keyword{ multivariate }
\concept{ independence }
\concept{ Gini distance correlation }
\concept{ Gini distance covariance }
\concept{ Gini statistics }
