\name{RcppKGmd}
\title{ Kernel Gini Mean Difference Statistics}
\alias{RcppKGmd}
\description{
  Computes Kernel Gini mean difference of X, sigma is the kernel parameter and returns the kernel Gini mean difference.
}
\usage{
  RcppKGmd(x, sigma)
}
\arguments{
  \item{x}{ data }
  \item{sigma}{ kernel parameter for Gaussian kernel}
}
\details{
  \code{RcppKGmd} compute kernel Gini mean difference for data 
  It is Rcpp version of \link{Kgmd}. 
}
\value{
  \code{RcppKGmd} returns the sample Kernel Gini distance
}
\seealso{
  \code{\link{gmd}}  \code{\link{Kgmd}} 
}
\examples{
  x<-iris[,1]
  RcppKGmd(x, sigma=1)
}
\keyword{ multivariate }
\concept{ independence }
\concept{ Gini distance correlation }
\concept{ Gini distance covariance }
\concept{ Gini statistics }
