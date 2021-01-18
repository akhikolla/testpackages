\name{RcppGmd}
\title{ Gini Mean Difference Statistics}
\alias{RcppGmd}
\description{
  Computes Gini mean difference of x, where alpha is an exponent on the Euclidean distance and return the Gini mean difference. The default value for alpha is 1. 
}

\usage{
  RcppGmd(x, alpha)
}
\arguments{
  \item{x}{ data }
  \item{alpha}{ exponent on Euclidean distance, in (0,2]}
}
\details{
  \code{RcppGmd} compute Gini mean difference statistics for data.
  It is a Rcpp version of \link{gmd}. 
}
\value{
  \code{RcppGmd} returns the sample Gini mean difference of \code{x}.
}
\seealso{
  \code{\link{RcppKgCov}}  \code{\link{RcppgCor}}   \code{\link{gCov}}  \code{\link{gCor}}
}
\examples{
  n=1000
  x<-runif(n)
  RcppGmd(x, alpha=1)
}
\keyword{ multivariate }
\concept{ independence }
\concept{ Gini distance correlation }
\concept{ Gini distance covariance }
\concept{ Gini statistics }
