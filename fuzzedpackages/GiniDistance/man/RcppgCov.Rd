\name{RcppgCov}
\title{ Gini Distance Covariance Statistics}
\alias{RcppgCov}
\description{
  Computes Gini distance covariance statistics,
  in which Xs are quantitative, Y are categorical, alpha is an exponent on Euclidean distance and returns the measures of dependence.
}
\usage{
  RcppgCov(x, y, alpha)
}
\arguments{
  \item{x}{ data }
  \item{y}{ label of data or univariate response variable}
  \item{alpha}{ exponent on Euclidean distance, in (0,2]}
}
\details{
  \code{RcppgCov} compute Gini distance covariance statistics.
  It is Rcpp version of \link{gCov}. 
}
\value{
  \code{RcppgCov} returns the sample Gini distance covariance
}
\seealso{
  \code{\link{RcppgCor}}  \code{\link{RcppKgCov}}  \code{\link{RcppKgCor}} 
}
\examples{
  x<-iris[,1:4]
  y<-unclass(iris[,5])
  RcppgCov(x, y, alpha=2)
}
\keyword{ multivariate }
\concept{ independence }
\concept{ Gini distance correlation }
\concept{ Gini distance covariance }
\concept{ Gini statistics }
