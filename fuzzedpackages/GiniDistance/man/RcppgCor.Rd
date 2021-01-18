\name{RcppgCor}
\title{Gini Distance Correlation Statistics}
\alias{RcppgCor}
\description{
  Computes Gini distance correlation statistics,
  in which Xs are quantitative, Y are categorical, alpha is exponent on the Euclidean distance and returns the measures of dependence.
}
\usage{
  RcppgCor(x, y, alpha)
}
\arguments{
  \item{x}{ data }
  \item{y}{ label of data or univariate response variable}
  \item{alpha}{ exponent on Euclidean distance, in (0,2]}
}
\details{
  \code{RcppgCor} compute Gini distance correlation statistic between \code{x} and \code{y}.
  It is a Rcpp version of \link{gCor}.  
}
\value{
  \code{RcppgCor} returns the sample Gini distance correlation
}
\seealso{
  \code{\link{RcppKgCov}}  \code{\link{RcppKgCor}}  \code{\link{RcppgCov}}
}
\examples{
  x<-iris[,1:4]
  y<-unclass(iris[,5])
  RcppgCor(x, y, alpha=2)
}
\keyword{ multivariate }
\concept{ independence }
\concept{ Gini distance correlation }
\concept{ Gini distance covariance }
\concept{ Gini statistics }
