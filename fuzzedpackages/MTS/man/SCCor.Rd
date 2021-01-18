\name{SCCor}
\alias{SCCor}
\title{Sample Constrained Correlations
}
\description{Compute the sample constrained correlation matrices 
}
\usage{
SCCor(rt,end,span,grp)
}
\arguments{
  \item{rt}{A T-by-k data matrix of a k-dimensional time series
}
 \item{end}{The time index of the last data point to be used in computing the sample 
correlations.}
 \item{span}{The size of the data span used to compute the correlations.}
 \item{grp}{A vector of group sizes. The time series in the same group are pooled to compute the 
 correlation matrix.}
}
\value{
\item{unconCor}{Un-constrained sample correlation matrix}
\item{conCor}{Constrained sample correlation matrix}
}
\author{Ruey S. Tsay
}
\note{This is an internal function, not intended to be a general command
}
\examples{
rt=matrix(rnorm(1000),200,5)
grp=c(3,2)
m1=SCCor(rt,200,200,grp)
m1$unconCor
m1$conCor
}




