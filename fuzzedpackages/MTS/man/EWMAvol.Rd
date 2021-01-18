\name{EWMAvol}
\alias{EWMAvol}
\title{Exponentially Weighted Moving-Average Volatility
}
\description{Use exponentially weighted moving-average method to 
compute the volatility matrix
}
\usage{
EWMAvol(rtn, lambda = 0.96)
}
\arguments{
  \item{rtn}{A T-by-k data matrix of k-dimensional asset returns, assuming the 
mean is zero
}
  \item{lambda}{Smoothing parameter. The default is 0.96. If lambda is negative, then the multivariate Gaussian likelihood is used to estimate the smoothing parameter.
}
}
\value{
\item{Sigma.t}{The volatility matrix with each row representing a volatility matrix}
\item{return}{The data}
\item{lambda}{The smoothing parameter lambda used}
}
\references{Tsay (2014, Chapter 7). Multivariate Time Series Analysis with 
R and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\examples{
data("mts-examples",package="MTS")
rtn=log(ibmspko[,2:4]+1)
m1=EWMAvol(rtn)
}




