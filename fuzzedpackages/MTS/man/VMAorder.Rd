\name{VMAorder}
\alias{VMAorder}
\title{VMA Order Specification
}
\description{Performs multivariate Ljung-Box tests to specify the 
order of a VMA process
}
\usage{
VMAorder(x, lag = 20)
}
\arguments{
  \item{x}{Data matrix of the observed k-dimensional time series. 
Each column represents a time series.
}
  \item{lag}{The maximum VMA order entertained. Default is 20.
}
}
\details{For a given lag, the command computes the Ljung-Box 
statistic for testing rho_j = ... = rho_lag = 0, where 
j = 1, 2, ..., lag. For a VMA(q) process, the 
Ljung-Box statistics should be significant for the first q lags, 
and insignificant thereafter.
}
\value{The Q-statistics and p-value plot
}
\references{Tsay (2014). Multivariate Time Series Analysis with R and 
Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\examples{
zt=matrix(rnorm(600),200,3)
VMAorder(zt)
}




