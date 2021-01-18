\name{MCholV}
\alias{MCholV}
\title{Multivariate Cholesky Volatility Model
}
\description{Use Cholesky decomposition to obtain multivariate volatility models
}
\usage{
MCholV(rtn, size = 36, lambda = 0.96, p = 0)
}
\arguments{
  \item{rtn}{A T-by-k data matrix of a k-dimensional asset return series.
}
  \item{size}{The initial sample size used to start recursive least squares 
estimation
}
  \item{lambda}{The exponential smoothing parameter. Default is 0.96.
}
  \item{p}{VAR order for the mean equation. Default is 0.
}
}
\details{Use recursive least squares to perform the time-varying 
Cholesky decomposition. The least squares estimates are then smoothed 
via the exponentially weighted moving-average method with decaying rate 0.96.
University GARCH(1,1) model is used for the innovations of each 
linear regression. 
}
\value{
\item{betat}{Recursive least squares estimates of the linear transformations in 
Cholesky decomposition}
\item{bt}{The transformation residual series}
\item{Vol}{The volatility series of individual innovations}
\item{Sigma.t}{Volatility matrices}
}
\references{Tsay (2014, Chapter 7)
}
\author{Ruey S. Tsay
}
\seealso{fGarch
}




