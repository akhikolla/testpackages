\name{BEKK11}
\alias{BEKK11}
\title{BEKK Model 
}
\description{Estimation of a BEKK(1,1) Model for a k-dimensional time series. 
Only k = 2 or 3 is available
}
\usage{
BEKK11(rt, include.mean = T, cond.dist = "normal", ini.estimates = NULL)
}
\arguments{
  \item{rt}{A T-by-k data matrix of k-dimensional asset returns
}
  \item{include.mean}{A logical switch to include a constant vector in the 
mean equation.  Default is with a constant vector.
}
  \item{cond.dist}{Conditional innovation distribution. 
Only Gaussian innovations are used in the current version.
}
  \item{ini.estimates}{Optional initial estimates.
}
}
\value{
\item{estimates}{Parameter estimates}
\item{HessianMtx}{Hessian matrix of the estimates}
\item{Sigma.t}{The multivariate volatilities, each row contains k-by-k elements of 
the volatility matrix Sigma(t)}
}
\references{Tsay (2014, Chapter 7)
}
\author{Ruey S. Tsay
}
\examples{
 #data("mts-examples",package="MTS")
 #da=ibmspko
 #rtn=log(da[,2:3]+1)
 #m1=BEKK11(rtn)
}




