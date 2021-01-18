\name{SCMid}
\alias{SCMid}
\title{Scalar Component Identification
}
\description{Find the overall order of a VARMA process via the 
scalar component model approach
}
\usage{
SCMid(zt, maxp = 5, maxq = 5, h = 0, crit = 0.05, output = FALSE)
}
\arguments{
  \item{zt}{The T-by-k data matrix of a k-dimensional time series
}
  \item{maxp}{Maximum AR order entertained. Default is 5.
}
  \item{maxq}{Maximum MA order entertained. Default is 5.
}
  \item{h}{The additional past lags used in canonical correlation analysis.
Default is 0.
}
  \item{crit}{Type-I error of the chi-square tests used.
}
  \item{output}{A logical switch to control the output.
}
}
\value{
\item{Nmtx}{The table of the numbers of zero canonical correlations}
\item{DDmtx}{The diagonal difference table of the number of zero canonical correlations}
}
\references{Tsay (2014, Chapter 4). Multivariate Time Series Analysis with R and Financial 
Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\examples{
phi=matrix(c(0.2,-0.6,0.3,1.1),2,2); sigma=diag(2)
m1=VARMAsim(300,arlags=c(1),phi=phi,sigma=sigma)
zt=m1$series
m2=SCMid(zt)
}




