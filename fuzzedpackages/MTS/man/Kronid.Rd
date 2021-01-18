\name{Kronid}
\alias{Kronid}
\title{Kronecker Index Identification
}
\description{Find the Kronecker indices of a k-dimensional time series
}
\usage{
Kronid(x, plag = 5, crit = 0.05)
}
\arguments{
  \item{x}{Data matrix (T-by-k) of a k-dimensional time series
}
  \item{plag}{The number of lags used to represent the past vector. Default 
is 5.
}
  \item{crit}{Type-I error used in testing for zero canonical correlations. 
Deafult is 0.05.
}
}
\value{
\item{index}{Kronecker indices}
\item{tests}{Chi-square test statistics}
}
\references{Tsay (2014, Chapter 4). Multivariate Time Series Analysis with 
R and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\examples{
phi=matrix(c(0.2,-0.6,.3,1.1),2,2); sigma=diag(2); theta=-0.5*sigma
m1=VARMAsim(300,arlags=c(1),malags=c(1),phi=phi,theta=theta,sigma=sigma)
zt=m1$series
Kronid(zt)
}




