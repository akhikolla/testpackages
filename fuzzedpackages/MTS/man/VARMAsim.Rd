\name{VARMAsim}
\alias{VARMAsim}
\title{Generating a VARMA Process
}
\description{Performs simulation of a given VARMA model
}
\usage{
VARMAsim(nobs, arlags = NULL, malags = NULL, 
  cnst = NULL, phi = NULL, theta = NULL, 
  skip = 200, sigma)
}
\arguments{
  \item{nobs}{Sample size
}
  \item{arlags}{The exact lags of the VAR matrix polynomial.
}
  \item{malags}{The exact lags of the VMA matrix polynomial.
}
  \item{cnst}{Constant vector, Phi0
}
  \item{phi}{Matrix of VAR coefficient matrices in the order of the given arlags.
}
  \item{theta}{Matrix of VMA coefficient matrices in the order of the given 
malags.
}
  \item{skip}{The number of initial data to be omitted. Default is 200.
}
  \item{sigma}{Covariance matrix (k-by-k, positive definite) 
of the innovations
}
}
\details{Use multivariate Gaussian distribution to generate random shocks. 
Then, generate a given VARMA model. The first skip data points were discarded.
}
\value{
\item{series}{Generated series}
\item{noises}{The noise series}
}
\references{Tsay (2014, Chapter 3). Multivariate Time Series Analysis with R and 
Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\examples{
p1=matrix(c(0.2,-0.6,0.3,1,1),2,2)
sig=matrix(c(4,0.8,0.8,1),2,2)
th1=matrix(c(-0.5,0,0,-0.6),2,2)
m1=VARMAsim(300,arlags=c(1),malags=c(1),phi=p1,theta=th1,sigma=sig)
zt=m1$series
}




