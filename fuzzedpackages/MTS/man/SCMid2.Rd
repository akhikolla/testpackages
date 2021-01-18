\name{SCMid2}
\alias{SCMid2}
\title{Scalar Component Model Specification II
}
\description{Provides detailed analysis of scalar component models 
for a specified VARMA model. The overall model is specified by 
SCMid.
}
\usage{
SCMid2(zt, maxp = 2, maxq = 2, h = 0, crit = 0.05, sseq = NULL)
}
\arguments{
  \item{zt}{The T-by-k data matrix of a k-dimensional time series
}
  \item{maxp}{Maximum AR order specified. Default is 2.
}
  \item{maxq}{Maximum MA order specified. Default is 2.
}
  \item{h}{The additional past lags used in canonical correlation analysis. 
Default is zero.
}
  \item{crit}{Type-I error used in testing. Default is 0.05.
}
  \item{sseq}{The search sequence for SCM components. Default sequence 
starts with AR order.
}
}
\value{
\item{Tmatrix}{The transformation matrix T}
\item{SCMorder}{The orders of SCM components}
}
\references{Tsay (2014, Chapter 4). Multivariate Time Series 
Analysis with R and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\seealso{SCMid
}
\examples{
phi=matrix(c(0.2,-0.6,0.3,1.1),2,2); sigma=diag(2)
m1=VARMAsim(300,arlags=c(1),phi=phi,sigma=sigma)
zt=m1$series
m2=SCMid2(zt)
names(m2)
}




