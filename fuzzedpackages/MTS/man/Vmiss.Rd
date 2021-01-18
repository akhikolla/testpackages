\name{Vmiss}
\alias{Vmiss}
\title{VARMA Model with Missing Value
}
\description{Assuming that the model is known, this program estimates 
the value of a missing data point. The whole data point is missing.
}
\usage{
Vmiss(zt, piwgt, sigma, tmiss, cnst = NULL, output = T)
}
\arguments{
  \item{zt}{A T-by-k data matrix of a k-dimensional time series
}
  \item{piwgt}{The pi-weights of a VARMA model defined as
piwgt=[pi0, pi1, pi2, ....]
}
  \item{sigma}{Positive definite covariance matrix of the innovations
}
  \item{tmiss}{Time index of the missing data point
}
  \item{cnst}{Constant term of the model
}
  \item{output}{A logical switch to control output
}
}
\details{Use the least squares method to estimate a missing data point. 
The missing is random.
}
\value{Estimates of the missing values
}
\references{Tsay (2014, Chapter 6). Multivariate Time Series Analysis with 
R and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\seealso{Vpmiss
}
\examples{
	data("mts-examples",package="MTS")
	gdp=log(qgdp[,3:5])
	m1=VAR(gdp,3)
	piwgt=m1$Phi; Sig=m1$Sigma; cnst=m1$Ph0
	m2=Vmiss(gdp,piwgt,Sig,50,cnst)
}



