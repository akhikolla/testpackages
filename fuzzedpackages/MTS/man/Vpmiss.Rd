\name{Vpmiss}
\alias{Vpmiss}
\title{Partial Missing Value of a VARMA Series
}
\description{Assuming that the data is only partially missing, 
this program estimates those missing values. The model is assumed to be 
known.
}
\usage{
Vpmiss(zt, piwgt, sigma, tmiss, mdx, cnst = NULL, output = T)
}
\arguments{
  \item{zt}{A T-by-k data matrix of a k-dimensional time series
}
  \item{piwgt}{pi-weights of the model in the form 
piwgt[pi0, pi1, pi2, ....]
}
  \item{sigma}{Residual covariance matrix
}
  \item{tmiss}{Time index of the partially missing data point
}
  \item{mdx}{A k-dimensional indicator with "0" denoting missing 
component and ""1" denoting observed value.
}
  \item{cnst}{Constant term of the model
}
  \item{output}{values of the partially missing data
}
}
\value{Estimates of the missing values
}
\references{Tsay (2014, Chapter 6). Multivariate Time Series 
Analysis with R and Financial Applications. 
John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\seealso{Vmiss
}
\examples{
#data("mts-examples",package="MTS")
#gdp=log(qgdp[,3:5])
#m1=VAR(gdp,1)
#piwgt=m1$Phi; cnst=m1$Ph0; Sig=m1$Sigma
#mdx=c(0,1,1)
#m2=Vpmiss(gdp,piwgt,Sig,50,mdx,cnst)
}




