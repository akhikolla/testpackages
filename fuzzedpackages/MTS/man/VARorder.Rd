\name{VARorder}
\alias{VARorder}
\title{VAR Order Specification
}
\description{Computes information criteria and the sequential Chi-square 
statistics for a vector autoregressive process
}
\usage{
VARorder(x, maxp = 13, output = T)
}
\arguments{
  \item{x}{Data matrix of dimension T-by-k with T being the sample size and 
k the number of time series
}
  \item{maxp}{The maximum VAR order entertained. Default is 13.
}
  \item{output}{A logical switch to control the output. Default is to provide output
}
}
\details{For a given maxp, the command computes Akaike, Bayesian and 
Hannan-Quinn information criteria for various VAR models using the 
data from t=maxp+1 to T. It also computes the Tiao-Box sequential Chi-square 
statistics and their p-values.
}
\value{\item{aic}{Akaike information criterion}
\item{bic}{Bayesian information criterion}
\item{hq}{Hannan and Quinn information criterion}
\item{aicor, bicor, hqor}{Orders selected by various criteria}
\item{Mstat}{Chi-square test statistics}
\item{Mpv}{p-values of the Mstat}
}
\references{Tsay (2014). Multivariate Time Series Analysis with R 
and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\seealso{VARorderI
}
\examples{
data("mts-examples",package="MTS")
zt=diffM(log(qgdp[,3:5]))
VARorder(zt,maxp=8)
}




