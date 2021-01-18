\name{Eccm}
\alias{Eccm}
\title{Extended Cross-Correlation Matrices 
}
\description{Compute the extended cross-correlation matrices and the 
associated two-way table of p-values of multivariate Ljung-Box 
statistics of a vector time series.
}
\usage{
Eccm(zt, maxp = 5, maxq = 6, include.mean = FALSE, rev = TRUE)
}
\arguments{
  \item{zt}{Data matrix (T-by-k) of a vector time series, where T is the sample size and 
k is the dimension. 
}
  \item{maxp}{Maximum AR order entertained. Default is 5.
}
  \item{maxq}{Maximum MA order entertained. Default is 6.
}
  \item{include.mean}{A logical switch controlling the mean vector in estimation. 
Default assumes zero mean.
}
  \item{rev}{A logical switch to control the cross-correlation matrices used 
to compute the multivariate Ljung-Box statistics. Traditional way is to 
compute test statistics from lag-1 to lag-m. If rev = TRUE, then the 
test statistics are compute from lag-(m-1) to lag-m, from lag-(m-2) to lag-m, etc.
}
}
\value{
\item{pEccm}{A two-way table of the p-values of extended cross-correlation 
matrices}
\item{vEccm}{The sample extended cross-correlation matrices}
\item{ARcoef}{AR coefficient matrices of iterated VAR fitting}
}
\references{Tsay (2014, Chapter 3). Multivariate Time Series Analysis with 
R and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\examples{
zt=matrix(rnorm(900),300,3)
m1=Eccm(zt)
}




