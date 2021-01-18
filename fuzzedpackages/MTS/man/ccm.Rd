\name{ccm}
\alias{ccm}
\title{Cross-Correlation Matrices
}
\description{Computes sample cross-correlation matrices of 
a multivariate time series, including simplified ccm matrix 
and p-value plot of Ljung-Box statistics.
}
\usage{
ccm(x, lags = 12, level = FALSE, output = T)
}
\arguments{
  \item{x}{A matrix of vector time series, each column represents a series.
}
  \item{lags}{
The number of lags of CCM to be computed. Default is 12.
}
  \item{level}{
A logical switch. When level=T, numerical values of CCM is printed. 
Default is no printing of CCM.
}
  \item{output}{
A logical switch. If ouput=F, no output is given. Default is with output.
}
}
\details{
The p-value of Ljung-Box statistics does not include any adjustment 
in degrees of freedom.
}
\value{
\item{ccm}{Sample cross-correlation matrices}
\item{pvalue}{p-values for each lag of CCM being a zero matrix}
}
\references{Tsay (2014, Chapter 1). Multivariate Time Series Analysis with R and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\examples{
xt=matrix(rnorm(1500),500,3)
ccm(xt)
ccm(xt,lag=20)
}




