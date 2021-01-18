\name{VARs}
\alias{VARs}
\title{VAR Model with Selected Lags
}
\description{This is a modified version of VAR command by allowing the 
users to specify which AR lags to be included in the model. 
}
\usage{
VARs(x, lags, include.mean = T, output = T, fixed = NULL)
}
\arguments{
  \item{x}{A T-by-k data matrix of k-dimensional time series with T observations
}
  \item{lags}{A vector of non-zero AR lags. For instance, lags=c(1,3) 
denotes a VAR(3) model with Phi2 = 0.
}
  \item{include.mean}{A logical switch to include the mean vector
}
  \item{output}{A logical switch to control output
}
  \item{fixed}{A logical matrix to fix parameters to zero.
}
}
\details{Performs VAR estimation by allowing certain lag 
coefficient matrices being zero.
}
\value{
\item{data}{Observed time series data}
\item{lags}{The selected VAR lags}
\item{order}{The VAR order}
\item{cnst}{A logical switch to include the mean vector}
\item{coef}{Parameter estimates}
\item{aic,bic}{Information criteria of the fitted model}
\item{residuals}{Residual series}
\item{secoef}{Standard errors of the estimates}
\item{Sigma}{Residual covariance matrix}
\item{Phi}{VAR coefficient matrix}
\item{Ph0}{A constant vector}
}
\references{Tsay (2014, Chapter 2). Multivariate Time Series Analysis with 
R and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\seealso{VAR command
}
\examples{
data("mts-examples",package="MTS")
zt=log(qgdp[,3:5])
m1=VARs(zt,lags=c(1,2,4))
}




