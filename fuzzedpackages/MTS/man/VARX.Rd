\name{VARX}
\alias{VARX}
\title{VAR Model with Exogenous Variables
}
\description{Estimation of a VARX model
}
\usage{
VARX(zt, p, xt = NULL, m = 0, include.mean = T, fixed = NULL, output = T)
}
\arguments{
  \item{zt}{A T-by-k data matrix of a k-dimensional time series
}
  \item{p}{The VAR order
}
  \item{xt}{A T-by-kx data matrix of kx exogenous variables
}
  \item{m}{The number of lags of exogenous variables
}
  \item{include.mean}{A logical switch to include the constant vector. 
Default is to include the constant.
}
  \item{fixed}{A logical matrix for setting parameters to zero.
}
  \item{output}{A logical switch to control output
}
}
\details{Performs least squares estimation of a VARX(p,s) model
}
\value{
\item{data}{The observed time series}
\item{xt}{The data matrix of explanatory variables}
\item{aror}{VAR order}
\item{m}{The number of lags of explanatory variables used}
\item{Ph0}{The constant vector}
\item{Phi}{VAR coefficient matrix}
\item{beta}{The regression coefficient matrix}
\item{residuals}{Residual series}
\item{coef}{The parameter estimates to be used in model simplification}
\item{se.coef}{Standard errors of the parameter estimates}
\item{include.mean}{A logical switch to include the mean vector}
}
\references{Tsay (2014, Chapter 6). Multivariate Time Series Analysis with 
R and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}




