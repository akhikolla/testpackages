\name{VARXorder}
\alias{VARXorder}
\title{VARX Order Specification
}
\description{Specifies the orders of a VARX model, including AR order 
and the number of lags of exogenous variables
}
\usage{
VARXorder(x, exog, maxp = 13, maxm = 3, output = T)
}
\arguments{
  \item{x}{A T-by-k data matrix of a k-dimensional time series
}
  \item{exog}{A T-by-v data matrix of exogenous variables
}
  \item{maxp}{The maximum VAR order entertained
}
  \item{maxm}{The maximum lags of exogenous variables entertained
}
  \item{output}{A logical switch to control output
}
}
\details{Computes the information criteria of a VARX process
}
\value{
\item{aic}{Akaike information criterion}
\item{aicor}{Order selected by AIC}
\item{bic}{Bayesian information criterion}
\item{bicor}{Order selected by BIC}
\item{hq}{Hannan and Quinn information criterion}
\item{hqor}{Order selected by hq}
}
\references{Tsay (2014, Chapter 6). Multivariate Time Series Analysis with R 
and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}




