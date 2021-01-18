\name{sVARMA}
\alias{sVARMA}
\title{Seasonal VARMA Model Estimation
}
\description{Performs conditional maximum likelihood estimation of 
a seasonal VARMA model
}
\usage{
sVARMA(da, order, sorder, s, include.mean = T, fixed = NULL, details = F, switch = F)
}
\arguments{
  \item{da}{A T-by-k data matrix of a k-dimensional seasonal time series
}
  \item{order}{Regular order (p,d,q) of the model
}
  \item{sorder}{Seasonal order (P,D,Q) of the model
}
  \item{s}{Seasonality. s=4 for quarterly data and s=12 for monthly series
}
  \item{include.mean}{A logical switch to include the mean vector. Default is 
to include the mean
}
  \item{fixed}{A logical matrix to set zero parameter constraints
}
  \item{details}{A logical switch for output
}
  \item{switch}{A logical switch to exchange the ordering of the regular 
and seasonal VMA factors. Default is theta(B)*Theta(B).
}
}
\details{Estimation of a seasonal VARMA model
}
\value{
\item{data}{The data matrix of the observed k-dimensional time series}
\item{order}{The regular order (p,d,q)}
\item{sorder}{The seasonal order (P,D,Q)}
\item{period}{Seasonality}
\item{cnst}{A logical switch for the constant term}
\item{ceof}{Parameter estimates for use in model simplification}
\item{secoef}{Standard errors of the parameter estimates}
\item{residuals}{Residual series}
\item{Sigma}{Residual covariance matrix}
\item{aic,bic}{Information criteria of the fitted model}
\item{regPhi}{Regular AR coefficients, if any}
\item{seaPhi}{Seasonal AR coefficients}
\item{regTheta}{Regular MA coefficients}
\item{seaTheta}{Seasonal MA coefficients}
\item{Ph0}{The constant vector, if any}
\item{switch}{The logical switch to change the ordering of matrix product}
}
\references{Tsay (2014, Chapter 6). Multivariate Time Series Analysis with R 
and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}




