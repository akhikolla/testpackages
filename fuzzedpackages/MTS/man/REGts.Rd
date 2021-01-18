\name{REGts}
\alias{REGts}
\title{Regression Model with Time Series Errors
}
\description{Perform the maximum likelihood estimation of a multivariate 
linear regression model with time-series errors
}
\usage{
REGts(zt, p, xt, include.mean = T, fixed = NULL, par = NULL, se.par = NULL, details = F)
}
\arguments{
  \item{zt}{A T-by-k data matrix of a k-dimensional time series
}
  \item{p}{The VAR order
}
  \item{xt}{A T-by-v data matrix of independent variables, where 
v denotes the number of independent variables (excluding constant 1).
}
  \item{include.mean}{A logical switch to include the constant term. Default is to include the constant term.
}
  \item{fixed}{A logical matrix used to set parameters to zero
}
  \item{par}{Initial parameter estimates of the beta coefficients, if any.
}
  \item{se.par}{Standard errors of the parameters in par, if any.
}
  \item{details}{A logical switch to control the output
}
}
\details{Perform the maximum likelihood estimation of a multivariate linear 
regression model with time series errors. Use multivariate linear regression 
to obtain initial estimates of regression coefficients if not provided
}
\value{
\item{data}{The observed k-dimensional time series}
\item{xt}{The data matrix of independent variables}
\item{aror}{VAR order}
\item{include.mean}{Logical switch for the constant vector}
\item{Phi}{The VAR coefficients}
\item{se.Phi}{The standard errors of Phi coefficients}
\item{beta}{The regression coefficients}
\item{se.beta}{The standard errors of beta}
\item{residuals}{The residual series}
\item{Sigma}{Residual covariance matrix}
\item{coef}{Parameter estimates, to be used in model simplification.}
\item{se.coef}{Standard errors of parameter estimates}
}
\references{Tsay (2014, Chapter 6). Multivariate Time Series Analysis with R and 
Financial Applications. John Wiley. Hoboken NJ.
}
\author{Ruey S. Tsay
}




