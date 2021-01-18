\name{VMAe}
\alias{VMAe}
\title{VMA Estimation with Exact likelihood
}
\description{Estimation of a VMA(q) model using the exact likelihood method. 
Multivariate Gaussian likelihood function is used.
}
\usage{
VMAe(da, q = 1, include.mean = T, coef0 = NULL, 
    secoef0 = NULL, fixed = NULL, prelim = F, 
    details = F, thres = 2)
}
\arguments{
  \item{da}{Data matrix (T-by-k) for a k-dimensional VMA process
}
  \item{q}{The order of a VMA model
}
  \item{include.mean}{A logical switch to include the mean vector in estimation. Default is to include the mean vector.
}
  \item{coef0}{Initial estimates of the coefficients used mainly in model 
refinement
}
  \item{secoef0}{Standard errors of the initial estimates
}
  \item{fixed}{A logical matrix to put zero parameter constraints
}
  \item{prelim}{A logical switch for preliminary estimation
}
  \item{details}{A logical switch to control output in estimation
}
  \item{thres}{The threshold value for zero parameter constraints
}
}
\value{
\item{data}{The observed time series}
\item{MAorder}{The VMA order}
\item{cnst}{A logical switch to include the mean vector}
\item{coef}{Parameter estimates}
\item{secoef}{Standard errors of parameter estimates}
\item{residuals}{Residual series}
\item{Sigma}{Residual covariance matrix}
\item{Theta}{VMA coefficient matrix}
\item{mu}{The mean vector}
\item{aic,bic}{The information criteria of the fitted model}
}
\references{Tsay (2014). Multivariate Time Series Analysis with R 
and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\seealso{VMA
}



