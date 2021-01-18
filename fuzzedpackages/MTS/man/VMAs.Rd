\name{VMAs}
\alias{VMAs}
\title{VMA Model with Selected Lags
}
\description{Performs the conditional maximum likelihood estimation 
of a VMA model with selected lags in the model
}
\usage{
VMAs(da, malags, include.mean = T, fixed = NULL, prelim = F, details = F, thres = 2)
}
\arguments{
  \item{da}{A T-by-k matrix of a k-dimensional time series with T observations
}
  \item{malags}{A vector consisting of non-zero MA lags
}
  \item{include.mean}{A logical switch to include the mean vector
}
  \item{fixed}{A logical matrix to fix coefficients to zero
}
  \item{prelim}{A logical switch concerning initial estimation
}
  \item{details}{A logical switch to control output level
}
  \item{thres}{A threshold value for setting coefficient estimates to zero
}
}
\details{A modified version of VMA model by allowing the user to select 
non-zero MA lags 
}
\value{
\item{data}{The observed time series}
\item{MAlags}{The VMA lags}
\item{cnst}{A logical switch to include the mean vector}
\item{coef}{The parameter estimates}
\item{secoef}{The standard errors of the estimates}
\item{residuals}{Residual series}
\item{aic,bic}{The information criteria of the fitted model}
\item{Sigma}{Residual covariance matrix}
\item{Theta}{The VMA matrix polynomial}
\item{mu}{The mean vector}
\item{MAorder}{The VMA order}
}
\references{Tsay (2014, Chapter 3). Multivariate Time Series 
Analysis with R and Financial Applications. 
John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\seealso{VMA
}




