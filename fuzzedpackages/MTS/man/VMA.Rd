\name{VMA}
\alias{VMA}
\title{Vector Moving Average Model
}
\description{Performs VMA estimation using the conditional multivariate 
Gaussian likelihood function
}
\usage{
VMA(da, q = 1, include.mean = T, fixed = NULL, 
    beta=NULL, sebeta=NULL, prelim = F, 
    details = F, thres = 2)
}
\arguments{
  \item{da}{Data matrix of a k-dimensional VMA process with each column 
containing one time series
}
  \item{q}{The order of VMA model
}
  \item{include.mean}{A logical switch to include the mean vector. 
The default is to include the mean vector in estimation.
}
  \item{fixed}{A logical matrix used to fix parameter to zero
}
  \item{beta}{Parameter estimates for use in model simplification
}
 \item{sebeta}{Standard errors of parameter estimates for use in model 
simplification
}
  \item{prelim}{A logical switch to select parameters to be included in estimation
}
  \item{details}{A logical switch to control the amount of output
}
  \item{thres}{Threshold for t-ratio used to fix parameter to zero. Default is 2.
}
}
\value{
\item{data}{The data of the observed time series}
\item{MAorder}{The VMA order}
\item{cnst}{A logical switch to include the mean vector}
\item{coef}{Parameter estimates}
\item{secoef}{Standard errors of the parameter estimates}
\item{residuals}{Residual series}
\item{Sigma}{Residual covariance matrix}
\item{Theta}{The VAR coefficient matrix}
\item{mu}{The constant vector}
\item{aic,bic}{The information criteria of the fitted model}
}
\references{Tsay (2014, Chapter 3).
}
\author{Ruey S. Tsay
}
\examples{
theta=matrix(c(0.5,0.4,0,0.6),2,2); sigma=diag(2)
m1=VARMAsim(200,malags=c(1),theta=theta,sigma=sigma)
zt=m1$series
m2=VMA(zt,q=1,include.mean=FALSE)
}




