\name{VARMA}
\alias{VARMA}
\title{Vector Autoregressive Moving-Average Models
}
\description{Performs conditional maximum likelihood estimation of 
a VARMA model. Multivariate Gaussian likelihood function is used.
}
\usage{
VARMA(da, p = 0, q = 0, include.mean = T, 
    fixed = NULL, beta=NULL, sebeta=NULL, 
    prelim = F, details = F, thres = 2)
}
\arguments{
  \item{da}{Data matrix (T-by-k) of a k-dimensional time series with sample 
size T.
}
  \item{p}{AR order
}
  \item{q}{MA order
}
  \item{include.mean}{A logical switch to control estimation of the mean vector.
Default is to include the mean in estimation.
}
  \item{fixed}{A logical matrix to control zero coefficients in estimation. 
It is mainly used by the command refVARMA.
}
 \item{beta}{Parameter estimates to be used in model simplification, if needed
}
 \item{sebeta}{Standard errors of parameter estimates for use in model simplification
}
  \item{prelim}{A logical switch to control preliminary estimation. 
Default is none.
}
  \item{details}{A logical switch to control the amount of output.
}
  \item{thres}{A threshold used to set zero parameter constraints based on 
individual t-ratio.
Default is 2.
}
}
\details{The fixed command is used for model refinement
}
\value{
\item{data}{Observed data matrix}
\item{ARorder}{VAR order}
\item{MAorder}{VMA order}
\item{cnst}{A logical switch to include the mean vector}
\item{coef}{Parameter estimates}
\item{secoef}{Standard errors of the estimates}
\item{residuals}{Residual matrix}
\item{Sigma}{Residual covariance matrix}
\item{aic,bic}{Information criteria of the fitted model}
\item{Phi}{VAR coefficients}
\item{Theta}{VMA coefficients}
\item{Ph0}{The constant vector}
}
\references{Tsay (2014, Chapter 3). Multivariate Time Series Analysis with 
R and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\seealso{refVARMA
}
\examples{
phi=matrix(c(0.2,-0.6,0.3,1.1),2,2); theta=matrix(c(-0.5,0,0,-0.5),2,2)
sigma=diag(2)
m1=VARMAsim(300,arlags=c(1),malags=c(1),phi=phi,theta=theta,sigma=sigma)
zt=m1$series
m2=VARMA(zt,p=1,q=1,include.mean=FALSE)
}




