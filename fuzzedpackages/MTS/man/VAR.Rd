\name{VAR}
\alias{VAR}
\title{Vector Autoregressive Model
}
\description{Perform least squares estimation of a VAR model
}
\usage{
VAR(x, p = 1, output = T, include.mean = T, fixed = NULL)
}
\arguments{
  \item{x}{A T-by-k matrix of k-dimensional time series
}
  \item{p}{Order of VAR model. Default is 1.
}
  \item{output}{A logical switch to control output. Default is with output.
}
  \item{include.mean}{A logical switch. It is true if mean vector is estimated.
}
  \item{fixed}{A logical matrix used in constrained estimation. 
It is used mainly in model simplification, e.g., removing insignificant 
estimates.
}
}
\details{To remove insignificant estimates, one specifies a threshold 
for individual t-ratio. The fixed matrix is then defined automatically to 
identify those parameters for removal. 
}
\value{
\item{data }{Observed data}
\item{cnst }{A logical switch to include the mean constant vector}
\item{order }{VAR order}
\item{coef }{Coefficient matrix}
\item{aic,bic,hq}{Information criteria of the fitted model}
\item{residuals}{Residuals}
\item{secoef}{Standard errors of the coefficients to be used in model refinement}
\item{Sigma }{Residual covariance matrix}
\item{Phi }{AR coefficient polynomial}
\item{Ph0 }{The constant vector}
}
\references{Tsay (2014, Chapter 3). Multivariate Time Series Analysis with 
R and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\seealso{refVAR command
}
\examples{
data("mts-examples",package="MTS")
gdp=log(qgdp[,3:5])
zt=diffM(gdp)
m1=VAR(zt,p=2)
}




