\name{SCMfit}
\alias{SCMfit}
\title{Scalar Component Model Fitting
}
\description{Perform estimation of a VARMA model specified via 
the SCM approach
}
\usage{
SCMfit(da, scms, Tdx, include.mean = T, fixed = NULL, 
    prelim = F, details = F, thres = 1, ref = 0, 
    SCMpar=NULL, seSCMpar=NULL)
}
\arguments{
  \item{da}{The T-by-k data matrix of a k-dimensional time series
}
  \item{scms}{A k-by-2 matrix of the orders of SCMs
}
  \item{Tdx}{A k-dimensional vector for locating "1" of each row in the 
transformation matrix.
}
  \item{include.mean}{A logical switch to include the mean vector. 
Default is to include mean vector.
}
  \item{fixed}{A logical matrix to set parameters to zero
}
  \item{prelim}{A logical switch for preliminary estimation. 
Default is false. 
}
  \item{details}{A logical switch to control details of output
}
  \item{thres}{Threshold for individual t-ratio when setting parameters 
to zero. Default is 1.
}
  \item{ref}{A switch to use SCMmod in model specification.
}
 \item{SCMpar}{Parameter estimates of the SCM model, to be used in model refinement}
 \item{seSCMpar}{Standard errors of the parameter estimates in SCMpar}
}
\details{Perform conditional maximum likelihood estimation of a VARMA model 
specified by the scalar component model approach, including the transformation matrix.
}
\value{
\item{data}{Observed time series}
\item{SCMs}{The specified SCMs}
\item{Tdx}{Indicator vector for the transformation matrix. The length of Tdx is k.}
\item{locTmtx}{Specification of estimable parameters of the transformation matrix}
\item{locAR}{Locators for the estimable parameters of the VAR coefficients}
\item{locMA}{Locators for the estimable parameters of the VMA coefficients}
\item{cnst}{A logical switch to include the constant vector in the model}
\item{coef}{The parameter estimates}
\item{secoef}{Standard errors of the parameter estimates}
\item{residuals}{Residual series}
\item{Sigma}{Residual covariance matrix}
\item{aic,bic}{Information criteria of the fitted model}
\item{Ph0}{Estimates of the constant vector, if any}
\item{Phi}{Estimates of the VAR coefficients}
\item{Theta}{Estimates of the VMA coefficients}
}
\references{Tsay (2014, Chapter 4). Multivariate Time Series Analysis with R and Financial 
Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}




